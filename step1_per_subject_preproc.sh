#!/bin/bash
###############################################################################
# STEP 1 of 3 — Per-session Pass 1 (run FIRST)
#
# For each session this does the per-subject work that must happen BEFORE the
# cohort-level responsemean step (step2).
#
# CHANGE FROM PREVIOUS VERSION: Phase A and Phase B now run CONCURRENTLY.
#   Previously the script ran ALL of Phase A serially for every session, and
#   only then started Phase B. Since recon-all consumes nothing that Phase A
#   produces (it needs only T1w/T2w), that ordering wasted hours of wall-clock
#   before the long pole even started. Now each session's recon-all is launched
#   into the background immediately, and that session's Phase A runs in the
#   foreground alongside the recon-all jobs already churning.
#
#   Phase A (serial, foreground):  mrconvert -> dwidenoise -> dwibiascorrect
#                                  -> dwi2response (per subject, pre-upsample)
#                                  -> upsample DWI + mask to 0.8 mm
#   Phase B (parallel, background): recon-all (needed later by 5ttgen in step3)
#
# Phase A remains SERIAL across sessions on purpose: dwibiascorrect ants has no
# -nthreads control and ANTs contention under parallel load is still a live
# hypothesis for the recurring dwibiascorrect failures. Do not parallelise
# Phase A until that root cause is resolved.
#
# NOTE on grad_dev.nii.gz: it is NOT used in this FBA/CSD arm. It belongs to the
# separate DTI arm (dtifit --gradnonlin). It is intentionally not referenced here.
#
# After this completes for ALL sessions, run step2 (responsemean), then step3.
###############################################################################

set -u

# =========================  EDIT THESE — must match across all 3 scripts  ====
DATA_ROOT="/home/amandaa/data_dell/thesis60s/subjects"
OUT_ROOT="/home/amandaa/data_dell/thesis60s/processed_VX"
FS_ROOT="/home/amandaa/data_dell/freesurfer"

NTHREADS=8              # threads for MRtrix commands (Phase A, one session at a time)
MAX_PARALLEL_RECON=20    # concurrent recon-all jobs. 12 is the suggested next step
                        # up; 16 was confirmed feasible by resource checks but
                        # hold there until dwibiascorrect root cause is settled.
OPENMP_RECON=2          # -openmp threads per recon-all job (2-4 is the sweet spot)

MIN_FREE_GB=200         # abort if the output filesystem has less than this free

#If running only recon-all
SKIP_PHASE_A="${SKIP_PHASE_A:-0}"   # set to 1 to jump straight to recon-all

# Sessions are taken from the command line. Example:
#   ./step1_per_subject_preproc.sh HCA6005242_V1_MR HCA6174770_V1_MR
if [ $# -eq 0 ]; then
    echo "Usage: $0 <SESSION_ID> [SESSION_ID ...]"
    echo "Example: $0 HCA6005242_V1_MR HCA6174770_V1_MR"
    exit 1
fi
SESSIONS=("$@")
# =============================================================================

mkdir -p "$OUT_ROOT" "$FS_ROOT"
STATUS_LOG="${OUT_ROOT}/step1_status_$(date +%Y%m%d_%H%M%S).log"
touch "$STATUS_LOG"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

###############################################################################
# Preflight: disk headroom. A prior R7525 outage was caused by a full disk, so
# fail loudly up front rather than part-way through a 10-hour recon-all batch.
###############################################################################
FREE_GB=$(df -BG --output=avail "$OUT_ROOT" 2>/dev/null | tail -1 | tr -dc '0-9')
if [ -n "${FREE_GB:-}" ] && [ "$FREE_GB" -lt "$MIN_FREE_GB" ]; then
    log "ABORT: only ${FREE_GB}GB free on $(df -h --output=target "$OUT_ROOT" | tail -1), need >= ${MIN_FREE_GB}GB"
    log "       Free space or lower MIN_FREE_GB if you are certain this batch fits."
    exit 1
fi
log "Preflight OK — ${FREE_GB:-unknown}GB free on output filesystem"
log "Config: NTHREADS=$NTHREADS  MAX_PARALLEL_RECON=$MAX_PARALLEL_RECON  OPENMP_RECON=$OPENMP_RECON"
log "Peak concurrent threads ~ $(( MAX_PARALLEL_RECON * OPENMP_RECON + NTHREADS )) (recon-all + one Phase A)"
log "Status log: $STATUS_LOG"
log ""

# Runs a command, checks its exit code, and records PASS/FAIL to STATUS_LOG.
# Usage: run_step "step_name" cmd arg1 arg2 ...
# Returns the command's exit status so callers can `|| continue`.
run_step() {
    local step_name="$1"; shift
    "$@"
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "FAIL  ${sess}  ${step_name}  exit=${rc}" >> "$STATUS_LOG"
        log "  FAILED at: $step_name (exit $rc) — skipping rest of Phase A for $sess"
    fi
    return $rc
}

# Blocks until fewer than MAX_PARALLEL_RECON background jobs are running.
# Counting live jobs directly (rather than tracking a counter) stays accurate
# even though Phase A now runs between launches and jobs finish out of order.
throttle_recon() {
    while [ "$(jobs -rp | wc -l)" -ge "$MAX_PARALLEL_RECON" ]; do
        wait -n
    done
}

###############################################################################
# MAIN LOOP — Phase B launched into background, Phase A run in foreground
###############################################################################
for sess in "${SESSIONS[@]}"; do

    log "=========================================================="
    log "SESSION: $sess"
    log "=========================================================="

    T1="${DATA_ROOT}/${sess}/T1w/T1w_acpc_dc_restore.nii.gz"
    T2="${DATA_ROOT}/${sess}/T1w/T2w_acpc_dc_restore.nii.gz"

    ###########################################################################
    # PHASE B — launch recon-all FIRST, in background. Depends only on T1w/T2w.
    ###########################################################################
    if [ -f "${FS_ROOT}/${sess}/scripts/recon-all.done" ]; then
        # NOTE: .done presence alone does not prove success — the true signal is
        # "finished without error" in recon-all.log. Verify separately if unsure.
        log "PHASE B: recon-all already done — skipping"
        echo "SKIP  ${sess}  recon-all  already_done" >> "$STATUS_LOG"
    elif [ ! -f "$T1" ]; then
        log "PHASE B: MISSING $T1 — cannot run recon-all"
        echo "FAIL  ${sess}  recon-all  missing_T1" >> "$STATUS_LOG"
    else
        throttle_recon
        log "PHASE B: launching recon-all in background ($(jobs -rp | wc -l) already running)"
        (
            if [ -f "$T2" ]; then
                recon-all -s "$sess" -i "$T1" -T2 "$T2" -T2pial \
                          -sd "$FS_ROOT" -all -openmp $OPENMP_RECON
            else
                recon-all -s "$sess" -i "$T1" \
                          -sd "$FS_ROOT" -all -openmp $OPENMP_RECON
            fi
            rc=$?
            # recon-all's own .done/.error signal is more reliable than its exit
            # code, since some non-fatal warnings can still return 0.
            if [ -f "${FS_ROOT}/${sess}/scripts/recon-all.done" ]; then
                echo "PASS  ${sess}  recon-all" >> "$STATUS_LOG"
            else
                echo "FAIL  ${sess}  recon-all  exit=${rc}" >> "$STATUS_LOG"
            fi
        ) &
    fi

    ###########################################################################
    # PHASE A — diffusion Pass 1, foreground, while recon-all jobs run.
    # A `continue` here only abandons THIS session's Phase A; the recon-all job
    # launched above is already detached and keeps running regardless.
    ###########################################################################
    log "PHASE A: starting diffusion preproc"

    DIFF="${DATA_ROOT}/${sess}/T1w/Diffusion"
    DWI="${DIFF}/data.nii.gz"
    BVEC="${DIFF}/bvecs"
    BVAL="${DIFF}/bvals"
    MASK_IN="${DIFF}/nodif_brain_mask.nii.gz" # does mask sit right 

    # Sanity check inputs exist
    missing=0
    for f in "$DWI" "$BVEC" "$BVAL" "$MASK_IN"; do
        if [ ! -f "$f" ]; then
            log "  MISSING $f — skipping Phase A for $sess"
            echo "FAIL  ${sess}  missing_input  file=${f}" >> "$STATUS_LOG"
            missing=1
            break
        fi
    done
    [ "$missing" -eq 1 ] && continue

    WORK="${OUT_ROOT}/${sess}"
    mkdir -p "${WORK}/response"
    cd "$WORK" || { log "  cannot cd $WORK"; echo "FAIL  ${sess}  cd_work" >> "$STATUS_LOG"; continue; }

    # Idempotency: skip if Pass-1 output already present
    if [ -f "${WORK}/dwi_upsampled.mif" ] && [ -f "${WORK}/response/wm.txt" ]; then
        log "  Pass-1 already done — skipping"
        echo "SKIP  ${sess}  phase_a  already_done" >> "$STATUS_LOG"
        continue
    fi

    # Clear any stale dwibiascorrect scratch dir from a previous failed run.
    # Stale dwibiascorrect-tmp-* dirs are the diagnostic indicator for the
    # recurring ANTs failures — leaving them around muddies the signal.
    for stale in "${WORK}"/dwibiascorrect-tmp-*; do
        [ -d "$stale" ] && { log "  removing stale scratch: $stale"; rm -rf "$stale"; }
    done

    # 1) Convert + attach gradient table
    run_step "mrconvert" mrconvert "$DWI" raw_dwi.mif -fslgrad "$BVEC" "$BVAL" -force -nthreads $NTHREADS || continue

    # 2) Denoise (MP-PCA). Cannot be dropped for later cohorts — doing so would
    #    introduce a confound aligned with the primary variable of interest.
    run_step "dwidenoise" dwidenoise raw_dwi.mif dwi_den.mif -noise noise.mif -force -nthreads $NTHREADS || continue

    # 3) Bias field correction (topup/eddy already done by HCP)
    run_step "dwibiascorrect" dwibiascorrect ants dwi_den.mif dwi_unbiased.mif -bias bias.mif -force || continue

    # 4) Bring in the HCP mask (preferred over dwi2mask for cross-subject consistency)
    run_step "mrconvert_mask" mrconvert "$MASK_IN" mask.mif -force || continue

    # 5) Response function — on NATIVE (non-upsampled) data, per subject.
    #    These are averaged into a single cohort response set in step2.
    run_step "dwi2response" dwi2response dhollander dwi_unbiased.mif \
        response/wm.txt response/gm.txt response/csf.txt \
        -voxels response/voxels.mif -mask mask.mif -force -nthreads $NTHREADS || continue

    # 6) Upsample DWI to 0.8 mm  
    run_step "mrgrid_dwi" mrgrid dwi_unbiased.mif regrid -voxel 0.8 dwi_upsampled.mif -force -nthreads $NTHREADS || continue

    # 7) Regrid the mask to match the upsampled grid (nearest-neighbour keeps it binary)
    run_step "mrgrid_mask" mrgrid mask.mif regrid -template dwi_upsampled.mif -interp nearest \
        mask_upsampled.mif -force -nthreads $NTHREADS || continue #check community page - concensus, updated, DO MASK AGAIN 

    echo "PASS  ${sess}  phase_a" >> "$STATUS_LOG"
    log "PHASE A: done for $sess"
done

###############################################################################
# Wait for any recon-all jobs still running after the last session's Phase A
###############################################################################
STILL_RUNNING=$(jobs -rp | wc -l)
if [ "$STILL_RUNNING" -gt 0 ]; then
    log ""
    log "All Phase A work finished. Waiting on $STILL_RUNNING remaining recon-all job(s)..."
    log "(Monitor with: grep -E '^#@#' ${FS_ROOT}/<SESSION>/scripts/recon-all.log | tail)"
fi
wait

###############################################################################
# SUMMARY
###############################################################################
log ""
log "=== SUMMARY (full detail in $STATUS_LOG) ==="
log "Phase A diffusion — PASS: $(grep -c '^PASS.*phase_a'    "$STATUS_LOG")   FAIL: $(grep '^FAIL' "$STATUS_LOG" | grep -vc 'recon-all')   SKIP(cached): $(grep -c '^SKIP.*phase_a' "$STATUS_LOG")"
log "Phase B recon-all  — PASS: $(grep -c '^PASS.*recon-all' "$STATUS_LOG")   FAIL: $(grep -c '^FAIL.*recon-all' "$STATUS_LOG")   SKIP(cached): $(grep -c '^SKIP.*recon-all' "$STATUS_LOG")"

if grep -q '^FAIL' "$STATUS_LOG"; then
    log ""
    log "Failed sessions/steps:"
    grep '^FAIL' "$STATUS_LOG" | while read -r line; do log "  $line"; done
    log ""
    log "STEP 1 completed WITH FAILURES. Do NOT run step2 until every session"
    log "in this cohort has passed — responsemean needs the full set."
    exit 1
fi

log ""
log "STEP 1 complete for all sessions, no failures."
log "Next: run step2_group_responsemean.sh"