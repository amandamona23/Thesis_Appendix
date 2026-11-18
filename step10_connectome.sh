#!/bin/bash
###############################################################################
# step10_connectome.sh  -  structural connectome per session (Aim 2)
#
# Per session, in native diffusion space (HCP T1w space: no registration step, because
# the HCP diffusion data and the T1 that recon-all ran on are already aligned):
#   1. 5ttgen hsvs     5-tissue image from the FreeSurfer surfaces + aseg   (HSVS_OPTS)
#                      checked with 5ttcheck
#   2. 5tt2gmwmi       grey/white matter interface, used as the seed region
#   3. labelconvert    aparc+aseg -> 84 Desikan-Killiany nodes (MRtrix fs_default)
#   4. tckgen          iFOD2 + ACT on wmfod_norm.mif: -backtrack -crop_at_gmwmi
#                      -seed_gmwmi, -select N_TRACKS, -maxlength 250, -cutoff 0.06
#   5. tcksift2        streamline weights, -act, -out_mu
#   6. tck2connectome  sum of SIFT2 weights per node pair, -symmetric -zero_diagonal
#                      (multiply by mu at extraction: step10_extract.py)
#   7. connectome2tck  streamlines + weights for the hippocampus <-> PFC edges only
#   8. tckedit         200,000-streamline subset for figures
#   9. QC              streamline count, mu, % of streamlines reaching 2 nodes, key node sizes,
#                      time and peak memory per stage
# The 10-million-streamline tractogram, weights and assignments stay on LOCAL scratch and
# are deleted once the session is finished.
#
#   READS   ${OUT_ROOT}/<SESS>/wmfod_norm.mif (step 3a), ${FS_ROOT}/<SESS>/ (recon-all)
#   WRITES  ${CONN_DIR}/<SESS>/  5tt.mif.gz gmwmi.mif.gz nodes.mif connectome.csv sift_mu.txt
#                                edges/ tracks_200k.tck qc.tsv log.txt .done
#           ${CONN_DIR}/connectome_qc.tsv  (--summary)
#
# Usage (run in tmux):
#   ./step10_connectome.sh --list pilot.txt               # 2 sessions first: time, memory, disk
#   ./step10_connectome.sh --list sessions.txt            # all; finished sessions are skipped
#   N_TRACKS=5000000 ./step10_connectome.sh --list ...    # fallback if 10M is too slow / heavy
#   ./step10_connectome.sh --summary
# Settings: N_TRACKS (10000000)  MAX_PARALLEL (7)  NTHREADS (27)  SIFT_SLOTS (= MAX_PARALLEL)
#           HSVS_OPTS ("-hippocampi aseg -thalami aseg")  SUBSET (200000)
#           SCRATCH_ROOT ($HOME/fba10_work)  CONN_DIR ($HOME/data_dell/fba/connectome)
###############################################################################
set -u
SELF="$(readlink -f "$0")"
source "$(dirname "$SELF")/pipeline_config_v2.sh"

N_TRACKS="${N_TRACKS:-10000000}"; MAX_PAR="${MAX_PARALLEL:-7}"; NTHREADS="${NTHREADS:-27}"
SIFT_SLOTS="${SIFT_SLOTS:-$MAX_PAR}"; SUBSET="${SUBSET:-200000}"
HSVS_OPTS="${HSVS_OPTS:--hippocampi aseg -thalami aseg}"
FOD_CUTOFF="${FOD_CUTOFF:-0.06}"; MAXLEN="${MAXLEN:-250}"
FS_ROOT="${FS_ROOT:-$HOME/data_dell/freesurfer}"
SCRATCH_ROOT="${SCRATCH_ROOT:-$HOME/fba10_work}"; CONN_DIR="${CONN_DIR:-$HOME/data_dell/fba/connectome}"
export FREESURFER_HOME="${FREESURFER_HOME:-/opt/freesurfer/freesurfer-7.1.1}"
LUT="${FS_LUT:-${FREESURFER_HOME}/FreeSurferColorLUT.txt}"
FSDEF="${FS_DEFAULT:-$(dirname "$(readlink -f "$(command -v labelconvert)")")/../share/mrtrix3/labelconvert/fs_default.txt}"
# hippocampus and ipsilateral PFC nodes (fs_default, 1-based): edges kept as streamline files
EDGE_NODES="${EDGE_NODES:-11,13,25,40,47,60,62,74}"
GB_PER_M=2   # local disk per million streamlines (tractogram + weights + assignments), with margin
LOCKDIR="/tmp/${USER}_step10_locks"; LOCK="/tmp/${USER}_step10.lock"
export N_TRACKS NTHREADS SIFT_SLOTS SUBSET HSVS_OPTS FOD_CUTOFF MAXLEN FS_ROOT SCRATCH_ROOT CONN_DIR LUT FSDEF EDGE_NODES GB_PER_M LOCKDIR OUT_ROOT

hrs()      { awk -v s="$1" 'BEGIN{printf "%.2f", s/3600}'; }
avail_gb() { df -BG --output=avail "$1" 2>/dev/null | tail -1 | tr -dc '0-9'; }

###############################################################################
# One session
###############################################################################
run_session() {
    local s=$1 W="${SCRATCH_ROOT}/$1" O="${CONN_DIR}/$1"
    local FOD="${OUT_ROOT}/$1/wmfod_norm.mif" FSD="${FS_ROOT}/$1" T5="${CONN_DIR}/$1/5tt.mif.gz"
    mkdir -p "$O" "$W"; exec >> "${O}/log.txt" 2>&1
    log "start ${s}: ${N_TRACKS} streamlines, ${NTHREADS} threads"
    declare -A T M
    timed() {  # timed <label> <command...>: records hours and peak memory (GB)
        local lab=$1; shift; local t0=$(date +%s) tf="${W}/.time_${lab}" rc
        if [ -x /usr/bin/time ]; then /usr/bin/time -v -o "$tf" "$@"; rc=$?; else "$@"; rc=$?; fi
        T[$lab]=$(hrs $(( $(date +%s) - t0 )))
        M[$lab]=$(awk -F: '/Maximum resident set size/ {printf "%.1f", $2/1048576}' "$tf" 2>/dev/null)
        log "  ${lab}: exit ${rc}, ${T[$lab]} h, peak ${M[$lab]:-NA} GB"; return $rc
    }
    fail() { log "FAIL ${s}: $*"; echo "FAIL ${s}: $*"; return 1; }

    # 1-2. tissue image and seed interface
    if [ ! -f "$T5" ]; then
        rm -rf "${W}/hsvs_tmp"
        timed 5ttgen 5ttgen hsvs "$FSD" "${W}/5tt.mif" $HSVS_OPTS -scratch "${W}/hsvs_tmp" -nthreads "$NTHREADS" -force -quiet \
            || { fail "5ttgen hsvs"; return 1; }
        5ttcheck "${W}/5tt.mif" > "${O}/5ttcheck.txt" 2>&1 || log "  WARNING: 5ttcheck reported problems (see 5ttcheck.txt)"
        mrconvert "${W}/5tt.mif" "${T5}.tmp.mif.gz" -force -quiet && mv "${T5}.tmp.mif.gz" "$T5" || { fail "saving 5tt"; return 1; }
    fi
    if [ ! -f "${O}/gmwmi.mif.gz" ]; then
        5tt2gmwmi "$T5" "${O}/.gmwmi.mif.gz" -force -quiet && mv "${O}/.gmwmi.mif.gz" "${O}/gmwmi.mif.gz" || { fail "5tt2gmwmi"; return 1; }
    fi
    # 3. nodes
    if [ ! -f "${O}/nodes.mif" ]; then
        labelconvert "${FSD}/mri/aparc+aseg.mgz" "$LUT" "$FSDEF" "${O}/.nodes.mif" -force -quiet \
            && mv "${O}/.nodes.mif" "${O}/nodes.mif" || { fail "labelconvert"; return 1; }
        local mx; mx=$(mrstats "${O}/nodes.mif" -output max | tr -d ' ')
        [ "${mx%.*}" = "84" ] || log "  WARNING: highest node label is ${mx}, expected 84"
    fi

    if [ ! -f "${O}/connectome.csv" ]; then
        # 4. tractography (local scratch only); wait for disk if other sessions are using it
        if [ ! -f "${W}/.tracks.done" ]; then
            local need=$(( N_TRACKS / 1000000 * GB_PER_M + 10 ))
            until [ "$(avail_gb "$SCRATCH_ROOT")" -ge "$need" ]; do log "  waiting for ${need} GB of local disk"; sleep 300; done
            timed tckgen tckgen "$FOD" "${W}/tracks.tck" -algorithm iFOD2 -act "$T5" -backtrack -crop_at_gmwmi \
                -seed_gmwmi "${O}/gmwmi.mif.gz" -select "$N_TRACKS" -maxlength "$MAXLEN" -cutoff "$FOD_CUTOFF" \
                -nthreads "$NTHREADS" -force -quiet || { fail "tckgen"; return 1; }
            local n; n=$(tckinfo "${W}/tracks.tck" 2>/dev/null | awk -F: '/^ *count/ {gsub(/ /,"",$2); print $2; exit}')
            [ "${n:-0}" -ge "$N_TRACKS" ] 2>/dev/null || { fail "tckgen produced ${n:-0} of ${N_TRACKS}"; return 1; }
            touch "${W}/.tracks.done"
        fi
        # 5. SIFT2, at most SIFT_SLOTS at once (memory)
        if [ ! -f "${W}/.sift.done" ]; then
            local slot="" i
            mkdir -p "$LOCKDIR"
            while [ -z "$slot" ]; do
                for i in $(seq 1 "$SIFT_SLOTS"); do
                    exec {fd}>"${LOCKDIR}/sift${i}.lock"
                    if flock -n "$fd"; then slot=$i; break; fi
                    exec {fd}>&-
                done
                [ -n "$slot" ] || { log "  waiting for a SIFT2 slot"; sleep 120; }
            done
            timed tcksift2 tcksift2 "${W}/tracks.tck" "$FOD" "${W}/weights.txt" -act "$T5" -out_mu "${O}/sift_mu.txt" \
                -nthreads "$NTHREADS" -force -quiet
            local rc=$?; exec {fd}>&-
            [ "$rc" = 0 ] || { fail "tcksift2"; return 1; }
            touch "${W}/.sift.done"
        fi
        # 6. connectome (sum of SIFT2 weights; mu applied at extraction)
        timed tck2connectome tck2connectome "${W}/tracks.tck" "${O}/nodes.mif" "${O}/.connectome.csv" \
            -tck_weights_in "${W}/weights.txt" -symmetric -zero_diagonal -out_assignments "${W}/assignments.txt" \
            -nthreads "$NTHREADS" -force -quiet || { fail "tck2connectome"; return 1; }
        # 7. hippocampus <-> PFC edges as streamline files (with their weights)
        rm -rf "${O}/edges"; mkdir -p "${O}/edges"
        connectome2tck "${W}/tracks.tck" "${W}/assignments.txt" "${O}/edges/edge-" -nodes "$EDGE_NODES" -exclusive \
            -files per_edge -tck_weights_in "${W}/weights.txt" -prefix_tck_weights_out "${O}/edges/weights-" \
            -nthreads "$NTHREADS" -force -quiet || log "  WARNING: connectome2tck failed (edges not saved)"
        # 8. subset for figures
        tckedit "${W}/tracks.tck" "${O}/tracks_200k.tck" -number "$SUBSET" -force -quiet || log "  WARNING: tckedit subset failed"
        # 9. QC
        local total assigned pct mu nv k nodevols=""
        total=$(wc -l < "${W}/assignments.txt")
        assigned=$(awk '$1!=0 && $2!=0' "${W}/assignments.txt" | wc -l)
        pct=$(awk -v a="$assigned" -v t="$total" 'BEGIN{printf "%.1f", (t>0?100*a/t:0)}')
        mu=$(tr -d ' \n' < "${O}/sift_mu.txt")
        for k in 40 47 11 13 25 60 62 74; do
            mrcalc "${O}/nodes.mif" "$k" -eq "${W}/n.mif" -datatype bit -force -quiet
            nv=$(mrstats "${W}/n.mif" -mask "${W}/n.mif" -output count 2>/dev/null | tr -d ' '); nodevols+="${nv:-0}\t"
        done
        printf 'session\tn_tracks\tmu\tpct_assigned\tvox_L_hippo\tvox_R_hippo\tvox_L_LOF\tvox_L_MOF\tvox_L_RAC\tvox_R_LOF\tvox_R_MOF\tvox_R_RAC\th_5ttgen\th_tckgen\tGB_tckgen\th_sift2\tGB_sift2\th_connectome\n' > "${O}/qc.tsv"
        printf "%s\t%s\t%s\t%s\t${nodevols}%s\t%s\t%s\t%s\t%s\t%s\n" "$s" "$N_TRACKS" "$mu" "$pct" \
            "${T[5ttgen]:-done_before}" "${T[tckgen]:-done_before}" "${M[tckgen]:-NA}" "${T[tcksift2]:-done_before}" \
            "${M[tcksift2]:-NA}" "${T[tck2connectome]:-NA}" >> "${O}/qc.tsv"
        mv "${O}/.connectome.csv" "${O}/connectome.csv"
    fi
    rm -rf "$W"; touch "${O}/.done"
    log "OK ${s}"; echo "OK   ${s}"
}
if [ "${1:-}" = "--worker" ]; then run_session "$2"; exit $?; fi

summary() {
    local f; f="${CONN_DIR}/connectome_qc.tsv"
    { head -1 "$(ls "${CONN_DIR}"/HCA*/qc.tsv | head -1)"; for q in "${CONN_DIR}"/HCA*/qc.tsv; do tail -n +2 "$q"; done; } > "$f" 2>/dev/null
    echo "Finished sessions: $(ls "${CONN_DIR}"/HCA*/.done 2>/dev/null | wc -l)"
    column -t -s$'\t' "$f"
    echo; echo "Saved: ${f}"
}

###############################################################################
# Controller
###############################################################################
LIST=""; MODE=run
while [ $# -gt 0 ]; do
    case "$1" in --list) LIST="$(readlink -f "$2")"; shift ;; --summary) MODE=summary ;; *) die "Unknown option: $1" ;; esac; shift
done
[ "$MODE" = summary ] && { summary; exit 0; }
[ -f "$LIST" ] || die "Usage: $0 --list sessions.txt   (or --summary)"
exec 9>"$LOCK"; flock -n 9 || die "Step 10 is already running."
for c in 5ttgen 5tt2gmwmi 5ttcheck labelconvert tckgen tcksift2 tck2connectome connectome2tck tckedit tckinfo mrstats mrcalc mrconvert; do
    command -v "$c" >/dev/null || die "Command not found: $c"
done
[ -x /usr/bin/time ] || log "WARNING: /usr/bin/time not found; times are recorded but not peak memory"
command -v mri_convert >/dev/null || log "WARNING: FreeSurfer tools not on PATH; 5ttgen hsvs needs them (source \$FREESURFER_HOME/SetUpFreeSurfer.sh)"
[ -f "$LUT" ] || die "FreeSurfer lookup table not found: $LUT (set FS_LUT)"
[ -f "$FSDEF" ] || die "MRtrix fs_default.txt not found: $FSDEF (set FS_DEFAULT)"
mkdir -p "$CONN_DIR" "$SCRATCH_ROOT" "$LOCKDIR"

mapfile -t SESSIONS < <(grep -o 'HCA[0-9]*_V[0-9]*_MR' "$LIST" | awk '!seen[$0]++')
problems=(); QUEUE=()
for s in "${SESSIONS[@]}"; do
    [ -f "${OUT_ROOT}/${s}/wmfod_norm.mif" ] || problems+=("${s}: no wmfod_norm.mif")
    for f in mri/aparc+aseg.mgz mri/aseg.mgz surf/lh.white surf/rh.white surf/lh.pial surf/rh.pial; do
        [ -f "${FS_ROOT}/${s}/${f}" ] || problems+=("${s}: no FreeSurfer ${f}")
    done
done
if [ "${#problems[@]}" -gt 0 ]; then printf '  %s\n' "${problems[@]}"; die "${#problems[@]} problem(s) above"; fi
for s in "${SESSIONS[@]}"; do [ -f "${CONN_DIR}/${s}/.done" ] || QUEUE+=("$s"); done
log "Step 10: ${#SESSIONS[@]} sessions, ${#QUEUE[@]} to do, ${MAX_PAR} at once (${NTHREADS} threads each), ${N_TRACKS} streamlines, SIFT2 slots ${SIFT_SLOTS}"
log "Local scratch: ${SCRATCH_ROOT} ($(avail_gb "$SCRATCH_ROOT") GB free); about $(( N_TRACKS / 1000000 * GB_PER_M )) GB per running session"
[ "${#QUEUE[@]}" -gt 0 ] && printf '%s\n' "${QUEUE[@]}" | xargs -P "$MAX_PAR" -I{} bash "$SELF" --worker {}
summary
