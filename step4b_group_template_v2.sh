#!/bin/bash
###############################################################################
# step4b_group_template_v2.sh  -  Build the group FOD template (FBA arm)
#
# Runs population_template ONCE on the 40 template subjects' halfway
# templates (from step 4a), with their halfway-space masks.
#
#   INPUT   ${TEMPLATE_DIR}/subject_templates/<SUBJ>/halfway_template.mif
#           ${TEMPLATE_DIR}/subject_templates/<SUBJ>/halfway_mask.mif
#           for every subject in template40_subjects.txt
#   OUTPUT  ${TEMPLATE_DIR}/group/wmfod_template.mif        the group template
#           ${TEMPLATE_DIR}/group/template_mask_40.mif      intersection of the
#                                                            40 masks (for checking
#                                                            coverage; the final
#                                                            analysis mask comes
#                                                            from all 80 in step 4c)
#           ${TEMPLATE_DIR}/group/step4b_run.tsv            timing, peak scratch,
#                                                            peak RAM
#
# Step 4c (next) registers ALL 80 halfway templates to this template, so the 40
# used to build it are treated exactly like the other 40. No per-subject warps
# are kept from this step for that reason.
#
# SCRATCH is the big unknown. Two inputs used ~45-53 GB; if that scales with
# the number of inputs, 40 inputs need roughly 1 TB. You must choose where it
# goes (SCRATCH_ROOT). A watchdog stops the run cleanly if free space on that
# disk falls below MIN_FREE_GB, rather than letting a shared disk fill up.
# The peak scratch it records is the real number for any future run.
#
# PROTECTION: population_template runs in its own process group, so Ctrl+C in
# this window does NOT stop it. To stop it on purpose, use the kill command
# printed at the start (also saved to ${TEMPLATE_DIR}/group/HOW_TO_STOP.txt).
#
# Usage:
#   SCRATCH_ROOT=/path/with/space NTHREADS=48 ./step4b_group_template_v2.sh --list template40_subjects.txt
###############################################################################

set -u
SELF="$(readlink -f "$0")"
: "${NTHREADS:=48}"; export NTHREADS
source "$(dirname "$SELF")/pipeline_config_v2.sh"

SUBJ_TEMPLATE_DIR="${TEMPLATE_DIR}/subject_templates"
GROUP_DIR="${TEMPLATE_DIR}/group"
LOGF="${GROUP_DIR}/step4b_population_template.log"
RUNREC="${GROUP_DIR}/step4b_run.tsv"
MIN_FREE_GB="${MIN_FREE_GB:-30}"          # stop the run if scratch disk gets this low
EST_PER_INPUT_GB="${EST_PER_INPUT_GB:-26}" # from 4a: ~45-53 GB for 2 inputs
POLL_SECS="${POLL_SECS:-60}"
DU_EVERY="${DU_EVERY:-5}"                  # measure scratch size every N polls
LOCK="/tmp/${USER}_step4b.pid"

avail_gb() { df -BG --output=avail "$1" 2>/dev/null | tail -1 | tr -dc '0-9'; }
fs_of()    { df --output=source "$1" 2>/dev/null | tail -1; }
alive()    { [ -e "/proc/$1" ] && ! grep -q 'Z (zombie)' "/proc/$1/status" 2>/dev/null; }
hours()    { awk -v a="$1" -v b="$2" 'BEGIN{printf "%.2f", (b-a)/3600}'; }

# ---------------------------------------------------------------------------
# Arguments and one-at-a-time lock
# ---------------------------------------------------------------------------
ASSUME_YES=0; ACCEPT_SPACE_RISK=0; LIST=""
while [ $# -gt 0 ]; do
    case "$1" in
        --list) LIST="$(readlink -f "$2")"; shift ;;
        --yes) ASSUME_YES=1 ;;
        --accept-space-risk) ACCEPT_SPACE_RISK=1 ;;
        *) die "Unknown option: $1" ;;
    esac
    shift
done
[ -n "$LIST" ] && [ -f "$LIST" ] || die "Usage: SCRATCH_ROOT=/path NTHREADS=48 $0 --list template40_subjects.txt"
[ -n "${SCRATCH_ROOT:-}" ] || die "Set SCRATCH_ROOT to a folder on a disk with ~1 TB free (see header)."
mkdir -p "$SCRATCH_ROOT" "$GROUP_DIR" || die "Cannot create ${SCRATCH_ROOT} or ${GROUP_DIR}"
SCRATCH_ROOT="$(readlink -f "$SCRATCH_ROOT")"

if [ -f "$LOCK" ] && alive "$(cat "$LOCK")"; then
    die "Step 4b is already running (PID $(cat "$LOCK"))."
fi
echo $$ > "$LOCK"

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
mapfile -t SUBJS < <(grep -o 'HCA[0-9]*' "$LIST" | awk '!seen[$0]++')
[ "${#SUBJS[@]}" -ge 2 ] || die "Template list has ${#SUBJS[@]} subjects"

problems=()
for s in "${SUBJS[@]}"; do
    sd="${SUBJ_TEMPLATE_DIR}/${s}"
    [ -f "$sd/halfway_template.mif" ]     || problems+=("${s}: no halfway_template.mif")
    [ -f "$sd/halfway_mask.mif" ]         || problems+=("${s}: no halfway_mask.mif")
    [ -f "$sd/mask_halfway_space.done" ]  || problems+=("${s}: mask not confirmed in halfway space")
done
QCM="${TEMPLATE_DIR}/qc4a/qc4a_metrics.tsv"
qc_warn=()
if [ -f "$QCM" ]; then
    for s in "${SUBJS[@]}"; do
        line=$(awk -F'\t' -v s="$s" '$1==s' "$QCM")
        if [ -z "$line" ]; then qc_warn+=("${s}: not in QC table yet (run ./qc4a.sh)")
        elif [ "$(echo "$line" | cut -f3)" != "OK" ]; then problems+=("${s}: QC integrity $(echo "$line" | cut -f3): $(echo "$line" | cut -f4)")
        fi
    done
else
    qc_warn+=("no step 4a QC table found; run ./qc4a.sh first if you can")
fi
if [ "${#problems[@]}" -gt 0 ]; then
    printf '  %s\n' "${problems[@]}"
    die "Fix the problems above before building the template."
fi

others=$(pgrep -u "$(id -u)" -f population_template | wc -l)
runner_alive=0
[ -f "/tmp/${USER}_run4a_controller.pid" ] && alive "$(cat "/tmp/${USER}_run4a_controller.pid")" && runner_alive=1

# WHERE THE INPUTS COME FROM (INPUTS=...). population_template rereads every
# input at every stage, so they must not be read over the plain network mount.
#   copy            (default) copy them next to the scratch; costs ~100 GB there
#   shm             copy them into /dev/shm (RAM); keeps ~100 GB of RAM in use
#                   for the whole run; deleted automatically at the end
#   cached:<root>   link them through a second, caching NAS mount at <root>
#                   (a mount of the same NAS folder as ~/data_dell); no copy
INPUTS="${INPUTS:-copy}"
size_inputs() { for s in "${SUBJS[@]}"; do du -bL "${SUBJ_TEMPLATE_DIR}/$s/halfway_template.mif" "${SUBJ_TEMPLATE_DIR}/$s/halfway_mask.mif"; done | awk '{t+=$1} END{printf "%d", t/1073741824 + 1}'; }
stage_gb=0; inputs_gb=$(size_inputs); STAGE=0; IN_SRC_DIR="$SUBJ_TEMPLATE_DIR"
case "$INPUTS" in
    copy)
        IN_ROOT_4B="${SCRATCH_ROOT}/step4b_inputs"
        if [ "$(fs_of "$SCRATCH_ROOT")" != "$(fs_of "$SUBJ_TEMPLATE_DIR")" ]; then STAGE=1; stage_gb=$inputs_gb; fi
        IN_DESC=$([ "$STAGE" -eq 1 ] && echo "copied next to scratch (~${inputs_gb} GB there)" || echo "linked (scratch is on the NAS)") ;;
    shm)
        IN_ROOT_4B="/dev/shm/${USER}_step4b_inputs"; STAGE=1
        shm_free=$(avail_gb /dev/shm); mem_free=$(awk '/MemAvailable/ {printf "%d", $2/1048576}' /proc/meminfo)
        [ "${shm_free:-0}" -ge $((inputs_gb + 5)) ] || die "/dev/shm has ${shm_free} GB free; inputs need ~${inputs_gb} GB"
        [ "${mem_free:-0}" -ge $((inputs_gb + 60)) ] || die "Only ${mem_free} GB of memory available; RAM inputs need ~$((inputs_gb + 60)) GB (inputs + running jobs)"
        IN_DESC="copied into RAM (/dev/shm, ~${inputs_gb} GB; ${mem_free} GB memory available now)" ;;
    cached:*)
        CROOT="$(readlink -f "${INPUTS#cached:}")"
        case "$SUBJ_TEMPLATE_DIR" in "$HOME/data_dell"/*) ;; *) die "cached mode expects the subject templates under ~/data_dell" ;; esac
        IN_SRC_DIR="${CROOT}/${SUBJ_TEMPLATE_DIR#"$HOME/data_dell/"}"
        [ -f "${IN_SRC_DIR}/${SUBJS[0]}/halfway_template.mif" ] || die "Cannot see ${IN_SRC_DIR}/${SUBJS[0]}: is the caching mount at ${CROOT} up?"
        IN_ROOT_4B="${SCRATCH_ROOT}/step4b_inputs"
        IN_DESC="read through the caching mount ${CROOT} (no copy)" ;;
    *) die "INPUTS must be copy, shm or cached:<mount root>" ;;
esac
est_gb=$(( ${#SUBJS[@]} * EST_PER_INPUT_GB ))
need_gb=$(( est_gb + stage_gb + MIN_FREE_GB ))
free_gb=$(avail_gb "$SCRATCH_ROOT")

cat <<EOF

-------------------------------------------------------------------
Step 4b: group FOD template
  subjects:            ${#SUBJS[@]} (from $(basename "$LIST"))
  threads:             ${NTHREADS}
  scratch folder:      ${SCRATCH_ROOT}
  inputs:              ${IN_DESC}
  estimated scratch:   ~${est_gb} GB  (${EST_PER_INPUT_GB} GB per input; the real peak gets recorded)
  space needed there:  ~${need_gb} GB including a ${MIN_FREE_GB} GB safety margin
  space free:          ${free_gb} GB
  other population_template processes running: ${others}
  step 4a runner still running: $([ "$runner_alive" -eq 1 ] && echo yes || echo no)
  output:              ${GROUP_DIR}/wmfod_template.mif
-------------------------------------------------------------------
EOF
[ "${#qc_warn[@]}" -gt 0 ] && { echo "Notes:"; printf '  %s\n' "${qc_warn[@]}"; echo; }
[ "$others" -gt 0 ] && echo "WARNING: other template jobs are running and will compete for CPU, disk and network."
if [ "${free_gb:-0}" -lt "$need_gb" ]; then
    echo "WARNING: less free space than the estimate. The watchdog will stop the run if free space"
    echo "         falls below ${MIN_FREE_GB} GB, which could waste many hours of compute."
    [ "$ACCEPT_SPACE_RISK" -eq 1 ] || die "Not starting. Use a bigger disk, or add --accept-space-risk to try anyway."
fi
if [ "$ASSUME_YES" -eq 0 ]; then
    read -r -p "Start the group template? This runs for a long time. [y/N] " ok
    [[ "$ok" =~ ^[Yy]$ ]] || die "Aborted."
fi

# ---------------------------------------------------------------------------
# Inputs: links (NAS scratch) or copies (other disk), same names for FOD and mask
# ---------------------------------------------------------------------------
# RAM copies must never outlive the run, whatever happens
trap 'rm -f "$LOCK"; [ "$INPUTS" = shm ] && [ -n "${IN_ROOT_4B:-}" ] && rm -rf "$IN_ROOT_4B"' EXIT
rm -rf "$IN_ROOT_4B"; mkdir -p "$IN_ROOT_4B/fod_input" "$IN_ROOT_4B/mask_input"
t0=$(date +%s)
for s in "${SUBJS[@]}"; do
    sd="${IN_SRC_DIR}/${s}"
    if [ "$STAGE" -eq 1 ]; then
        cp "$sd/halfway_template.mif" "$IN_ROOT_4B/fod_input/${s}.mif"  || die "copy failed: ${s}"
        cp "$sd/halfway_mask.mif"     "$IN_ROOT_4B/mask_input/${s}.mif" || die "copy failed: ${s}"
    else
        ln -sfn "$sd/halfway_template.mif" "$IN_ROOT_4B/fod_input/${s}.mif"
        ln -sfn "$sd/halfway_mask.mif"     "$IN_ROOT_4B/mask_input/${s}.mif"
    fi
done
t1=$(date +%s)
log "Inputs ready in $(hours "$t0" "$t1") h"

# population_template's own template-mask output, if this version has it
TMASK_OPT=()
population_template -help 2>/dev/null | grep -q -- '-template_mask' && \
    TMASK_OPT=(-template_mask "${GROUP_DIR}/template_mask_40.mif")

TIMEF="${GROUP_DIR}/step4b_time.txt"
TIMECMD=()
[ -x /usr/bin/time ] && /usr/bin/time -v true >/dev/null 2>&1 && TIMECMD=(/usr/bin/time -v -o "$TIMEF")

WORK_SCRATCH="${SCRATCH_ROOT}/step4b_scratch"
mkdir -p "$WORK_SCRATCH"

# ---------------------------------------------------------------------------
# Launch in its own process group (immune to Ctrl+C here), then monitor
# ---------------------------------------------------------------------------
setsid "${TIMECMD[@]}" population_template "$IN_ROOT_4B/fod_input" "${GROUP_DIR}/wmfod_template.mif" \
    -mask_dir "$IN_ROOT_4B/mask_input" \
    "${TMASK_OPT[@]}" \
    -voxel_size "$UPSAMPLE_VOXEL" \
    -scratch "$WORK_SCRATCH" \
    -nthreads "$NTHREADS" -force < /dev/null >> "$LOGF" 2>&1 &
PT=$!
STOP_CMD="kill -TERM -- -${PT}"
echo "$STOP_CMD" > "${GROUP_DIR}/HOW_TO_STOP.txt"
log "population_template started (PID ${PT}). Log: ${LOGF}"
log "Ctrl+C here will NOT stop it. To stop it on purpose:  ${STOP_CMD}"

trap 'log "Ctrl+C ignored: the template build keeps running. To stop it on purpose:  ${STOP_CMD}"' INT

peak_b=0; n=0; status="RUNNING"; last_report=$(date +%s)
while alive "$PT"; do
    sleep "$POLL_SECS"
    n=$((n+1))
    free_now=$(avail_gb "$SCRATCH_ROOT")
    if [ -n "$free_now" ] && [ "$free_now" -lt "$MIN_FREE_GB" ]; then
        log "Free space on scratch disk down to ${free_now} GB: stopping the run to protect the disk."
        kill -TERM -- "-${PT}" 2>/dev/null; sleep 20; kill -KILL -- "-${PT}" 2>/dev/null
        status="STOPPED_LOW_SPACE"; break
    fi
    if [ $((n % DU_EVERY)) -eq 0 ]; then
        cur=$(du -sb "$WORK_SCRATCH" 2>/dev/null | cut -f1)
        [ -n "$cur" ] && [ "$cur" -gt "$peak_b" ] && peak_b=$cur
    fi
    now=$(date +%s)
    if [ $((now - last_report)) -ge 1800 ]; then
        stage=$(grep -o 'Optimising template[^)]*)' "$LOGF" 2>/dev/null | tail -1)
        log "running $(hours "$t1" "$now") h | scratch peak $(awk -v b=$peak_b 'BEGIN{printf "%.0f", b/1073741824}') GB | free ${free_now} GB | ${stage:-starting}"
        last_report=$now
    fi
done
wait "$PT" 2>/dev/null; rc=$?
t2=$(date +%s)
[ "$status" = "RUNNING" ] && { [ "$rc" -eq 0 ] && status="OK" || status="FAILED_rc${rc}"; }

# ---------------------------------------------------------------------------
# Check the result and record the run
# ---------------------------------------------------------------------------
if [ "$status" = "OK" ]; then
    sz=$(mrinfo "${GROUP_DIR}/wmfod_template.mif" -size 2>/dev/null)
    [ "$(echo "$sz" | awk '{print $4}')" = "45" ] || { status="FAILED_OUTPUT_CHECK"; log "Template unreadable or not 45 volumes: ${sz}"; }
fi
rss_gb=0
[ -f "$TIMEF" ] && rss_gb=$(grep 'Maximum resident set size' "$TIMEF" | awk '{printf "%.1f", $NF/1048576}')
peak_gb=$(awk -v b=$peak_b 'BEGIN{printf "%.1f", b/1073741824}')
[ -f "$RUNREC" ] || printf 'start\tend\thours\tn_subjects\tthreads\tscratch_root\tpeak_scratch_GB\tpeak_rss_GB\tstatus\n' > "$RUNREC"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$(date -d "@$t0" '+%F %T')" "$(date -d "@$t2" '+%F %T')" \
    "$(hours "$t0" "$t2")" "${#SUBJS[@]}" "$NTHREADS" "$SCRATCH_ROOT" "$peak_gb" "$rss_gb" "$status" >> "$RUNREC"

log "Step 4b ${status} after $(hours "$t0" "$t2") h; peak scratch ${peak_gb} GB; peak RAM ${rss_gb} GB"
if [ "$status" = "OK" ]; then
    rm -rf "$IN_ROOT_4B"
    rm -rf "$WORK_SCRATCH"
    cat <<EOF

Group template: ${GROUP_DIR}/wmfod_template.mif
Check it before step 4c:
  mrview ${GROUP_DIR}/wmfod_template.mif -odf.load_sh ${GROUP_DIR}/wmfod_template.mif$( [ -f "${GROUP_DIR}/template_mask_40.mif" ] && echo " -overlay.load ${GROUP_DIR}/template_mask_40.mif" )
Look for sharp anatomy (not blurred), and for the mask covering the temporal poles
and medial temporal lobe, where the uncinate and fornix run.
Next: step 4c registers all 80 halfway templates to this template.
EOF
elif [ "$status" = "STOPPED_LOW_SPACE" ]; then
    # The point of stopping was to free the disk, so clear everything this run wrote
    rm -rf "$WORK_SCRATCH" "$IN_ROOT_4B"
    log "Scratch and staged inputs deleted to free the disk. Peak scratch reached before stopping: ${peak_gb} GB."
    log "Rerun with SCRATCH_ROOT on a disk with more room."
else
    log "Scratch kept for inspection: ${WORK_SCRATCH} (it may be very large)"
    log "Delete it when done:  rm -rf ${WORK_SCRATCH}$([ "$INPUTS" != shm ] && echo " ${IN_ROOT_4B}")"
    log "Details: ${LOGF}"
fi
