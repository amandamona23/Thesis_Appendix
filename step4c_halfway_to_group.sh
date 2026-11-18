#!/bin/bash
###############################################################################
# step4c_halfway_to_group.sh  -  Register ALL 80 halfway templates to the group
#                                template, and build the final analysis mask
#
# For every subject (the 40 that built the template AND the other 40, so all 80
# are treated identically), it runs the standard FBA registration step:
#
#   mrregister halfway_template.mif  group/wmfod_template.mif
#              -mask1 halfway_mask.mif  -nl_warp halfway2group.mif group2halfway.mif
#
# Default mrregister settings (affine then non-linear, FOD-based), as in the
# MRtrix FBA tutorial and the two-step method of Lebrun et al.
#
#   READS   ${TEMPLATE_DIR}/group/wmfod_template.mif             (from step 4b)
#           ${TEMPLATE_DIR}/subject_templates/<SUBJ>/halfway_template.mif, halfway_mask.mif
#   WRITES  per subject, in ${TEMPLATE_DIR}/subject_templates/<SUBJ>/ :
#             halfway2group_warp.mif    halfway -> group   (4D deformation field)
#             group2halfway_warp.mif    group -> halfway
#             halfway_mask_in_group.mif the halfway mask warped into group space
#             halfway2group.done        marker: this subject is finished
#           ${TEMPLATE_DIR}/group/template_mask.mif   the FINAL analysis mask:
#             voxels that are brain in ALL 80 subjects (intersection)
#           ${TEMPLATE_DIR}/run4c/run4c_results.tsv   one row per subject
#
# WHY THESE WARPS ARE EASIER THAN THE 4a ONES: mrregister -nl_warp writes plain
# 4D deformation fields, which chain directly in step 5. (The 4a warps from
# population_template are 5D "full" warps and need converting first.)
#
# SAFETY
#   - Each subject runs as its own detached worker: Ctrl+C in the controller
#     window only stops NEW subjects starting.
#   - Results are written to a temporary folder and moved into place only after
#     they pass checks (warp is 4D, on the group grid; mask on the group grid),
#     so a half-finished run is never mistaken for a finished one.
#   - Refuses to start if any subject's 4a result is missing or unconfirmed
#     (this catches a session whose halfway template is still being remade).
#   - Only one controller at a time. Graceful stop:  touch ${TEMPLATE_DIR}/STOP_4C
#   - The final mask is built only when ALL subjects are done.
#   - Never use pkill -f with this script's name: it would kill the workers.
#
# Usage:
#   NTHREADS=16 MAX_PARALLEL_4C=4 ./step4c_halfway_to_group.sh --list sessions.txt
#   ./step4c_halfway_to_group.sh --list sessions.txt --mask-only     # just rebuild the final mask
###############################################################################

set -u
SELF="$(readlink -f "$0")"
: "${NTHREADS:=16}"; export NTHREADS
source "$(dirname "$SELF")/pipeline_config_v2.sh"

SUBJ_TEMPLATE_DIR="${TEMPLATE_DIR}/subject_templates"
GROUP_DIR="${TEMPLATE_DIR}/group"
GT="${GROUP_DIR}/wmfod_template.mif"
RUN_DIR="${TEMPLATE_DIR}/run4c"
RESULTS="${RUN_DIR}/run4c_results.tsv"
STOP_FILE="${TEMPLATE_DIR}/STOP_4C"
MAX_PAR="${MAX_PARALLEL_4C:-4}"
POLL_SECS="${POLL_SECS:-30}"
STAGGER_SECS="${STAGGER_SECS:-60}"
RESULT_LOCK="/tmp/${USER}_run4c_results.lock"

minutes() { awk -v a="$1" -v b="$2" 'BEGIN{printf "%.1f", (b-a)/60}'; }
alive() { [ -e "/proc/$1" ] && ! grep -q 'Z (zombie)' "/proc/$1/status" 2>/dev/null; }
grid3() { mrinfo "$1" -size 2>/dev/null | awk '{print $1,$2,$3}'; }
is_done() { [ -f "$SUBJ_TEMPLATE_DIR/$1/halfway2group.done" ]; }
is_running() {
    local u; u=$(id -u)
    pgrep -u "$u" -f "step4c_halfway_to_group.sh --worker $1\$" >/dev/null && return 0
    pgrep -u "$u" -f "subject_templates/$1/" >/dev/null && return 0
    return 1
}
record_row() {  # subj t0 t1 status
    local line
    line=$(printf '%s\t%s\t%s\t%s\t%s\t%s' "$1" "$(date -d "@$2" '+%F %T')" "$(date -d "@$3" '+%F %T')" \
           "$(minutes "$2" "$3")" "$NTHREADS" "$4")
    ( flock 9; printf '%s\n' "$line" >> "$RESULTS" ) 9>"$RESULT_LOCK"
}

###############################################################################
# WORKER: one subject
###############################################################################
run_worker() {
    local subj=$1 SD="${SUBJ_TEMPLATE_DIR}/$1" tmp t0 t1 rc
    tmp="${SD}/tmp_4c"
    mkdir -p "${RUN_DIR}/logs"; exec >> "${RUN_DIR}/logs/${subj}.log" 2>&1
    log "worker start: ${subj}, ${NTHREADS} threads"
    t0=$(date +%s)
    rm -rf "$tmp" "$SD/halfway2group_warp.mif" "$SD/group2halfway_warp.mif" "$SD/halfway_mask_in_group.mif"
    mkdir -p "$tmp"

    mrregister "$SD/halfway_template.mif" "$GT" \
        -mask1 "$SD/halfway_mask.mif" \
        -nl_warp "$tmp/h2g.mif" "$tmp/g2h.mif" \
        -nthreads "$NTHREADS" -force
    rc=$?
    if [ "$rc" -ne 0 ] || [ ! -f "$tmp/h2g.mif" ] || [ ! -f "$tmp/g2h.mif" ]; then
        t1=$(date +%s); record_row "$subj" "$t0" "$t1" "FAIL_REGISTER"; rm -rf "$tmp"; log "FAIL_REGISTER"; return 1
    fi

    # The warp must be a plain 4D deformation field on the group template's grid
    if [ "$(mrinfo "$tmp/h2g.mif" -ndim 2>/dev/null)" != "4" ] || [ "$(grid3 "$tmp/h2g.mif")" != "$(grid3 "$GT")" ]; then
        t1=$(date +%s); record_row "$subj" "$t0" "$t1" "FAIL_WARP_CHECK"
        log "warp is not a 4D field on the group grid: $(mrinfo "$tmp/h2g.mif" -ndim 2>&1) / $(grid3 "$tmp/h2g.mif") vs $(grid3 "$GT")"
        rm -rf "$tmp"; return 1
    fi

    mrtransform "$SD/halfway_mask.mif" -warp "$tmp/h2g.mif" -interp nearest -datatype bit \
        "$tmp/mask_g.mif" -force -quiet
    if [ ! -f "$tmp/mask_g.mif" ] || [ "$(grid3 "$tmp/mask_g.mif")" != "$(grid3 "$GT")" ]; then
        t1=$(date +%s); record_row "$subj" "$t0" "$t1" "FAIL_MASK"; log "mask not on the group grid"; rm -rf "$tmp"; return 1
    fi

    mv "$tmp/h2g.mif"    "$SD/halfway2group_warp.mif"
    mv "$tmp/g2h.mif"    "$SD/group2halfway_warp.mif"
    mv "$tmp/mask_g.mif" "$SD/halfway_mask_in_group.mif"
    rm -rf "$tmp"
    touch "$SD/halfway2group.done"
    t1=$(date +%s); record_row "$subj" "$t0" "$t1" "OK"
    log "OK in $(minutes "$t0" "$t1") min"
    return 0
}

###############################################################################
# Final analysis mask (intersection of all subjects' masks in group space)
###############################################################################
build_final_mask() {
    local s masks=() out="${GROUP_DIR}/template_mask.mif"
    for s in "${SUBJS[@]}"; do
        is_done "$s" || { log "Final mask not built: ${s} is not finished"; return 1; }
        masks+=("${SUBJ_TEMPLATE_DIR}/${s}/halfway_mask_in_group.mif")
    done
    mrmath "${masks[@]}" min "$out" -datatype bit -force || { log "mrmath failed"; return 1; }
    [ "$(grid3 "$out")" = "$(grid3 "$GT")" ] || { log "Final mask is not on the group grid"; return 1; }

    local n_final n40="-"
    n_final=$(mrstats "$out" -mask "$out" -output count 2>/dev/null)
    [ -f "${GROUP_DIR}/template_mask_40.mif" ] && n40=$(mrstats "${GROUP_DIR}/template_mask_40.mif" -mask "${GROUP_DIR}/template_mask_40.mif" -output count 2>/dev/null)
    log "Final analysis mask: ${out}"
    log "  voxels in mask, all ${#SUBJS[@]} subjects: ${n_final}    (40-subject mask: ${n40})"

    # Which subjects cover the least brain? Useful if the mask looks too small.
    for s in "${SUBJS[@]}"; do
        f="${SUBJ_TEMPLATE_DIR}/${s}/halfway_mask_in_group.mif"
        echo "$(mrstats "$f" -mask "$f" -output count 2>/dev/null) ${s}"
    done | sort -n | head -3 | awk '{printf "  smallest coverage in group space: %s (%s voxels)\n", $2, $1}' | while read -r l; do log "$l"; done
    echo "Look at it before step 5: mrview ${GT} -overlay.load ${out}"
    echo "Check the temporal poles and medial temporal lobes are inside it (uncinate and fornix)."
}

if [ "${1:-}" = "--worker" ]; then run_worker "$2"; exit $?; fi

###############################################################################
# CONTROLLER
###############################################################################
LIST=""; MASK_ONLY=0; ASSUME_YES=0
while [ $# -gt 0 ]; do
    case "$1" in
        --list) LIST="$(readlink -f "$2")"; shift ;;
        --mask-only) MASK_ONLY=1 ;;
        --yes) ASSUME_YES=1 ;;
        *) die "Unknown option: $1" ;;
    esac; shift
done
[ -n "$LIST" ] && [ -f "$LIST" ] || die "Usage: NTHREADS=16 MAX_PARALLEL_4C=4 $0 --list sessions.txt"
mapfile -t SUBJS < <(cut -d_ -f1 "$LIST" | grep -o 'HCA[0-9]*' | awk '!seen[$0]++')

[ -f "$GT" ] && [ "$(mrinfo "$GT" -size 2>/dev/null | awk '{print $4}')" = "45" ] \
    || die "Group template missing or not a readable 45-volume FOD: ${GT} (run step 4b first)"
problems=()
for s in "${SUBJS[@]}"; do
    sd="${SUBJ_TEMPLATE_DIR}/${s}"
    [ -f "$sd/halfway_template.mif" ]    || problems+=("${s}: no halfway_template.mif")
    [ -f "$sd/halfway_mask.mif" ]        || problems+=("${s}: no halfway_mask.mif")
    [ -f "$sd/mask_halfway_space.done" ] || problems+=("${s}: 4a not finished or mask not confirmed in halfway space")
done
if [ "${#problems[@]}" -gt 0 ]; then printf '  %s\n' "${problems[@]}"; die "Step 4a is not complete for ${#problems[@]} item(s) above."; fi

mkdir -p "$RUN_DIR/logs"
[ -f "$RESULTS" ] || printf 'subject\tstart\tend\tminutes\tthreads\tstatus\n' > "$RESULTS"

if [ "$MASK_ONLY" -eq 1 ]; then build_final_mask; exit $?; fi

CTRL_LOCK="/tmp/${USER}_step4c_controller.pid"
if [ -f "$CTRL_LOCK" ] && alive "$(cat "$CTRL_LOCK")"; then die "A step 4c controller is already running (PID $(cat "$CTRL_LOCK"))."; fi
echo $$ > "$CTRL_LOCK"
[ -f "$STOP_FILE" ] && { rm -f "$CTRL_LOCK"; die "Stop file exists (${STOP_FILE}). Delete it first."; }

declare -A LANE=()
on_exit() { rm -f "$CTRL_LOCK"; [ "${#LANE[@]}" -gt 0 ] && log "Controller stopped; ${#LANE[@]} subject(s) keep running and will finish on their own."; }
trap on_exit EXIT
trap 'log "Interrupted: no new subjects will start."; exit 130' INT TERM

QUEUE=(); n_done=0; n_busy=0
for s in "${SUBJS[@]}"; do
    if is_done "$s"; then n_done=$((n_done+1)); elif is_running "$s"; then n_busy=$((n_busy+1)); else QUEUE+=("$s"); fi
done
running_now=$(pgrep -u "$(id -u)" -fc 'bin/population_template|bin/mrregister')
cat <<EOF

-------------------------------------------------------------------
Step 4c: halfway templates to group template
  subjects:            ${#SUBJS[@]}   already finished: ${n_done}   running elsewhere: ${n_busy}   to do: ${#QUEUE[@]}
  at once:             ${MAX_PAR}, ${NTHREADS} threads each
  group template:      ${GT}
  results:             ${RESULTS}
  other template or registration processes running now: ${running_now}
  stop gracefully:     touch ${STOP_FILE}
-------------------------------------------------------------------
EOF
[ "${#QUEUE[@]}" -gt 0 ] || { log "Nothing to register."; build_final_mask; exit 0; }
if [ "$ASSUME_YES" -eq 0 ]; then read -r -p "Start? [y/N] " ok; [[ "$ok" =~ ^[Yy]$ ]] || die "Aborted."; fi

reap() { local p; for p in "${!LANE[@]}"; do alive "$p" || { log "FINISHED ${LANE[$p]}: $(tail -1 "${RUN_DIR}/logs/${LANE[$p]}.log" 2>/dev/null)"; unset "LANE[$p]"; }; done; }
i=0; stopping=0
while true; do
    reap
    if [ -f "$STOP_FILE" ] && [ "$stopping" -eq 0 ]; then log "Stop file found: nothing new will start."; stopping=1; fi
    if [ "$stopping" -eq 0 ]; then
        while [ "$i" -lt "${#QUEUE[@]}" ] && [ "${#LANE[@]}" -lt "$MAX_PAR" ]; do
            s="${QUEUE[$i]}"; i=$((i+1))
            if is_done "$s" || is_running "$s"; then log "skip ${s}: already done or in progress"; continue; fi
            setsid "$SELF" --worker "$s" < /dev/null > /dev/null 2>&1 &
            LANE[$!]="$s"; log "START ${s} (running: ${#LANE[@]} of ${MAX_PAR})"
            sleep "$STAGGER_SECS"; reap
        done
    fi
    if [ "${#LANE[@]}" -eq 0 ] && { [ "$i" -ge "${#QUEUE[@]}" ] || [ "$stopping" -eq 1 ]; }; then break; fi
    sleep "$POLL_SECS"
done

fin=0; for s in "${SUBJS[@]}"; do is_done "$s" && fin=$((fin+1)); done
log "Finished: ${fin} of ${#SUBJS[@]} subjects registered to the group template."
if [ "$fin" -eq "${#SUBJS[@]}" ]; then build_final_mask
else log "Not all done. Rerun the same command to retry the rest (finished subjects are skipped). Failures: grep -v OK ${RESULTS}"; fi
