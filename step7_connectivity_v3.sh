#!/bin/bash
###############################################################################
# step7_connectivity_v3.sh  -  Template tractography, fixel connectivity, smoothing
#
# MRtrix multi-tissue FBA tutorial steps 19 to 22, restricted to your tracts:
#   1. analysis mask  union of all tract fixel masks from step 6 (cutoff MASK_CUTOFF,
#                     default 1 = the most inclusive, so every stricter mask used in
#                     the statistics is inside it)
#   2. tckgen         whole-brain tractography on the group FOD template (tutorial 19)
#                     tutorial settings: -angle 22.5 -maxlen 250 -minlen 10 -power 1.0
#                     -cutoff 0.06, seeded and masked by the template mask
#   3. tcksift        reduce to N_SIFT streamlines so densities match the FODs (tutorial 20)
#   4. fixelconnectivity  fixel-to-fixel connectivity, ONLY for the tract fixels
#                     (-mask), so the matrix stays small (tutorial 21)
#   5. fixelfilter    smooth FD, log FC and FDC along that connectivity, again only
#                     within the tract fixels (tutorial 22). Also smooths the 0.06
#                     versions (fd_fmls0p06 ...) if step 5 made them.
# Statistics (tutorial 23) are step 9: fixelcfestats ... -mask <tract mask>.
#
#   READS   group/wmfod_template.mif, group/template_mask.mif, fixel_mask/,
#           tracts/fixels/*/*_ge<MASK_CUTOFF>.mif (step 6), fd/ log_fc/ fdc/ (step 5)
#   WRITES  ${TEMPLATE_DIR}/step7/   analysis_mask/, tracks_sift.tck, matrix/, record files
#           ${TEMPLATE_DIR}/fd_smooth/ log_fc_smooth/ fdc_smooth/  (and *_fmls0p06_smooth/)
#
# The tractography (tckgen) and SIFT are the heavy parts: hours, and SIFT needs a lot
# of memory. RUN THE TEST FIRST: it does everything at 1/20 size in step7_test/ and
# records time and peak memory, so the full run can be planned.
#
# SAFETY: stages are resumable (a finished stage is skipped on rerun); results are
# written under hidden names and renamed when complete; the 20-million-streamline
# file is kept on LOCAL disk and deleted after SIFT; one run at a time.
#
# Usage:
#   ./step7_connectivity_v3.sh --test          # 1 million -> 100,000 streamlines, in step7_test/
#   ./step7_connectivity_v3.sh                 # full: 20 million -> 2 million (tutorial)
#   ./step7_connectivity_v3.sh --smooth-only   # redo only the smoothing (e.g. new step 5 data)
#   ./step7_connectivity_v3.sh --summary
# Settings (environment): N_TRACKS  N_SIFT  NTHREADS (default 48)  MASK_CUTOFF (1)
#                         SCRATCH_ROOT ($HOME/fba7_work)  KEEP_TRACKS=1 keeps the big file
###############################################################################

set -u
SELF="$(readlink -f "$0")"
: "${NTHREADS:=48}"; export NTHREADS
source "$(dirname "$SELF")/pipeline_config_v2.sh"

TGROUP="${TEMPLATE_DIR}/group"; GT="${TGROUP}/wmfod_template.mif"; TMASK="${TGROUP}/template_mask.mif"
FIXEL_MASK="${TEMPLATE_DIR}/fixel_mask"; TRACTS_DIR="${TEMPLATE_DIR}/tracts"
MASK_CUTOFF="${MASK_CUTOFF:-1}"
SCRATCH_ROOT="${SCRATCH_ROOT:-$HOME/fba7_work}"
LOCK="/tmp/${USER}_step7.lock"

MODE=full
while [ $# -gt 0 ]; do
    case "$1" in
        --test) MODE=test ;; --smooth-only) MODE=smooth ;; --summary) MODE=summary ;;
        *) die "Unknown option: $1" ;;
    esac; shift
done
if [ "$MODE" = "test" ]; then
    N_TRACKS="${N_TRACKS:-1000000}"; N_SIFT="${N_SIFT:-100000}"; OUT="${TEMPLATE_DIR}/step7_test"; SMOOTH_ROOT="$OUT"
else
    N_TRACKS="${N_TRACKS:-20000000}"; N_SIFT="${N_SIFT:-2000000}"; OUT="${TEMPLATE_DIR}/step7"; SMOOTH_ROOT="$TEMPLATE_DIR"
fi
AMASK_DIR="${OUT}/analysis_mask"; AMASK="${AMASK_DIR}/analysis_mask.mif"
WORK="${SCRATCH_ROOT}/$(basename "$OUT")"; TRACKS="${WORK}/tracks_all.tck"
SIFT="${OUT}/tracks_sift.tck"; MATRIX="${OUT}/matrix"; REC="${OUT}/step7_record.txt"

nfix()   { mrinfo "$1" -size 2>/dev/null | awk '{print $1}'; }
count()  { mrstats "$1" -mask "$1" -output count 2>/dev/null | tr -d ' '; }
avail_gb() { df -BG --output=avail "$1" 2>/dev/null | tail -1 | tr -dc '0-9'; }
mem_gb() { awk '/MemAvailable/ {printf "%d", $2/1048576}' /proc/meminfo; }
hrs()    { awk -v s="$1" 'BEGIN{printf "%.2f", s/3600}'; }
tck_count() { tckinfo "$1" 2>/dev/null | awk -F: '/^ *count/ {gsub(/ /,"",$2); print $2; exit}'; }
record() { echo "$*" >> "$REC"; log "$*"; }
timed() {  # timed <label> <command...> : runs it, records hours and peak memory
    local label=$1; shift; local t0 rc tf="${OUT}/.time_${label}.txt" peak="NA"
    t0=$(date +%s)
    if [ -x /usr/bin/time ]; then /usr/bin/time -v -o "$tf" "$@"; rc=$?
        peak=$(awk -F: '/Maximum resident set size/ {printf "%.1f", $2/1048576}' "$tf" 2>/dev/null)
    else "$@"; rc=$?; fi
    record "${label}: exit ${rc}, $(hrs $(( $(date +%s) - t0 ))) h, peak memory ${peak:-NA} GB"
    return $rc
}
summary() {
    local d any=0
    for d in "${TEMPLATE_DIR}/step7_test" "${TEMPLATE_DIR}/step7"; do
        [ -f "${d}/step7_record.txt" ] || continue
        any=1; echo "== ${d##*/}"; cat "${d}/step7_record.txt"
    done
    [ "$any" = 1 ] || echo "Nothing recorded yet."
}

[ "$MODE" = "summary" ] && { summary; exit 0; }      # read-only, allowed while a run is going
exec 9>"$LOCK"; flock -n 9 || die "Step 7 is already running."

###############################################################################
# Inputs
###############################################################################
[ -f "$GT" ] && [ -f "$TMASK" ] || die "Group template or template mask missing (steps 4b, 4c)"
[ -f "${FIXEL_MASK}/index.mif" ] || die "Fixel mask missing (step 5 --build-mask)"
[ -f "${TRACTS_DIR}/tracts.done" ] || die "Tract masks missing (step 6)"
mapfile -t TMASKS < <(ls "${TRACTS_DIR}"/fixels/*/*_ge"${MASK_CUTOFF}".mif 2>/dev/null | grep -v '_both_')
[ "${#TMASKS[@]}" -gt 0 ] || die "No step 6 fixel masks at cutoff ${MASK_CUTOFF}"
for c in tckgen tcksift tckinfo fixelconnectivity fixelfilter mrcalc mrstats; do
    command -v "$c" >/dev/null || die "Command not found: ${c}"
done
mkdir -p "$OUT" "$WORK"

cat <<EOF

-------------------------------------------------------------------
Step 7 ($( [ "$MODE" = test ] && echo "TEST, 1/20 size" || echo "$MODE" ))
  streamlines:     ${N_TRACKS} generated, SIFT to ${N_SIFT}
  threads:         ${NTHREADS}
  tract masks:     ${#TMASKS[@]} at cutoff ${MASK_CUTOFF}: $(for f in "${TMASKS[@]}"; do basename "$f" .mif; done | tr '\n' ' ')
  output:          ${OUT}
  local scratch:   ${WORK} ($(avail_gb "$SCRATCH_ROOT") GB free)
  memory free now: $(mem_gb) GB
-------------------------------------------------------------------
EOF

###############################################################################
# 1. Analysis mask: union of the tract fixel masks
###############################################################################
if [ ! -f "${AMASK_DIR}/.done" ]; then
    rm -rf "${AMASK_DIR}" "${AMASK_DIR}.building"; mkdir -p "${AMASK_DIR}.building"
    cp "${FIXEL_MASK}/index.mif" "${FIXEL_MASK}/directions.mif" "${AMASK_DIR}.building/"
    cp "${TMASKS[0]}" "${AMASK_DIR}.building/analysis_mask.mif"
    for f in "${TMASKS[@]:1}"; do
        mrcalc "${AMASK_DIR}.building/analysis_mask.mif" "$f" -max "${AMASK_DIR}.building/u.mif" -datatype bit -force -quiet \
            && mv "${AMASK_DIR}.building/u.mif" "${AMASK_DIR}.building/analysis_mask.mif" || die "mask union failed"
    done
    [ "$(nfix "${AMASK_DIR}.building/analysis_mask.mif")" = "$(nfix "${FIXEL_MASK}/directions.mif")" ] || die "analysis mask has the wrong size"
    mv "${AMASK_DIR}.building" "${AMASK_DIR}"; touch "${AMASK_DIR}/.done"
    record "analysis mask: $(count "$AMASK") fixels (union of ${#TMASKS[@]} tract masks at cutoff ${MASK_CUTOFF})"
fi

if [ "$MODE" != "smooth" ]; then
###############################################################################
# 2. Whole-brain tractography on the template (tutorial step 19)
###############################################################################
if [ ! -f "${OUT}/.tracks.done" ] && [ ! -f "${OUT}/.sift.done" ]; then
    need=$(( N_TRACKS / 500000 + ${DISK_MARGIN_GB:-30} ))   # about 2 GB per million streamlines, plus margin (estimate)
    [ "$(avail_gb "$SCRATCH_ROOT")" -ge "$need" ] || die "Local disk: need about ${need} GB free in ${SCRATCH_ROOT}"
    rm -f "$TRACKS"
    timed tckgen tckgen -angle 22.5 -maxlen 250 -minlen 10 -power 1.0 "$GT" -seed_image "$TMASK" -mask "$TMASK" \
        -select "$N_TRACKS" -cutoff 0.06 "$TRACKS" -nthreads "$NTHREADS" -force || die "tckgen failed"
    n=$(tck_count "$TRACKS"); record "tckgen: ${n} streamlines, file $(du -h "$TRACKS" | cut -f1)"
    [ "${n:-0}" -ge "$N_TRACKS" ] 2>/dev/null || die "tckgen produced ${n} of ${N_TRACKS} streamlines"
    touch "${OUT}/.tracks.done"
fi

###############################################################################
# 3. SIFT (tutorial step 20)
###############################################################################
if [ ! -f "${OUT}/.sift.done" ]; then
    [ -f "$TRACKS" ] || die "The tractogram is missing (${TRACKS}); delete ${OUT}/.tracks.done and rerun"
    record "SIFT starting with $(mem_gb) GB memory free"
    timed tcksift tcksift "$TRACKS" "$GT" "${SIFT}.building.tck" -term_number "$N_SIFT" -nthreads "$NTHREADS" -force \
        || die "tcksift failed"
    n=$(tck_count "${SIFT}.building.tck"); record "tcksift: ${n} streamlines kept"
    [ "${n:-0}" -gt 0 ] 2>/dev/null || die "SIFT output is empty"
    mv "${SIFT}.building.tck" "$SIFT"; touch "${OUT}/.sift.done"
    [ "${KEEP_TRACKS:-0}" = "1" ] || { rm -f "$TRACKS"; record "deleted the full tractogram from local disk"; }
fi

###############################################################################
# 4. Fixel-fixel connectivity, tract fixels only (tutorial step 21)
###############################################################################
if [ ! -f "${OUT}/.matrix.done" ]; then
    rm -rf "${MATRIX}" "${MATRIX}.building"
    timed fixelconnectivity fixelconnectivity "$FIXEL_MASK" "$SIFT" "${MATRIX}.building" -mask "$AMASK" -nthreads "$NTHREADS" \
        || die "fixelconnectivity failed"
    [ -n "$(ls -A "${MATRIX}.building" 2>/dev/null)" ] || die "connectivity matrix is empty"
    mv "${MATRIX}.building" "$MATRIX"; touch "${OUT}/.matrix.done"
    record "connectivity matrix: $(du -sh "$MATRIX" | cut -f1)"
fi
fi

###############################################################################
# 5. Smoothing (tutorial step 22)
###############################################################################
[ -f "${OUT}/.matrix.done" ] || die "No connectivity matrix yet: run the full step first"
METRICS=(fd log_fc fdc)
for t in "${TEMPLATE_DIR}"/fd_fmls*; do
    case "$t" in *_smooth) continue ;; esac          # never smooth an already-smoothed folder
    [ -d "$t" ] && { tag=${t##*/fd}; METRICS+=("fd${tag}" "log_fc${tag}" "fdc${tag}"); }
done
[ "$MODE" = "test" ] && METRICS=(fd)      # the test only times one
for m in "${METRICS[@]}"; do
    src="${TEMPLATE_DIR}/${m}"; dst="${SMOOTH_ROOT}/${m}_smooth"
    [ -d "$src" ] || { log "skip ${m}: no folder"; continue; }
    nin=$(ls "$src" | grep -c '^HCA.*\.mif$')
    if [ -f "${dst}/.done" ] && [ "$MODE" != "smooth" ]; then log "${m}_smooth already done"; continue; fi
    rm -rf "${dst}.building"
    timed "smooth_${m}" fixelfilter "$src" smooth "${dst}.building" -matrix "$MATRIX" -mask "$AMASK" -nthreads "$NTHREADS" \
        || die "fixelfilter failed for ${m}"
    nout=$(ls "${dst}.building" | grep -c '^HCA.*\.mif$')
    [ "$nout" = "$nin" ] || die "${m}: smoothed ${nout} of ${nin} files"
    rm -rf "$dst"; mv "${dst}.building" "$dst"; touch "${dst}/.done"
    record "smoothed ${m}: ${nout} sessions -> ${dst}"
done

log "Step 7 ${MODE} finished."
echo; summary
if [ "$MODE" = "test" ]; then cat <<EOF

For the full run, scale from the test: tckgen and SIFT time grow about 20-fold; SIFT memory
grows with the number of streamlines, so check the SIFT peak above against the memory free
(about $(mem_gb) GB now) before starting. Paste these lines and I will estimate it with you.
EOF
fi
exit 0
