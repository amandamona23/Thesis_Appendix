#!/bin/bash
###############################################################################
# step5_fixels_v3.sh  -  Fixel mask and per-session FBA metrics (FD, FC, FDC)
#
# This is steps 12 to 18 of the MRtrix multi-tissue FBA tutorial.
#
# PART 1, once (--build-mask)                          tutorial step 12
#   Segment the group template into fixels. The result, fixel_mask/, is the
#   fixed list of "slots" (one per fibre population per voxel) that every
#   session is measured at, and where every statistical test happens.
#   CHECK IT before part 2: changing it later means redoing all sessions.
#
# PART 2, per session, in parallel (--list sessions.txt)
#   a) build ONE warp from the session straight to the group template, by
#      chaining  session -> halfway (step 4a)  with  halfway -> group (step 4c)
#   b) warp the FOD with that warp, no reorientation yet     tutorial step 13
#   c) segment it into fixels, FD = fixel integral (-afd)     tutorial step 14
#   d) rotate the fixel directions by the warp                tutorial step 15
#   e) match the session's fixels to the template slots       tutorial step 16
#   f) FC from the warp, then log FC                          tutorial step 17
#   g) FDC = FD x FC                                          tutorial step 18
#
# WHY CHAIN THE WARPS: the 4a warp carries each visit's change relative to the
# subject's own halfway space, so FC differs between a subject's visits (the
# longitudinal signal). One combined warp means the FOD is interpolated ONCE.
# The 4a warps are population_template's 5D "full" warps, so they are first
# converted to a plain deformation field; the 4c warps already are one.
#
#   READS   group/wmfod_template.mif, group/template_mask.mif     (4b, 4c)
#           subject_templates/<SUBJ>/warps/<SESS>.mif, halfway_template.mif (4a)
#           subject_templates/<SUBJ>/halfway2group_warp.mif               (4c)
#           ${OUT_ROOT}/<SESS>/wmfod_norm.mif                             (3a)
#   WRITES  ${TEMPLATE_DIR}/fixel_mask/            template fixels
#           ${TEMPLATE_DIR}/fd/ fc/ log_fc/ fdc/   one <SESS>.mif per session
#           ${TEMPLATE_DIR}/session_warps/<SESS>.mif   combined warp (step 8 reuses it for DTI)
#           ${TEMPLATE_DIR}/step5/step5_results.tsv    one row per session, with QC numbers
#
# QC NUMBERS PER SESSION (in the results table)
#   r_compose   one-step warp vs applying the two warps in turn (should be ~0.98+;
#               below 0.9 the session FAILS: the chained warp is wrong)
#   r_template  warped session vs group template (alignment; compare across sessions)
#   fd_zero_pct template fixels with no matching session fixel (FD = 0)
#   fc_mean     mean FC (near 1; below 1 = smaller than the template)
#
# SAFETY
#   - Heavy temporary files go on LOCAL disk (SCRATCH_ROOT); only finished,
#     checked results are copied to the NAS, under a hidden name and then
#     renamed, so the fd/ fc/ ... folders only ever hold complete files.
#   - A session counts as finished only when its done marker exists.
#   - Workers are detached: Ctrl+C in the controller only stops new launches.
#   - Graceful stop: touch ${TEMPLATE_DIR}/STOP_5
#   - New launches pause while local free space is below MIN_FREE_GB.
#   - Do not edit or replace this file while it is running.
#
# Usage:
#   ./step5_fixels_v3.sh --build-mask                       # part 1, then check it
#   QC_KEEP=1 ./step5_fixels_v3.sh --list test5.txt         # one session, keep QC images
#   NTHREADS=8 MAX_PARALLEL_5=8 ./step5_fixels_v3.sh --list sessions.txt
#   ./step5_fixels_v3.sh --summary                          # QC table
#   SESSION_FMLS=0.06 NTHREADS=8 MAX_PARALLEL_5=8 ./step5_fixels_v3.sh --list sessions.txt
#       Same pipeline with a different peak cutoff for the session fixels. Results go to
#       fd_fmls0p06/ fc_fmls0p06/ log_fc_fmls0p06/ fdc_fmls0p06/ and step5_fmls0p06/, and
#       nothing from the standard run is touched. Use the same SESSION_FMLS with --summary.
#       Only one step 5 controller can run at a time.
###############################################################################

set -u
SELF="$(readlink -f "$0")"
: "${NTHREADS:=8}"; export NTHREADS
source "$(dirname "$SELF")/pipeline_config_v2.sh"

SUBJ_TEMPLATE_DIR="${TEMPLATE_DIR}/subject_templates"
TGROUP="${TEMPLATE_DIR}/group"
GT="${TGROUP}/wmfod_template.mif"
TMASK="${TGROUP}/template_mask.mif"
FIXEL_MASK="${TEMPLATE_DIR}/fixel_mask"
# SESSION_FMLS: peak cutoff for segmenting each SESSION's FOD into fixels (fod2fixel).
# Unset = fod2fixel's default (0.1), as in the MRtrix tutorial. If set (e.g. 0.06), all
# session results go to separate folders tagged with the value (fd_fmls0p06/ ...), so the
# standard results are never touched and the two versions can be compared.
SESSION_FMLS="${SESSION_FMLS:-}"; export SESSION_FMLS
TAG=""; [ -n "$SESSION_FMLS" ] && TAG="_fmls${SESSION_FMLS//./p}"
FD_DIR="${TEMPLATE_DIR}/fd${TAG}"; FC_DIR="${TEMPLATE_DIR}/fc${TAG}"
LOGFC_DIR="${TEMPLATE_DIR}/log_fc${TAG}"; FDC_DIR="${TEMPLATE_DIR}/fdc${TAG}"
SWARP_DIR="${TEMPLATE_DIR}/session_warps"          # shared: the warp does not depend on the cutoff
S5_MAIN="${TEMPLATE_DIR}/step5"                    # holds the shared fixel mask record
S5="${TEMPLATE_DIR}/step5${TAG}"
RESULTS="${S5}/step5_results.tsv"; DONE_DIR="${S5}/done"; LOG_DIR="${S5}/logs"; QC_DIR="${S5}/qc"
MASK_MARK="${S5_MAIN}/fixel_mask.done"; TL0="${S5_MAIN}/template_l0.mif"
STOP_FILE="${TEMPLATE_DIR}/STOP_5"

SCRATCH_ROOT="${SCRATCH_ROOT:-$HOME/fba5_work${TAG}}"
MAX_PAR="${MAX_PARALLEL_5:-8}"
NTHREADS_MASK="${NTHREADS_MASK:-32}"
FMLS="${FMLS_PEAK_VALUE:-0.06}"
MIN_R_COMPOSE="${MIN_R_COMPOSE:-0.9}"
MIN_FREE_GB="${MIN_FREE_GB:-30}"
POLL_SECS="${POLL_SECS:-20}"
STAGGER_SECS="${STAGGER_SECS:-30}"
QC_KEEP="${QC_KEEP:-0}"; export QC_KEEP
RESULT_LOCK="/tmp/${USER}_step5_results.lock"

minutes() { awk -v a="$1" -v b="$2" 'BEGIN{printf "%.1f", (b-a)/60}'; }
alive()   { [ -e "/proc/$1" ] && ! grep -q 'Z (zombie)' "/proc/$1/status" 2>/dev/null; }
grid3()   { mrinfo "$1" -size 2>/dev/null | awk '{print $1,$2,$3}'; }
nvol()    { mrinfo "$1" -size 2>/dev/null | awk '{print $4}'; }
ndim()    { mrinfo "$1" -ndim 2>/dev/null; }
nfix()    { mrinfo "$1" -size 2>/dev/null | awk '{print $1}'; }
avail_gb(){ df -BG --output=avail "$1" 2>/dev/null | tail -1 | tr -dc '0-9'; }
is_done() { [ -f "${DONE_DIR}/$1.done" ]; }
is_running() {
    local u; u=$(id -u)
    pgrep -u "$u" -f "step5_fixels_v3.sh --worker $1\$" >/dev/null && return 0
    pgrep -u "$u" -f "${SCRATCH_ROOT}/$1/" >/dev/null && return 0
    return 1
}
pearson() {  # a b mask -> r within mask (same method as qc4a.sh)
    local a=$1 b=$2 m=$3 am as bm bs abm
    am=$(mrstats "$a" -mask "$m" -output mean) || return 1
    as=$(mrstats "$a" -mask "$m" -output std)  || return 1
    bm=$(mrstats "$b" -mask "$m" -output mean) || return 1
    bs=$(mrstats "$b" -mask "$m" -output std)  || return 1
    abm=$(mrcalc "$a" "$b" -mult - -quiet | mrstats - -mask "$m" -output mean) || return 1
    awk -v am="$am" -v as="$as" -v bm="$bm" -v bs="$bs" -v abm="$abm" \
        'BEGIN{ if (as*bs==0) exit 1; printf "%.4f", (abm-am*bm)/(as*bs) }'
}
outputs_of() {  # every file part 2 publishes for one session
    printf '%s\n' "${FD_DIR}/$1.mif" "${FC_DIR}/$1.mif" "${LOGFC_DIR}/$1.mif" "${FDC_DIR}/$1.mif"
    [ -z "$TAG" ] && printf '%s\n' "${SWARP_DIR}/$1.mif"
}
unpublish() { local f; while read -r f; do rm -f "$f" "$(dirname "$f")/.$(basename "$f").tmp"; done < <(outputs_of "$1"); rm -f "${DONE_DIR}/$1.done"; }
publish() {  # src dest: copy under a hidden name, then rename (rename is atomic)
    local tmp; tmp="$(dirname "$2")/.$(basename "$2").tmp"
    cp "$1" "$tmp" && mv -f "$tmp" "$2"
}
NF=NA; FDM=NA; FDZ=NA; FCM=NA; RT=NA; RC=NA
record_row() {  # sess t0 t1 status
    local line
    line=$(printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' "$1" "$(date -d "@$2" '+%F %T')" "$(date -d "@$3" '+%F %T')" \
        "$(minutes "$2" "$3")" "$NTHREADS" "$NF" "$FDM" "$FDZ" "$FCM" "$RT" "$RC" "$4")
    ( flock 9; printf '%s\n' "$line" >> "$RESULTS" ) 9>"$RESULT_LOCK"
}

###############################################################################
# PART 1: template fixel mask (tutorial step 12)
###############################################################################
check_template() {
    [ -f "$GT" ] && [ "$(nvol "$GT")" = "45" ] || die "Group template missing or not 45 volumes: ${GT} (step 4b)"
    [ -f "$TMASK" ] || die "Final analysis mask missing: ${TMASK} (step 4c)"
    [ "$(grid3 "$TMASK")" = "$(grid3 "$GT")" ] || die "Template mask is not on the template grid"
}
mask_is_current() {
    [ -f "$MASK_MARK" ] || return 1
    [ "$GT" -nt "$MASK_MARK" ] && return 1
    [ "$TMASK" -nt "$MASK_MARK" ] && return 1
    return 0
}
view_help() {
    cat <<EOF

Check the fixel mask before running part 2 (changing it later means redoing every session):
  1. Copy to your Mac:  ${TL0}  and the folder  ${FIXEL_MASK}
  2. mrview template_l0.mif, then Tools > Fixel plot > open fixel_mask/index.mif
  3. Look for missing fixels, especially where bundles cross, and along the
     fornix and uncinate. Missing fixels: rebuild with a lower value, e.g.
       FMLS_PEAK_VALUE=0.05 ./step5_fixels_v3.sh --rebuild-mask
EOF
}
build_mask() {
    check_template
    if mask_is_current; then log "Fixel mask already built: $(cat "$MASK_MARK")"; view_help; return 0; fi
    if [ -f "$MASK_MARK" ]; then die "The template or template mask changed after the fixel mask was built. Use --rebuild-mask."; fi
    [ -e "$FIXEL_MASK" ] && die "${FIXEL_MASK} exists but was not made by this script. Move it aside first."
    local d
    for d in "$FD_DIR" "$FC_DIR" "$LOGFC_DIR" "$FDC_DIR"; do
        [ -e "$d" ] && die "${d} already exists from an earlier run. Move it aside first."
    done
    mkdir -p "$S5" "$DONE_DIR" "$LOG_DIR" "$SWARP_DIR"
    local tmp="${TEMPLATE_DIR}/.fixel_mask_building"; rm -rf "$tmp"
    log "Segmenting the group template into fixels (-fmls_peak_value ${FMLS}, ${NTHREADS_MASK} threads)..."
    fod2fixel -mask "$TMASK" -fmls_peak_value "$FMLS" "$GT" "$tmp" -nthreads "$NTHREADS_MASK" -force \
        || { rm -rf "$tmp"; die "fod2fixel failed on the template"; }
    [ -f "$tmp/index.mif" ] && [ -f "$tmp/directions.mif" ] || { rm -rf "$tmp"; die "fixel mask incomplete"; }
    mv "$tmp" "$FIXEL_MASK"
    for d in "$FD_DIR" "$FC_DIR" "$LOGFC_DIR" "$FDC_DIR"; do
        mkdir -p "$d"; cp "${FIXEL_MASK}/index.mif" "${FIXEL_MASK}/directions.mif" "$d/"
    done
    mrconvert "$GT" -coord 3 0 "$TL0" -force -quiet || die "could not make ${TL0}"
    local n; n=$(nfix "${FIXEL_MASK}/directions.mif")
    printf '%s fixels | -fmls_peak_value %s | built %s\n' "$n" "$FMLS" "$(date '+%F %T')" > "$MASK_MARK"
    log "Template fixel mask: ${n} fixels"
    log "  (the tutorial expects several hundred thousand at 1.25 mm; 0.8 mm has about 3.8x as many voxels, so expect more)"
    view_help
}
rebuild_mask() {
    local old="${S5}/old_$(date '+%Y%m%d_%H%M%S')" d n=0
    n=$(ls "$DONE_DIR"/*.done 2>/dev/null | wc -l)
    echo "This moves the fixel mask, all session results (${n} finished sessions) and their warps to:"
    echo "  ${old}"
    if [ "$ASSUME_YES" -eq 0 ]; then read -r -p "Type REBUILD to continue: " ok; [ "$ok" = "REBUILD" ] || die "Aborted."; fi
    mkdir -p "$old"
    for d in "$FIXEL_MASK" "$FD_DIR" "$FC_DIR" "$LOGFC_DIR" "$FDC_DIR" "$SWARP_DIR" "$DONE_DIR" "$RESULTS" "$MASK_MARK"; do
        [ -e "$d" ] && mv "$d" "$old/"
    done
    build_mask
}

###############################################################################
# QC SUMMARY
###############################################################################
summary() {
    [ -f "$RESULTS" ] || { echo "No results yet: ${RESULTS}"; return 0; }
    local latest; latest=$(mktemp)
    tail -n +2 "$RESULTS" | awk -F'\t' '{last[$1]=$0} END{for (s in last) print last[s]}' | sort > "$latest"
    echo "Sessions in results: $(wc -l < "$latest")"
    cut -f12 "$latest" | sort | uniq -c | awk '{printf "  %-22s %s\n", $2, $1}'
    awk -F'\t' '$12=="OK"' "$latest" > "${latest}.ok"
    if [ -s "${latest}.ok" ]; then
        echo; echo "Lowest r_template (alignment to the group template):"
        sort -t$'\t' -k10,10g "${latest}.ok" | head -5 | awk -F'\t' '{printf "  %-18s %s\n", $1, $10}'
        echo "Highest fd_zero_pct (template fixels with no matching fixel):"
        sort -t$'\t' -k8,8gr "${latest}.ok" | head -5 | awk -F'\t' '{printf "  %-18s %s%%\n", $1, $8}'
        echo "Lowest r_compose (one-step vs two-step warp; expect about 0.98 or more):"
        sort -t$'\t' -k11,11g "${latest}.ok" | head -3 | awk -F'\t' '{printf "  %-18s %s\n", $1, $11}'
        awk -F'\t' 'NR==1{mn=mx=$9} {if($9<mn)mn=$9; if($9>mx)mx=$9; s+=$9} END{printf "fc_mean across sessions: min %s  mean %.3f  max %s\n", mn, s/NR, mx}' "${latest}.ok"
    fi
    if awk -F'\t' '$12!="OK"' "$latest" | grep -q .; then
        echo; echo "Not OK (latest attempt):"; awk -F'\t' '$12!="OK"{printf "  %-18s %s\n", $1, $12}' "$latest"
    fi
    rm -f "$latest" "${latest}.ok"
}

###############################################################################
# PART 2 WORKER: one session (tutorial steps 13 to 18)
###############################################################################
run_worker() {
    local sess=$1 subj=${1%%_*}
    local W="${SCRATCH_ROOT}/$1" SD FOD W4A H2G HT t0 nt
    SD="${SUBJ_TEMPLATE_DIR}/${subj}"; FOD="${OUT_ROOT}/${sess}/wmfod_norm.mif"
    W4A="${SD}/warps/${sess}.mif"; H2G="${SD}/halfway2group_warp.mif"; HT="${SD}/halfway_template.mif"
    mkdir -p "$LOG_DIR"; exec >> "${LOG_DIR}/${sess}.log" 2>&1
    log "worker start: ${sess}, ${NTHREADS} threads"
    t0=$(date +%s)
    fail() { log "$1: $2"; record_row "$sess" "$t0" "$(date +%s)" "$1"; [ "${KEEP_FAILED:-0}" = 1 ] || rm -rf "$W"; exit 1; }

    unpublish "$sess"
    rm -rf "$W"; mkdir -p "$W" || fail FAIL_SCRATCH "cannot create ${W}"
    nt=$(nfix "${FIXEL_MASK}/directions.mif")

    # a) session -> halfway (4a, 5D full warp) as a plain deformation field
    local nd; nd=$(ndim "$W4A")
    if [ "$nd" = "5" ]; then
        warpconvert "$W4A" warpfull2deformation "$W/s2h.mif" -template "$HT" -force -quiet || fail FAIL_CONVERT "warpconvert failed"
    else
        cp "$W4A" "$W/s2h.mif"
    fi
    [ "$(ndim "$W/s2h.mif")" = "4" ] && [ "$(grid3 "$W/s2h.mif")" = "$(grid3 "$HT")" ] \
        || fail FAIL_CONVERT "session->halfway field is not 4D on the halfway grid"
    #    ... chained with halfway -> group (4c), applied in that order
    transformcompose "$W/s2h.mif" "$H2G" "$W/s2g.mif" -template "$GT" -force -quiet || fail FAIL_COMPOSE "transformcompose failed"
    [ "$(ndim "$W/s2g.mif")" = "4" ] && [ "$(grid3 "$W/s2g.mif")" = "$(grid3 "$GT")" ] \
        || fail FAIL_COMPOSE "combined warp is not 4D on the group grid"

    # b) tutorial step 13: warp the FOD, no reorientation yet
    mrtransform "$FOD" -warp "$W/s2g.mif" -reorient_fod no "$W/fod_nr.mif" -nthreads "$NTHREADS" -force -quiet \
        || fail FAIL_WARP_FOD "mrtransform failed"
    [ "$(nvol "$W/fod_nr.mif")" = "45" ] && [ "$(grid3 "$W/fod_nr.mif")" = "$(grid3 "$GT")" ] \
        || fail FAIL_WARP_FOD "warped FOD is not 45 volumes on the group grid"

    # QC: the one-step result must match applying the two warps in turn
    mrconvert "$W/fod_nr.mif" -coord 3 0 "$W/l0_onestep.mif" -force -quiet || fail FAIL_QC "l0 extract"
    mrconvert "$FOD" -coord 3 0 "$W/l0_native.mif" -force -quiet || fail FAIL_QC "l0 extract"
    if [ "$nd" = "5" ]; then
        mrtransform "$W/l0_native.mif" -warp_full "$W4A" -template "$HT" "$W/l0_half.mif" -force -quiet || fail FAIL_QC "two-step 1"
    else
        mrtransform "$W/l0_native.mif" -warp "$W4A" "$W/l0_half.mif" -force -quiet || fail FAIL_QC "two-step 1"
    fi
    mrtransform "$W/l0_half.mif" -warp "$H2G" "$W/l0_twostep.mif" -force -quiet || fail FAIL_QC "two-step 2"
    RC=$(pearson "$W/l0_onestep.mif" "$W/l0_twostep.mif" "$TMASK") || RC=NA
    RT=$(pearson "$W/l0_onestep.mif" "$TL0" "$TMASK") || RT=NA
    if [ "$RC" = "NA" ] || awk -v r="$RC" -v m="$MIN_R_COMPOSE" 'BEGIN{exit !(r<m)}'; then
        fail FAIL_COMPOSE_CHECK "one-step vs two-step r=${RC} (minimum ${MIN_R_COMPOSE}): the combined warp is wrong"
    fi
    log "r_compose=${RC}  r_template=${RT}"

    # c) tutorial step 14: segment into fixels; FD is the fixel integral
    fod2fixel -mask "$TMASK" ${SESSION_FMLS:+-fmls_peak_value "$SESSION_FMLS"} "$W/fod_nr.mif" "$W/fixel_nr" -afd fd.mif -nthreads "$NTHREADS" -force -quiet \
        || fail FAIL_FOD2FIXEL "fod2fixel failed"
    [ -f "$W/fixel_nr/fd.mif" ] && [ -f "$W/fixel_nr/directions.mif" ] || fail FAIL_FOD2FIXEL "fixel output incomplete"
    NF=$(nfix "$W/fixel_nr/directions.mif")
    [ "${NF:-0}" -gt 0 ] 2>/dev/null || fail FAIL_FOD2FIXEL "no fixels found"
    rm -f "$W/fod_nr.mif"     # about 2 GB, nothing later needs it

    # d) tutorial step 15: rotate fixel directions by the warp
    fixelreorient "$W/fixel_nr" "$W/s2g.mif" "$W/fixel" -nthreads "$NTHREADS" -force -quiet \
        || fail FAIL_REORIENT "fixelreorient failed"
    [ -f "$W/fixel/fd.mif" ] || fail FAIL_REORIENT "reoriented fixels incomplete"

    # e) tutorial step 16: match session fixels to the template slots (no match = FD 0)
    fixelcorrespondence "$W/fixel/fd.mif" "$FIXEL_MASK" "$W/out_fd" "${sess}.mif" -nthreads "$NTHREADS" -quiet \
        || fail FAIL_CORRESPOND "fixelcorrespondence failed"
    [ "$(nfix "$W/out_fd/${sess}.mif")" = "$nt" ] || fail FAIL_CORRESPOND "FD file does not have ${nt} fixels"

    # f) tutorial step 17: FC from the warp, then log FC (log FC is what goes into stats)
    warp2metric "$W/s2g.mif" -fc "$FIXEL_MASK" "$W/out_fc" "${sess}.mif" -nthreads "$NTHREADS" -quiet \
        || fail FAIL_FC "warp2metric failed"
    [ "$(nfix "$W/out_fc/${sess}.mif")" = "$nt" ] || fail FAIL_FC "FC file does not have ${nt} fixels"
    mrcalc "$W/out_fc/${sess}.mif" -log "$W/logfc.mif" -force -quiet || fail FAIL_LOGFC "mrcalc -log failed"

    # g) tutorial step 18: FDC = FD x FC
    mrcalc "$W/out_fd/${sess}.mif" "$W/out_fc/${sess}.mif" -mult "$W/fdc.mif" -force -quiet || fail FAIL_FDC "mrcalc -mult failed"
    [ "$(nfix "$W/logfc.mif")" = "$nt" ] && [ "$(nfix "$W/fdc.mif")" = "$nt" ] || fail FAIL_FDC "log FC or FDC has the wrong size"

    # QC numbers
    FDM=$(mrstats "$W/out_fd/${sess}.mif" -output mean 2>/dev/null | tr -d ' ')
    FDZ=$(mrcalc "$W/out_fd/${sess}.mif" 0 -eq - -quiet | mrstats - -output mean 2>/dev/null | awk '{printf "%.2f", $1*100}')
    FCM=$(mrstats "$W/out_fc/${sess}.mif" -output mean 2>/dev/null | tr -d ' ')

    # Publish (only complete, checked files reach the NAS), done marker last
    publish "$W/out_fd/${sess}.mif" "${FD_DIR}/${sess}.mif"  || fail FAIL_PUBLISH "fd"
    publish "$W/out_fc/${sess}.mif" "${FC_DIR}/${sess}.mif"  || fail FAIL_PUBLISH "fc"
    publish "$W/logfc.mif"          "${LOGFC_DIR}/${sess}.mif" || fail FAIL_PUBLISH "log_fc"
    publish "$W/fdc.mif"            "${FDC_DIR}/${sess}.mif" || fail FAIL_PUBLISH "fdc"
    [ -z "$TAG" ] && { publish "$W/s2g.mif" "${SWARP_DIR}/${sess}.mif" || fail FAIL_PUBLISH "warp"; }
    if [ "$QC_KEEP" = "1" ]; then
        rm -rf "${QC_DIR}/${sess}"; mkdir -p "${QC_DIR}/${sess}"
        cp "$W/l0_onestep.mif" "$W/l0_twostep.mif" "${QC_DIR}/${sess}/"
        cp -r "$W/fixel" "${QC_DIR}/${sess}/session_fixels"
        log "QC images kept in ${QC_DIR}/${sess}"
    fi
    touch "${DONE_DIR}/${sess}.done"
    record_row "$sess" "$t0" "$(date +%s)" "OK"
    rm -rf "$W"
    log "OK in $(minutes "$t0" "$(date +%s)") min"
    return 0
}

if [ "${1:-}" = "--worker" ]; then run_worker "$2"; exit $?; fi

###############################################################################
# CONTROLLER
###############################################################################
LIST=""; MODE=""; ASSUME_YES=0
while [ $# -gt 0 ]; do
    case "$1" in
        --list) LIST="$(readlink -f "$2")"; MODE=run; shift ;;
        --build-mask) MODE=mask ;;
        --rebuild-mask) MODE=rebuild ;;
        --summary) MODE=summary ;;
        --yes) ASSUME_YES=1 ;;
        *) die "Unknown option: $1" ;;
    esac; shift
done
case "$MODE" in
    mask|rebuild) [ -z "$TAG" ] || die "The fixel mask is shared by all versions. Run this without SESSION_FMLS." ;;
esac
case "$MODE" in
    mask) build_mask; exit $? ;;
    rebuild) check_template; rebuild_mask; exit $? ;;
    summary) summary; exit 0 ;;
    run) ;;
    *) die "Usage: $0 --build-mask | --list sessions.txt | --summary | --rebuild-mask" ;;
esac

[ -f "$LIST" ] || die "Session list not found: ${LIST}"
check_template
if [ -f "$MASK_MARK" ] && ! mask_is_current; then die "The template or template mask changed after the fixel mask was built. Use --rebuild-mask."; fi
mask_is_current || die "No fixel mask yet. Run: $0 --build-mask (and check it) first."
mapfile -t SESSIONS < <(grep -o 'HCA[0-9]*_V[0-9]*_MR' "$LIST" | awk '!seen[$0]++')
[ "${#SESSIONS[@]}" -gt 0 ] || die "No sessions in ${LIST}"

problems=()
for s in "${SESSIONS[@]}"; do
    sd="${SUBJ_TEMPLATE_DIR}/${s%%_*}"
    [ -f "${OUT_ROOT}/${s}/wmfod_norm.mif" ] || problems+=("${s}: no FOD (step 3a)")
    [ -f "${sd}/warps/${s}.mif" ]             || problems+=("${s}: no 4a warp")
    [ -f "${sd}/halfway_template.mif" ]       || problems+=("${s}: no halfway template (4a)")
    [ -f "${sd}/halfway2group.done" ] && [ -f "${sd}/halfway2group_warp.mif" ] || problems+=("${s}: 4c not finished")
done
if [ "${#problems[@]}" -gt 0 ]; then printf '  %s\n' "${problems[@]}"; die "${#problems[@]} problem(s) above."; fi

mkdir -p "$S5" "$DONE_DIR" "$LOG_DIR" "$SWARP_DIR" "$SCRATCH_ROOT"
for d in "$FD_DIR" "$FC_DIR" "$LOGFC_DIR" "$FDC_DIR"; do
    mkdir -p "$d"; [ -f "$d/index.mif" ] || cp "${FIXEL_MASK}/index.mif" "${FIXEL_MASK}/directions.mif" "$d/"
done
[ -f "$RESULTS" ] || printf 'session\tstart\tend\tminutes\tthreads\tn_fixels\tfd_mean\tfd_zero_pct\tfc_mean\tr_template\tr_compose\tstatus\n' > "$RESULTS"

CTRL_LOCK="/tmp/${USER}_step5_controller.pid"
if [ -f "$CTRL_LOCK" ] && alive "$(cat "$CTRL_LOCK")"; then die "A step 5 controller is already running (PID $(cat "$CTRL_LOCK"))."; fi
echo $$ > "$CTRL_LOCK"
[ -f "$STOP_FILE" ] && { rm -f "$CTRL_LOCK"; die "Stop file exists (${STOP_FILE}). Delete it first."; }

declare -A LANE=()
on_exit() { rm -f "$CTRL_LOCK"; [ "${#LANE[@]}" -gt 0 ] && log "Controller stopped; ${#LANE[@]} session(s) keep running and will finish on their own."; }
trap on_exit EXIT
trap 'log "Interrupted: no new sessions will start."; exit 130' INT TERM

QUEUE=(); n_done=0; n_busy=0; n_orphan=0
for s in "${SESSIONS[@]}"; do
    if is_done "$s"; then n_done=$((n_done+1))
    elif is_running "$s"; then n_busy=$((n_busy+1))
    else
        while read -r f; do [ -e "$f" ] && n_orphan=$((n_orphan+1)); done < <(outputs_of "$s")
        unpublish "$s"; QUEUE+=("$s")
    fi
done
cat <<EOF

-------------------------------------------------------------------
Step 5: fixels, FD, FC, log FC, FDC     $(cat "$MASK_MARK")
  sessions:        ${#SESSIONS[@]}   finished: ${n_done}   running elsewhere: ${n_busy}   to do: ${#QUEUE[@]}
  at once:         ${MAX_PAR}, ${NTHREADS} threads each
  local scratch:   ${SCRATCH_ROOT} ($(avail_gb "$SCRATCH_ROOT") GB free; launches pause below ${MIN_FREE_GB} GB)
  session cutoff:  ${SESSION_FMLS:-0.1 (fod2fixel default, as in the tutorial)}
  outputs go to:   ${FD_DIR}  (and fc, log_fc, fdc with the same ending)
  keep QC images:  $([ "$QC_KEEP" = 1 ] && echo yes || echo no)
  results:         ${RESULTS}
  stop gracefully: touch ${STOP_FILE}
EOF
[ "$n_orphan" -gt 0 ] && echo "  removed ${n_orphan} leftover output file(s) from unfinished sessions"
echo "-------------------------------------------------------------------"
[ "${#QUEUE[@]}" -gt 0 ] || { log "Nothing to do."; summary; exit 0; }
if [ "$ASSUME_YES" -eq 0 ]; then read -r -p "Start? [y/N] " ok; [[ "$ok" =~ ^[Yy]$ ]] || die "Aborted."; fi

reap() {
    local p m
    for p in "${!LANE[@]}"; do
        alive "$p" && continue
        m=$(tail -1 "${LOG_DIR}/${LANE[$p]}.log" 2>/dev/null)
        [ -n "$m" ] || m="no log written: the worker did not start"
        log "FINISHED ${LANE[$p]}: ${m}"; unset "LANE[$p]"
    done
}
i=0; stopping=0; low_space_noted=0
while true; do
    reap
    if [ -f "$STOP_FILE" ] && [ "$stopping" -eq 0 ]; then log "Stop file found: nothing new will start."; stopping=1; fi
    if [ "$stopping" -eq 0 ]; then
        while [ "$i" -lt "${#QUEUE[@]}" ] && [ "${#LANE[@]}" -lt "$MAX_PAR" ]; do
            if [ "$(avail_gb "$SCRATCH_ROOT")" -lt "$MIN_FREE_GB" ]; then
                [ "$low_space_noted" -eq 0 ] && log "Local free space below ${MIN_FREE_GB} GB: waiting before launching more."
                low_space_noted=1; break
            fi
            low_space_noted=0
            s="${QUEUE[$i]}"; i=$((i+1))
            if is_done "$s" || is_running "$s"; then log "skip ${s}: already done or in progress"; continue; fi
            setsid bash "$SELF" --worker "$s" < /dev/null > /dev/null 2>&1 &
            LANE[$!]="$s"; log "START ${s} (running: ${#LANE[@]} of ${MAX_PAR})"
            sleep "$STAGGER_SECS"; reap
        done
    fi
    if [ "${#LANE[@]}" -eq 0 ] && { [ "$i" -ge "${#QUEUE[@]}" ] || [ "$stopping" -eq 1 ]; }; then break; fi
    sleep "$POLL_SECS"
done

fin=0; for s in "${SESSIONS[@]}"; do is_done "$s" && fin=$((fin+1)); done
log "Finished: ${fin} of ${#SESSIONS[@]} sessions."
[ "$fin" -lt "${#SESSIONS[@]}" ] && log "Rerun the same command to retry the rest (finished sessions are skipped)."
echo; summary
