#!/bin/bash
###############################################################################
# step8_dti_tracts_v3.sh  -  DTI values (FA, MD, RD, AD) in each tract, per session
#
# The DTI arm, using the SAME tracts as the FBA arm, on the group template:
#   1. masks (once)  For each tract (FX_left, FX_right, UF_left, UF_right, and the
#                    left+right unions FX_both, UF_both), two kinds of voxel mask:
#                      <tract>_ts.mif     TractSeg's own voxel mask (step 6), limited to the
#                                         80-subject analysis mask (template_mask)
#                      <tract>_ge<N>.mif  voxels holding at least one tract fixel at
#                                         cutoff N (1, 3, 10): the FBA tract definition
#                                         turned into voxels (fixel2voxel ... max).
#                    ge3 matches the FBA primary mask; the others are for sensitivity.
#   2. per session   warp fa/md/rd/ad (step 3c) onto the template with the session's
#                    combined warp from step 5 (session_warps/), linear interpolation
#                    (no overshoot at edges), then mean and median of each metric in
#                    every mask, plus the % of mask voxels where FA is 0 (= outside
#                    that session's brain mask after warping; should be about 0).
#
#   READS   ${OUT_ROOT}/<SESS>/dti/{fa,md,rd,ad}.mif         (step 3c)
#           ${TEMPLATE_DIR}/session_warps/<SESS>.mif         (step 5)
#           ${TEMPLATE_DIR}/tracts/voxel_masks, tracts/fixels (step 6)
#   WRITES  ${TEMPLATE_DIR}/dti_tracts/
#             masks/            the 24 voxel masks, and mask_voxels.tsv (their sizes)
#             rows/<SESS>.tsv   one file per finished session (written whole, then renamed)
#             dti_by_session.tsv  all sessions together (rebuilt at the end of every run)
#             maps/<SESS>/      warped maps, only with KEEP_MAPS=1 (about 170 MB a session)
#
# Resumable: a session with rows/<SESS>.tsv is skipped. Ctrl+C stops it; rerun the
# same command to continue. Run it in tmux.
#
# Usage:
#   ./step8_dti_tracts_v3.sh --masks-only                 # (re)build and list the masks; any earlier
#                                                         # session results are moved aside, not mixed
#   ./step8_dti_tracts_v3.sh --list test8.txt             # a few sessions first
#   MAX_PARALLEL_8=8 ./step8_dti_tracts_v3.sh --list sessions.txt
#   ./step8_dti_tracts_v3.sh --summary
#   ./step8_dti_tracts_v3.sh --list sessions.txt --redo   # recompute every listed session
###############################################################################

set -u
SELF="$(readlink -f "$0")"
: "${NTHREADS:=4}"; export NTHREADS
source "$(dirname "$SELF")/pipeline_config_v2.sh"

TGROUP="${TEMPLATE_DIR}/group"; TMASK="${TGROUP}/template_mask.mif"
TRACTS_DIR="${TEMPLATE_DIR}/tracts"; SWARP_DIR="${TEMPLATE_DIR}/session_warps"
OUT8="${TEMPLATE_DIR}/dti_tracts"; MASKS="${OUT8}/masks"; ROWS="${OUT8}/rows"; LOGS="${OUT8}/logs"
RESULTS="${OUT8}/dti_by_session.tsv"; MASK_TABLE="${MASKS}/mask_voxels.tsv"
CUTOFFS="${TDI_CUTOFFS:-1 3 10}"
SCRATCH_ROOT="${SCRATCH_ROOT:-$HOME/fba8_work}"; MAX_PAR="${MAX_PARALLEL_8:-8}"; KEEP_MAPS="${KEEP_MAPS:-0}"
export SCRATCH_ROOT KEEP_MAPS
LOCK="/tmp/${USER}_step8.lock"
HEADER='session\ttract\tmask\tvoxels\tfa_mean\tmd_mean\trd_mean\tad_mean\tfa_median\tmd_median\trd_median\tad_median\tfa_zero_pct'

grid3() { mrinfo "$1" -size 2>/dev/null | awk '{print $1,$2,$3}'; }
count() { mrstats "$1" -mask "$1" -output count 2>/dev/null | tr -d ' '; }
four()  { tr -s ' \t\n' ' ' | sed 's/^ //; s/ $//' | awk 'NF==4{printf "%s\t%s\t%s\t%s", $1,$2,$3,$4; ok=1} END{if(!ok) printf "NA\tNA\tNA\tNA"}'; }

###############################################################################
# 1. Masks
###############################################################################
build_masks() {
    local t b c tmp="${OUT8}/.masks_building" n
    [ -f "${TRACTS_DIR}/tracts.done" ] || die "Step 6 tracts missing"
    rm -rf "$tmp"; mkdir -p "$tmp"
    local n0 n1
    for t in FX_left FX_right UF_left UF_right; do
        [ -f "${TRACTS_DIR}/voxel_masks/${t}.mif" ] || die "No TractSeg voxel mask for ${t}"
        # limit TractSeg's mask to the 80-subject analysis mask, as the fixel-based masks already are,
        # so no voxel lies outside any session's brain after warping (those read as FA = 0)
        mrcalc "${TRACTS_DIR}/voxel_masks/${t}.mif" "$TMASK" -mult "${tmp}/${t}_ts.mif" -datatype bit -force -quiet \
            || die "could not restrict ${t} to the analysis mask"
        n0=$(count "${TRACTS_DIR}/voxel_masks/${t}.mif"); n1=$(count "${tmp}/${t}_ts.mif")
        [ "${n1:-0}" -gt 0 ] 2>/dev/null || die "${t}: nothing left of the TractSeg mask inside the analysis mask"
        log "  ${t} ts: ${n0} voxels -> ${n1} inside the analysis mask ($(echo "$n0 $n1" | awk '{printf "%.1f", 100*($1-$2)/$1}')% dropped)"
    done
    for b in FX UF; do
        mrcalc "${tmp}/${b}_left_ts.mif" "${tmp}/${b}_right_ts.mif" -max "${tmp}/${b}_both_ts.mif" -datatype bit -force -quiet || die "union failed"
    done
    for t in FX_left FX_right FX_both UF_left UF_right UF_both; do
        for c in $CUTOFFS; do
            fixel2voxel "${TRACTS_DIR}/fixels/${t}/${t}_ge${c}.mif" max "${tmp}/f2v.mif" -force -quiet || die "fixel2voxel failed for ${t} ge${c}"
            mrcalc "${tmp}/f2v.mif" 0 -gt "${tmp}/${t}_ge${c}.mif" -datatype bit -force -quiet || die "mask failed for ${t} ge${c}"
        done
    done
    rm -f "${tmp}/f2v.mif"
    printf 'mask\tvoxels\n' > "${tmp}/mask_voxels.tsv"
    for m in "${tmp}"/*.mif; do
        [ "$(grid3 "$m")" = "$(grid3 "$TMASK")" ] || die "$(basename "$m") is not on the template grid"
        n=$(count "$m"); [ "${n:-0}" -gt 0 ] 2>/dev/null || die "$(basename "$m") is empty"
        printf '%s\t%s\n' "$(basename "$m" .mif)" "$n" >> "${tmp}/mask_voxels.tsv"
    done
    rm -rf "$MASKS"; mv "$tmp" "$MASKS"; touch "${MASKS}/.done"
    if ls "${ROWS}"/HCA*.tsv >/dev/null 2>&1; then        # results made with other masks must not be mixed with new ones
        local old="${OUT8}/old_$(date '+%Y%m%d_%H%M%S')"; mkdir -p "$old"
        mv "$ROWS" "${old}/rows"; mkdir -p "$ROWS"; [ -f "$RESULTS" ] && mv "$RESULTS" "$old/"
        log "Masks were rebuilt: earlier session results moved to ${old}. Run --list again to recompute them."
    fi
    log "Masks built: $(ls "$MASKS"/*.mif | wc -l) in ${MASKS}"
    column -t -s$'\t' "$MASK_TABLE"
}

###############################################################################
# 2. One session
###############################################################################
run_session() {
    local s=$1 W="${SCRATCH_ROOT}/$1" m M name tract type n means meds z pct row
    mkdir -p "$LOGS"; exec >> "${LOGS}/${s}.log" 2>&1
    log "start ${s}"
    rm -rf "$W"; mkdir -p "$W"
    for m in fa md rd ad; do
        mrtransform "${OUT_ROOT}/${s}/dti/${m}.mif" -warp "${SWARP_DIR}/${s}.mif" -interp linear "${W}/${m}.mif" \
            -nthreads "$NTHREADS" -force -quiet || { log "FAIL warp ${m}"; rm -rf "$W"; return 1; }
    done
    [ "$(grid3 "${W}/fa.mif")" = "$(grid3 "$TMASK")" ] || { log "FAIL: warped FA not on the template grid"; rm -rf "$W"; return 1; }
    mrcat "${W}/fa.mif" "${W}/md.mif" "${W}/rd.mif" "${W}/ad.mif" -axis 3 "${W}/dti4.mif" -force -quiet \
        || { log "FAIL mrcat"; rm -rf "$W"; return 1; }
    : > "${W}/rows.tsv"
    while IFS=$'\t' read -r name n; do
        [ "$name" = "mask" ] && continue
        M="${MASKS}/${name}.mif"; tract=${name%_*}; type=${name##*_}
        means=$(mrstats "${W}/dti4.mif" -mask "$M" -output mean 2>/dev/null | four)
        meds=$(mrstats "${W}/dti4.mif" -mask "$M" -output median 2>/dev/null | four)
        mrcalc "${W}/fa.mif" 0 -eq "$M" -mult "${W}/z.mif" -datatype bit -force -quiet
        z=$(count "${W}/z.mif"); z=${z:-0}
        pct=$(echo "$z $n" | awk '{printf "%.2f", 100*$1/$2}')
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$s" "$tract" "$type" "$n" "$means" "$meds" "$pct" >> "${W}/rows.tsv"
    done < "$MASK_TABLE"
    grep -q $'\tNA\t' "${W}/rows.tsv" && log "WARNING: some statistics came back empty (NA)"
    if [ "$KEEP_MAPS" = "1" ]; then mkdir -p "${OUT8}/maps/${s}"; cp "${W}"/{fa,md,rd,ad}.mif "${OUT8}/maps/${s}/"; fi
    cp "${W}/rows.tsv" "${ROWS}/.${s}.tsv.tmp" && mv "${ROWS}/.${s}.tsv.tmp" "${ROWS}/${s}.tsv"
    rm -rf "$W"
    log "OK ${s}: $(wc -l < "${ROWS}/${s}.tsv") rows"
}
if [ "${1:-}" = "--worker" ]; then run_session "$2"; exit $?; fi

###############################################################################
# Summary and combining
###############################################################################
combine() {
    { printf "${HEADER}\n"; cat "${ROWS}"/HCA*.tsv 2>/dev/null; } > "${RESULTS}.tmp" && mv "${RESULTS}.tmp" "$RESULTS"
}
summary() {
    [ -f "$RESULTS" ] || { echo "No results yet: ${RESULTS}"; return 0; }
    echo "Sessions: $(tail -n +2 "$RESULTS" | cut -f1 | sort -u | wc -l)    rows: $(tail -n +2 "$RESULTS" | wc -l)"
    echo "Mean across sessions, per tract and mask (FA, MD in mm2/s) and the worst FA-zero %:"
    local na; na=$(tail -n +2 "$RESULTS" | awk -F'\t' '$5=="NA"' | wc -l)
    [ "$na" -gt 0 ] && echo "WARNING: ${na} rows have missing statistics (NA); they are left out of the means below. Check those sessions' logs."
    tail -n +2 "$RESULTS" | awk -F'\t' '$5!="NA" {k=$2" "$3; n[k]++; fa[k]+=$5; md[k]+=$6; if($13+0>mz[k]) mz[k]=$13+0}
        END{for(k in n) printf "  %-16s FA %.3f  MD %.6f  max FA-zero %.1f%%\n", k, fa[k]/n[k], md[k]/n[k], mz[k]}' | sort
}

###############################################################################
# Controller
###############################################################################
LIST=""; MODE=run; REDO=0
while [ $# -gt 0 ]; do
    case "$1" in
        --list) LIST="$(readlink -f "$2")"; shift ;;
        --masks-only) MODE=masks ;; --summary) MODE=summary ;; --redo) REDO=1 ;;
        *) die "Unknown option: $1" ;;
    esac; shift
done
[ "$MODE" = "summary" ] && { summary; exit 0; }
exec 9>"$LOCK"; flock -n 9 || die "Step 8 is already running."
mkdir -p "$OUT8" "$ROWS" "$LOGS" "$SCRATCH_ROOT"
if [ "$MODE" = "masks" ] || [ ! -f "${MASKS}/.done" ]; then build_masks; fi
[ "$MODE" = "masks" ] && exit 0
[ -f "$LIST" ] || die "Usage: $0 --list sessions.txt   (or --masks-only, --summary)"

mapfile -t SESSIONS < <(grep -o 'HCA[0-9]*_V[0-9]*_MR' "$LIST" | awk '!seen[$0]++')
problems=(); QUEUE=()
for s in "${SESSIONS[@]}"; do
    for m in fa md rd ad; do [ -f "${OUT_ROOT}/${s}/dti/${m}.mif" ] || problems+=("${s}: no ${m}.mif (step 3c)"); done
    [ -f "${SWARP_DIR}/${s}.mif" ] || problems+=("${s}: no session warp (step 5)")
done
if [ "${#problems[@]}" -gt 0 ]; then printf '  %s\n' "${problems[@]}"; die "${#problems[@]} problem(s) above."; fi
for s in "${SESSIONS[@]}"; do      # only after every check has passed
    [ "$REDO" = "1" ] && rm -f "${ROWS}/${s}.tsv"
    [ -f "${ROWS}/${s}.tsv" ] || QUEUE+=("$s")
done
log "Step 8: ${#SESSIONS[@]} sessions, $(( ${#SESSIONS[@]} - ${#QUEUE[@]} )) already done, ${#QUEUE[@]} to do, ${MAX_PAR} at once (${NTHREADS} threads each)"
if [ "${#QUEUE[@]}" -gt 0 ]; then
    printf '%s\n' "${QUEUE[@]}" | xargs -P "$MAX_PAR" -I{} bash "$SELF" --worker {}
fi
fail=0; for s in "${SESSIONS[@]}"; do [ -f "${ROWS}/${s}.tsv" ] || { fail=$((fail+1)); log "not finished: ${s} (see ${LOGS}/${s}.log)"; }; done
combine
log "Finished: $(( ${#SESSIONS[@]} - fail )) of ${#SESSIONS[@]} sessions. Results: ${RESULTS}"
[ "$fail" -gt 0 ] && log "Rerun the same command to retry the unfinished ones."
echo; summary
exit 0
