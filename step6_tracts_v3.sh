#!/bin/bash
###############################################################################
# step6_tracts_v3.sh  -  Define the fornix and uncinate on the group template
#
# Runs ONCE, on the group template. The result is shared by all 160 sessions and
# by both arms (FBA and DTI), so both use exactly the same definition of each tract.
#
# What it does (the same commands you tested by hand in step6_test):
#   1. peaks       sh2peaks -num 3 on the template, saved as NIfTI in the HCP/FSL
#                  orientation TractSeg expects (-strides -1,2,3,4)
#   2. TractSeg    tract_segmentation, endings_segmentation, TOM (TractSeg 2.x;
#                  peaks are its default input)
#   3. Tracking    only the tracts you need, NR_FIBERS streamlines each, saved as .tck
#   4. checks      streamline count and mean length per tract; left/right check:
#                  each _left tract must sit at a lower x than its _right partner
#                  (MRtrix x runs towards the right). A swap STOPS the script.
#   5. DTI masks   each TractSeg bundle mask on the template grid (voxel_masks/)
#   6. FBA masks   streamlines -> fixels (tck2fixel), then masks at each cutoff in
#                  TDI_CUTOFFS (fixels with at least that many streamlines), plus
#                  left+right unions (e.g. FX_both)
#
#   READS   group/wmfod_template.mif, group/template_mask.mif   (4b, 4c)
#           fixel_mask/                                        (step 5)
#   WRITES  ${TEMPLATE_DIR}/tracts/
#             tract_summary.tsv          one row per tract: streamlines, length, x, voxels, fixels
#             tractseg/                  TractSeg output, incl. TOM_trackings/<tract>.tck
#             voxel_masks/<tract>.mif    for the DTI arm (step 8)
#             fixels/<tract>/<tract>_ge<N>.mif   fixel masks for the stats (step 9)
#             tracts.done                settings used (for your methods)
#
# SAFETY: builds in a hidden folder and renames it only when every check passes,
# so a half-finished run never appears as 'tracts/'. Refuses to run twice at once.
# An existing result is never overwritten: --redo moves it aside first.
#
# Usage:
#   ./step6_tracts_v3.sh                         # build (about 15 to 30 minutes, mostly TractSeg)
#   ./step6_tracts_v3.sh --summary               # show the summary table
#   ./step6_tracts_v3.sh --redo                  # rebuild; if it succeeds, the old result is moved to tracts_old_<date>
#   ./step6_tracts_v3.sh --zero-table --list sessions.txt
#   SESSION_FMLS=0.06 ./step6_tracts_v3.sh --zero-table --list sessions.txt
#       per session and tract: % of tract fixels with FD = 0, mean FD, mean log FC
#
# Settings (environment): NR_FIBERS=5000  TDI_CUTOFFS="1 3 10"  PRIMARY_CUTOFF=3
#                         NTHREADS_PEAKS=32  NR_CPUS=32
###############################################################################

set -u
SELF="$(readlink -f "$0")"
source "$(dirname "$SELF")/pipeline_config_v2.sh"

TGROUP="${TEMPLATE_DIR}/group"; GT="${TGROUP}/wmfod_template.mif"; TMASK="${TGROUP}/template_mask.mif"
FIXEL_MASK="${TEMPLATE_DIR}/fixel_mask"; MASK_MARK="${TEMPLATE_DIR}/step5/fixel_mask.done"
OUT="${TEMPLATE_DIR}/tracts"; BUILD="${TEMPLATE_DIR}/.tracts_building"; MARK="${OUT}/tracts.done"
if declare -p TRACTS >/dev/null 2>&1; then TRACT_LIST=("${TRACTS[@]}"); else TRACT_LIST=(FX_left FX_right UF_left UF_right); fi
NR_FIBERS="${NR_FIBERS:-5000}"; CUTOFFS="${TDI_CUTOFFS:-1 3 10}"; PRIMARY="${PRIMARY_CUTOFF:-3}"
NTHREADS_PEAKS="${NTHREADS_PEAKS:-32}"; NR_CPUS="${NR_CPUS:-32}"
LOCK="/tmp/${USER}_step6.lock"

grid3() { mrinfo "$1" -size 2>/dev/null | awk '{print $1,$2,$3}'; }
nvol()  { mrinfo "$1" -size 2>/dev/null | awk '{print $4}'; }
cx()    { mrcentroid "$1" 2>/dev/null | awk '{print $1}'; }
count() { mrstats "$1" -mask "$1" -output count 2>/dev/null | tr -d ' '; }
pairs() {  # print the base name of every tract that has both _left and _right
    local t b
    for t in "${TRACT_LIST[@]}"; do
        b=${t%_left}; [ "$b" != "$t" ] || continue
        printf '%s\n' "${TRACT_LIST[@]}" | grep -qx "${b}_right" && echo "$b"
    done
}
summary() {
    [ -f "${OUT}/tract_summary.tsv" ] || { echo "No result yet: ${OUT}"; return 1; }
    cat "$MARK"; echo; column -t -s$'\t' "${OUT}/tract_summary.tsv"
}

###############################################################################
# Per-session zero table (read-only on step 5 results)
###############################################################################
zero_table() {
    local list=$1 TAG="" s t c M n z fd lf pct R
    [ -f "$MARK" ] || die "Build the tracts first."
    [ -n "${SESSION_FMLS:-}" ] && TAG="_fmls${SESSION_FMLS//./p}"
    [ -d "${TEMPLATE_DIR}/fd${TAG}" ] || die "No step 5 results in ${TEMPLATE_DIR}/fd${TAG}"
    R="${OUT}/zeros_by_session${TAG}.tsv"
    export MRTRIX_NTHREADS=4
    local -A N
    local all=("${TRACT_LIST[@]}"); for b in $(pairs); do all+=("${b}_both"); done
    for t in "${all[@]}"; do for c in $CUTOFFS; do N[$t$c]=$(count "${OUT}/fixels/$t/${t}_ge$c.mif"); done; done
    printf 'session\ttract\tcutoff\tfixels\tfd_zero_pct\tmean_fd\tmean_logfc\n' > "$R"
    log "Zero table for fd${TAG} -> ${R}"
    for s in $(grep -o 'HCA[0-9]*_V[0-9]*_MR' "$list" | sort -u); do
        [ -f "${TEMPLATE_DIR}/fd${TAG}/$s.mif" ] || { log "skip ${s}: no FD file"; continue; }
        for t in "${all[@]}"; do for c in $CUTOFFS; do
            M="${OUT}/fixels/$t/${t}_ge$c.mif"; n=${N[$t$c]}
            mrcalc "${TEMPLATE_DIR}/fd${TAG}/$s.mif" 0 -eq "$M" -mult "${OUT}/.z${TAG}.mif" -datatype bit -force -quiet
            z=$(count "${OUT}/.z${TAG}.mif"); z=${z:-0}
            fd=$(mrstats "${TEMPLATE_DIR}/fd${TAG}/$s.mif" -mask "$M" -output mean 2>/dev/null | tr -d ' '); fd=${fd:-NA}
            lf=$(mrstats "${TEMPLATE_DIR}/log_fc${TAG}/$s.mif" -mask "$M" -output mean 2>/dev/null | tr -d ' '); lf=${lf:-NA}
            pct=$(echo "$z $n" | awk '{printf "%.2f", 100*$1/$2}')
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$s" "$t" "$c" "$n" "$pct" "$fd" "$lf" >> "$R"
        done; done
    done
    rm -f "${OUT}/.z${TAG}.mif"
    log "Done: $(tail -n +2 "$R" | wc -l) rows in ${R}"
}

###############################################################################
# Build
###############################################################################
build() {
    local t b c n len x vox ok
    [ -f "$GT" ] && [ "$(nvol "$GT")" = "45" ] || die "Group template missing or not 45 volumes: ${GT}"
    [ -f "$TMASK" ] || die "Template mask missing: ${TMASK}"
    [ -f "${FIXEL_MASK}/index.mif" ] && [ -f "${FIXEL_MASK}/directions.mif" ] && [ -f "$MASK_MARK" ] \
        || die "Fixel mask missing: run step 5 --build-mask first"
    for c in sh2peaks mrconvert TractSeg Tracking tckstats mrcentroid tck2fixel mrthreshold mrstats mrcalc; do
        command -v "$c" >/dev/null || die "Command not found: ${c}"
    done
    command -v nvidia-smi >/dev/null && log "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)" \
        || log "No GPU found: TractSeg will run on the CPU (slower)"
    rm -rf "$BUILD"; mkdir -p "$BUILD"
    local TS="${BUILD}/tractseg" BUNDLES; BUNDLES=$(IFS=,; echo "${TRACT_LIST[*]}")

    log "1/6 Peaks from the template (${NTHREADS_PEAKS} threads)"
    sh2peaks -num 3 "$GT" "${BUILD}/peaks.mif" -nthreads "$NTHREADS_PEAKS" -force -quiet || die "sh2peaks failed"
    mrconvert "${BUILD}/peaks.mif" "${BUILD}/peaks.nii.gz" -strides -1,2,3,4 -force -quiet || die "mrconvert failed"
    [ "$(nvol "${BUILD}/peaks.nii.gz")" = "9" ] || die "peaks.nii.gz does not have 9 volumes"
    rm -f "${BUILD}/peaks.mif"

    log "2/6 TractSeg (3 runs; logs in ${OUT}/tractseg_*.log afterwards)"
    for c in tract_segmentation endings_segmentation TOM; do
        TractSeg -i "${BUILD}/peaks.nii.gz" -o "$TS" --output_type "$c" > "${BUILD}/tractseg_${c}.log" 2>&1 \
            || die "TractSeg ${c} failed. Log: ${BUILD}/tractseg_${c}.log"
        if grep -qi 'flip' "${BUILD}/tractseg_${c}.log"; then
            log "WARNING: TractSeg ${c} mentioned flipping the peaks:"; grep -i 'flip' "${BUILD}/tractseg_${c}.log" | head -3
        fi
    done
    for t in "${TRACT_LIST[@]}"; do
        [ -f "${TS}/bundle_segmentations/${t}.nii.gz" ] || die "No TractSeg mask for ${t}"
    done

    log "3/6 Tracking ${BUNDLES}, ${NR_FIBERS} streamlines each (${NR_CPUS} CPUs)"
    Tracking -i "${BUILD}/peaks.nii.gz" -o "$TS" --nr_fibers "$NR_FIBERS" --tracking_format tck \
        --bundles "$BUNDLES" --nr_cpus "$NR_CPUS" > "${BUILD}/tracking.log" 2>&1 \
        || die "Tracking failed. Log: ${BUILD}/tracking.log"
    for t in "${TRACT_LIST[@]}"; do
        [ -f "${TS}/TOM_trackings/${t}.tck" ] || die "No streamlines for ${t}"
    done

    log "4/6 Checks: streamlines, lengths, left/right"
    local -A X
    local brain_x; brain_x=$(cx "$TMASK")
    printf 'tract\tstreamlines\tmean_length_mm\tcentroid_x_mm\tvoxels' > "${BUILD}/tract_summary.tsv"
    for c in $CUTOFFS; do printf '\tfixels_ge%s' "$c" >> "${BUILD}/tract_summary.tsv"; done
    echo >> "${BUILD}/tract_summary.tsv"
    for t in "${TRACT_LIST[@]}"; do
        n=$(tckstats "${TS}/TOM_trackings/${t}.tck" -output count 2>/dev/null | tr -d ' ')
        len=$(tckstats "${TS}/TOM_trackings/${t}.tck" -output mean 2>/dev/null | tr -d ' ')
        X[$t]=$(cx "${TS}/bundle_segmentations/${t}.nii.gz")
        [ "${n:-0}" -ge $((NR_FIBERS / 10)) ] 2>/dev/null || die "${t}: only ${n:-0} streamlines (asked for ${NR_FIBERS})"
        log "  ${t}: ${n} streamlines, mean length ${len} mm, centroid x ${X[$t]} (brain centre ${brain_x})"
    done
    for b in $(pairs); do
        ok=$(awk -v l="${X[${b}_left]}" -v r="${X[${b}_right]}" 'BEGIN{print (l<r)?1:0}')
        [ "$ok" = "1" ] || die "${b}: left tract (x=${X[${b}_left]}) is not left of the right tract (x=${X[${b}_right]}). Left/right swapped: flip the peaks and rerun."
        awk -v l="${X[${b}_left]}" -v r="${X[${b}_right]}" -v m="$brain_x" 'BEGIN{exit !(l<m && r>m)}' \
            || log "  WARNING ${b}: the pair is ordered correctly but does not straddle the brain centre (x=${brain_x}). Check the picture."
        log "  ${b}: left/right OK"
    done

    log "5/6 Voxel masks for the DTI arm"
    mkdir -p "${BUILD}/voxel_masks"
    for t in "${TRACT_LIST[@]}"; do
        mrconvert "${TS}/bundle_segmentations/${t}.nii.gz" "${BUILD}/voxel_masks/${t}.mif" -strides "$TMASK" -datatype bit -force -quiet \
            || die "mask conversion failed for ${t}"
        [ "$(grid3 "${BUILD}/voxel_masks/${t}.mif")" = "$(grid3 "$TMASK")" ] || die "${t} voxel mask is not on the template grid"
    done

    log "6/6 Tract fixels and fixel masks at cutoffs: ${CUTOFFS}"
    mkdir -p "${BUILD}/fixels"
    local -A F
    for t in "${TRACT_LIST[@]}"; do
        tck2fixel "${TS}/TOM_trackings/${t}.tck" "$FIXEL_MASK" "${BUILD}/fixels/${t}" "${t}_tdi.mif" -force -quiet \
            || die "tck2fixel failed for ${t}"
        for c in $CUTOFFS; do
            mrthreshold "${BUILD}/fixels/${t}/${t}_tdi.mif" -abs "$c" "${BUILD}/fixels/${t}/${t}_ge${c}.mif" -force -quiet \
                || die "mrthreshold failed for ${t} at ${c}"
            F[$t$c]=$(count "${BUILD}/fixels/${t}/${t}_ge${c}.mif")
            [ "${F[$t$c]:-0}" -gt 0 ] 2>/dev/null || die "${t}: no fixels at cutoff ${c}"
        done
    done
    for b in $(pairs); do
        mkdir -p "${BUILD}/fixels/${b}_both"; cp "${FIXEL_MASK}/index.mif" "${FIXEL_MASK}/directions.mif" "${BUILD}/fixels/${b}_both/"
        for c in $CUTOFFS; do
            mrcalc "${BUILD}/fixels/${b}_left/${b}_left_ge${c}.mif" "${BUILD}/fixels/${b}_right/${b}_right_ge${c}.mif" -max \
                "${BUILD}/fixels/${b}_both/${b}_both_ge${c}.mif" -datatype bit -force -quiet || die "union failed for ${b}"
            F[${b}_both$c]=$(count "${BUILD}/fixels/${b}_both/${b}_both_ge${c}.mif")
        done
    done

    for t in "${TRACT_LIST[@]}"; do
        vox=$(count "${BUILD}/voxel_masks/${t}.mif")
        printf '%s\t%s\t%s\t%s\t%s' "$t" "$(tckstats "${TS}/TOM_trackings/${t}.tck" -output count 2>/dev/null | tr -d ' ')" \
            "$(tckstats "${TS}/TOM_trackings/${t}.tck" -output mean 2>/dev/null | tr -d ' ')" "${X[$t]}" "$vox" >> "${BUILD}/tract_summary.tsv"
        for c in $CUTOFFS; do printf '\t%s' "${F[$t$c]}" >> "${BUILD}/tract_summary.tsv"; done
        echo >> "${BUILD}/tract_summary.tsv"
    done
    for b in $(pairs); do
        printf '%s\t-\t-\t-\t-' "${b}_both" >> "${BUILD}/tract_summary.tsv"
        for c in $CUTOFFS; do printf '\t%s' "${F[${b}_both$c]}" >> "${BUILD}/tract_summary.tsv"; done
        echo >> "${BUILD}/tract_summary.tsv"
    done

    {
        echo "built $(date '+%F %T')"
        echo "TractSeg $(pip show TractSeg 2>/dev/null | awk '/^Version/{print $2}')  peaks: sh2peaks -num 3, strides -1,2,3,4"
        echo "tracts: ${TRACT_LIST[*]}   streamlines per tract: ${NR_FIBERS}"
        echo "fixel mask cutoffs (streamlines per fixel): ${CUTOFFS}   primary: ${PRIMARY}"
    } > "${BUILD}/tracts.done"
    if [ -d "$OUT" ]; then local old="${OUT}_old_$(date '+%Y%m%d_%H%M%S')"; mv "$OUT" "$old"; log "Previous result moved to ${old}"; fi
    mv "$BUILD" "$OUT"
    log "Step 6 finished: ${OUT}"
    summary
    cat <<EOF

Check before using it:
  - the picture: copy ${OUT}/tractseg/TOM_trackings/*.tck and look at them over the template
  - the counts above should be close to your step6_test numbers (tracking has a random element)
The primary fixel mask for statistics is <tract>_ge${PRIMARY}.mif; the others are for sensitivity checks.
EOF
}

###############################################################################
# Main
###############################################################################
MODE=build; LIST=""
while [ $# -gt 0 ]; do
    case "$1" in
        --summary) MODE=summary ;;
        --redo) MODE=redo ;;
        --zero-table) MODE=zero ;;
        --list) LIST="$(readlink -f "$2")"; shift ;;
        *) die "Unknown option: $1" ;;
    esac; shift
done
exec 9>"$LOCK"; flock -n 9 || die "Step 6 is already running."
case "$MODE" in
    summary) summary; exit $? ;;
    zero) [ -f "$LIST" ] || die "Use --zero-table --list sessions.txt"; zero_table "$LIST"; exit 0 ;;
    redo) build ;;   # the old result is moved aside only if the new build succeeds
    build) if [ -f "$MARK" ]; then log "Already built. Use --redo to rebuild."; summary; exit 0; fi
           [ -d "$OUT" ] && die "${OUT} exists but is incomplete. Move it aside, or use --redo."
           build ;;
esac
