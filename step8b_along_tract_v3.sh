#!/bin/bash
###############################################################################
# step8b_along_tract_v3.sh  -  Along-tract profiles (thesis plan section 1.7.4)
#
# Instead of one mean per tract, measure each metric at N positions along the tract
# (default 20), so you can see WHERE along the fornix and uncinate change happens.
#
# Everything is defined ONCE on the group template, so the segments are identical in
# every session and in both arms:
#   1. centreline    The step 6 streamlines of each tract (5000, from TractSeg's Tracking)
#                    are resampled to 100 points, aligned to run the same way, and
#                    averaged into one centreline. Only streamlines that span at least 75% of
#                    the median length are used (partial ones would blur it), and the average
#                    is repeated over the 80% of them closest to it. The centreline is then
#                    extended straight beyond both ends (EXTEND_MM, default 20) so the parts of
#                    a tract that fan out past its ends are spread over several segments.
#   2. orientation   Every centreline runs POSTERIOR -> ANTERIOR, so segment 1 is the
#                    posterior end and segment N the anterior end, for left and right.
#                    Fornix: segment 1 = hippocampal end, N = anterior (columns).
#                    Uncinate: segment 1 = temporal end, N = frontal end.
#   3. segments      Each voxel of the template is given the number of the centreline
#                    point nearest to it, 1..N in equal steps of length along the part of
#                    the centreline that the tract mask covers (the analysis mask can trim
#                    the ends of a tract, so the segments span what is actually measurable;
#                    "covered_mm" in the centreline table gives that length).
#                    The tract's voxel mask and fixel mask (cutoff ALONG_CUTOFF,
#                    default 3 = the FBA primary mask) are then cut into N pieces.
#                    Fixels take the segment of the voxel they sit in (voxel2fixel).
#   4. per session   DTI: FA, MD, RD, AD (step 3c maps, warped with the step 5 warp)
#                         mean in each voxel segment.
#                    FBA: FD, log FC, FDC (unsmoothed, from step 5) mean in each fixel
#                         segment, plus % of the segment's fixels with FD = 0.
#
#   READS   tracts/tractseg/TOM_trackings/<tract>.tck, tracts/fixels/<tract>/*_ge<N>.mif (step 6)
#           dti_tracts/masks/<tract>_ge<N>.mif (step 8), fd/ log_fc/ fdc/ (step 5),
#           session_warps/, ${OUT_ROOT}/<SESS>/dti/ (steps 5, 3c)
#   WRITES  ${TEMPLATE_DIR}/alongtract/
#             segments/                 labels, per-segment masks, centrelines, segments.tsv
#             alongtract_dti_by_session.tsv     session tract segment voxels fa md rd ad
#             alongtract_fixel_by_session.tsv   session tract segment fixels fd log_fc fdc fd_zero_pct
#             (SESSION_FMLS=0.06 gives alongtract_fixel_by_session_fmls0p06.tsv; DTI is unaffected)
#
# NOT DONE HERE: the pre-/post-commissural split of the fornix (thesis plan 3.3.2).
# Segment sizes are not equal in voxels: the tract is thicker in some places than others.
# Needs python3 with nibabel, numpy and scipy (TractSeg installs them).
#
# Usage:
#   ./step8b_along_tract_v3.sh --build-segments          # once; prints the segment table
#   ./step8b_along_tract_v3.sh --list test8.txt          # a few sessions first
#   MAX_PARALLEL_8B=8 ./step8b_along_tract_v3.sh --list sessions.txt
#   SESSION_FMLS=0.06 MAX_PARALLEL_8B=8 ./step8b_along_tract_v3.sh --list sessions.txt
#   ./step8b_along_tract_v3.sh --summary
#   ./step8b_along_tract_v3.sh --rebuild-segments        # earlier results are set aside
#   ./step8b_along_tract_v3.sh --list sessions.txt --redo
# Settings (environment): N_SEGMENTS=20  ALONG_CUTOFF=3  MIN_PER_SEGMENT=10  EXTEND_MM=20
###############################################################################

set -u
SELF="$(readlink -f "$0")"
: "${NTHREADS:=4}"; export NTHREADS
source "$(dirname "$SELF")/pipeline_config_v2.sh"

TGROUP="${TEMPLATE_DIR}/group"; TMASK="${TGROUP}/template_mask.mif"
TRACTS_DIR="${TEMPLATE_DIR}/tracts"; SWARP_DIR="${TEMPLATE_DIR}/session_warps"; FIXEL_MASK="${TEMPLATE_DIR}/fixel_mask"
STEP8_MASKS="${TEMPLATE_DIR}/dti_tracts/masks"
AT="${TEMPLATE_DIR}/alongtract"; SEG="${AT}/segments"; LOGS="${AT}/logs"
N_SEG="${N_SEGMENTS:-20}"; CUT="${ALONG_CUTOFF:-3}"; MIN_PER_SEG="${MIN_PER_SEGMENT:-10}"
ALONG_TRACTS=(FX_left FX_right UF_left UF_right)
SESSION_FMLS="${SESSION_FMLS:-}"; export SESSION_FMLS
TAG=""; [ -n "$SESSION_FMLS" ] && TAG="_fmls${SESSION_FMLS//./p}"
ROWS_DTI="${AT}/rows_dti"; ROWS_FIX="${AT}/rows_fixel${TAG}"
RES_DTI="${AT}/alongtract_dti_by_session.tsv"; RES_FIX="${AT}/alongtract_fixel_by_session${TAG}.tsv"
SCRATCH_ROOT="${SCRATCH_ROOT:-$HOME/fba8b_work}"; MAX_PAR="${MAX_PARALLEL_8B:-8}"
export N_SEG CUT SCRATCH_ROOT
LOCK="/tmp/${USER}_step8b.lock"

grid3() { mrinfo "$1" -size 2>/dev/null | awk '{print $1,$2,$3}'; }
count() { mrstats "$1" -mask "$1" -output count 2>/dev/null | tr -d ' '; }
four()  { tr -s ' \t\n' ' ' | sed 's/^ //; s/ $//' | awk 'NF==4{printf "%s\t%s\t%s\t%s", $1,$2,$3,$4; ok=1} END{if(!ok) printf "NA\tNA\tNA\tNA"}'; }
pad()   { printf '%02d' "$1"; }

###############################################################################
# Python: centreline and segment labels for one tract
###############################################################################
write_python() {
cat > "$1" << 'PYEOF'
import sys
import numpy as np
import nibabel as nib
from scipy.spatial import cKDTree

tck_path, tmask_path, out_path, n_seg, cl_path, tract_path = sys.argv[1:7]
n_seg = int(n_seg)
ext_mm = float(sys.argv[7]) if len(sys.argv) > 7 else 20.0
M = 100

def arclen(s):
    return np.r_[0, np.cumsum(np.linalg.norm(np.diff(s, axis=0), axis=1))]

def resample(s, m):
    d = arclen(s)
    t = np.linspace(0, d[-1], m)
    return np.stack([np.interp(t, d, s[:, k]) for k in range(3)], axis=1)

def unit(v):
    return v / np.linalg.norm(v)

tf = nib.streamlines.load(tck_path)
streams = [np.asarray(x, dtype=float) for x in tf.streamlines if len(x) >= 2]
lens = np.array([arclen(x)[-1] for x in streams])
good = lens > 0
streams = [x for x, g in zip(streams, good) if g]; lens = lens[good]
n_total = len(streams)
# the centreline is estimated from streamlines that span most of the tract (partial ones would blur it)
full = lens >= 0.75 * np.median(lens)
R = np.array([resample(x, M) for x, f in zip(streams, full) if f])
if len(R) < 50:
    sys.exit("fewer than 50 usable streamlines in " + tck_path)

# reference: the streamline whose midpoint is closest to the median midpoint
mid = R[:, M // 2]
ref = R[np.argmin(np.linalg.norm(mid - np.median(mid, axis=0), axis=1))].copy()
parity = np.zeros(len(R), bool)
for _ in range(4):                       # align every streamline to run the same way as the reference
    d_same = np.linalg.norm(R - ref, axis=2).mean(axis=1)
    d_flip = np.linalg.norm(R[:, ::-1] - ref, axis=2).mean(axis=1)
    f = d_flip < d_same
    R[f] = R[f, ::-1]
    parity ^= f
    ref = R.mean(axis=0)
for _ in range(2):                       # then average only the 80% of streamlines closest to the centreline
    dist = np.linalg.norm(R - ref, axis=2).mean(axis=1)
    ref = R[dist <= np.percentile(dist, 80)].mean(axis=0)
if ref[0, 1] > ref[-1, 1]:               # world y points anterior: start at the posterior end
    ref = ref[::-1]; R = R[:, ::-1]
dist = np.linalg.norm(R - ref, axis=2).mean(axis=1)
length = arclen(ref)[-1]

# extend the centreline straight beyond both ends, so that parts of the tract mask lying beyond the
# ends (where tracts fan out) are spread over several segments instead of all landing in the end ones
k = (np.arange(1, int(ext_mm / 0.5) + 1) * 0.5)[:, None]
ext0 = ref[0] + k * unit(ref[0] - ref[5])
ext1 = ref[-1] + k * unit(ref[-1] - ref[-6])
poly = np.vstack([ext0[::-1], ref, ext1])
s_mm = arclen(poly)
s_ref0 = s_mm[len(ext0)]; s_ref1 = s_mm[len(ext0) + M - 1]
tree = cKDTree(poly)

img = nib.load(tmask_path)
mask = np.squeeze(np.asarray(img.dataobj)) > 0
idx = np.argwhere(mask)
xyz = nib.affines.apply_affine(img.affine, idx)
_, j = tree.query(xyz)

# The segments span the part of the axis that the tract mask covers (0.5th to 99.5th percentile of the
# positions of its voxels). The tract mask lies inside the 80-subject analysis mask, which can trim tract ends.
timg = nib.load(tract_path)
tm = np.squeeze(np.asarray(timg.dataobj)) > 0
txyz = nib.affines.apply_affine(timg.affine, np.argwhere(tm))
_, jt = tree.query(txyz)
lo, hi = np.percentile(s_mm[jt], [0.5, 99.5])
if hi - lo < 10:
    sys.exit("the tract mask covers less than 10 mm of the centreline axis in %s" % tract_path)
pos = np.clip((s_mm[j] - lo) / (hi - lo), 0, 1)
seg = np.minimum((pos * n_seg).astype(int), n_seg - 1) + 1
lab = np.zeros(mask.shape, np.uint8)
lab[tuple(idx.T)] = seg
nib.save(nib.Nifti1Image(lab, img.affine), out_path)

try:      # the axis the segments were cut along, for checking in mrview
    axis = poly[(s_mm >= lo) & (s_mm <= hi)].astype(np.float32)
    nib.streamlines.save(nib.streamlines.Tractogram([axis], affine_to_rasmm=np.eye(4)), cl_path)
except Exception as e:
    print("note: centreline not saved (%s)" % e, file=sys.stderr)

print("\t".join(str(v) for v in [
    len(R), "%.1f" % (100 * parity.mean()), "%.1f" % length,
    "%.1f" % ref[0, 0], "%.1f" % ref[0, 1], "%.1f" % ref[0, 2],
    "%.1f" % ref[-1, 0], "%.1f" % ref[-1, 1], "%.1f" % ref[-1, 2],
    "%.2f" % np.median(dist), "%.2f" % np.percentile(dist, 95), "%.1f" % (hi - lo),
    n_total, "%.1f" % max(0.0, s_ref0 - lo), "%.1f" % max(0.0, hi - s_ref1)]))
PYEOF
}

###############################################################################
# 1. Segments (once)
###############################################################################
segments_current() {
    [ -f "${SEG}/.done" ] || return 1
    [ "${TRACTS_DIR}/tracts.done" -nt "${SEG}/.done" ] && return 1
    [ "${STEP8_MASKS}/.done" -nt "${SEG}/.done" ] && return 1
    return 0
}
build_segments() {
    local t k kk tmp="${AT}/.segments_building" py="${AT}/.centreline.py" out n_v n_f
    command -v python3 >/dev/null && python3 -c 'import nibabel, numpy, scipy' 2>/dev/null \
        || die "python3 with nibabel, numpy and scipy is needed (try: pip install --user nibabel scipy)"
    for c in mrconvert mrcalc mrstats mrcat mrtransform voxel2fixel; do command -v "$c" >/dev/null || die "Command not found: ${c}"; done
    [ -f "${TRACTS_DIR}/tracts.done" ] || die "Step 6 tracts missing"
    [ -f "${STEP8_MASKS}/.done" ] || die "Step 8 masks missing: run step8_dti_tracts_v3.sh --masks-only first"
    [ -f "${FIXEL_MASK}/index.mif" ] || die "Fixel mask missing (step 5)"
    for t in "${ALONG_TRACTS[@]}"; do
        for f in "${TRACTS_DIR}/tractseg/TOM_trackings/${t}.tck" "${TRACTS_DIR}/fixels/${t}/${t}_ge${CUT}.mif" "${STEP8_MASKS}/${t}_ge${CUT}.mif"; do
            [ -f "$f" ] || die "Missing: ${f}"
        done
    done
    mkdir -p "$AT"; rm -rf "$tmp"; mkdir -p "$tmp"; write_python "$py"
    mrconvert "$TMASK" "${tmp}/tmask.nii.gz" -force -quiet || die "could not convert the template mask"
    printf 'tract\tstreamlines_used\tpct_flipped\tlength_mm\tstart_x\tstart_y\tstart_z\tend_x\tend_y\tend_z\tmedian_dist_mm\tp95_dist_mm\tcovered_mm\tstreamlines_total\tbeyond_start_mm\tbeyond_end_mm\n' > "${tmp}/centrelines.tsv"
    printf 'tract\tsegment\tvoxels\tfixels\n' > "${tmp}/segments.tsv"
    for t in "${ALONG_TRACTS[@]}"; do
        mkdir -p "${tmp}/${t}"
        log "Centreline and labels: ${t}"
        mrconvert "${STEP8_MASKS}/${t}_ge${CUT}.mif" "${tmp}/${t}/tractmask.nii.gz" -force -quiet || die "could not convert the ${t} voxel mask"
        out=$(python3 "$py" "${TRACTS_DIR}/tractseg/TOM_trackings/${t}.tck" "${tmp}/tmask.nii.gz" "${tmp}/${t}/labels.nii.gz" "$N_SEG" "${tmp}/${t}/centreline.tck" "${tmp}/${t}/tractmask.nii.gz" "${EXTEND_MM:-20}") \
            || die "centreline failed for ${t}"
        rm -f "${tmp}/${t}/tractmask.nii.gz"
        printf '%s\t%s\n' "$t" "$out" >> "${tmp}/centrelines.tsv"
        mrconvert "${tmp}/${t}/labels.nii.gz" "${tmp}/${t}/labels.mif" -strides "$TMASK" -datatype uint8 -force -quiet || die "label conversion failed"
        [ "$(grid3 "${tmp}/${t}/labels.mif")" = "$(grid3 "$TMASK")" ] || die "${t}: labels are not on the template grid"
        rm -f "${tmp}/${t}/labels.nii.gz"
        voxel2fixel "${tmp}/${t}/labels.mif" "$FIXEL_MASK" "${tmp}/${t}/fixlab" labels.mif -force -quiet || die "voxel2fixel failed for ${t}"
        for k in $(seq 1 "$N_SEG"); do
            kk=$(pad "$k")
            mrcalc "${tmp}/${t}/labels.mif" "$k" -eq "${STEP8_MASKS}/${t}_ge${CUT}.mif" -mult "${tmp}/${t}/seg${kk}_vox.mif" -datatype bit -force -quiet || die "voxel segment failed"
            mrcalc "${tmp}/${t}/fixlab/labels.mif" "$k" -eq "${TRACTS_DIR}/fixels/${t}/${t}_ge${CUT}.mif" -mult "${tmp}/${t}/seg${kk}_fix.mif" -datatype bit -force -quiet || die "fixel segment failed"
            n_v=$(count "${tmp}/${t}/seg${kk}_vox.mif"); n_f=$(count "${tmp}/${t}/seg${kk}_fix.mif")
            printf '%s\t%s\t%s\t%s\n' "$t" "$k" "${n_v:-0}" "${n_f:-0}" >> "${tmp}/segments.tsv"
            [ "${n_v:-0}" -gt 0 ] 2>/dev/null && [ "${n_f:-0}" -gt 0 ] 2>/dev/null \
                || die "${t} segment ${k} is empty (voxels ${n_v:-0}, fixels ${n_f:-0}). Use fewer segments: N_SEGMENTS=10"
        done
        rm -rf "${tmp}/${t}/fixlab"
    done
    local bad
    bad=$(awk -F'\t' -v m="$MIN_PER_SEG" 'NR>1 && ($3<m || $4<m) {printf "%s seg%s (%s voxels, %s fixels)  ", $1, $2, $3, $4}' "${tmp}/segments.tsv")
    [ -n "$bad" ] && log "WARNING: segments with fewer than ${MIN_PER_SEG} voxels or fixels: ${bad}"
    # left and right should start at about the same place (same anatomical end)
    awk -F'\t' 'NR>1 {y[$1]=$6; ln[$1]=$4} END{
        if (y["FX_left"]!="" && y["FX_right"]!="" && (y["FX_left"]-y["FX_right"] > 15 || y["FX_right"]-y["FX_left"] > 15)) print "WARNING: left and right fornix start more than 15 mm apart in y"
        if (y["UF_left"]!="" && y["UF_right"]!="" && (y["UF_left"]-y["UF_right"] > 15 || y["UF_right"]-y["UF_left"] > 15)) print "WARNING: left and right uncinate start more than 15 mm apart in y" }' \
        "${tmp}/centrelines.tsv" | while read -r l; do log "$l"; done
    rm -rf "$SEG"; mv "$tmp" "$SEG"; touch "${SEG}/.done"
    if ls "${ROWS_DTI}"/HCA*.tsv "${AT}"/rows_fixel*/HCA*.tsv >/dev/null 2>&1; then
        local old="${AT}/old_$(date '+%Y%m%d_%H%M%S')"; mkdir -p "$old"
        mv "${AT}"/rows_dti "${AT}"/rows_fixel* "${AT}"/alongtract_*.tsv "$old/" 2>/dev/null
        log "Segments were rebuilt: earlier session results moved to ${old}. Run --list again."
    fi
    log "Segments built: ${SEG}"
    show_segments
}
show_segments() {
    echo; echo "Centrelines (start = posterior end; mm). streamlines_used = full-length ones the centreline is built from; median/p95_dist = how far they lie from it;"
    echo "covered_mm = length the segments span; beyond_start/end_mm = how far the tract mask reaches past the centreline ends:"
    column -t -s$'\t' "${SEG}/centrelines.tsv"
    echo; echo "Voxels and fixels per segment (cutoff ge${CUT}):"
    local t
    for t in "${ALONG_TRACTS[@]}"; do
        awk -F'\t' -v t="$t" '$1==t {v=v" "$3; f=f" "$4} END{printf "%-9s voxels:%s\n%-9s fixels:%s\n", t, v, t, f}' "${SEG}/segments.tsv"
    done
    echo; echo "To look at them: mrview template_l0.mif -overlay.load ${SEG}/<tract>/labels.mif -overlay.colourmap 3 -tractography.load ${SEG}/<tract>/centreline.tck"
}

###############################################################################
# 2. One session
###############################################################################
run_session() {
    local s=$1 W="${SCRATCH_ROOT}/$1" t k kk vals nv nf m
    mkdir -p "$LOGS"; exec >> "${LOGS}/${s}.log" 2>&1
    log "start ${s} (tag '${TAG}')"
    rm -rf "$W"; mkdir -p "$W"
    local -A NV NF
    while IFS=$'\t' read -r t k nv nf; do NV[$t.$k]=$nv; NF[$t.$k]=$nf; done < <(tail -n +2 "${SEG}/segments.tsv")

    if [ -z "$TAG" ] && [ ! -f "${ROWS_DTI}/${s}.tsv" ]; then
        for m in fa md rd ad; do
            mrtransform "${OUT_ROOT}/${s}/dti/${m}.mif" -warp "${SWARP_DIR}/${s}.mif" -interp linear "${W}/${m}.mif" \
                -nthreads "$NTHREADS" -force -quiet || { log "FAIL warp ${m}"; rm -rf "$W"; return 1; }
        done
        [ "$(grid3 "${W}/fa.mif")" = "$(grid3 "$TMASK")" ] || { log "FAIL: warped FA not on the template grid"; rm -rf "$W"; return 1; }
        mrcat "${W}/fa.mif" "${W}/md.mif" "${W}/rd.mif" "${W}/ad.mif" -axis 3 "${W}/dti4.mif" -force -quiet || { log "FAIL mrcat dti"; rm -rf "$W"; return 1; }
        : > "${W}/dti_rows.tsv"
        for t in "${ALONG_TRACTS[@]}"; do for k in $(seq 1 "$N_SEG"); do
            kk=$(pad "$k")
            vals=$(mrstats "${W}/dti4.mif" -mask "${SEG}/${t}/seg${kk}_vox.mif" -output mean 2>/dev/null | four)
            printf '%s\t%s\t%s\t%s\t%s\n' "$s" "$t" "$k" "${NV[$t.$k]}" "$vals" >> "${W}/dti_rows.tsv"
        done; done
        grep -q $'\tNA\t' "${W}/dti_rows.tsv" && log "WARNING: empty DTI statistics (NA)"
        mkdir -p "$ROWS_DTI"; cp "${W}/dti_rows.tsv" "${ROWS_DTI}/.${s}.tmp" && mv "${ROWS_DTI}/.${s}.tmp" "${ROWS_DTI}/${s}.tsv"
    fi

    if [ ! -f "${ROWS_FIX}/${s}.tsv" ]; then
        local FD="${TEMPLATE_DIR}/fd${TAG}/${s}.mif" LF="${TEMPLATE_DIR}/log_fc${TAG}/${s}.mif" FC="${TEMPLATE_DIR}/fdc${TAG}/${s}.mif"
        mrcalc "$FD" 0 -eq "${W}/z.mif" -datatype float32 -force -quiet || { log "FAIL zero indicator"; rm -rf "$W"; return 1; }
        mrcat "$FD" "$LF" "$FC" "${W}/z.mif" -axis 3 "${W}/fix4.mif" -force -quiet || { log "FAIL mrcat fixel"; rm -rf "$W"; return 1; }
        : > "${W}/fix_rows.tsv"
        for t in "${ALONG_TRACTS[@]}"; do for k in $(seq 1 "$N_SEG"); do
            kk=$(pad "$k")
            vals=$(mrstats "${W}/fix4.mif" -mask "${SEG}/${t}/seg${kk}_fix.mif" -output mean 2>/dev/null | four \
                   | awk -F'\t' 'BEGIN{OFS="\t"} {if ($4!="NA") $4=sprintf("%.2f", 100*$4); print}')
            printf '%s\t%s\t%s\t%s\t%s\n' "$s" "$t" "$k" "${NF[$t.$k]}" "$vals" >> "${W}/fix_rows.tsv"
        done; done
        grep -q $'\tNA\t' "${W}/fix_rows.tsv" && log "WARNING: empty fixel statistics (NA)"
        mkdir -p "$ROWS_FIX"; cp "${W}/fix_rows.tsv" "${ROWS_FIX}/.${s}.tmp" && mv "${ROWS_FIX}/.${s}.tmp" "${ROWS_FIX}/${s}.tsv"
    fi
    rm -rf "$W"; log "OK ${s}"
}
if [ "${1:-}" = "--worker" ]; then run_session "$2"; exit $?; fi

###############################################################################
# Combine and summarise
###############################################################################
combine() {
    local f
    if [ -z "$TAG" ] && ls "${ROWS_DTI}"/HCA*.tsv >/dev/null 2>&1; then
        { printf 'session\ttract\tsegment\tvoxels\tfa\tmd\trd\tad\n'; cat "${ROWS_DTI}"/HCA*.tsv; } > "${RES_DTI}.tmp" && mv "${RES_DTI}.tmp" "$RES_DTI"
    fi
    if ls "${ROWS_FIX}"/HCA*.tsv >/dev/null 2>&1; then
        { printf 'session\ttract\tsegment\tfixels\tfd\tlog_fc\tfdc\tfd_zero_pct\n'; cat "${ROWS_FIX}"/HCA*.tsv; } > "${RES_FIX}.tmp" && mv "${RES_FIX}.tmp" "$RES_FIX"
    fi
}
profile() {  # file column scale label decimals : mean across sessions, per segment, one line per tract
    local f=$1 col=$2 scale=$3 lab=$4 dec=${5:-3} t
    echo "  ${lab}, segments 1 to ${N_SEG}, posterior to anterior:"
    for t in "${ALONG_TRACTS[@]}"; do
        printf '    %-9s' "$t"
        awk -F'\t' -v t="$t" -v c="$col" -v sc="$scale" -v N="$N_SEG" -v d="$dec" 'NR>1 && $2==t && $c!="NA" {n[$3]++; s[$3]+=$c}
            END{for(k=1;k<=N;k++) printf " %." d "f", (n[k] ? sc*s[k]/n[k] : 0); printf "\n"}' "$f"
    done
}
summary() {
    [ -f "${SEG}/.done" ] && show_segments
    if [ -f "$RES_DTI" ]; then
        echo; echo "DTI: $(tail -n +2 "$RES_DTI" | cut -f1 | sort -u | wc -l) sessions, $(tail -n +2 "$RES_DTI" | wc -l) rows"
        local na; na=$(tail -n +2 "$RES_DTI" | awk -F'\t' '$5=="NA"' | wc -l); [ "$na" -gt 0 ] && echo "  WARNING: ${na} rows with NA"
        profile "$RES_DTI" 5 1 "FA"; profile "$RES_DTI" 6 1000 "MD x 1000 (10^-3 mm2/s)"
    fi
    if [ -f "$RES_FIX" ]; then
        echo; echo "FBA${TAG}: $(tail -n +2 "$RES_FIX" | cut -f1 | sort -u | wc -l) sessions, $(tail -n +2 "$RES_FIX" | wc -l) rows"
        local na; na=$(tail -n +2 "$RES_FIX" | awk -F'\t' '$5=="NA"' | wc -l); [ "$na" -gt 0 ] && echo "  WARNING: ${na} rows with NA"
        profile "$RES_FIX" 5 1 "FD"; profile "$RES_FIX" 8 1 "FD = 0 (% of fixels)" 1
    fi
    [ -f "$RES_DTI" ] || [ -f "$RES_FIX" ] || echo "No session results yet."
}

###############################################################################
# Controller
###############################################################################
LIST=""; MODE=run; REDO=0
while [ $# -gt 0 ]; do
    case "$1" in
        --list) LIST="$(readlink -f "$2")"; shift ;;
        --build-segments) MODE=build ;; --rebuild-segments) MODE=rebuild ;; --summary) MODE=summary ;; --redo) REDO=1 ;;
        *) die "Unknown option: $1" ;;
    esac; shift
done
[ "$MODE" = "summary" ] && { summary; exit 0; }
exec 9>"$LOCK"; flock -n 9 || die "Step 8b is already running."
mkdir -p "$AT" "$LOGS" "$SCRATCH_ROOT"
case "$MODE" in
    build)   if segments_current; then log "Segments already built and current."; show_segments; exit 0; fi
             [ -f "${SEG}/.done" ] && die "Segments exist but the step 6 tracts or step 8 masks are newer. Use --rebuild-segments."
             build_segments; exit 0 ;;
    rebuild) build_segments; exit 0 ;;
esac
[ -f "$LIST" ] || die "Usage: $0 --list sessions.txt   (or --build-segments, --summary)"
[ -f "${SEG}/.done" ] || die "No segments yet. Run --build-segments first."
segments_current || die "The step 6 tracts or step 8 masks changed after the segments were built. Use --rebuild-segments."

mapfile -t SESSIONS < <(grep -o 'HCA[0-9]*_V[0-9]*_MR' "$LIST" | awk '!seen[$0]++')
problems=(); QUEUE=()
for s in "${SESSIONS[@]}"; do
    [ -f "${SWARP_DIR}/${s}.mif" ] || problems+=("${s}: no session warp (step 5)")
    for d in fd log_fc fdc; do [ -f "${TEMPLATE_DIR}/${d}${TAG}/${s}.mif" ] || problems+=("${s}: no ${d}${TAG} file (step 5)"); done
    if [ -z "$TAG" ]; then for m in fa md rd ad; do [ -f "${OUT_ROOT}/${s}/dti/${m}.mif" ] || problems+=("${s}: no ${m}.mif (step 3c)"); done; fi
done
if [ "${#problems[@]}" -gt 0 ]; then printf '  %s\n' "${problems[@]}"; die "${#problems[@]} problem(s) above."; fi
for s in "${SESSIONS[@]}"; do      # only after every check has passed
    if [ "$REDO" = "1" ]; then rm -f "${ROWS_FIX}/${s}.tsv"; [ -z "$TAG" ] && rm -f "${ROWS_DTI}/${s}.tsv"; fi
    if [ -f "${ROWS_FIX}/${s}.tsv" ] && { [ -n "$TAG" ] || [ -f "${ROWS_DTI}/${s}.tsv" ]; }; then continue; fi
    QUEUE+=("$s")
done
log "Step 8b: ${#SESSIONS[@]} sessions, $(( ${#SESSIONS[@]} - ${#QUEUE[@]} )) already done, ${#QUEUE[@]} to do, ${MAX_PAR} at once, ${N_SEG} segments, $([ -n "$TAG" ] && echo "fixel part only, cutoff ${SESSION_FMLS}" || echo "DTI and FBA")"
if [ "${#QUEUE[@]}" -gt 0 ]; then
    printf '%s\n' "${QUEUE[@]}" | xargs -P "$MAX_PAR" -I{} bash "$SELF" --worker {}
fi
fail=0
for s in "${SESSIONS[@]}"; do
    ok=1; [ -f "${ROWS_FIX}/${s}.tsv" ] || ok=0; [ -z "$TAG" ] && { [ -f "${ROWS_DTI}/${s}.tsv" ] || ok=0; }
    [ "$ok" = 1 ] || { fail=$((fail+1)); log "not finished: ${s} (see ${LOGS}/${s}.log)"; }
done
combine
log "Finished: $(( ${#SESSIONS[@]} - fail )) of ${#SESSIONS[@]} sessions."
[ "$fail" -gt 0 ] && log "Rerun the same command to retry the unfinished ones."
echo; summary
exit 0
