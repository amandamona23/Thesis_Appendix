#!/bin/bash
###############################################################################
# STEP 3a v2 - FOD estimation and intensity normalisation
#
# READS  from ${IN_ROOT}/<sess>/   (symlinks to RDS, from build_session_index.sh)
#          dwi_upsampled.mif, mask_upsampled.mif
# WRITES to ${OUT_ROOT}/<sess>/    (local disk)
#          wmfod_norm.mif, gmfod_norm.mif, csffod_norm.mif,
#          mask_upsampled.mif (local copy), vf_qc.mif
#
# Reads over sshfs are far safer than writes: a dropped mount makes the read
# fail loudly, the session aborts, and a rerun picks it up. Writes stay local.
#
# DISK: wmfod.mif (un-normalised) is the same size as wmfod_norm.mif and is
# never used again after mtnormalise, so it is deleted straight away. Without
# that, step 3a would need roughly double the space.
#
# CONCURRENCY: each session reads ~11 GB over the network. Default is 6 at
# once. Test with 2 or 3 first and watch the mount survive.
#
# Usage: ./step3a_fod_normalise_v2.sh --list sessions.txt
###############################################################################

set -u
source "$(dirname "$0")/pipeline_config_v2.sh"
parse_sessions "$@"

MAX_PARALLEL_3A=${MAX_PARALLEL_3A:-6}

GROUP_WM="${GROUP_DIR}/group_wm.txt"
GROUP_GM="${GROUP_DIR}/group_gm.txt"
GROUP_CSF="${GROUP_DIR}/group_csf.txt"
for f in "$GROUP_WM" "$GROUP_GM" "$GROUP_CSF"; do
    [ -s "$f" ] || die "Missing $f. Run step 2 first."
done

mkdir -p "$OUT_ROOT"
log "Reading from ${IN_ROOT}, writing to ${OUT_ROOT}"
log "Up to ${MAX_PARALLEL_3A} sessions at once"
df -h "$OUT_ROOT" | tail -1 | awk '{print "[local disk] "$4" free of "$2}'

for sess in "${SESSIONS[@]}"; do
    SRC="${IN_ROOT}/${sess}"
    DST="${OUT_ROOT}/${sess}"

    [ -f "${DST}/wmfod_norm.mif" ] && { log "SKIP ${sess}: already done"; continue; }
    [ -f "${SRC}/dwi_upsampled.mif" ] || { log "SKIP ${sess}: no dwi_upsampled on RDS"; continue; }
    [ -f "${SRC}/mask_upsampled.mif" ] || { log "SKIP ${sess}: no mask_upsampled on RDS"; continue; }

    (
        set -e
        mkdir -p "$DST"
        cd "$DST"

        # Resolution guard: every session must be 0.8 mm
        sp=$(mrinfo "${SRC}/dwi_upsampled.mif" -spacing | awk '{printf "%.2f", $1}')
        if [ "$sp" != "$(printf '%.2f' "$UPSAMPLE_VOXEL")" ]; then
            echo "ABORT ${sess}: ${sp} mm, expected ${UPSAMPLE_VOXEL} mm"
            exit 1
        fi

        # Local copy of the mask. Small, and every later FBA step needs it.
        mrconvert "${SRC}/mask_upsampled.mif" mask_upsampled.mif -force -quiet

        # 1) Multi-shell multi-tissue CSD against the GROUP response
        dwi2fod msmt_csd "${SRC}/dwi_upsampled.mif" -mask mask_upsampled.mif \
            "$GROUP_WM"  wmfod.mif \
            "$GROUP_GM"  gmfod.mif \
            "$GROUP_CSF" csffod.mif \
            -force -nthreads $NTHREADS

        # 2) Joint intensity normalisation
        mtnormalise wmfod.mif  wmfod_norm.mif \
                    gmfod.mif  gmfod_norm.mif \
                    csffod.mif csffod_norm.mif \
                    -mask mask_upsampled.mif -force

        # Un-normalised versions are never used again
        rm -f wmfod.mif gmfod.mif csffod.mif

        # 3) Small QC image: CSF / GM / WM as three volumes
        mrconvert wmfod_norm.mif -coord 3 0 - -quiet | \
            mrcat csffod_norm.mif gmfod_norm.mif - vf_qc.mif -force -quiet

        echo "[$(date '+%H:%M:%S')] done ${sess}"
    ) || log "FAILED ${sess} (rerun the script to retry)" &

    throttle $MAX_PARALLEL_3A
done
wait

n_done=$(ls "${OUT_ROOT}"/*/wmfod_norm.mif 2>/dev/null | wc -l)
log "Step 3a finished. ${n_done} sessions have wmfod_norm.mif"
log "Rerun the same command to retry any that failed; finished sessions are skipped."
