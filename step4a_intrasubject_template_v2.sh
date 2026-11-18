#!/bin/bash
###############################################################################
# STEP 4a v2 - Intra-subject halfway templates (FBA ARM, level 1 of 2)
#
# Per subject: register V1 and V2 symmetrically to their own midpoint.
# Embarrassingly parallel across subjects.
#
# READS  ${OUT_ROOT}/<SESS>/wmfod_norm.mif, mask_upsampled.mif   (NAS, from 3a)
# WRITES ${TEMPLATE_DIR}/subject_templates/<SUBJ>/                (NAS)
#          halfway_template.mif, halfway_mask.mif,
#          warps/<SESS>.mif  (visit -> halfway, needed in step 5)
#
# CHANGES FROM v1
#   - sources pipeline_config_v2.sh (v1 loaded the old config and old paths)
#   - each subject gets its own -scratch folder; population_template otherwise
#     writes multi-GB temp dirs into whatever directory you launched from
#   - failed subjects are logged; rerunning retries only unfinished subjects
#   - MAX_PARALLEL_TEMPLATE can be set on the command line
#
# Usage: MAX_PARALLEL_TEMPLATE=4 ./step4a_intrasubject_template_v2.sh --list sessions.txt
###############################################################################

set -u
source "$(dirname "$0")/pipeline_config_v2.sh"
parse_sessions "$@"

MAX_PARALLEL_TEMPLATE=${MAX_PARALLEL_TEMPLATE:-4}
SUBJ_TEMPLATE_DIR="${TEMPLATE_DIR}/subject_templates"
mkdir -p "$SUBJ_TEMPLATE_DIR"

log "Reading FODs from ${OUT_ROOT}"
log "Writing templates to ${SUBJ_TEMPLATE_DIR}"
df -h "$SUBJ_TEMPLATE_DIR" | tail -1 | awk '{print "[output disk] "$4" free of "$2}'

# ---------------------------------------------------------------------------
# Group sessions by subject
# ---------------------------------------------------------------------------
declare -A subj_sessions
for sess in "${SESSIONS[@]}"; do
    subj="${sess%%_*}"
    subj_sessions[$subj]="${subj_sessions[$subj]:-} ${sess}"
done
log "${#subj_sessions[@]} subjects, up to ${MAX_PARALLEL_TEMPLATE} at once, ${NTHREADS} threads each"

for subj in "${!subj_sessions[@]}"; do
    read -r -a visits <<< "${subj_sessions[$subj]}"
    SUBJ_DIR="${SUBJ_TEMPLATE_DIR}/${subj}"

    if [ -f "${SUBJ_DIR}/halfway_template.mif" ] && [ -f "${SUBJ_DIR}/halfway_mask.mif" ]; then
        log "SKIP ${subj}: already done"
        continue
    fi

    ok=1
    for v in "${visits[@]}"; do
        [ -f "${OUT_ROOT}/${v}/wmfod_norm.mif" ]     || { log "  ${subj}: ${v} missing wmfod_norm"; ok=0; }
        [ -f "${OUT_ROOT}/${v}/mask_upsampled.mif" ] || { log "  ${subj}: ${v} missing mask_upsampled"; ok=0; }
    done
    [ "$ok" -eq 1 ] || { log "SKIP ${subj}: incomplete inputs"; continue; }

    if [ "${#visits[@]}" -ne 2 ]; then
        log "SKIP ${subj}: expected 2 visits, found ${#visits[@]} (${visits[*]})"
        continue
    fi

    (
        set -e
        # Clear any half-finished attempt so population_template starts clean
        rm -rf "${SUBJ_DIR}/warps" "${SUBJ_DIR}/scratch" "${SUBJ_DIR}/halfway_template.mif"
        mkdir -p "${SUBJ_DIR}/fod_input" "${SUBJ_DIR}/mask_input" "${SUBJ_DIR}/scratch"

        for v in "${visits[@]}"; do
            ln -sfn "${OUT_ROOT}/${v}/wmfod_norm.mif"     "${SUBJ_DIR}/fod_input/${v}.mif"
            ln -sfn "${OUT_ROOT}/${v}/mask_upsampled.mif" "${SUBJ_DIR}/mask_input/${v}.mif"
        done
        printf '%s\n' "${visits[@]}" > "${SUBJ_DIR}/visits.txt"

        # Symmetric registration of this subject's two visits to their midpoint
        population_template "${SUBJ_DIR}/fod_input" "${SUBJ_DIR}/halfway_template.mif" \
            -mask_dir "${SUBJ_DIR}/mask_input" \
            -warp_dir "${SUBJ_DIR}/warps" \
            -voxel_size $UPSAMPLE_VOXEL \
            -scratch "${SUBJ_DIR}/scratch" \
            -nthreads $NTHREADS -force

        # Intersection mask: voxels present in both visits
        mrmath "${SUBJ_DIR}/mask_input/${visits[0]}.mif" \
               "${SUBJ_DIR}/mask_input/${visits[1]}.mif" \
               min "${SUBJ_DIR}/halfway_mask.mif" -datatype bit -force

        rm -rf "${SUBJ_DIR}/scratch"
        echo "[$(date '+%H:%M:%S')] done ${subj}"
    ) || log "FAILED ${subj} (rerun to retry)" &

    throttle $MAX_PARALLEL_TEMPLATE
done
wait

n_done=$(ls "${SUBJ_TEMPLATE_DIR}"/*/halfway_template.mif 2>/dev/null | wc -l)
log "Step 4a finished. ${n_done} of ${#subj_sessions[@]} subjects have a halfway template."
log "Rerun the same command to retry any that failed; finished subjects are skipped."
