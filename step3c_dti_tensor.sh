#!/bin/bash
###############################################################################
# STEP 3c - DTI tensor fitting (DTI ARM)
#
# Fully independent. Needs only step 1 output: no group response, no FODs, no
# recon-all. Can run at any point after step 1, including right now, in parallel
# with everything else.
#
#   dwiextract -shell 0,1000   low b-value shell only
#   dwi2tensor                 one diffusion ellipsoid per voxel
#   tensor2metric              -> FA, MD, RD, AD
#
# SHELL EXTRACTION IS NOT OPTIONAL
# The tensor model assumes Gaussian diffusion, which only holds at low b. Fitting
# the full multi-shell data biases FA and MD. This is a silent error: dwi2tensor
# will happily accept multi-shell input and return wrong numbers.
#
# NATIVE RESOLUTION
# DTI runs on dwi_unbiased.mif, not the 0.8 mm upsampled version. Upsampling adds
# no real information for tensor fitting, and DTI's spatial resolution needs are
# lower than FBA's.
#
# Usage: ./step3c_dti_tensor.sh --list sessions.txt
###############################################################################

set -u
source "$(dirname "$0")/pipeline_config.sh"
parse_sessions "$@"

log "DTI shells: ${DTI_SHELLS}"
log "Running tensor fitting, up to ${MAX_PARALLEL_LIGHT} concurrent"

for sess in "${SESSIONS[@]}"; do
    WORK="${OUT_ROOT}/${sess}"

    [ -f "${WORK}/dwi_unbiased.mif" ] || { log "SKIP ${sess}: no dwi_unbiased"; continue; }
    [ -f "${WORK}/dti/fa.mif" ]       && { log "SKIP ${sess}: already done"; continue; }

    (
        cd "$WORK" || exit 1
        set -e
        mkdir -p dti

        # Confirm the requested shells actually exist in this data
        shells=$(mrinfo dwi_unbiased.mif -shell_bvalues)
        echo "  ${sess} available shells: ${shells}"

        # 1) Extract the low b-value shell
        dwiextract dwi_unbiased.mif -shells $DTI_SHELLS \
            dti/dti_shell.mif -force -nthreads $NTHREADS

        # 2) Fit the tensor
        dwi2tensor dti/dti_shell.mif dti/tensor.mif \
            -mask mask.mif -force -nthreads $NTHREADS

        # 3) Derive the scalar metrics from the eigenvalues
        #    FA  ellipsoid elongation
        #    MD  mean of all three eigenvalues
        #    AD  lambda1, diffusion along the primary axis
        #    RD  mean of lambda2 and lambda3, diffusion across it
        tensor2metric dti/tensor.mif \
            -fa  dti/fa.mif \
            -adc dti/md.mif \
            -rd  dti/rd.mif \
            -ad  dti/ad.mif \
            -mask mask.mif -force -nthreads $NTHREADS

        echo "[$(date '+%H:%M:%S')] done ${sess}"
    ) &

    throttle $MAX_PARALLEL_LIGHT
done
wait

log "Step 3c complete. FA/MD/RD/AD in ${OUT_ROOT}/<sess>/dti/ (native space)."
log ""
log "Sanity check FA range. Values should sit roughly 0 to 1, with WM around"
log "0.4 to 0.8. Values above 1 or widespread near-zero WM indicate a bad fit:"
log "  mrstats ${OUT_ROOT}/<sess>/dti/fa.mif -mask ${OUT_ROOT}/<sess>/mask.mif"
log ""
log "These maps stay in native space until step 8, which warps them into"
log "template space using the FBA warps and extracts per-tract values."
