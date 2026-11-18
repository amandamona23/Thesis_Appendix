#!/bin/bash
###############################################################################
# step9_make_change.sh  -  annual rate of change images for the change analysis
#
# For each subject: rate = (last visit - V1) / years between them, for the SMOOTHED
# fd, log_fc and fdc (step 7). Smoothing is linear, so this equals smoothing the change.
#   READS   covariates_by_subject.tsv (build_covariates.py), <metric>_smooth/ (step 7)
#   WRITES  ${TEMPLATE_DIR}/change/<metric>_rate/<subject>.mif  (+ index.mif, directions.mif)
# Usage:
#   ./step9_make_change.sh                     # fd, log_fc, fdc (standard)
#   TAG=_fmls0p06 ./step9_make_change.sh       # the 0.06 versions
###############################################################################
set -u
source "$(dirname "$(readlink -f "$0")")/pipeline_config_v2.sh"
COV="${COV:-$HOME/data_dell/fba/covariates/covariates_by_subject.tsv}"
TAG="${TAG:-}"
[ -f "$COV" ] || die "No subject table: $COV (run build_covariates.py first)"
for m in fd log_fc fdc; do
    src="${TEMPLATE_DIR}/${m}${TAG}_smooth"; dst="${TEMPLATE_DIR}/change/${m}${TAG}_rate"
    [ -d "$src" ] || die "Missing ${src}"
    mkdir -p "$dst"; cp "${src}/index.mif" "${src}/directions.mif" "$dst/"
    n=0
    while IFS=$'\t' read -r sid years first last; do
        [ -f "${dst}/${sid}.mif" ] && { n=$((n+1)); continue; }
        mrcalc "${src}/${last}.mif" "${src}/${first}.mif" -sub "$years" -div "${dst}/.${sid}.mif" -force -quiet \
            && mv "${dst}/.${sid}.mif" "${dst}/${sid}.mif" || die "mrcalc failed for ${sid} ${m}"
        n=$((n+1))
    done < <(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} {print $c["subject"]"\t"$c["years_v1_to_last"]"\t"$c["first_session"]"\t"$c["last_session"]}' "$COV")
    log "${m}${TAG}: ${n} rate images in ${dst}"
done
