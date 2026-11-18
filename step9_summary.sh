#!/bin/bash
###############################################################################
# step9_summary.sh  -  summarise the fixelcfestats results and make figure files
#
# For every finished run in ${TEMPLATE_DIR}/stats/ and every hypothesis in it:
#   fixels_sig   fixels with FWE-corrected p < 0.05 (fwe_1mpvalue > 0.95) inside the tract mask
#   pct_sig      as a % of the tract mask
#   min_p_fwe    smallest FWE-corrected p in the mask (= 1 - max of fwe_1mpvalue)
#   abs_effect   mean absolute effect in the significant fixels (metric units per decade for
#                age contrasts; per year for mean change); std_effect the same, standardised
# Writes stats/step9_summary.tsv and prints it.
#
# Figures (--figures): for each run and hypothesis with significant fixels, maps the FWE
# result onto that tract's TractSeg streamlines (fixel2tsf) so it can be shown as coloured
# streamlines in mrview: stats/figures/<run>_<hypothesis>_<side>.tsf
#
# Usage:   ./step9_summary.sh            ./step9_summary.sh --figures
###############################################################################
set -u
source "$(dirname "$(readlink -f "$0")")/pipeline_config_v2.sh"
STATS="${TEMPLATE_DIR}/stats"; FIX="${TEMPLATE_DIR}/tracts/fixels"; TCK="${TEMPLATE_DIR}/tracts/tractseg/TOM_trackings"
OUT="${STATS}/step9_summary.tsv"; FIG=0; [ "${1:-}" = "--figures" ] && FIG=1
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
printf 'run\tmodel\tmetric\ttract\thypothesis\tmask_fixels\tfixels_sig\tpct_sig\tmin_p_fwe\tabs_effect\tstd_effect\n' > "$OUT"

for d in "${STATS}"/*/; do
    d=${d%/}; run=$(basename "$d"); [ -f "$d/.done" ] || continue
    model=${run%%_*}; tract=$(echo "$run" | grep -o 'FX_both\|UF_both\|FX_left\|FX_right\|UF_left\|UF_right')
    cut=$(echo "$run" | grep -o '_ge[0-9]*' | head -1); metric=${run#${model}_}; metric=${metric%%_${tract}*}
    [ "$model" = "v1" ] && [[ "$run" == v1_vent_* ]] && { model="v1_vent"; metric=${run#v1_vent_}; metric=${metric%%_${tract}*}; }
    M="${FIX}/${tract}/${tract}${cut}.mif"; [ -f "$M" ] || { echo "no mask for ${run}"; continue; }
    nmask=$(mrstats "$M" -mask "$M" -output count | tr -d ' ')
    read -r -a HYP < "$d/hypotheses.txt"
    # one FWE file per hypothesis: fwe_1mpvalue_t0.mif, _t1.mif ... (or fwe_1mpvalue.mif for a single one)
    mapfile -t FWE < <(ls "$d"/fwe_1mpvalue*.mif 2>/dev/null | sort -V)
    [ "${#FWE[@]}" -gt 0 ] || { echo "no fwe_1mpvalue files in ${run}"; continue; }
    for i in "${!FWE[@]}"; do
        f=${FWE[$i]}; suf=$(basename "$f" .mif); suf=${suf#fwe_1mpvalue}
        h=${HYP[$i]:-h$i}
        mrcalc "$f" 0.95 -gt "$M" -mult "$W/sig.mif" -datatype bit -force -quiet
        nsig=$(mrstats "$W/sig.mif" -mask "$W/sig.mif" -output count 2>/dev/null | tr -d ' '); nsig=${nsig:-0}
        pmax=$(mrstats "$f" -mask "$M" -output max | tr -d ' ')
        ae="NA"; se="NA"
        if [ "${nsig}" -gt 0 ] 2>/dev/null; then
            [ -f "$d/abs_effect${suf}.mif" ] && ae=$(mrstats "$d/abs_effect${suf}.mif" -mask "$W/sig.mif" -output mean | tr -d ' ')
            [ -f "$d/std_effect${suf}.mif" ] && se=$(mrstats "$d/std_effect${suf}.mif" -mask "$W/sig.mif" -output mean | tr -d ' ')
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$run" "$model" "$metric" "$tract" "$h" "$nmask" "$nsig" \
            "$(awk -v a="$nsig" -v b="$nmask" 'BEGIN{printf "%.1f", 100*a/b}')" \
            "$(awk -v p="$pmax" 'BEGIN{printf "%.4f", 1-p}')" "$ae" "$se" >> "$OUT"
        if [ "$FIG" = 1 ] && [ "${nsig}" -gt 0 ] 2>/dev/null; then
            mkdir -p "${STATS}/figures"; base=${tract%_both}
            for side in left right; do
                t="${TCK}/${base}_${side}.tck"; [ -f "$t" ] || continue
                fixel2tsf "$f" "$t" "${STATS}/figures/${run}_${h}_${side}.tsf" -force -quiet \
                    && echo "figure: ${STATS}/figures/${run}_${h}_${side}.tsf"
            done
        fi
    done
done
echo; column -t -s$'\t' "$OUT"
echo; echo "Saved: ${OUT}"
