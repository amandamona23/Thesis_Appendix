#!/bin/bash
###############################################################################
# step9_fixelcfestats.sh  -  fixel-wise statistics (CFE + permutation testing)
#
# One fixelcfestats run per model x metric x tract, all contrasts of that model tested
# together in the run (they share the permutations):
#   v1      V1 sessions, smoothed data (<metric>_smooth/)
#           contrasts: age_neg (lower with age), age_pos (higher with age)
#   change  annual rate images (change/<metric>_rate/), -errors ise (sign-flipping, needed
#           to test the intercept = average change)
#           contrasts: mean_decline, mean_increase, age_faster_decline, age_slower_decline
#   v1_vent v1 plus ventricle size (fornix sensitivity analysis; add it to MODELS)
# Design files: ~/data_dell/fba/covariates/design/<model>/ (build_covariates.py).
#
#   WRITES  ${TEMPLATE_DIR}/stats/<model>_<metric><TAG>_<tract>_ge<cutoff>/
#           (written as .building, renamed when fixelcfestats finishes; a finished
#            run is skipped, so the same command resumes after a stop)
#           ${TEMPLATE_DIR}/stats/run_log.tsv    one line per run: time and exit code
#
# Usage (run in tmux; about 4 to 8 hours for the default 12 runs at 5000 permutations):
#   ./step9_fixelcfestats.sh
#   NSHUFFLES=100 ./step9_fixelcfestats.sh                 # quick dry run (separate folders)
#   MODELS="v1_vent" TRACTS="FX_both" ./step9_fixelcfestats.sh
#   SESSION_FMLS=0.06 ./step9_fixelcfestats.sh             # the 0.06 data
# Settings: MODELS (v1 change)  METRICS (fd log_fc fdc)  TRACTS (FX_both UF_both)
#           CUTOFF (3)  NSHUFFLES (5000)  NTHREADS (48)  MAX_PARALLEL (4)
###############################################################################
set -u
source "$(dirname "$(readlink -f "$0")")/pipeline_config_v2.sh"
MODELS="${MODELS:-v1 change}"; METRICS="${METRICS:-fd log_fc fdc}"; TRACTS="${TRACTS:-FX_both UF_both}"
CUTOFF="${CUTOFF:-3}"; NSHUFFLES="${NSHUFFLES:-5000}"; NTHREADS="${NTHREADS:-48}"; MAX_PAR="${MAX_PARALLEL:-4}"
TAG=""; [ -n "${SESSION_FMLS:-}" ] && TAG="_fmls${SESSION_FMLS//./p}"
DES="${DESIGN_DIR:-$HOME/data_dell/fba/covariates/design}"
STATS="${TEMPLATE_DIR}/stats"; MATRIX="${TEMPLATE_DIR}/step7/matrix"; LOGF="${STATS}/run_log.tsv"
SUFFIX=""; [ "$NSHUFFLES" != "5000" ] && SUFFIX="_n${NSHUFFLES}"
mkdir -p "$STATS"; [ -f "$LOGF" ] || printf 'run\tstart\tminutes\texit\n' > "$LOGF"

# combined contrast files: one row per hypothesis
for model in $MODELS; do
    d="${DES}/${model}"; [ -f "${d}/design.txt" ] || die "No design: ${d}/design.txt (run build_covariates.py)"
    case "$model" in
        change) rows="mean_decline mean_increase age_faster_decline age_slower_decline" ;;
        *)      rows="age_neg age_pos" ;;
    esac
    : > "${d}/contrasts_all.txt"
    for r in $rows; do cat "${d}/contrast_${r}.txt" >> "${d}/contrasts_all.txt"; done
    echo "$rows" > "${d}/contrasts_all_order.txt"
done
[ -d "$MATRIX" ] || die "No connectivity matrix: ${MATRIX} (step 7)"

# list of runs
JOBS=()
for model in $MODELS; do for m in $METRICS; do for t in $TRACTS; do
    case "$model" in change) in="${TEMPLATE_DIR}/change/${m}${TAG}_rate" ;; *) in="${TEMPLATE_DIR}/${m}${TAG}_smooth" ;; esac
    mask="${TEMPLATE_DIR}/tracts/fixels/${t}/${t}_ge${CUTOFF}.mif"
    [ -d "$in" ] || { log "skip ${model} ${m} ${t}: no input folder ${in}"; continue; }
    [ -f "$mask" ] || die "No mask ${mask}"
    n=$(ls "$in"/HCA*.mif 2>/dev/null | wc -l); need=$(wc -l < "${DES}/${model}/files.txt")
    [ "$n" -ge "$need" ] || { log "skip ${model} ${m} ${t}: ${n} of ${need} input files in ${in}"; continue; }
    name="${model}_${m}${TAG}_${t}_ge${CUTOFF}${SUFFIX}"
    [ -f "${STATS}/${name}/.done" ] && { log "done already: ${name}"; continue; }
    JOBS+=("${name}|${model}|${in}|${mask}")
done; done; done
log "Step 9: ${#JOBS[@]} runs to do, ${MAX_PAR} at once, ${NTHREADS} threads each, ${NSHUFFLES} permutations"
[ "${#JOBS[@]}" -gt 0 ] || exit 0

run_one() {
    IFS='|' read -r name model in mask <<< "$1"
    local out="${STATS}/${name}" b="${STATS}/${name}.building" t0 rc extra=""
    [ "$model" = "change" ] && extra="-errors ise"
    rm -rf "$b"; t0=$(date +%s)
    fixelcfestats "$in" "${DES}/${model}/files.txt" "${DES}/${model}/design.txt" "${DES}/${model}/contrasts_all.txt" \
        "$MATRIX" "$b" -mask "$mask" -nshuffles "$NSHUFFLES" -nthreads "$NTHREADS" $extra > "${b}.log" 2>&1
    rc=$?
    printf '%s\t%s\t%s\t%s\n' "$name" "$(date -d @$t0 '+%F %T')" "$(( ($(date +%s) - t0) / 60 ))" "$rc" >> "$LOGF"
    if [ "$rc" = 0 ]; then
        cp "${DES}/${model}/contrasts_all_order.txt" "$b/hypotheses.txt"; mv "${b}.log" "$b/fixelcfestats.log"
        rm -rf "$out"; mv "$b" "$out"; touch "$out/.done"; echo "OK   ${name}"
    else
        echo "FAIL ${name} (exit ${rc}, see ${b}.log)"
    fi
}
export -f run_one; export STATS DES MATRIX NSHUFFLES NTHREADS LOGF
printf '%s\n' "${JOBS[@]}" | xargs -P "$MAX_PAR" -I{} bash -c 'run_one "$@"' _ {}
log "Finished. Times and exit codes: ${LOGF}"
column -t -s$'\t' "$LOGF" | tail -n 30
