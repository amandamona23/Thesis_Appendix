#!/bin/bash
###############################################################################
# STEP 2 - Group response functions (HARD BARRIER)
#
# Averages per-subject response functions into one shared set. Every subject's
# FODs are then estimated against the same reference, so differences that show
# up downstream are real data differences rather than each subject having been
# calibrated to their own yardstick.
#
# This MUST run after step 1 has completed AND passed QC for every session, and
# BEFORE any dwi2fod call. It is a single coordinated command, not parallel.
#
# ONE SESSION PER SUBJECT (default): using both visits per subject would weight
# 2-visit subjects twice in the average. Response functions are robust enough
# that this is minor, but one-per-subject costs nothing. Pass --all-sessions to
# override.
#
# Reads responses from ${IN_ROOT} (symlinks to RDS, see build_session_index.sh)
# Usage: ./step2_group_responsemean_v2.sh --list sessions.txt
#        ./step2_group_responsemean.sh --all-sessions --list sessions.txt
###############################################################################

set -u
source "$(dirname "$0")/pipeline_config.sh"

ONE_PER_SUBJECT=1
if [ "${1:-}" == "--all-sessions" ]; then
    ONE_PER_SUBJECT=0
    shift
    log "Using ALL sessions (2-visit subjects will be double-weighted)"
fi

parse_sessions "$@"
mkdir -p "$GROUP_DIR"

# ---------------------------------------------------------------------------
# Reduce to one session per subject unless overridden.
# Session IDs look like HCA6002236_V1_MR, so the subject is field 1.
# Prefers V1; falls back to the first available visit.
# ---------------------------------------------------------------------------
if [ "$ONE_PER_SUBJECT" -eq 1 ]; then
    declare -A picked
    SELECTED=()
    for sess in "${SESSIONS[@]}"; do
        subj="${sess%%_*}"
        if [ -z "${picked[$subj]:-}" ]; then
            picked[$subj]="$sess"
        elif [[ "$sess" == *_V1_MR ]]; then
            picked[$subj]="$sess"   # prefer V1 if we see it later
        fi
    done
    for subj in "${!picked[@]}"; do
        SELECTED+=("${picked[$subj]}")
    done
    SESSIONS=("${SELECTED[@]}")
    log "Reduced to ${#SESSIONS[@]} sessions (one per subject)"
fi

# ---------------------------------------------------------------------------
# Precondition check. responsemean silently averages whatever it is given, so
# a session that failed QC but still has output files would be included without
# any error. Check explicitly.
# ---------------------------------------------------------------------------
WM_FILES=(); GM_FILES=(); CSF_FILES=()
missing=0

for sess in "${SESSIONS[@]}"; do
    wm="${IN_ROOT}/${sess}/response/wm.txt"
    gm="${IN_ROOT}/${sess}/response/gm.txt"
    csf="${IN_ROOT}/${sess}/response/csf.txt"

    for f in "$wm" "$gm" "$csf"; do
        if [ ! -s "$f" ]; then
            log "MISSING or EMPTY: $f"
            missing=1
        fi
    done

    WM_FILES+=("$wm"); GM_FILES+=("$gm"); CSF_FILES+=("$csf")
done

if [ "$missing" -eq 1 ]; then
    die "One or more response files missing. Fix step 1 before running responsemean."
fi

# ---------------------------------------------------------------------------
# Manual QC confirmation. This is the point of no return for the whole FBA and
# DTI arms: every downstream FOD depends on these three files.
# ---------------------------------------------------------------------------
cat <<EOF

-------------------------------------------------------------------
About to average responses across ${#SESSIONS[@]} sessions.

Confirm before proceeding:
  - Step 1 QC is 100 percent clean (no unresolved flagged sessions)
  - No known-bad sessions remain in the list
  - Cohort composition is final

Sessions:
$(printf '  %s\n' "${SESSIONS[@]}")
-------------------------------------------------------------------
EOF

read -r -p "Proceed? [y/N] " confirm
[[ "$confirm" =~ ^[Yy]$ ]] || die "Aborted by user."

# Record exactly which sessions built the group response, for the methods section
printf '%s\n' "${SESSIONS[@]}" > "${GROUP_DIR}/sessions_used.txt"

log "Averaging responses..."
responsemean "${WM_FILES[@]}"  "${GROUP_DIR}/group_wm.txt"  -force
responsemean "${GM_FILES[@]}"  "${GROUP_DIR}/group_gm.txt"  -force
responsemean "${CSF_FILES[@]}" "${GROUP_DIR}/group_csf.txt" -force

log "Step 2 complete. Group responses in ${GROUP_DIR}"
log "Session list recorded in ${GROUP_DIR}/sessions_used.txt"
log "Next: step3a_fod_normalise.sh"
