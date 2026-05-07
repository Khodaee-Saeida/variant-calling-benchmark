#!/usr/bin/env bash
# ==============================================================================
# run_all_callers.sh
# ------------------------------------------------------------------------------
# Master runner: invoke all four germline variant callers in sequence and
# write a summary table with wall time and variant counts.
#
# Order of execution:
#   1. GATK4 HaplotypeCaller   (04_call_gatk.sh)
#   2. Strelka2 germline        (06_call_strelka2.sh)
#   3. bcftools mpileup+call    (07_call_bcftools.sh)
#   4. FreeBayes                (08_call_freebayes.sh)
#
# Each caller runs to completion before the next starts; failures don't abort
# the run, so a crash in one caller doesn't waste the time spent on the others.
#
# Output:
#   ${OUTDIR}/logs/caller_summary.tsv     wall time, total/PASS variants per caller
#   plus per-caller logs in ${OUTDIR}/logs/
#
# Usage:
#   bash scripts/run_all_callers.sh
#
# Conda envs: managed by each individual caller script. This wrapper does not
# need to activate one.
# ==============================================================================
# Note: -e is intentionally NOT set, so a single caller failure doesn't abort
# the whole run. Each caller's exit status is captured and reported.
set -uo pipefail

# ---- Configuration ---------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${SCRIPT_DIR}/../config/config.sh"

if [[ ! -f "${CONFIG}" ]]; then
    cat >&2 <<EOF
[ERROR] config/config.sh not found.

Setup:
  cp config/config.example.sh config/config.sh
  \$EDITOR config/config.sh    # edit paths for your system

See README.md for details.
EOF
    exit 1
fi

# shellcheck disable=SC1090
source "${CONFIG}"

for v in SM OUTDIR VCF_DIR; do
    [[ -n "${!v:-}" ]] || {
        echo "[ERROR] ${v} not set in ${CONFIG}" >&2
        exit 1
    }
done

mkdir -p "${OUTDIR}/logs"
SUMMARY="${OUTDIR}/logs/caller_summary.tsv"
printf 'caller\twall_time_sec\ttotal_variants\tpass_variants\tstatus\n' > "${SUMMARY}"

# ---- Helpers ---------------------------------------------------------------
hr() {
    echo "============================================================"
}

run_caller() {
    local name="$1" script="$2"
    local start end vcf n_total n_pass

    echo
    hr
    echo "[$(date '+%F %T')] Starting ${name}"
    hr

    start=$(date +%s)
    if bash "${script}"; then
        end=$(($(date +%s) - start))
        vcf="${VCF_DIR}/${SM}.${name}.vcf.gz"
        n_total=0
        n_pass=0
        if [[ -f "${vcf}" ]]; then
            n_total=$(bcftools view -H "${vcf}" 2>/dev/null | wc -l)
            n_pass=$(bcftools view -H -f PASS,. "${vcf}" 2>/dev/null | wc -l)
        fi
        printf '%s\t%d\t%d\t%d\tOK\n' "${name}" "${end}" "${n_total}" "${n_pass}" >> "${SUMMARY}"
        echo "[$(date '+%F %T')] ${name} OK -- ${end}s, ${n_total} variants (${n_pass} PASS)"
    else
        end=$(($(date +%s) - start))
        printf '%s\t%d\t-\t-\tFAIL\n' "${name}" "${end}" >> "${SUMMARY}"
        echo "[$(date '+%F %T')] ${name} FAILED after ${end}s -- check ${OUTDIR}/logs/"
        FAIL_COUNT=$((FAIL_COUNT+1))
    fi
}

# ---- Run all callers -------------------------------------------------------
FAIL_COUNT=0
TOTAL_START=$(date +%s)

run_caller gatk      "${SCRIPT_DIR}/04_call_gatk.sh"
run_caller strelka2  "${SCRIPT_DIR}/06_call_strelka2.sh"
run_caller bcftools  "${SCRIPT_DIR}/07_call_bcftools.sh"
run_caller freebayes "${SCRIPT_DIR}/08_call_freebayes.sh"

TOTAL_TIME=$(($(date +%s) - TOTAL_START))

# ---- Final summary ---------------------------------------------------------
echo
hr
if [[ "${FAIL_COUNT}" -eq 0 ]]; then
    echo "ALL DONE -- 4/4 callers OK (total wall time: ${TOTAL_TIME}s)"
else
    echo "DONE WITH FAILURES -- ${FAIL_COUNT}/4 callers failed (total wall time: ${TOTAL_TIME}s)"
fi
hr
echo
column -t -s $'\t' "${SUMMARY}"
echo
echo "Summary file: ${SUMMARY}"

exit "${FAIL_COUNT}"
