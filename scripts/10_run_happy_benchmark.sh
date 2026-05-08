#!/usr/bin/env bash
# ==============================================================================
# 10_run_happy_benchmark.sh
# ------------------------------------------------------------------------------
# Benchmark the four germline callers (GATK, Strelka2, bcftools, FreeBayes)
# against the GIAB v4.2.1 truth set using hap.py + RTG vcfeval.
#
# This script assumes:
#   - 03_prepare_inputs.sh has run (truth files downloaded, RTG SDF built)
#   - 09_prepare_vcfs.sh has run (per-caller VCFs normalized, PASS-filtered,
#     evaluable region BED computed)
#
# Pipeline per caller:
#
#   PASS-filtered VCF
#         │
#         ▼  hap.py --engine=vcfeval
#         │  truth = ${TRUTH_VCF}, high-conf = ${TRUTH_BED}
#         │  target = ${EVAL_BED}, reference = ${REF}
#         │
#         ▼  ${SM}.${caller}.summary.csv     per-caller TP/FP/FN, recall, precision, F1
#            ${SM}.${caller}.extended.csv    per-stratification breakdown
#            ${SM}.${caller}.vcf.gz          annotated comparison VCF
#
# Then aggregated:
#   ${SM}.${caller}.summary.csv x 4   ->   all_callers_summary.csv
#
# Inputs (all from config/config.sh):
#   ${REF}, ${TRUTH_VCF}, ${TRUTH_BED}, ${EVAL_BED}, ${SDF_DIR}
#   ${VCF_DIR}/${SM}.<caller>.pass.vcf.gz   one per caller
#
# Outputs (in ${HAPPY_OUT_DIR}):
#   per-caller hap.py outputs + aggregated all_callers_summary.csv
#
# Conda env: `happy` (hap.py + RTG vcfeval)
# ==============================================================================
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

for v in SM REF VCF_DIR TRUTH_VCF TRUTH_BED EVAL_BED SDF_DIR HAPPY_OUT_DIR THREADS_CALL; do
    [[ -n "${!v:-}" ]] || {
        echo "[ERROR] ${v} not set in ${CONFIG}" >&2
        exit 1
    }
done

# ---- Required inputs -------------------------------------------------------
require_file "${REF}"
require_file "${TRUTH_VCF}"
require_file "${TRUTH_BED}"
require_file "${EVAL_BED}"
[[ -d "${SDF_DIR}" ]] || {
    echo "[ERROR] RTG SDF not found at ${SDF_DIR}" >&2
    echo "        Run scripts/03_prepare_inputs.sh first." >&2
    exit 1
}

# ---- Conda init ------------------------------------------------------------
if [[ -z "${CONDA_EXE:-}" ]]; then
    for base in "$HOME/miniconda3" "$HOME/anaconda3" "/opt/conda" "/opt/miniconda3"; do
        if [[ -f "${base}/etc/profile.d/conda.sh" ]]; then
            # shellcheck disable=SC1091
            source "${base}/etc/profile.d/conda.sh"
            break
        fi
    done
else
    # shellcheck disable=SC1091
    source "$(dirname "$(dirname "${CONDA_EXE}")")/etc/profile.d/conda.sh"
fi

CALLERS=(gatk strelka2 bcftools freebayes)

mkdir -p "${HAPPY_OUT_DIR}"
log() { echo "[$(date '+%F %T')] $*"; }
hr()  { echo "============================================================"; }

# ============================================================================
# Run hap.py per caller
# ============================================================================
conda activate happy

FAIL_COUNT=0

for caller in "${CALLERS[@]}"; do
    QUERY_VCF="${VCF_DIR}/${SM}.${caller}.pass.vcf.gz"
    OUT_PREFIX="${HAPPY_OUT_DIR}/${SM}.${caller}"
    HAPPY_LOG="${HAPPY_OUT_DIR}/${SM}.${caller}.happy.log"

    if [[ ! -f "${QUERY_VCF}" ]]; then
        log "SKIP ${caller}: ${QUERY_VCF} not found (run 09_prepare_vcfs.sh first)"
        FAIL_COUNT=$((FAIL_COUNT+1))
        continue
    fi

    echo
    hr
    log "Benchmarking: ${caller}"
    hr

    if hap.py \
            "${TRUTH_VCF}" \
            "${QUERY_VCF}" \
            -f "${TRUTH_BED}" \
            -T "${EVAL_BED}" \
            -r "${REF}" \
            --engine=vcfeval \
            --engine-vcfeval-template "${SDF_DIR}" \
            --threads "${THREADS_CALL}" \
            -o "${OUT_PREFIX}" \
            > "${HAPPY_LOG}" 2>&1
    then
        if [[ -f "${OUT_PREFIX}.summary.csv" ]]; then
            log "${caller} done -> ${OUT_PREFIX}.summary.csv"
        else
            log "${caller} ran but no summary produced -- check ${HAPPY_LOG}"
            FAIL_COUNT=$((FAIL_COUNT+1))
        fi
    else
        log "${caller} FAILED -- check ${HAPPY_LOG}"
        FAIL_COUNT=$((FAIL_COUNT+1))
    fi
done

conda deactivate

# ============================================================================
# Aggregate per-caller summaries into one CSV
# ============================================================================
SUMMARY="${HAPPY_OUT_DIR}/all_callers_summary.csv"
log "Aggregating per-caller summaries -> ${SUMMARY}"

first=1
for caller in "${CALLERS[@]}"; do
    f="${HAPPY_OUT_DIR}/${SM}.${caller}.summary.csv"
    if [[ ! -f "${f}" ]]; then
        log "  (skip ${caller}, no summary)"
        continue
    fi
    if (( first )); then
        head -1 "${f}" | sed 's/^/caller,/' > "${SUMMARY}"
        first=0
    fi
    tail -n +2 "${f}" | sed "s/^/${caller},/" >> "${SUMMARY}"
done

# ============================================================================
# Final report
# ============================================================================
echo
hr
if [[ "${FAIL_COUNT}" -eq 0 ]]; then
    echo "ALL DONE -- 4/4 callers benchmarked"
else
    echo "DONE WITH FAILURES -- ${FAIL_COUNT}/4 callers failed"
fi
hr
echo

if [[ -f "${SUMMARY}" ]]; then
    column -t -s, "${SUMMARY}"
    echo
    echo "Full results    : ${HAPPY_OUT_DIR}"
    echo "Combined table  : ${SUMMARY}"
fi

exit "${FAIL_COUNT}"
