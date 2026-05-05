#!/usr/bin/env bash
# ==============================================================================
# 09_evaluate_happy.sh
# ------------------------------------------------------------------------------
# Evaluate every caller's VCF against the GIAB v4.2.1 truth set with hap.py.
#
# Evaluation region = AgilentV5 capture ∩ GIAB high-confidence BED.
# Restricting to this intersection is essential for a fair WES benchmark:
#   - the GIAB BED defines where calls are trusted at all;
#   - the capture BED defines where reads are expected.
# Any region outside both is excluded from TP/FP/FN counting.
#
# Outputs: ${EVAL_DIR}/happy/<caller>/<caller>.summary.csv  (and full hap.py output tree)
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/config.sh"

require_file "${REF}"
require_file "${TRUTH_VCF}"
require_file "${TRUTH_BED}"
require_file "${CAPTURE_BED}"

mkdir -p "${EVAL_DIR}/happy" "$(dirname "${EVAL_BED}")"

# Build the evaluation region (capture ∩ high-confidence) once.
if [[ ! -f "${EVAL_BED}" ]]; then
  echo "[INFO] Building evaluation BED = capture ∩ GIAB-HC -> ${EVAL_BED}"
  bedtools intersect -a "${CAPTURE_BED}" -b "${TRUTH_BED}" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - > "${EVAL_BED}"
fi

CALLERS=(gatk deepvariant strelka2 bcftools freebayes)

for c in "${CALLERS[@]}"; do
  qvcf="${VCF_DIR}/${SM}.${c}.vcf.gz"
  if [[ ! -f "${qvcf}" ]]; then
    echo "[WARN] Skipping ${c}: ${qvcf} not found"
    continue
  fi

  out="${EVAL_DIR}/happy/${c}"
  mkdir -p "${out}"

  echo "[INFO] hap.py ${c}"
  hap.py \
    "${TRUTH_VCF}" \
    "${qvcf}" \
    -r "${REF}" \
    -f "${EVAL_BED}" \
    -o "${out}/${c}" \
    --threads "${THREADS_CALL}" \
    --engine vcfeval \
    --engine-vcfeval-template "${REF%.*}.sdf" \
    ${STRAT_TSV:+--stratification "${STRAT_TSV}"}

  echo "[OK] ${out}/${c}.summary.csv"
done

echo "[DONE] hap.py evaluation complete"
