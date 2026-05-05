#!/usr/bin/env bash
# ==============================================================================
# 10_evaluate_vcfeval.sh
# ------------------------------------------------------------------------------
# Independent cross-check of hap.py results using RTG vcfeval. The two tools
# use different match algorithms (hap.py: xcmp/vcfeval; RTG: ROC-aware
# pathfinding), so agreement between them is a strong sanity check.
#
# Outputs: ${EVAL_DIR}/vcfeval/<caller>/{summary.txt,weighted_roc.tsv.gz,...}
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/config.sh"

require_file "${TRUTH_VCF}"
require_file "${EVAL_BED}"

SDF="${REF%.*}.sdf"
[[ -d "${SDF}" ]] || { echo "[ERROR] RTG SDF not found: ${SDF} (run 01_index_reference.sh)"; exit 1; }

mkdir -p "${EVAL_DIR}/vcfeval"

CALLERS=(gatk deepvariant strelka2 bcftools freebayes)

for c in "${CALLERS[@]}"; do
  qvcf="${VCF_DIR}/${SM}.${c}.vcf.gz"
  if [[ ! -f "${qvcf}" ]]; then
    echo "[WARN] Skipping ${c}: ${qvcf} not found"
    continue
  fi

  out="${EVAL_DIR}/vcfeval/${c}"
  rm -rf "${out}"   # vcfeval insists the dir not exist

  echo "[INFO] rtg vcfeval ${c}"
  rtg vcfeval \
    -b "${TRUTH_VCF}" \
    -c "${qvcf}" \
    -t "${SDF}" \
    --bed-regions "${EVAL_BED}" \
    -o "${out}" \
    --threads "${THREADS_CALL}" \
    --output-mode annotate

  echo "[OK] ${out}/summary.txt"
done

echo "[DONE] vcfeval evaluation complete"
