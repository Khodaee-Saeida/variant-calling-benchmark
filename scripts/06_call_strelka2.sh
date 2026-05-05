#!/usr/bin/env bash
# ==============================================================================
# 06_call_strelka2.sh
# ------------------------------------------------------------------------------
# Strelka2 germline workflow on a single WES sample.
#
# Strelka2 expects a bgzipped+tabixed BED for --callRegions. We build that
# from the AgilentV5 capture BED on first run.
#
# Output: ${VCF_DIR}/${SM}.strelka2.vcf.gz
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/config.sh"

mkdir -p "${VCF_DIR}"

BAM="${OUTDIR}/${SM}.merged.dedup.bqsr.bam"
[[ -f "${BAM}" ]] || BAM="${OUTDIR}/${SM}.merged.dedup.bam"
require_file "${BAM}"
require_file "${REF}"
require_file "${CAPTURE_BED}"

# Strelka requires a bgzipped + tabixed BED.
CALL_REGIONS="${CAPTURE_BED%.bed}.bed.gz"
if [[ ! -f "${CALL_REGIONS}" ]]; then
  echo "[INFO] bgzip + tabix capture BED for Strelka2"
  bgzip -kf "${CAPTURE_BED}"
  tabix -p bed -f "${CALL_REGIONS}"
fi

RUN_DIR="${OUTDIR}/strelka2_${SM}"
rm -rf "${RUN_DIR}"

echo "[INFO] Configuring Strelka2 germline workflow"
configureStrelkaGermlineWorkflow.py \
  --bam "${BAM}" \
  --referenceFasta "${REF}" \
  --callRegions "${CALL_REGIONS}" \
  --exome \
  --runDir "${RUN_DIR}"

echo "[INFO] Running Strelka2 (-m local -j ${THREADS_CALL})"
"${RUN_DIR}/runWorkflow.py" -m local -j "${THREADS_CALL}"

# Canonical output location for downstream evaluation scripts.
cp -f "${RUN_DIR}/results/variants/variants.vcf.gz"     "${VCF_DIR}/${SM}.strelka2.vcf.gz"
cp -f "${RUN_DIR}/results/variants/variants.vcf.gz.tbi" "${VCF_DIR}/${SM}.strelka2.vcf.gz.tbi"

echo "[DONE] Strelka2 -> ${VCF_DIR}/${SM}.strelka2.vcf.gz"
