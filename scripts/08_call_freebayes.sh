#!/usr/bin/env bash
# ==============================================================================
# 08_call_freebayes.sh
# ------------------------------------------------------------------------------
# FreeBayes (Bayesian haplotype-based caller). We split the capture BED into
# regions and run freebayes-parallel for speed; outputs are concatenated and
# normalised with bcftools.
#
# Output: ${VCF_DIR}/${SM}.freebayes.vcf.gz
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

REGIONS="${OUTDIR}/intervals/freebayes_regions.txt"
mkdir -p "$(dirname "${REGIONS}")"
awk 'BEGIN{OFS=""}{print $1,":",$2,"-",$3}' "${CAPTURE_BED}" > "${REGIONS}"

RAW="${VCF_DIR}/${SM}.freebayes.raw.vcf"
VCF="${VCF_DIR}/${SM}.freebayes.vcf.gz"

echo "[INFO] freebayes-parallel ($(wc -l < "${REGIONS}") regions, ${THREADS_CALL} threads)"
freebayes-parallel "${REGIONS}" "${THREADS_CALL}" \
  -f "${REF}" \
  --min-mapping-quality 20 \
  --min-base-quality 20 \
  --min-coverage 5 \
  --use-best-n-alleles 4 \
  "${BAM}" > "${RAW}"

echo "[INFO] Normalize, decompose multi-allelics, bgzip+tabix -> ${VCF}"
bcftools norm -f "${REF}" -m -any -O u "${RAW}" \
  | bcftools filter -e 'QUAL<20' -s LowQual -O z -o "${VCF}"
tabix -p vcf -f "${VCF}"
rm -f "${RAW}"

echo "[DONE] FreeBayes -> ${VCF}"
