#!/usr/bin/env bash
# ==============================================================================
# 07_call_bcftools.sh
# ------------------------------------------------------------------------------
# bcftools mpileup | bcftools call (multiallelic caller, -mv).
#
# Conservative defaults appropriate for a WES single-sample benchmark:
#   -q 20  : minimum mapping quality
#   -Q 20  : minimum base quality
#   --max-depth 10000 : avoid mpileup's default cap masking deep exome regions
#
# Output: ${VCF_DIR}/${SM}.bcftools.vcf.gz
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

VCF="${VCF_DIR}/${SM}.bcftools.vcf.gz"

echo "[INFO] bcftools mpileup | bcftools call -> ${VCF}"
bcftools mpileup \
    --threads "${THREADS_CALL}" \
    -f "${REF}" \
    -R "${CAPTURE_BED}" \
    --max-depth 10000 \
    -q 20 -Q 20 \
    -a FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD \
    -O u \
    "${BAM}" \
  | bcftools call \
    --threads "${THREADS_CALL}" \
    -m -v \
    -f GQ,GP \
    -O u \
  | bcftools norm \
    --threads "${THREADS_CALL}" \
    -f "${REF}" \
    -m -any \
    -O z -o "${VCF}"

tabix -p vcf -f "${VCF}"

# Light filtering — keep PASS but mark low-quality records.
FILT="${VCF_DIR}/${SM}.bcftools.filtered.vcf.gz"
bcftools filter -e 'QUAL<20 || INFO/DP<10' -s LowQual -O z -o "${FILT}" "${VCF}"
tabix -p vcf -f "${FILT}"

cp -f "${FILT}"     "${VCF}"
cp -f "${FILT}.tbi" "${VCF}.tbi"

echo "[DONE] bcftools -> ${VCF}"
