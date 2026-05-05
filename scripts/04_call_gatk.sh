#!/usr/bin/env bash
# ==============================================================================
# 04_call_gatk.sh
# ------------------------------------------------------------------------------
# GATK HaplotypeCaller in GVCF mode -> GenotypeGVCFs -> single-sample VCF.
# Confined to the AgilentV5 capture intervals with a 100 bp padding on each
# side (-ip 100), which matches the recommended Best-Practices behaviour for
# exome data.
#
# Output: ${VCF_DIR}/${SM}.gatk.vcf.gz
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/config.sh"

mkdir -p "${VCF_DIR}"

# Use BQSR'd BAM if it exists, else the dedup BAM.
BAM="${OUTDIR}/${SM}.merged.dedup.bqsr.bam"
[[ -f "${BAM}" ]] || BAM="${OUTDIR}/${SM}.merged.dedup.bam"
require_file "${BAM}"
require_file "${REF}"
require_file "${CAPTURE_BED}"

GVCF="${VCF_DIR}/${SM}.gatk.g.vcf.gz"
VCF="${VCF_DIR}/${SM}.gatk.vcf.gz"

echo "[INFO] HaplotypeCaller -> ${GVCF}"
gatk HaplotypeCaller \
  -R "${REF}" \
  -I "${BAM}" \
  -L "${CAPTURE_BED}" -ip 100 \
  -ERC GVCF \
  --native-pair-hmm-threads "${THREADS_CALL}" \
  -O "${GVCF}"

echo "[INFO] GenotypeGVCFs -> ${VCF}"
gatk GenotypeGVCFs \
  -R "${REF}" \
  -V "${GVCF}" \
  -L "${CAPTURE_BED}" -ip 100 \
  -O "${VCF}"

# Hard-filter recommended for single-sample exome (VQSR needs cohort-scale data).
SNP_VCF="${VCF_DIR}/${SM}.gatk.snps.filtered.vcf.gz"
INDEL_VCF="${VCF_DIR}/${SM}.gatk.indels.filtered.vcf.gz"
FINAL="${VCF_DIR}/${SM}.gatk.filtered.vcf.gz"

echo "[INFO] Hard filtering (GATK recommended thresholds)"
gatk SelectVariants -R "${REF}" -V "${VCF}" --select-type-to-include SNP \
  -O "${VCF_DIR}/${SM}.gatk.snps.vcf.gz"
gatk SelectVariants -R "${REF}" -V "${VCF}" --select-type-to-include INDEL \
  -O "${VCF_DIR}/${SM}.gatk.indels.vcf.gz"

gatk VariantFiltration -R "${REF}" -V "${VCF_DIR}/${SM}.gatk.snps.vcf.gz" \
  --filter-expression "QD < 2.0"              --filter-name "QD2"   \
  --filter-expression "FS > 60.0"             --filter-name "FS60"  \
  --filter-expression "MQ < 40.0"             --filter-name "MQ40"  \
  --filter-expression "MQRankSum < -12.5"     --filter-name "MQRS"  \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "RPRS"  \
  --filter-expression "SOR > 3.0"             --filter-name "SOR3"  \
  -O "${SNP_VCF}"

gatk VariantFiltration -R "${REF}" -V "${VCF_DIR}/${SM}.gatk.indels.vcf.gz" \
  --filter-expression "QD < 2.0"              --filter-name "QD2"   \
  --filter-expression "FS > 200.0"            --filter-name "FS200" \
  --filter-expression "ReadPosRankSum < -20.0" --filter-name "RPRS" \
  --filter-expression "SOR > 10.0"            --filter-name "SOR10" \
  -O "${INDEL_VCF}"

echo "[INFO] Merging filtered SNPs + INDELs -> ${FINAL}"
gatk MergeVcfs -I "${SNP_VCF}" -I "${INDEL_VCF}" -O "${FINAL}"

# Convenience: a copy under the canonical name expected by evaluation scripts.
cp -f "${FINAL}"      "${VCF_DIR}/${SM}.gatk.vcf.gz"
cp -f "${FINAL}.tbi"  "${VCF_DIR}/${SM}.gatk.vcf.gz.tbi" 2>/dev/null || \
  tabix -p vcf "${VCF_DIR}/${SM}.gatk.vcf.gz"

echo "[DONE] GATK -> ${VCF_DIR}/${SM}.gatk.vcf.gz"
