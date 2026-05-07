#!/usr/bin/env bash
# ==============================================================================
# 04_call_gatk.sh
# ------------------------------------------------------------------------------
# Call germline variants on the analysis-ready BAM using GATK4
# (HaplotypeCaller → GenotypeGVCFs → hard-filter), restricted to the AgilentV5
# capture region with a 100 bp interval padding.
#
#
# Filters used are the standard GATK Best Practices thresholds for germline
# short variants. They were tuned on WGS and are known to be over-aggressive
# on exome capture data (especially MQ < 40 and ReadPosRankSum < -8)
#
# Inputs:
#   ${OUTDIR}/${SM}.merged.dedup.bam   analysis-ready BAM (from 02_validate_bam.sh)
#   ${REF}                              reference FASTA + .fai + .dict
#   ${CAPTURE_BED}                      capture intervals
#
# Outputs (in ${VCF_DIR}):
#   ${SM}.gatk.g.vcf.gz                  per-sample GVCF (intermediate)
#   ${SM}.gatk.raw.vcf.gz                genotyped, unfiltered
#   ${SM}.gatk.{snps,indels}.vcf.gz      split by variant type
#   ${SM}.gatk.{snps,indels}.filt.vcf.gz hard-filtered
#   ${SM}.gatk.vcf.gz                    merged final VCF (benchmark this one)
#
# Conda env: `variant_benchmark` (gatk4 + bcftools + tabix)
# ==============================================================================
set -euo pipefail

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

for v in SM REF OUTDIR VCF_DIR CAPTURE_BED THREADS_CALL; do
    [[ -n "${!v:-}" ]] || {
        echo "[ERROR] ${v} not set in ${CONFIG}" >&2
        exit 1
    }
done

# ---- Conda env -------------------------------------------------------------
if [[ "${CONDA_DEFAULT_ENV:-}" != "variant_benchmark" ]]; then
    if [[ -z "${CONDA_EXE:-}" ]]; then
        for base in "$HOME/miniconda3" "$HOME/anaconda3" "/opt/conda" "/opt/miniconda3"; do
            if [[ -f "${base}/etc/profile.d/conda.sh" ]]; then
                source "${base}/etc/profile.d/conda.sh"
                break
            fi
        done
    else
        source "$(dirname "$(dirname "${CONDA_EXE}")")/etc/profile.d/conda.sh"
    fi
    conda activate variant_benchmark
fi

# ---- Inputs and output paths -----------------------------------------------
mkdir -p "${VCF_DIR}" "${OUTDIR}/logs"

BAM="${OUTDIR}/${SM}.merged.dedup.bam"
require_file "${BAM}"
require_file "${REF}"
require_file "${CAPTURE_BED}"

GVCF="${VCF_DIR}/${SM}.gatk.g.vcf.gz"
RAW="${VCF_DIR}/${SM}.gatk.raw.vcf.gz"
SNPS="${VCF_DIR}/${SM}.gatk.snps.vcf.gz"
INDELS="${VCF_DIR}/${SM}.gatk.indels.vcf.gz"
SNP_F="${VCF_DIR}/${SM}.gatk.snps.filt.vcf.gz"
INDEL_F="${VCF_DIR}/${SM}.gatk.indels.filt.vcf.gz"
FINAL="${VCF_DIR}/${SM}.gatk.vcf.gz"
LOG="${OUTDIR}/logs/04_call_gatk.log"

log() { echo "[$(date '+%F %T')] $*" | tee -a "${LOG}"; }
: > "${LOG}"

# ---- 1. HaplotypeCaller (per-sample GVCF) ----------------------------------
log "HaplotypeCaller -> ${GVCF##*/}"
gatk HaplotypeCaller \
    -R "${REF}" \
    -I "${BAM}" \
    -L "${CAPTURE_BED}" \
    -ip 100 \
    -ERC GVCF \
    --native-pair-hmm-threads "${THREADS_CALL}" \
    -O "${GVCF}" 2>> "${LOG}"

# ---- 2. GenotypeGVCFs (single-sample joint genotyping) ---------------------
log "GenotypeGVCFs -> ${RAW##*/}"
gatk GenotypeGVCFs \
    -R "${REF}" \
    -V "${GVCF}" \
    -L "${CAPTURE_BED}" \
    -ip 100 \
    -O "${RAW}" 2>> "${LOG}"

# ---- 3. Split by variant type ----------------------------------------------
log "Split SNPs and INDELs"
gatk SelectVariants -R "${REF}" -V "${RAW}" \
    --select-type-to-include SNP   -O "${SNPS}"   2>> "${LOG}"
gatk SelectVariants -R "${REF}" -V "${RAW}" \
    --select-type-to-include INDEL -O "${INDELS}" 2>> "${LOG}"

# ---- 4. GATK Best Practices hard filters -----------------------------------
log "Hard filter SNPs"
gatk VariantFiltration -R "${REF}" -V "${SNPS}" \
    --filter-expression "QD < 2.0"              --filter-name "QD2"   \
    --filter-expression "FS > 60.0"             --filter-name "FS60"  \
    --filter-expression "MQ < 40.0"             --filter-name "MQ40"  \
    --filter-expression "MQRankSum < -12.5"     --filter-name "MQRS"  \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "RPRS"  \
    --filter-expression "SOR > 3.0"             --filter-name "SOR3"  \
    -O "${SNP_F}" 2>> "${LOG}"

log "Hard filter INDELs"
gatk VariantFiltration -R "${REF}" -V "${INDELS}" \
    --filter-expression "QD < 2.0"               --filter-name "QD2"   \
    --filter-expression "FS > 200.0"             --filter-name "FS200" \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "RPRS"  \
    --filter-expression "SOR > 10.0"             --filter-name "SOR10" \
    -O "${INDEL_F}" 2>> "${LOG}"

# ---- 5. Merge filtered SNPs + INDELs into the final VCF --------------------
log "Merge filtered SNPs + INDELs -> ${FINAL##*/}"
gatk MergeVcfs -I "${SNP_F}" -I "${INDEL_F}" -O "${FINAL}" 2>> "${LOG}"
[[ -f "${FINAL}.tbi" ]] || tabix -p vcf -f "${FINAL}"

n_total=$(bcftools view -H "${FINAL}" 2>/dev/null | wc -l)
n_pass=$(bcftools view -H -f PASS,. "${FINAL}" 2>/dev/null | wc -l)
log "Done: ${n_total} variants total, ${n_pass} PASS"
