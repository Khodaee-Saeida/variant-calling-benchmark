#!/usr/bin/env bash
# ==============================================================================
# 08_call_freebayes.sh
# ------------------------------------------------------------------------------
# Call germline variants on the analysis-ready BAM using FreeBayes, parallelized
# by splitting the AgilentV5 capture region into per-interval chunks.
#
# FreeBayes calling parameters:
#   --min-mapping-quality 20    skip reads with MQ < 20
#   --min-base-quality    20    skip bases with BQ < 20
#   --min-coverage         5    require >= 5 reads to consider a site
#   --use-best-n-alleles   4    cap alt alleles per site (perf vs. completeness)
#
# Inputs:
#   ${OUTDIR}/${SM}.merged.dedup.bam   analysis-ready BAM
#   ${REF}                              reference FASTA + .fai
#   ${CAPTURE_BED}                      capture intervals (split into chunks)
#
# Outputs (in ${VCF_DIR}):
#   ${SM}.freebayes.vcf.gz                final VCF
#   ${SM}.freebayes.vcf.gz.tbi            tabix index
#
# Conda env: `variant_benchmark` (freebayes + bcftools + tabix + GNU parallel)
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
mkdir -p "${VCF_DIR}" "${OUTDIR}/logs" "${OUTDIR}/intervals"

BAM="${OUTDIR}/${SM}.merged.dedup.bam"
require_file "${BAM}"
require_file "${REF}"
require_file "${CAPTURE_BED}"

REGIONS="${OUTDIR}/intervals/freebayes_regions.txt"
RAW="${VCF_DIR}/${SM}.freebayes.raw.vcf"
VCF="${VCF_DIR}/${SM}.freebayes.vcf.gz"
LOG="${OUTDIR}/logs/08_call_freebayes.log"

log() { echo "[$(date '+%F %T')] $*" | tee -a "${LOG}"; }
: > "${LOG}"

# ---- 1. Build region list from the capture BED -----------------------------
awk 'BEGIN{OFS=""}{print $1,":",$2,"-",$3}' "${CAPTURE_BED}" > "${REGIONS}"
n_regions=$(wc -l < "${REGIONS}")

# ---- 2. Run FreeBayes in parallel ------------------------------------------
log "freebayes-parallel (${n_regions} regions, ${THREADS_CALL} threads)"
freebayes-parallel "${REGIONS}" "${THREADS_CALL}" \
        -f "${REF}" \
        --min-mapping-quality 20 \
        --min-base-quality 20 \
        --min-coverage 5 \
        --use-best-n-alleles 4 \
        "${BAM}" > "${RAW}" 2>> "${LOG}"

# ---- 3. Normalize and filter -----------------------------------------------
log "Normalize + filter -> ${VCF##*/}"
bcftools norm -f "${REF}" -m -any -O u "${RAW}" 2>> "${LOG}" \
  | bcftools filter -e 'QUAL<20' -s LowQual -O z -o "${VCF}" 2>> "${LOG}"

tabix -p vcf -f "${VCF}"
rm -f "${RAW}"

# ---- Summary ---------------------------------------------------------------
n_total=$(bcftools view -H "${VCF}" 2>/dev/null | wc -l)
n_pass=$(bcftools view -H -f PASS,. "${VCF}" 2>/dev/null | wc -l)
log "Done: ${n_total} variants total, ${n_pass} PASS"
