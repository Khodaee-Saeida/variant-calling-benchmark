#!/usr/bin/env bash
# ==============================================================================
# 07_call_bcftools.sh
# ------------------------------------------------------------------------------
# Call germline variants on the analysis-ready BAM using bcftools, restricted
# to the AgilentV5 capture region.
#
# Filter design (matches GATK Best Practices spirit, tuned for bcftools):
#   - QUAL >= 20  : per-site phred-scaled quality threshold
#   - DP   >= 10  : minimum total depth at the site
# Records failing either are kept but tagged "LowQual"; downstream the
# benchmark uses `bcftools view -f PASS,.` to keep only PASS records.
#
# Inputs:
#   ${OUTDIR}/${SM}.merged.dedup.bam   analysis-ready BAM
#   ${REF}                              reference FASTA + .fai
#   ${CAPTURE_BED}                      capture intervals
#
# Outputs (in ${VCF_DIR}):
#   ${SM}.bcftools.vcf.gz                final VCF
#   ${SM}.bcftools.vcf.gz.tbi            tabix index
#
# Conda env: `variant_benchmark` (bcftools + tabix)
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

VCF="${VCF_DIR}/${SM}.bcftools.vcf.gz"
LOG="${OUTDIR}/logs/07_call_bcftools.log"

log() { echo "[$(date '+%F %T')] $*" | tee -a "${LOG}"; }
: > "${LOG}"

# ---- Variant calling pipeline ----------------------------------------------
log "bcftools mpileup | call | norm | filter -> ${VCF##*/}"

bcftools mpileup \
        --threads "${THREADS_CALL}" \
        -f "${REF}" \
        -R "${CAPTURE_BED}" \
        --max-depth 10000 \
        -q 20 -Q 20 \
        -a FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD \
        -O u "${BAM}" 2>> "${LOG}" \
  | bcftools call \
        --threads "${THREADS_CALL}" \
        -m -v -f GQ,GP \
        -O u 2>> "${LOG}" \
  | bcftools norm \
        --threads "${THREADS_CALL}" \
        -f "${REF}" \
        -m -any \
        -O u 2>> "${LOG}" \
  | bcftools filter \
        -e 'QUAL<20 || INFO/DP<10' \
        -s LowQual \
        -O z -o "${VCF}" 2>> "${LOG}"

tabix -p vcf -f "${VCF}"

# ---- Summary ---------------------------------------------------------------
n_total=$(bcftools view -H "${VCF}" 2>/dev/null | wc -l)
n_pass=$(bcftools view -H -f PASS,. "${VCF}" 2>/dev/null | wc -l)
log "Done: ${n_total} variants total, ${n_pass} PASS"
