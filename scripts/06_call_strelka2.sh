#!/usr/bin/env bash
# ==============================================================================
# 06_call_strelka2.sh
# ------------------------------------------------------------------------------
# Call germline variants on the analysis-ready BAM using Strelka2's exome
# workflow, restricted to the AgilentV5 capture region.

# Strelka2 requires the call-regions BED to be bgzipped and tabix-indexed;
# the script handles that automatically the first time it runs.
#
# Inputs:
#   ${OUTDIR}/${SM}.merged.dedup.bam   analysis-ready BAM
#   ${REF}                              reference FASTA + .fai
#   ${CAPTURE_BED}                      capture intervals (bgzipped on first run)
#
# Outputs (in ${VCF_DIR}):
#   ${SM}.strelka2.vcf.gz                final VCF
#   ${SM}.strelka2.vcf.gz.tbi            tabix index
#   plus a Strelka2 work directory at ${OUTDIR}/strelka2_${SM}/
#
# Conda envs (this script switches between them):
#   - `strelka`            for Strelka2 itself (Python 2.7)
#   - `variant_benchmark`  for the bcftools variant-count at the end
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

# ---- Conda init ------------------------------------------------------------
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

# Strelka2 lives in its own env (needs Python 2.7)
if [[ "${CONDA_DEFAULT_ENV:-}" != "strelka" ]]; then
    conda activate strelka
fi

# ---- Inputs and output paths -----------------------------------------------
mkdir -p "${VCF_DIR}" "${OUTDIR}/logs"

BAM="${OUTDIR}/${SM}.merged.dedup.bam"
require_file "${BAM}"
require_file "${REF}"
require_file "${CAPTURE_BED}"

RUN_DIR="${OUTDIR}/strelka2_${SM}"
OUT_VCF="${VCF_DIR}/${SM}.strelka2.vcf.gz"
LOG="${OUTDIR}/logs/06_call_strelka2.log"

log() { echo "[$(date '+%F %T')] $*" | tee -a "${LOG}"; }
: > "${LOG}"

# ---- 1. Prepare bgzipped + tabixed call-regions BED ------------------------
CALL_REGIONS="${CAPTURE_BED}.gz"
if [[ ! -f "${CALL_REGIONS}" || ! -f "${CALL_REGIONS}.tbi" ]]; then
    log "Preparing bgzipped call-regions BED"
    bgzip -kf "${CAPTURE_BED}"
    tabix -p bed -f "${CALL_REGIONS}"
fi

# ---- 2. Configure the Strelka2 workflow ------------------------------------
rm -rf "${RUN_DIR}"
log "Configure Strelka2 germline workflow (exome mode)"
configureStrelkaGermlineWorkflow.py \
    --bam "${BAM}" \
    --referenceFasta "${REF}" \
    --callRegions "${CALL_REGIONS}" \
    --exome \
    --runDir "${RUN_DIR}" 2>> "${LOG}"

# ---- 3. Run the workflow ---------------------------------------------------
log "Run Strelka2 workflow (-j ${THREADS_CALL})"
"${RUN_DIR}/runWorkflow.py" -m local -j "${THREADS_CALL}" 2>> "${LOG}"

# ---- 4. Publish the final VCF to the standard VCF_DIR ----------------------
cp -f "${RUN_DIR}/results/variants/variants.vcf.gz"     "${OUT_VCF}"
cp -f "${RUN_DIR}/results/variants/variants.vcf.gz.tbi" "${OUT_VCF}.tbi"

# ---- 5. Quick variant count (uses bcftools from variant_benchmark env) -----
conda deactivate
conda activate variant_benchmark

n_total=$(bcftools view -H "${OUT_VCF}" 2>/dev/null | wc -l)
n_pass=$(bcftools view -H -f PASS,. "${OUT_VCF}" 2>/dev/null | wc -l)
log "Done: ${n_total} variants total, ${n_pass} PASS"
