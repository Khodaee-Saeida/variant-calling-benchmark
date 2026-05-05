#!/usr/bin/env bash
# ==============================================================================
# 05_call_deepvariant.sh
# ------------------------------------------------------------------------------
# DeepVariant via the official Singularity (or Docker) image.
# Model type WES is selected — DeepVariant ships separate models for WGS,
# WES, PACBIO, and HYBRID; using the wrong one materially hurts accuracy.
#
# Output: ${VCF_DIR}/${SM}.deepvariant.vcf.gz
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

VCF="${VCF_DIR}/${SM}.deepvariant.vcf.gz"
GVCF="${VCF_DIR}/${SM}.deepvariant.g.vcf.gz"

# Folder mounts for the container — DeepVariant needs everything visible inside.
INPUT_DIR="$(dirname "${BAM}")"
REF_DIR="$(dirname "${REF}")"
BED_DIR="$(dirname "${CAPTURE_BED}")"

case "${CONTAINER_RUNTIME}" in
  singularity)
    RUN="singularity exec --bind ${REF_DIR}:${REF_DIR},${INPUT_DIR}:${INPUT_DIR},${BED_DIR}:${BED_DIR},${VCF_DIR}:${VCF_DIR} ${DEEPVARIANT_IMG}"
    ;;
  docker)
    RUN="docker run --rm \
      -v ${REF_DIR}:${REF_DIR} \
      -v ${INPUT_DIR}:${INPUT_DIR} \
      -v ${BED_DIR}:${BED_DIR} \
      -v ${VCF_DIR}:${VCF_DIR} \
      ${DEEPVARIANT_IMG#docker://}"
    ;;
  *)
    echo "[ERROR] CONTAINER_RUNTIME must be 'singularity' or 'docker' for DeepVariant"
    exit 1
    ;;
esac

echo "[INFO] DeepVariant (WES model) -> ${VCF}"
${RUN} /opt/deepvariant/bin/run_deepvariant \
  --model_type=WES \
  --ref="${REF}" \
  --reads="${BAM}" \
  --regions="${CAPTURE_BED}" \
  --output_vcf="${VCF}" \
  --output_gvcf="${GVCF}" \
  --num_shards="${THREADS_CALL}" \
  --intermediate_results_dir="${OUTDIR}/tmp/dv"

tabix -p vcf -f "${VCF}"
echo "[DONE] DeepVariant -> ${VCF}"
