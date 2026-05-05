#!/usr/bin/env bash
# ==============================================================================
# 03_bqsr.sh
# ------------------------------------------------------------------------------
# Optional: GATK Base Quality Score Recalibration (BQSR).
#
# BQSR uses a known-sites VCF (dbSNP, Mills indels) to learn systematic biases
# in reported base qualities and apply a corrective table. It is part of the
# GATK Best Practices germline pipeline. Skip with -DRY_RUN if you want to
# benchmark "raw" qualities — but most production pipelines apply BQSR.
#
# Inputs : ${OUTDIR}/${SM}.merged.dedup.bam
# Output : ${OUTDIR}/${SM}.merged.dedup.bqsr.bam (+ .bai, +.recal.table)
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/config.sh"

require_file "${REF}"
require_file "${OUTDIR}/${SM}.merged.dedup.bam"
require_file "${KNOWN_DBSNP}"
require_file "${KNOWN_INDELS}"

IN="${OUTDIR}/${SM}.merged.dedup.bam"
TABLE="${OUTDIR}/${SM}.recal.table"
OUT="${OUTDIR}/${SM}.merged.dedup.bqsr.bam"

echo "[INFO] BaseRecalibrator -> ${TABLE}"
gatk BaseRecalibrator \
  -R "${REF}" \
  -I "${IN}" \
  --known-sites "${KNOWN_DBSNP}" \
  --known-sites "${KNOWN_INDELS}" \
  -L "${CAPTURE_BED}" -ip 100 \
  -O "${TABLE}"

echo "[INFO] ApplyBQSR -> ${OUT}"
gatk ApplyBQSR \
  -R "${REF}" \
  -I "${IN}" \
  --bqsr-recal-file "${TABLE}" \
  -O "${OUT}"

echo "[DONE] BQSR complete: ${OUT}"
