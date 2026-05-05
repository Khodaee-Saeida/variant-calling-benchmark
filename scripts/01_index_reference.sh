#!/usr/bin/env bash
# ==============================================================================
# 01_index_reference.sh
# ------------------------------------------------------------------------------
# Build all indices needed by the downstream pipeline:
#   - bwa-mem2 index  (alignment)
#   - samtools faidx  (random access to the FASTA)
#   - GATK sequence dictionary (.dict)  (required by GATK / Picard tools)
#   - RTG SDF format  (required by RTG vcfeval)
#
# Idempotent: each step is skipped if its output already exists.
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/config.sh"

require_file "${REF}"

# bwa-mem2 index ---------------------------------------------------------------
if [[ ! -f "${REF}.bwt.2bit.64" && ! -f "${REF}.0123" ]]; then
  echo "[INFO] bwa-mem2 index..."
  bwa-mem2 index "${REF}"
else
  echo "[SKIP] bwa-mem2 index already present"
fi

# samtools faidx ---------------------------------------------------------------
if [[ ! -f "${REF}.fai" ]]; then
  echo "[INFO] samtools faidx..."
  samtools faidx "${REF}"
else
  echo "[SKIP] .fai already present"
fi

# GATK sequence dictionary -----------------------------------------------------
DICT="${REF%.*}.dict"
if [[ ! -f "${DICT}" ]]; then
  echo "[INFO] gatk CreateSequenceDictionary..."
  gatk CreateSequenceDictionary -R "${REF}" -O "${DICT}"
else
  echo "[SKIP] sequence dictionary already present"
fi

# RTG SDF for vcfeval ----------------------------------------------------------
SDF="${REF%.*}.sdf"
if [[ ! -d "${SDF}" ]]; then
  echo "[INFO] rtg format -> ${SDF}"
  rtg format -o "${SDF}" "${REF}"
else
  echo "[SKIP] RTG SDF already present"
fi

echo "[DONE] Reference indexing complete"
