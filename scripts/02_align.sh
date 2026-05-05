#!/usr/bin/env bash
# ==============================================================================
# 02_align.sh
# ------------------------------------------------------------------------------
# Per-lane alignment of paired-end FASTQs with BWA-MEM2, then merge across
# lanes, then mark duplicates with Picard.
#
# Read-group construction follows the SAM specification and GATK Best Practices:
#   - One @RG per FASTQ pair (per-lane / per-barcode).
#   - PU = FLOWCELL.LANE.BARCODE  (parsed from FASTQ header + filename).
#   - LB = barcode (so duplicates are marked correctly per library).
#   - ID = SM.PU                  (unique across the whole dataset).
#
# This is critical: GATK BQSR and Picard MarkDuplicates rely on accurate RG
# tags, and incorrect PU/LB silently degrades variant-calling accuracy.
#
# Inputs : FASTQs in ${FASTQ_DIR} matching *_R1_*.{fastq,fq}.gz (or *_1.*).
# Outputs:
#   ${OUTDIR}/lanes/<ID>.bam            per-lane sorted+indexed BAMs
#   ${OUTDIR}/${SM}.merged.bam          merge of lane BAMs
#   ${OUTDIR}/${SM}.merged.dedup.bam    duplicates marked  ← INPUT to callers
#   ${OUTDIR}/${SM}.markdup.metrics.txt
#   ${OUTDIR}/${SM}.dedup.flagstat.txt
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/config.sh"

require_file "${REF}"
[[ -d "${FASTQ_DIR}" ]] || { echo "[ERROR] FASTQ_DIR not found: ${FASTQ_DIR}"; exit 1; }

mkdir -p "${OUTDIR}/lanes" "${OUTDIR}/logs" "${OUTDIR}/tmp"

# Index reference if needed (delegates to script 01)
if [[ ! -f "${REF}.bwt.2bit.64" && ! -f "${REF}.0123" ]]; then
  bash "${SCRIPT_DIR}/01_index_reference.sh"
fi

# ---------- helpers ----------------------------------------------------------

# Filename convention assumed: SAMPLE_BARCODE_S?_L00X_R[12]_001.fastq.gz
get_barcode_from_fname() {
  basename "$1" | awk -F'_' '{print $2}'
}

# Lane number from filename (_L001_ -> 1) — fallback if FASTQ header lacks it.
get_lane_from_fname() {
  local f; f="$(basename "$1")"
  [[ "$f" =~ _L([0-9]{3})_ ]] && echo "$((10#${BASH_REMATCH[1]}))"
}

# Flowcell + lane from Illumina FASTQ header line:
#   @<instr>:<run>:<flowcell>:<lane>:<tile>:<x>:<y>...
get_fc_lane_from_header() {
  local fq="$1"
  local first
  first="$( (zcat -f -- "$fq" || cat "$fq") 2>/dev/null | head -n1 | tr -d '\r' )"
  IFS=':' read -r _ _ FLOWCELL LANE _ <<< "$first"
  [[ -z "${LANE:-}" ]] && LANE="$(get_lane_from_fname "$fq")"
  echo "${FLOWCELL:-NA}:${LANE:-NA}"
}

# PU = FLOWCELL.LANE.BARCODE — guaranteed unique per lane and library.
get_pu_from_fastq() {
  local fq="$1"
  local fc_lane; fc_lane="$(get_fc_lane_from_header "$fq")"
  local fc="${fc_lane%%:*}"
  local ln="${fc_lane##*:}"
  local bc; bc="$(get_barcode_from_fname "$fq")"
  echo "${fc}.${ln}.${bc}"
}

# Derive R2 from R1 path.
pair_r2() { echo "$1" | sed -E 's/_R1_/_R2_/; s/_1\.fq/_2.fq/; s/_1\.fastq/_2.fastq/'; }

# ---------- discover lanes ---------------------------------------------------

mapfile -t R1S < <(find "${FASTQ_DIR}" -type f \
  \( -name "*_R1_*.fastq.gz" -o -name "*_R1_*.fq.gz" \
     -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" \) \
  | sort)

echo "[INFO] Found ${#R1S[@]} R1 file(s)"
[[ ${#R1S[@]} -gt 0 ]] || { echo "[ERROR] No R1 FASTQs found in ${FASTQ_DIR}"; exit 1; }

# ---------- align each lane --------------------------------------------------

LANE_BAMS=()
for R1 in "${R1S[@]}"; do
  R2="$(pair_r2 "$R1")"
  [[ -f "$R2" ]] || { echo "[ERROR] Missing mate for $R1 (expected $R2)"; exit 1; }

  PU="$(get_pu_from_fastq "$R1")"      # H7AP8ADXX.1.TAAGGCGA
  BARCODE="${PU##*.}"
  ID="${SM}.${PU}"                     # HG001.H7AP8ADXX.1.TAAGGCGA
  LB="${BARCODE}"
  RG=$'@RG\tID:'"${ID}"$'\tSM:'"${SM}"$'\tLB:'"${LB}"$'\tPL:ILLUMINA\tPU:'"${PU}"

  out="${OUTDIR}/lanes/${ID}.bam"
  log="${OUTDIR}/logs/${ID}.align.log"
  tmp="${OUTDIR}/tmp/${ID}.tmp"

  if [[ -f "$out" && -f "${out}.bai" ]]; then
    echo "[SKIP] ${ID} (already aligned)"
    LANE_BAMS+=("${out}")
    continue
  fi

  echo "[INFO] Aligning ${ID}"
  bwa-mem2 mem -t "${THREADS_ALIGN}" -K 100000000 -M -Y \
    -R "${RG}" "${REF}" "$R1" "$R2" \
  | samtools sort -@ "${THREADS_SORT}" -m 2G -O BAM -T "${tmp}" -o "${out}" \
    2> "${log}"
  samtools index -@ "${THREADS_SORT}" "${out}"
  LANE_BAMS+=("${out}")
done

# ---------- merge lanes ------------------------------------------------------

MERGED="${OUTDIR}/${SM}.merged.bam"
echo "[INFO] Merging ${#LANE_BAMS[@]} lane BAM(s) -> ${MERGED}"
if [[ ${#LANE_BAMS[@]} -eq 1 ]]; then
  cp -f "${LANE_BAMS[0]}" "${MERGED}"
  cp -f "${LANE_BAMS[0]}.bai" "${MERGED}.bai"
else
  samtools merge -f -@ "${THREADS_SORT}" -O BAM -o "${MERGED}" "${LANE_BAMS[@]}"
  samtools index -@ "${THREADS_SORT}" "${MERGED}"
fi

# ---------- mark duplicates --------------------------------------------------

DEDUP="${OUTDIR}/${SM}.merged.dedup.bam"
METRICS="${OUTDIR}/${SM}.markdup.metrics.txt"
echo "[INFO] Marking duplicates (Picard)"
picard MarkDuplicates \
  I="${MERGED}" O="${DEDUP}" M="${METRICS}" \
  VALIDATION_STRINGENCY=LENIENT \
  CREATE_INDEX=true

# ---------- QC ---------------------------------------------------------------

echo "[INFO] flagstat -> ${OUTDIR}/${SM}.dedup.flagstat.txt"
samtools flagstat -@ "${THREADS_SORT}" "${DEDUP}" \
  | tee "${OUTDIR}/${SM}.dedup.flagstat.txt"

echo "[DONE] Alignment + dedup complete: ${DEDUP}"
