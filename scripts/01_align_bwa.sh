#!/usr/bin/env bash
# ==============================================================================
# 01_align_bwa.sh
# ------------------------------------------------------------------------------
# Per-lane alignment of HG001 exome FASTQs and produce an analysis-ready BAM.
#
#   FASTQ pairs (R1/R2)  →  bwa-mem2 mem  →  samtools sort   (per-lane BAM)
#                        →  samtools merge                   (one merged BAM)
#                        →  Picard MarkDuplicates            (deduped BAM)
#                        →  samtools flagstat                (QC summary)
#
# Read groups are built from the FASTQ header (FLOWCELL, LANE) and filename
# (BARCODE), giving a unique PU = FLOWCELL.LANE.BARCODE per lane. This is the
# minimum identifier downstream tools (Picard, GATK BQSR) need to recognize
# distinct sequencing batches.
#
# Inputs:
#   - Paired FASTQs in ${FASTQ_DIR}, named like *_R1_*.fastq.gz / *_R2_*.fastq.gz
#   - Reference FASTA at ${REF} (auto-indexed on first run)
#
# Outputs (under ${OUTDIR}):
#   lanes/${SM}.${PU}.bam            one BAM per lane (sorted + indexed)
#   ${SM}.merged.bam                 all lanes merged
#   ${SM}.merged.dedup.bam           analysis-ready BAM (post-MarkDuplicates)
#   ${SM}.markdup.metrics.txt        Picard duplication metrics
#   ${SM}.dedup.flagstat.txt         samtools flagstat summary
#
# Conda envs: `variant_benchmark` (bwa-mem2, samtools) and `picard`
# ==============================================================================
set -euo pipefail

# ---- Configuration ---------------------------------------------------------
# All paths and parameters come from config/config.sh, which is gitignored.
# To set up:
#   cp config/config.example.sh config/config.sh
#   $EDITOR config/config.sh
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

for v in SM REF FASTQ_DIR OUTDIR THREADS_ALIGN THREADS_SORT; do
    [[ -n "${!v:-}" ]] || {
        echo "[ERROR] ${v} not set in ${CONFIG}" >&2
        exit 1
    }
done

# ----------------------------------------------------------------------------

mkdir -p "${OUTDIR}"/{lanes,logs,tmp}

log() { echo "[$(date '+%F %T')] $*"; }

# Index the reference once if the bwa-mem2 index files are missing
if [[ ! -f "${REF}.bwt.2bit.64" && ! -f "${REF}.0123" ]]; then
    log "Indexing reference for bwa-mem2 (one-time, ~20 min)"
    bwa-mem2 index "${REF}"
fi
[[ -f "${REF}.fai" ]] || samtools faidx "${REF}"

# ---- Helpers ---------------------------------------------------------------

# Barcode is the second underscore-separated token in the filename.
# e.g. H7AP8ADXX_TAAGGCGA_L001_R1_001.fastq.gz -> TAAGGCGA
get_barcode_from_fname() {
    basename "$1" | awk -F'_' '{print $2}'
}

# Lane number parsed from filename (_L001_ -> 1); used as a fallback when
# the FASTQ header lacks a lane field.
get_lane_from_fname() {
    local f; f="$(basename "$1")"
    [[ "$f" =~ _L([0-9]{3})_ ]] && echo "$((10#${BASH_REMATCH[1]}))"
}

# Read FLOWCELL and LANE from the first FASTQ record header.
# Illumina format: @<instrument>:<run>:<flowcell>:<lane>:...
get_fc_lane_from_header() {
    local fq="$1"
    local first
    first="$( (zcat -f -- "$fq" || cat "$fq") 2>/dev/null | head -n1 | tr -d '\r' )"
    IFS=':' read -r _ _ FLOWCELL LANE _ <<< "$first"
    [[ -z "${LANE}" ]] && LANE="$(get_lane_from_fname "$fq")"
    echo "${FLOWCELL:-NA}:${LANE:-NA}"
}

# Build a unique platform-unit string: PU = FLOWCELL.LANE.BARCODE
get_pu_from_fastq() {
    local fc_lane flowcell lane barcode
    fc_lane="$(get_fc_lane_from_header "$1")"
    flowcell="${fc_lane%%:*}"
    lane="${fc_lane##*:}"
    barcode="$(get_barcode_from_fname "$1")"
    echo "${flowcell}.${lane}.${barcode}"
}

# Pair an R1 file with its R2 counterpart by string substitution.
pair_r2() { echo "$1" | sed -E 's/_R1_/_R2_/'; }

# ---- Discover input lanes --------------------------------------------------

mapfile -t R1S < <(find "${FASTQ_DIR}" -type f \
    \( -name "*_R1_*.fastq.gz" -o -name "*_R1_*.fq.gz" \
       -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" \) \
    | sort)

(( ${#R1S[@]} > 0 )) || { log "ERROR: no R1 FASTQs found in ${FASTQ_DIR}"; exit 1; }
log "Found ${#R1S[@]} lanes to align"

# ---- Align each lane -------------------------------------------------------

LANE_BAMS=()

for R1 in "${R1S[@]}"; do
    R2="$(pair_r2 "$R1")"
    PU="$(get_pu_from_fastq "$R1")"          # e.g. H7AP8ADXX.1.TAAGGCGA
    BARCODE="${PU##*.}"
    ID="${SM}.${PU}"                          # unique read-group ID
    LB="${BARCODE}"                           # one library per barcode

    RG="@RG\tID:${ID}\tSM:${SM}\tLB:${LB}\tPL:ILLUMINA\tPU:${PU}"

    out="${OUTDIR}/lanes/${ID}.bam"
    log_file="${OUTDIR}/logs/${ID}.align.log"
    tmp="${OUTDIR}/tmp/${ID}.tmp"

    log "Aligning ${ID}"
    bwa-mem2 mem -t "${THREADS_ALIGN}" -K 100000000 -M -Y \
            -R "${RG}" "${REF}" "${R1}" "${R2}" \
      | samtools sort -@ "${THREADS_SORT}" -m 2G -O BAM \
            -T "${tmp}" -o "${out}" 2> "${log_file}"
    samtools index -@ "${THREADS_SORT}" "${out}"

    LANE_BAMS+=("${out}")
done

# ---- Merge lanes -----------------------------------------------------------

MERGED="${OUTDIR}/${SM}.merged.bam"
log "Merging ${#LANE_BAMS[@]} lanes -> $(basename "${MERGED}")"
samtools merge -f -@ "${THREADS_SORT}" -O BAM -o "${MERGED}" "${LANE_BAMS[@]}"
samtools index -@ "${THREADS_SORT}" "${MERGED}"

# ---- Mark duplicates -------------------------------------------------------

DEDUP="${OUTDIR}/${SM}.merged.dedup.bam"
METRICS="${OUTDIR}/${SM}.markdup.metrics.txt"

log "Running Picard MarkDuplicates"
picard MarkDuplicates \
    I="${MERGED}" \
    O="${DEDUP}" \
    M="${METRICS}" \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=true

# ---- QC --------------------------------------------------------------------

samtools flagstat -@ "${THREADS_SORT}" "${DEDUP}" \
  | tee "${OUTDIR}/${SM}.dedup.flagstat.txt"

log "Done. Analysis-ready BAM: ${DEDUP}"
