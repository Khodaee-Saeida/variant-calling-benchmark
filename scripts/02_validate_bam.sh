#!/usr/bin/env bash
# ==============================================================================
# 02_validate_bam.sh
# ------------------------------------------------------------------------------
# Validate the merged + dedup BAM before variant calling. Catches the kinds
# of issues that silently corrupt every downstream caller.
#
#   1. File integrity         samtools quickcheck     BAM not truncated
#   2. BAM index              .bai present (auto-create if missing)
#   3. Sort order             @HD shows SO:coordinate
#   4. Reference compatibility  @SQ in BAM matches reference .fai
#   5. Spec compliance        Picard ValidateSamFile  (SAM/BAM spec, EOF, RG)
#   6. Read-group sanity      @RG ID/SM/LB/PL/PU populated
#   7. Alignment stats        samtools flagstat       (% mapped/paired/duplicate)
#   8. Detailed stats         samtools stats          (error rate, insert size)
#   9. Coverage on target     samtools depth on AgilentV5 BED
#
# Usage:
#   bash scripts/02_validate_bam.sh                       # uses default BAM
#   bash scripts/02_validate_bam.sh /path/to/sample.bam   # validate a specific file
#
# Exit codes:
#   0   PASS (possibly with warnings) — BAM is safe to use
#   1   FAIL — at least one hard failure; do NOT run callers
#
# Conda envs: `variant_benchmark` (samtools) and `picard`
# ==============================================================================
set -uo pipefail

# ---- Configuration ---------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${SCRIPT_DIR}/../config/config.sh" ]]; then
    source "${SCRIPT_DIR}/../config/config.sh"
fi
: "${SM:=HG001}"
: "${REF:=/mnt/vdb/WES_PPMI/benchmark/references/GRCh38_GIAB_noalt_masked.fa}"
: "${OUTDIR:=/mnt/vdb/WES_PPMI/benchmark/work/${SM}}"
: "${CAPTURE_BED:=/mnt/vdb/variants_benchmark/AgilentV5_GRCh38.bed}"
: "${THREADS_SORT:=8}"
if ! declare -F require_file >/dev/null; then
    require_file() {
        [[ -f "$1" ]] || { echo "[ERROR] required file missing: $1" >&2; exit 1; }
    }
fi
if ! declare -F require_file >/dev/null; then
    require_file() {
        [[ -f "$1" ]] || { echo "[ERROR] required file missing: $1" >&2; exit 1; }
    }
fi

BAM="${1:-${OUTDIR}/${SM}.merged.dedup.bam}"
[[ -f "${BAM}" ]] || { echo "[ERROR] BAM not found: ${BAM}"; exit 1; }
require_file "${REF}"

QCDIR="${OUTDIR}/qc"
mkdir -p "${QCDIR}"

REPORT="${QCDIR}/${SM}.validate.report.txt"
FULL_VAL="${QCDIR}/${SM}.validate.full.txt"
FLAGSTAT="${QCDIR}/${SM}.flagstat.txt"
STATS="${QCDIR}/${SM}.stats.txt"
DEPTH="${QCDIR}/${SM}.depth.summary.txt"

GREEN=$'\033[1;32m'; YELLOW=$'\033[1;33m'; RED=$'\033[1;31m'
BOLD=$'\033[1m'; OFF=$'\033[0m'

FAIL_COUNT=0
WARN_COUNT=0

status_ok()   { printf "  %sPASS%s  %s\n"  "${GREEN}"  "${OFF}" "$*"; echo "PASS  $*" >> "${REPORT}"; }
status_warn() { printf "  %sWARN%s  %s\n"  "${YELLOW}" "${OFF}" "$*"; echo "WARN  $*" >> "${REPORT}"; WARN_COUNT=$((WARN_COUNT+1)); }
status_fail() { printf "  %sFAIL%s  %s\n"  "${RED}"    "${OFF}" "$*"; echo "FAIL  $*" >> "${REPORT}"; FAIL_COUNT=$((FAIL_COUNT+1)); }

: > "${REPORT}"
{
    echo "BAM validation report"
    echo "  BAM      : ${BAM}"
    echo "  Reference: ${REF}"
    echo "  Date     : $(date)"
    echo "  Host     : $(hostname)"
    echo
} >> "${REPORT}"

echo "${BOLD}== 1. File integrity (samtools quickcheck) ==${OFF}"
if samtools quickcheck "${BAM}"; then
    status_ok "BAM not truncated; EOF block present"
else
    status_fail "samtools quickcheck FAILED — BAM appears truncated or corrupt"
fi
echo

echo "${BOLD}== 2. BAM index ==${OFF}"
if [[ -f "${BAM}.bai" || -f "${BAM%.bam}.bai" ]]; then
    status_ok ".bai index present"
else
    status_warn "No .bai index; creating one now"
    samtools index -@ "${THREADS_SORT}" "${BAM}"
    status_ok "Created index"
fi
echo

echo "${BOLD}== 3. Sort order ==${OFF}"
hd_line="$(samtools view -H "${BAM}" | grep '^@HD' || true)"
if [[ -z "${hd_line}" ]]; then
    status_warn "No @HD line in header; can't confirm sort order"
elif [[ "${hd_line}" == *"SO:coordinate"* ]]; then
    status_ok "BAM is coordinate-sorted"
else
    status_fail "BAM is NOT coordinate-sorted — sort with 'samtools sort' before calling"
    echo "    @HD line: ${hd_line}"
fi
echo

echo "${BOLD}== 4. Reference compatibility (BAM @SQ vs reference .fai) ==${OFF}"
if [[ ! -f "${REF}.fai" ]]; then
    status_warn "No ${REF}.fai found; can't cross-check sequence dictionary"
else
    bam_seqs="$(samtools view -H "${BAM}" \
                  | awk '$1=="@SQ"{
                          for(i=2;i<=NF;i++){
                            if($i~/^SN:/) sn=substr($i,4);
                            if($i~/^LN:/) ln=substr($i,4);
                          }
                          print sn"\t"ln
                        }' | sort)"
    ref_seqs="$(awk '{print $1"\t"$2}' "${REF}.fai" | sort)"
    n_bam=$(echo "${bam_seqs}" | grep -c .)
    n_ref=$(echo "${ref_seqs}" | grep -c .)
    mismatches="$(comm -23 <(echo "${bam_seqs}") <(echo "${ref_seqs}") | head -5)"
    if [[ -z "${mismatches}" ]]; then
        status_ok "All ${n_bam} BAM @SQ entries match the reference (${n_ref} total contigs)"
    else
        status_fail "BAM @SQ does NOT match the reference — likely aligned to a different FASTA"
        echo "    First mismatched entries (BAM contigs absent or different length in .fai):"
        echo "${mismatches}" | sed 's/^/      /'
    fi
fi
echo

echo "${BOLD}== 5. Picard ValidateSamFile ==${OFF}"
if conda run -n picard picard -Xmx16g ValidateSamFile \
        I="${BAM}" R="${REF}" \
        MODE=SUMMARY \
        IGNORE=MATE_NOT_FOUND \
        IGNORE=INVALID_TAG_NM \
        > "${FULL_VAL}" 2>&1; then
    if grep -q "No errors found" "${FULL_VAL}"; then
        status_ok "Picard ValidateSamFile: no errors"
    else
        nerr=$(grep -c "ERROR" "${FULL_VAL}" || true)
        nwarn=$(grep -c "WARNING" "${FULL_VAL}" || true)
        if [[ "${nerr}" -gt 0 ]]; then
            status_fail "Picard reported ${nerr} ERROR types (see ${FULL_VAL})"
        else
            status_warn "Picard reported ${nwarn} WARNING types (see ${FULL_VAL})"
        fi
    fi
else
    status_fail "Picard ValidateSamFile crashed (see ${FULL_VAL})"
fi
echo

echo "${BOLD}== 6. Read-group headers ==${OFF}"
RG_LINES="$(samtools view -H "${BAM}" | grep -c '^@RG' || true)"
if [[ "${RG_LINES}" -eq 0 ]]; then
    status_fail "No @RG lines in header — variant callers will reject this BAM"
else
    status_ok "${RG_LINES} @RG line(s) present"
    first_rg="$(samtools view -H "${BAM}" | grep '^@RG' | head -1)"
    for tag in ID SM LB PL PU; do
        if [[ "${first_rg}" == *$'\t'"${tag}:"* ]]; then
            status_ok "@RG has ${tag}"
        else
            status_fail "@RG missing tag ${tag}"
        fi
    done
fi
echo

echo "${BOLD}== 7. Alignment statistics (samtools flagstat) ==${OFF}"
samtools flagstat -@ "${THREADS_SORT}" "${BAM}" > "${FLAGSTAT}"
cat "${FLAGSTAT}"

mapped_pct=$(awk '/mapped \(/ && !/primary mapped/ {gsub(/[(%]/,""); print $5; exit}' "${FLAGSTAT}")
ppair_pct=$(awk '/properly paired/ {gsub(/[(%]/,""); print $6; exit}' "${FLAGSTAT}")
dup_pct=$(awk '/duplicates/ && !/primary duplicates/ {dup=$1}
               /^[0-9]+ \+ 0 in total/ {tot=$1}
               END {if (tot>0) printf "%.2f", 100*dup/tot}' "${FLAGSTAT}")
echo

[[ -n "${mapped_pct}" ]] && {
    awk -v p="${mapped_pct}" 'BEGIN{exit !(p>=98)}' \
        && status_ok   "Mapped reads: ${mapped_pct}% (>= 98% expected for WES)" \
        || status_warn "Mapped reads: ${mapped_pct}% (low; expect >= 98% for WES)"
}
[[ -n "${ppair_pct}" ]] && {
    awk -v p="${ppair_pct}" 'BEGIN{exit !(p>=90)}' \
        && status_ok   "Properly paired: ${ppair_pct}% (>= 90% expected)" \
        || status_warn "Properly paired: ${ppair_pct}% (low; expect >= 90%)"
}
[[ -n "${dup_pct}" ]] && {
    awk -v p="${dup_pct}" 'BEGIN{exit !(p>=2 && p<=40)}' \
        && status_ok   "Duplicates: ${dup_pct}% (typical 5-30% for AgilentV5)" \
        || status_warn "Duplicates: ${dup_pct}% (outside typical 5-30% range)"
}
echo

echo "${BOLD}== 8. Detailed stats (samtools stats) ==${OFF}"
samtools stats -@ "${THREADS_SORT}" --reference "${REF}" "${BAM}" > "${STATS}"
err_rate=$(awk '/^SN[[:space:]]+error rate:/         {print $4}' "${STATS}")
ins_size=$(awk '/^SN[[:space:]]+insert size average:/ {print $5}' "${STATS}")
read_len=$(awk '/^SN[[:space:]]+average length:/      {print $4}' "${STATS}")

status_ok "Average read length:  ${read_len}"
status_ok "Average insert size:  ${ins_size}"
[[ -n "${err_rate}" ]] && {
    awk -v p="${err_rate}" 'BEGIN{exit !(p<=0.01)}' \
        && status_ok   "Error rate: ${err_rate} (<= 1% expected)" \
        || status_warn "Error rate: ${err_rate} (high; > 1%)"
}
echo

echo "${BOLD}== 9. Coverage on AgilentV5 capture regions ==${OFF}"
if [[ -f "${CAPTURE_BED}" ]]; then
    samtools depth -a -b "${CAPTURE_BED}" "${BAM}" \
      | awk 'BEGIN{n=0; s=0; below10=0}
             {n++; s+=$3; if($3<10) below10++}
             END {if(n>0) printf "mean_depth=%.1f\nfraction_below_10x=%.3f\nbases_in_target=%d\n",
                                 s/n, below10/n, n}' > "${DEPTH}"
    cat "${DEPTH}"
    mean_depth=$(awk -F= '/mean_depth/{print $2}' "${DEPTH}")
    frac_low=$(  awk -F= '/fraction_below_10x/{print $2}' "${DEPTH}")
    [[ -n "${mean_depth}" ]] && {
        awk -v p="${mean_depth}" 'BEGIN{exit !(p>=30)}' \
            && status_ok   "Mean on-target depth: ${mean_depth}x (>= 30x for WES)" \
            || status_warn "Mean on-target depth: ${mean_depth}x (< 30x; calling will be limited)"
    }
    [[ -n "${frac_low}" ]] && {
        awk -v p="${frac_low}" 'BEGIN{exit !(p<=0.20)}' \
            && status_ok   "Fraction below 10x: ${frac_low} (<= 20% expected)" \
            || status_warn "Fraction below 10x: ${frac_low} (> 20%; uneven capture)"
    }
else
    status_warn "CAPTURE_BED not set or not found — skipping coverage check"
fi

echo
echo "${BOLD}== Verdict ==${OFF}"
if [[ "${FAIL_COUNT}" -gt 0 ]]; then
    printf "  %sFAIL%s — %d hard failure(s); fix before running callers.\n  Report: %s\n" \
        "${RED}" "${OFF}" "${FAIL_COUNT}" "${REPORT}"
    exit 1
elif [[ "${WARN_COUNT}" -gt 0 ]]; then
    printf "  %sPASS WITH WARNINGS%s — %d warning(s); BAM is usable but review %s\n" \
        "${YELLOW}" "${OFF}" "${WARN_COUNT}" "${REPORT}"
    exit 0
else
    printf "  %sPASS%s — BAM is ready for variant calling.\n  Report: %s\n" \
        "${GREEN}" "${OFF}" "${REPORT}"
    exit 0
fi
