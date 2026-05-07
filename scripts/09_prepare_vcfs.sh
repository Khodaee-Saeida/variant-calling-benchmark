#!/usr/bin/env bash
# ==============================================================================
# 09_prepare_vcfs.sh
# ------------------------------------------------------------------------------
# Pre-flight VCF preparation for hap.py benchmarking.
#
# Runs once, after all four callers have finished, and produces:
#   - The evaluable-region BED (capture × GIAB high-confidence)
#   - Normalized, PASS-filtered, indexed VCFs per caller — ready to feed hap.py
#   - A summary table with sample names + variant counts (total, in eval region,
#     by type) so you can sanity-check before benchmarking

# Sanity checks performed on each prepped VCF:
#   1. File is readable (bcftools view -h works)
#   2. Sample name matches expected ${SM}
#   3. ##contig lines are consistent with the reference .fai
#   4. Variant count > 0 (caller didn't silently produce an empty file)
#
# Inputs (all from config/config.sh):
#   ${REF}                              reference FASTA + .fai
#   ${CAPTURE_BED}                      AgilentV5 capture intervals
#   ${TRUTH_BED}                        GIAB high-confidence regions
#   ${VCF_DIR}/${SM}.<caller>.vcf.gz    one per caller
#
# Outputs:
#   ${EVAL_BED}                                       capture × GIAB-HC
#   ${VCF_DIR}/${SM}.<caller>.norm.vcf.gz             per caller
#   ${VCF_DIR}/${SM}.<caller>.pass.vcf.gz (+ .tbi)    per caller — feed hap.py
#
# Conda env: `variant_benchmark` (bcftools + bedtools)
# ==============================================================================
set -uo pipefail

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

for v in SM REF VCF_DIR CAPTURE_BED TRUTH_BED EVAL_BED; do
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
                # shellcheck disable=SC1091
                source "${base}/etc/profile.d/conda.sh"
                break
            fi
        done
    else
        # shellcheck disable=SC1091
        source "$(dirname "$(dirname "${CONDA_EXE}")")/etc/profile.d/conda.sh"
    fi
    conda activate variant_benchmark
fi

# ---- Required inputs -------------------------------------------------------
require_file "${REF}"
require_file "${REF}.fai"
require_file "${CAPTURE_BED}"
require_file "${TRUTH_BED}"

CALLERS=(gatk strelka2 bcftools freebayes)
EXPECTED_SAMPLE="${SM}"

# ---- Pretty status output --------------------------------------------------
GREEN=$'\033[1;32m'; YELLOW=$'\033[1;33m'; RED=$'\033[1;31m'
BOLD=$'\033[1m'; OFF=$'\033[0m'
log() { echo "[$(date '+%F %T')] $*"; }
hr()  { echo "------------------------------------------------------------"; }

FAIL_COUNT=0
WARN_COUNT=0
status_ok()   { printf "  %sOK%s    %s\n"   "${GREEN}"  "${OFF}" "$*"; }
status_warn() { printf "  %sWARN%s  %s\n"   "${YELLOW}" "${OFF}" "$*"; WARN_COUNT=$((WARN_COUNT+1)); }
status_fail() { printf "  %sFAIL%s  %s\n"   "${RED}"    "${OFF}" "$*"; FAIL_COUNT=$((FAIL_COUNT+1)); }

# ============================================================================
# STEP 1. Build the evaluable region BED  (capture × GIAB high-confidence)
# ----------------------------------------------------------------------------
# This is the region in which hap.py will compute precision/recall.  It
# excludes (a) regions outside the exome capture (where we don't expect
# variants) and (b) regions outside GIAB's high-confidence call set (where
# the truth itself is uncertain).
# ============================================================================
hr; log "STEP 1. Build evaluable region BED"; hr

if [[ -f "${EVAL_BED}" ]]; then
    log "  ${EVAL_BED} exists, skipping rebuild"
else
    bedtools intersect -a "${CAPTURE_BED}" -b "${TRUTH_BED}" \
      | sort -k1,1 -k2,2n \
      | bedtools merge -i - \
      > "${EVAL_BED}"
    log "  Wrote ${EVAL_BED}"
fi

n_capture=$(awk '{s+=$3-$2} END{print s}' "${CAPTURE_BED}")
n_eval=$(   awk '{s+=$3-$2} END{print s}' "${EVAL_BED}")
pct=$(awk -v e="${n_eval}" -v c="${n_capture}" 'BEGIN{printf "%.1f", 100*e/c}')
log "  Capture   : ${n_capture} bases (~$((n_capture/1000000)) Mb)"
log "  Evaluable : ${n_eval} bases (~$((n_eval/1000000)) Mb, ${pct}% of capture)"

# ============================================================================
# STEP 2. Per-caller VCF preparation
# ============================================================================
hr; log "STEP 2. Per-caller VCF normalize + PASS filter + checks"; hr

# Reference contig set used for compatibility checking
ref_contigs="$(cut -f1 "${REF}.fai" | sort -u)"

echo
printf '  %-10s %12s %12s %12s %12s   %s\n' \
    "caller" "raw" "normalized" "pass" "in_eval" "sample"
printf '  %-10s %12s %12s %12s %12s   %s\n' \
    "------" "---" "----------" "----" "-------" "------"

for caller in "${CALLERS[@]}"; do
    RAW_VCF="${VCF_DIR}/${SM}.${caller}.vcf.gz"
    NORM_VCF="${VCF_DIR}/${SM}.${caller}.norm.vcf.gz"
    PASS_VCF="${VCF_DIR}/${SM}.${caller}.pass.vcf.gz"

    # ----- Check 1: raw VCF exists and is readable --------------------------
    if [[ ! -f "${RAW_VCF}" ]]; then
        printf '  %-10s %s\n' "${caller}" "MISSING ${RAW_VCF}"
        FAIL_COUNT=$((FAIL_COUNT+1))
        continue
    fi
    if ! bcftools view -h "${RAW_VCF}" >/dev/null 2>&1; then
        printf '  %-10s %s\n' "${caller}" "UNREADABLE ${RAW_VCF}"
        FAIL_COUNT=$((FAIL_COUNT+1))
        continue
    fi

    # ----- Step: normalize --------------------------------------------------
    # -m -any            : split multi-allelic records, atomize MNPs
    # --check-ref s      : silently fix REF mismatches against the fasta
    if [[ ! -f "${NORM_VCF}" ]]; then
        bcftools norm \
            -f "${REF}" \
            -m -any \
            --check-ref s \
            -Oz -o "${NORM_VCF}" \
            "${RAW_VCF}" 2>/dev/null
        bcftools index -tf "${NORM_VCF}" 2>/dev/null
    fi

    # ----- Step: PASS filter ------------------------------------------------
    # PASS,. keeps records where FILTER is "PASS" or empty (callers that
    # don't populate FILTER, like bcftools/FreeBayes by default, would
    # otherwise be dropped).
    if [[ ! -f "${PASS_VCF}" ]]; then
        bcftools view -f PASS,. -Oz -o "${PASS_VCF}" "${NORM_VCF}" 2>/dev/null
        bcftools index -tf "${PASS_VCF}" 2>/dev/null
    fi

    # ----- Check 2: sample name matches expected ----------------------------
    sample_name=$(bcftools query -l "${PASS_VCF}" | head -1)
    if [[ "${sample_name}" != "${EXPECTED_SAMPLE}" ]]; then
        sample_name="${sample_name} (!= ${EXPECTED_SAMPLE})"
    fi

    # ----- Check 3: reference contig compatibility --------------------------
    # Compare ##contig lines in VCF header to reference .fai. We tolerate
    # the VCF having a *subset* of reference contigs (typical for an exome
    # restricted to autosomes + X/Y), but flag any contigs in the VCF that
    # are NOT in the reference.
    vcf_contigs="$(bcftools view -h "${PASS_VCF}" 2>/dev/null \
                    | awk -F'[=,]' '/^##contig=<ID=/{print $3}' | sort -u)"
    extra="$(comm -23 <(echo "${vcf_contigs}") <(echo "${ref_contigs}") | head -3)"
    ref_ok=true
    [[ -n "${extra}" ]] && ref_ok=false

    # ----- Counts -----------------------------------------------------------
    n_raw=$(bcftools view -H "${RAW_VCF}" 2>/dev/null | wc -l)
    n_norm=$(bcftools view -H "${NORM_VCF}" 2>/dev/null | wc -l)
    n_pass=$(bcftools view -H "${PASS_VCF}" 2>/dev/null | wc -l)
    n_eval=$(bcftools view -R "${EVAL_BED}" -H "${PASS_VCF}" 2>/dev/null | wc -l)

    printf '  %-10s %12d %12d %12d %12d   %s\n' \
        "${caller}" "${n_raw}" "${n_norm}" "${n_pass}" "${n_eval}" "${sample_name}"

    # Per-caller diagnostic notes (suppressed in the table; printed if not OK)
    if [[ "${n_pass}" -eq 0 ]]; then
        status_fail "${caller}: PASS-filtered file has 0 variants"
    fi
    if [[ "${ref_ok}" == false ]]; then
        status_warn "${caller}: VCF has contigs absent from reference: $(echo "${extra}" | tr '\n' ' ')"
    fi
    if [[ "${sample_name}" == *"!= "* ]]; then
        status_warn "${caller}: sample name mismatch — ${sample_name}"
    fi
done

# ============================================================================
# STEP 3. SNP / INDEL breakdown inside the evaluable region
# ----------------------------------------------------------------------------
# Quick sanity check: counts should be in roughly the same order of magnitude
# across the four callers. Wildly different numbers usually mean a caller
# silently failed or produced malformed output.
# ============================================================================
hr; log "STEP 3. SNP / INDEL counts in the evaluable region"; hr

echo
printf '  %-10s %12s %12s %12s\n' "caller" "all_eval" "snps_eval" "indels_eval"
printf '  %-10s %12s %12s %12s\n' "------" "--------" "---------" "-----------"

for caller in "${CALLERS[@]}"; do
    PASS_VCF="${VCF_DIR}/${SM}.${caller}.pass.vcf.gz"
    [[ -f "${PASS_VCF}" ]] || continue
    n_all=$(bcftools view -R "${EVAL_BED}" -H "${PASS_VCF}" 2>/dev/null | wc -l)
    n_snp=$(bcftools view -R "${EVAL_BED}" -v snps   -H "${PASS_VCF}" 2>/dev/null | wc -l)
    n_ind=$(bcftools view -R "${EVAL_BED}" -v indels -H "${PASS_VCF}" 2>/dev/null | wc -l)
    printf '  %-10s %12d %12d %12d\n' "${caller}" "${n_all}" "${n_snp}" "${n_ind}"
done

# ============================================================================
# Verdict
# ============================================================================
echo
hr
if [[ "${FAIL_COUNT}" -gt 0 ]]; then
    printf "  %sFAIL%s — %d hard failure(s); fix before running hap.py.\n" \
        "${RED}" "${OFF}" "${FAIL_COUNT}"
    exit 1
elif [[ "${WARN_COUNT}" -gt 0 ]]; then
    printf "  %sPASS WITH WARNINGS%s — %d warning(s); review above before benchmarking.\n" \
        "${YELLOW}" "${OFF}" "${WARN_COUNT}"
    exit 0
else
    printf "  %sPASS%s — all VCFs ready for hap.py benchmarking.\n" \
        "${GREEN}" "${OFF}"
    echo
    echo "  Ready to run: bash scripts/10_run_happy_benchmark.sh"
    exit 0
fi
