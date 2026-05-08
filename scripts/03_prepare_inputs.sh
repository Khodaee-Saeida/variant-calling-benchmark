#!/usr/bin/env bash
# ==============================================================================
# 03_prepare_inputs.sh
# ------------------------------------------------------------------------------
# One-time setup of all reference and target files needed before variant
# calling. Idempotent: each step skips if its output already exists.
#
#   A. CAPTURE TARGET BED
#      Post-process the raw Agilent SureSelect V5 (S04380110) BED:
#        - strip headers/tracks
#        - add 'chr' prefix if missing
#        - sort + merge
#      → ${CAPTURE_BED}
#
#   B. GIAB TRUTH SET
#      Download (if missing) the v4.2.1 truth VCF, its index, and the
#      high-confidence regions BED.
#      → ${TRUTH_VCF}, ${TRUTH_BED}
#
#   C. EVALUABLE REGION
#      Intersect capture × GIAB high-confidence — the region where
#      precision/recall is actually computed.
#      → ${EVAL_BED}
#
#   D. RTG SDF
#      Pre-format the reference for RTG vcfeval (one-time, ~5-10 min, ~3 GB).
#      → ${SDF_DIR}
#
# Per-caller VCF normalization, PASS filtering, and sample-name validation
# happen in scripts/09_prepare_vcfs.sh, which runs *after* the callers.
#
# Conda envs: `variant_benchmark` (bcftools/bedtools/wget) and `happy` (rtg)
# ==============================================================================
set -uo pipefail

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

for v in REF RAW_AGILENT_BED CAPTURE_BED TRUTH_VCF TRUTH_BED EVAL_BED SDF_DIR; do
    [[ -n "${!v:-}" ]] || {
        echo "[ERROR] ${v} not set in ${CONFIG}" >&2
        exit 1
    }
done

require_file "${REF}"

# ---- Conda init ------------------------------------------------------------
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

# ---- GIAB download URL -----------------------------------------------------
GIAB_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38"

log() { echo "[$(date '+%F %T')] $*"; }
hr()  { echo "------------------------------------------------------------"; }

# ============================================================================
# A. Build the capture target BED
# ============================================================================
hr; log "STEP A: Capture target BED"; hr

if [[ -f "${CAPTURE_BED}" ]]; then
    log "  ${CAPTURE_BED} exists, skipping"
else
    if [[ ! -f "${RAW_AGILENT_BED}" ]]; then
        log "  ERROR: raw Agilent BED not found at ${RAW_AGILENT_BED}"
        log "         Download S04380110 from the Agilent SureDesign portal first."
        exit 1
    fi
    conda activate variant_benchmark
    grep -v '^@\|^#\|^browser\|^track' "${RAW_AGILENT_BED}" \
      | awk 'BEGIN{OFS="\t"} {if ($1 !~ /^chr/) $1="chr"$1; print}' \
      | sort -k1,1 -k2,2n \
      | bedtools merge -i - \
      > "${CAPTURE_BED}"
    conda deactivate
    log "  Wrote ${CAPTURE_BED}"
fi
n_intervals=$(wc -l < "${CAPTURE_BED}")
n_bases=$(awk '{s+=$3-$2} END{print s}' "${CAPTURE_BED}")
log "  Capture: ${n_intervals} intervals, ${n_bases} bases (~$((n_bases/1000000)) Mb)"

# ============================================================================
# B. GIAB truth set
# ============================================================================
hr; log "STEP B: GIAB truth set (v4.2.1)"; hr

conda activate variant_benchmark
for f in "${TRUTH_VCF}" "${TRUTH_VCF}.tbi" "${TRUTH_BED}"; do
    if [[ -f "${f}" ]]; then
        log "  ${f} exists"
    else
        url="${GIAB_BASE}/$(basename "${f}")"
        log "  Downloading ${url}"
        wget -q -O "${f}" "${url}"
    fi
done
conda deactivate

n_truth=$(zcat "${TRUTH_VCF}" | grep -vc '^#')
n_hc=$(wc -l < "${TRUTH_BED}")
log "  Truth: ${n_truth} variants, ${n_hc} high-confidence intervals"

# ============================================================================
# C. Evaluable region (capture × GIAB high-conf)
# ============================================================================
hr; log "STEP C: Evaluable region (capture × GIAB HC)"; hr

if [[ -f "${EVAL_BED}" ]]; then
    log "  ${EVAL_BED} exists, skipping"
else
    conda activate variant_benchmark
    bedtools intersect -a "${CAPTURE_BED}" -b "${TRUTH_BED}" \
      | sort -k1,1 -k2,2n \
      | bedtools merge -i - \
      > "${EVAL_BED}"
    conda deactivate
    log "  Wrote ${EVAL_BED}"
fi
n_eval=$(awk '{s+=$3-$2} END{print s}' "${EVAL_BED}")
pct=$(awk -v e="${n_eval}" -v c="${n_bases}" 'BEGIN{printf "%.1f", 100*e/c}')
log "  Evaluable: ${n_eval} bases (~$((n_eval/1000000)) Mb, ${pct}% of capture)"

# ============================================================================
# D. RTG SDF (pre-formatted reference for vcfeval)
# ============================================================================
hr; log "STEP D: RTG SDF"; hr

if [[ -d "${SDF_DIR}" ]]; then
    log "  ${SDF_DIR} exists, skipping"
else
    conda activate happy
    rtg format -o "${SDF_DIR}" "${REF}"
    conda deactivate
    log "  Wrote ${SDF_DIR}"
fi

# ============================================================================
hr
log "Setup complete. Next steps:"
echo "  1. Run callers: scripts/04_call_gatk.sh, 06, 07, 08 (or run_all_callers.sh)"
echo "  2. Prepare VCFs for hap.py:    scripts/09_prepare_vcfs.sh"
echo "  3. Run benchmark:              scripts/10_run_happy_benchmark.sh"
hr
echo
echo "Generated/checked artefacts:"
echo "  Capture BED  : ${CAPTURE_BED}"
echo "  Truth VCF    : ${TRUTH_VCF}"
echo "  Truth BED    : ${TRUTH_BED}"
echo "  Evaluable    : ${EVAL_BED}"
echo "  RTG SDF      : ${SDF_DIR}"
