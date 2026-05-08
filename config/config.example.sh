#!/usr/bin/env bash
# ==============================================================================
# config.example.sh — template configuration for the variant-calling benchmark
# ------------------------------------------------------------------------------
# Setup:
#   1. cp config/config.example.sh config/config.sh
#   2. Edit config/config.sh with your real paths
#   3. config/config.sh is gitignored, so your local paths stay local
# ==============================================================================

# ---- Sample identity --------------------------------------------------------
SM="HG001"

# ---- Reference --------------------------------------------------------------
REF="/path/to/GRCh38_GIAB_noalt_masked.fa"

# ---- Project root (everything else derives from this) ----------------------
PROJECT_DIR="/path/to/variant-calling-benchmark"

# ---- Input / output roots --------------------------------------------------
FASTQ_DIR="${PROJECT_DIR}/fastq/${SM}"
OUTDIR="${PROJECT_DIR}/work/${SM}"
VCF_DIR="${PROJECT_DIR}/vcfs"
HAPPY_OUT_DIR="${PROJECT_DIR}/happy_out"

# ---- Capture target --------------------------------------------------------
CAPTURE_BED="${PROJECT_DIR}/AgilentV5_GRCh38.bed"

# ---- Truth set + evaluation files ------------------------------------------
TRUTH_VCF="${PROJECT_DIR}/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="${PROJECT_DIR}/HG001_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
EVAL_BED="${PROJECT_DIR}/AgilentV5_x_GIAB_HC.bed"
SDF_DIR="${PROJECT_DIR}/GRCh38_GIAB_noalt_masked.sdf"

# ---- Threads ---------------------------------------------------------------
THREADS_ALIGN=16
THREADS_SORT=8
THREADS_CALL=8

# ---- Helpers ---------------------------------------------------------------
require_file() {
    [[ -f "$1" ]] || { echo "[ERROR] required file missing: $1" >&2; exit 1; }
}
