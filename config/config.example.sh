#!/usr/bin/env bash
# ==============================================================================
# config.example.sh — template configuration for the variant-calling benchmark
# ------------------------------------------------------------------------------
# Single source of truth for all paths and parameters. Every script in
# scripts/ sources config/config.sh (NOT this file) and reads its variables
# from there.
#
# Setup:
#   1. cp config/config.example.sh config/config.sh
#   2. Edit config/config.sh — fill in REF, PROJECT_DIR, RAW_AGILENT_BED
#   3. config/config.sh is gitignored, so your local paths stay local.
#
# Variables that must be edited (placeholders below):
#   REF, PROJECT_DIR, RAW_AGILENT_BED
# Everything else derives automatically from PROJECT_DIR and SM.
# ==============================================================================

# ---- Sample identity --------------------------------------------------------
# Sample name — used in every output filename and as the expected SM tag in
# the BAM read groups.
SM="HG001"

# ---- Reference genome -------------------------------------------------------
# GRCh38 with no alt contigs and the GIAB-recommended decoy/PAR mask.
# Must have a sibling .fai (samtools faidx) and .dict (Picard CreateSequenceDictionary).
REF="/path/to/GRCh38_GIAB_noalt_masked.fa"

# ---- Project root -----------------------------------------------------------
# All other paths derive from this. Change PROJECT_DIR and the entire pipeline
# relocates — useful for moving between machines.
PROJECT_DIR="/path/to/variant-calling-benchmark"

# ---- Input data -------------------------------------------------------------
# Paired FASTQs to be aligned. File names should look like:
#   <sample>_<barcode>_L00N_R1_001.fastq.gz
#   <sample>_<barcode>_L00N_R2_001.fastq.gz
FASTQ_DIR="${PROJECT_DIR}/fastq/${SM}"

# Raw Agilent SureSelect V5 BED downloaded from the Agilent SureDesign portal.
# 03_prepare_inputs.sh post-processes this into ${CAPTURE_BED}.
RAW_AGILENT_BED="${PROJECT_DIR}/S04380110_Covered.bed"

# ---- Working / output directories -------------------------------------------
# Per-sample working dir (BAMs, intermediates, logs, QC reports go here).
OUTDIR="${PROJECT_DIR}/work/${SM}"

# Where final per-caller VCFs land (output of 04/06/07/08; input to 09/10).
VCF_DIR="${PROJECT_DIR}/vcfs"

# Where hap.py results land.
HAPPY_OUT_DIR="${PROJECT_DIR}/happy_out"

# ---- Capture target ---------------------------------------------------------
# Cleaned-up Agilent V5 BED produced by 03_prepare_inputs.sh:
#   - chr-prefixed
#   - sorted, merged
CAPTURE_BED="${PROJECT_DIR}/AgilentV5_GRCh38.bed"

# ---- GIAB truth set (downloaded by 03_prepare_inputs.sh if missing) --------
TRUTH_VCF="${PROJECT_DIR}/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="${PROJECT_DIR}/HG001_GRCh38_1_22_v4.2.1_benchmark.bed"

# ---- Evaluable region -------------------------------------------------------
# Intersection of CAPTURE_BED and TRUTH_BED — the region in which hap.py
# computes precision/recall. Built by 03_prepare_inputs.sh or 09_prepare_vcfs.sh.
EVAL_BED="${PROJECT_DIR}/AgilentV5_x_GIAB_HC.bed"

# ---- RTG SDF (built by 03_prepare_inputs.sh) -------------------------------
# Pre-formatted reference for vcfeval (the comparison engine inside hap.py).
SDF_DIR="${PROJECT_DIR}/GRCh38_GIAB_noalt_masked.sdf"

# ---- Threads ----------------------------------------------------------------
THREADS_ALIGN=16   # bwa-mem2
THREADS_SORT=8     # samtools sort / index / flagstat / stats / depth
THREADS_CALL=8     # GATK / Strelka2 / bcftools / FreeBayes / hap.py

# ---- Helpers ----------------------------------------------------------------
# Standardized "this file must exist or we exit" check used by every script.
require_file() {
    [[ -f "$1" ]] || { echo "[ERROR] required file missing: $1" >&2; exit 1; }
}
