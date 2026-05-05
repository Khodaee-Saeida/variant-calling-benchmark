#!/usr/bin/env bash
# ==============================================================================
# Centralized configuration for the WES variant-calling benchmark pipeline.
# All scripts source this file. Edit paths and parameters here only.
# ==============================================================================

# ---------- Sample / library ----------
export SM="HG001"                  # sample name (must match the sample in GIAB truth)
export PROJECT="WES_GIAB_bench"    # used in some output filenames

# ---------- Reference ----------
# GIAB-recommended GRCh38, alt contigs removed, decoys/false duplications masked.
export REF="/data/references/GRCh38_GIAB_noalt_masked.fa"

# ---------- Input / output roots ----------
export FASTQ_DIR="/data/HG001_fastq"
export OUTDIR="/data/work/${SM}"
export VCF_DIR="${OUTDIR}/vcfs"
export EVAL_DIR="${OUTDIR}/eval"

# ---------- Capture and truth set ----------
# Agilent SureSelect Human All Exon V5 capture intervals (sorted, GRCh38).
export CAPTURE_BED="/data/intervals/AgilentV5_GRCh38.bed"

# GIAB v4.2.1 truth VCF + high-confidence BED for HG001.
export TRUTH_VCF="/data/truth/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
export TRUTH_BED="/data/truth/HG001_GRCh38_1_22_v4.2.1_benchmark.bed"

# Evaluation region = capture ∩ high-confidence (built once by 09_evaluate_happy.sh).
export EVAL_BED="${OUTDIR}/intervals/AgilentV5_x_GIAB_HC.bed"

# Optional stratification BEDs from GIAB (low-complexity, segdup, etc.) — leave empty if unused.
export STRAT_TSV=""   # e.g. /data/stratification/v3.0/GRCh38/GRCh38-all-stratifications.tsv

# ---------- Known-sites VCFs (for BQSR) ----------
export KNOWN_DBSNP="/data/known/dbsnp_146.hg38.vcf.gz"
export KNOWN_INDELS="/data/known/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

# ---------- Resource controls ----------
export THREADS_ALIGN=16
export THREADS_SORT=8
export THREADS_CALL=16

# ---------- Containers (DeepVariant / Strelka2) ----------
# Set CONTAINER_RUNTIME to "singularity" or "docker"; leave blank to use native installs.
export CONTAINER_RUNTIME="singularity"
export DEEPVARIANT_IMG="docker://google/deepvariant:1.6.1"
export STRELKA2_IMG=""   # Strelka2 is python-based, native install recommended

# ---------- Sanity check ----------
# Bail out with a clear message if any required file is missing.
require_file() {
  if [[ ! -f "$1" ]]; then
    echo "[CONFIG ERROR] Required file not found: $1" >&2
    return 1
  fi
}
