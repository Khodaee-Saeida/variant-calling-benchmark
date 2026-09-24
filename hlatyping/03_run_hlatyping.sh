#!/usr/bin/env bash
# Step 4: Run nf-core/hlatyping on real infrastructure with container access
# (Tübingen HPC / de.NBI cloud — NOT a laptop sandbox without registry access).
#
# Before running for real: check `nextflow run nf-core/hlatyping --help` and
# https://nf-co.re/hlatyping/usage against the version you pin below —
# parameters do change between releases, this script should not be trusted
# blindly, verify against the current docs the day you run it.
#
# Usage: ./03_run_hlatyping.sh <samplesheet.csv> <output_dir>

set -euo pipefail
 
SAMPLESHEET="${1:?Usage: $0 <samplesheet.csv> <output_dir>}"
OUTDIR="${2:?Usage: $0 <samplesheet.csv> <output_dir>}"
 
# Pinned below the Nextflow 25.10 config-parser breaking change (see header).
# Override by exporting NXF_VER before running this script if needed.
export NXF_VER="${NXF_VER:-24.10.5}"
export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-$HOME/singularity_cache}"
mkdir -p "$NXF_SINGULARITY_CACHEDIR" "$OUTDIR"
 
echo "Using Nextflow version: $NXF_VER (pinned — see header comment for why)"
 
nextflow run nf-core/hlatyping \
  -r 2.0.0 \
  -profile singularity \
  --input "$SAMPLESHEET" \
  --outdir "$OUTDIR" \
  --genome GRCh37 \
  -with-report "$OUTDIR/execution_report.html" \
  -with-trace "$OUTDIR/execution_trace.txt"
 
echo "Done. OptiType result should be under: $OUTDIR/optitype/"