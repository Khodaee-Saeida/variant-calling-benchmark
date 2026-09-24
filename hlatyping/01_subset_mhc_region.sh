#!/usr/bin/env bash
# Step 1-2: Subset the existing NA12878 exome BAM to the MHC region.
#
# Usage: ./01_subset_mhc_region.sh <input_bam> <output_dir>

set -euo pipefail

INPUT_BAM="${1:?Usage: $0 <input_bam> <output_dir>}"
OUTDIR="${2:?Usage: $0 <input_bam> <output_dir>}"
mkdir -p "$OUTDIR"

# IMPORTANT: check the BAM's chromosome naming convention first.
# Ensembl-style GRCh38 references use "6", UCSC/NCBI-style use "chr6".
# This one line tells you which:
echo "Detecting chromosome naming convention..."
samtools view -H "$INPUT_BAM" | grep -E "^@SQ" | grep -E "SN:(chr)?6" || true

# MHC region, GRCh38 coordinates (~6:28,510,120-33,480,577).
# Edit the region string below to match the naming your BAM actually uses.
REGION="chr6:29000000-33000000"   

echo "Subsetting to MHC region: $REGION"
samtools view -b -h "$INPUT_BAM" "$REGION" > "$OUTDIR/na12878_mhc.bam"
samtools index "$OUTDIR/na12878_mhc.bam"

echo "Read count in subset:"
samtools view -c "$OUTDIR/na12878_mhc.bam"

echo "Done. Output: $OUTDIR/na12878_mhc.bam"
