#!/usr/bin/env bash
# Step 3: Convert the MHC-subset BAM to FASTQ and write the samplesheet
# nf-core/hlatyping expects.
#
# Usage: ./02_bam_to_fastq_and_samplesheet.sh <subset_bam> <output_dir>

set -euo pipefail

SUBSET_BAM="${1:?Usage: $0 <subset_bam> <output_dir>}"
OUTDIR="${2:?Usage: $0 <subset_bam> <output_dir>}"
mkdir -p "$OUTDIR"

# samtools fastq needs name-sorted (or collated) input for correct pairing
echo "Name-sorting..."
samtools sort -n -@ 2 "$SUBSET_BAM" -o "$OUTDIR/na12878_mhc.namesorted.bam"

echo "Converting to paired FASTQ..."
samtools fastq \
  -1 "$OUTDIR/NA12878_mhc_R1.fastq.gz" \
  -2 "$OUTDIR/NA12878_mhc_R2.fastq.gz" \
  -0 /dev/null -s /dev/null \
  -n "$OUTDIR/na12878_mhc.namesorted.bam"

echo "Read counts:"
echo "R1: $(($(zcat "$OUTDIR/NA12878_mhc_R1.fastq.gz" | wc -l) / 4))"
echo "R2: $(($(zcat "$OUTDIR/NA12878_mhc_R2.fastq.gz" | wc -l) / 4))"

# nf-core/hlatyping's current samplesheet format (verified against its usage
# docs): sample,fastq_1,fastq_2,seq_type
SAMPLESHEET="$OUTDIR/samplesheet.csv"
cat > "$SAMPLESHEET" <<EOF
sample,fastq_1,fastq_2,seq_type
NA12878,$(realpath "$OUTDIR/NA12878_mhc_R1.fastq.gz"),$(realpath "$OUTDIR/NA12878_mhc_R2.fastq.gz"),dna
EOF

echo "Samplesheet written to: $SAMPLESHEET"
cat "$SAMPLESHEET"
