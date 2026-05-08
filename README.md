# Germline Variant Calling Benchmark on GIAB HG001 (Exome)

End-to-end pipeline for short-read germline variant calling and benchmarking against the Genome in a Bottle (GIAB) HG001 truth set, with a side-by-side comparison of four widely used callers: **GATK HaplotypeCaller**, **Strelka2**, **bcftools**, and **FreeBayes**.

The goal is to produce a transparent, reproducible benchmark that quantifies the precision/recall tradeoffs between these callers on real exome data.

---

## Results

Benchmark of HG001 against GIAB v4.2.1 truth, restricted to the AgilentV5 × GIAB high-confidence evaluable region (~46.7 Mb).

| Caller     | SNP F1     | SNP Recall | SNP Precision | INDEL F1   | INDEL Recall | INDEL Precision |
|------------|------------|------------|---------------|------------|--------------|-----------------|
| FreeBayes  | **0.9114** | 0.8461     | 0.9878        | 0.7030     | 0.6587       | 0.7536          |
| Strelka2   | 0.9077     | 0.8449     | 0.9805        | 0.7389     | 0.6844       | 0.8028          |
| bcftools   | 0.8999     | 0.8257     | **0.9887**    | 0.6672     | 0.5883       | 0.7706          |
| GATK       | 0.8541     | 0.7631     | 0.9696        | **0.7434** | **0.7264**   | 0.7613          |


**Key findings**

- **FreeBayes and Strelka2** are the strongest SNP callers on this exome dataset (F1 ≈ 0.91), with Strelka2 marginally more conservative.
- **bcftools** has the highest SNP precision (0.989) but pays for it with lower recall — a useful profile when false positives are expensive.
- **GATK** has the best INDEL profile (F1 = 0.743, recall = 0.726), reflecting HaplotypeCaller's local reassembly strength on small insertions/deletions.
- **GATK SNP recall (0.763) is ~8 points below the other callers.** Diagnostic showed the default Best Practices hard filters (designed for WGS) discard 13% of raw exome SNPs. This is a well-documented limitation of the default thresholds on capture data; the report uses the unmodified Best Practices filters for reproducibility.
- **All four callers show much lower performance on INDELs than SNPs** — expected, and consistent across the field. INDELs in tandem repeats and homopolymers remain the hardest variant class.

---

## Data and References

**Sample.** HG001 / NA12878 (Coriell), Illumina HiSeq exome library prepared by Garvan Institute, downloaded from the GIAB FTP:

> `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_indexes/NA12878/sequence.index.NA12878_Illumina_HiSeq_Exome_Garvan_fastq_09252015`

**Reference genome.** GRCh38 with no alt contigs and the GIAB-recommended decoy/PAR mask (`GRCh38_GIAB_noalt_masked.fa`). This is the reference recommended by the GIAB consortium for benchmarking against their truth sets, eliminating common alt-contig and PAR ambiguities.

**Capture target.** Agilent SureSelect Human All Exon V5, design ID **S04380110**, downloaded from Agilent SureDesign in GRCh38 coordinates. The `S04380110_Covered.bed` (per-bait covered region — the correct file for variant calling) was post-processed (chr-prefix normalization, sort, merge) to produce `AgilentV5_GRCh38.bed` (230,719 intervals, ~50.4 Mb).

**Truth set.** GIAB HG001 v4.2.1 benchmark VCF and high-confidence regions BED:

> `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38`

**Evaluable region.** Intersection of the AgilentV5 capture and the GIAB high-confidence regions = 46,689,743 bases (92.6% of the capture). All recall/precision metrics use this as the denominator.

---

## Pipeline Overview

The pipeline runs in six stages, each implemented by one numbered script in `scripts/`.

**Step 1 — Alignment** (`01_align_bwa.sh`)
Per-lane bwa-mem2 alignment of the Garvan HG001 exome FASTQs to `GRCh38_GIAB_noalt_masked`, with read groups built from FASTQ header (FLOWCELL.LANE) and filename (BARCODE). Lanes are sorted, merged, and deduped with Picard MarkDuplicates to produce a single analysis-ready BAM.

**Step 2 — BAM validation** (`02_validate_bam.sh`)
Nine-section QC pass on the BAM: file integrity (`samtools quickcheck`), index, sort order, reference compatibility (`@SQ` vs reference `.fai`), Picard `ValidateSamFile`, read-group sanity (ID/SM/LB/PL/PU populated), `flagstat` thresholds, `samtools stats`, and on-target depth across the AgilentV5 capture. Hard-fails on anything that would silently corrupt downstream calling.

**Step 3 — Prepare inputs** (`03_prepare_inputs.sh`)
One-time setup of reference and target files: post-process the raw Agilent V5 BED (`S04380110_Covered.bed` → `AgilentV5_GRCh38.bed`), download the GIAB v4.2.1 truth VCF + high-confidence BED, build the evaluable region (capture × GIAB-HC), and pre-format the reference as an RTG SDF for `vcfeval`.

**Step 4 — Variant calling**
Four callers run on the same analysis-ready BAM, restricted to the AgilentV5 capture region. Either run them individually or via the wrapper `run_all_callers.sh`, which times each call and writes a summary table.

- **GATK4 HaplotypeCaller** (`04_call_gatk.sh`) — `HaplotypeCaller -ERC GVCF` → `GenotypeGVCFs` → split SNP/INDEL → `VariantFiltration` with Best Practices hard filters → merged final VCF.
- **Strelka2 germline** (`06_call_strelka2.sh`) — `configureStrelkaGermlineWorkflow.py --exome` followed by `runWorkflow.py`, with a bgzipped + tabixed call-regions BED.
- **bcftools** (`07_call_bcftools.sh`) — single piped `bcftools mpileup | call -m | norm | filter` chain, soft-filtering `QUAL<20 || INFO/DP<10` as `LowQual`.
- **FreeBayes** (`08_call_freebayes.sh`) — `freebayes-parallel` with the capture region split into per-interval chunks, then bcftools normalization and a `QUAL<20` `LowQual` filter.

**Step 5 — Prepare VCFs for benchmarking** (`09_prepare_vcfs.sh`)
Per-caller `bcftools norm -m -any --check-ref s` (left-align, decompose multi-allelics, fix REF mismatches), then `bcftools view -f PASS,.`, then sample-name + reference-contig + variant-count sanity checks. Builds the evaluable-region BED if not already present. Produces `${SM}.<caller>.pass.vcf.gz` ready to feed hap.py.

**Step 6 — Benchmark** (`10_run_happy_benchmark.sh`)
`hap.py --engine=vcfeval` per caller against the GIAB v4.2.1 truth set, restricted to the evaluable region. Aggregates the four per-caller `summary.csv` files into a single `all_callers_summary.csv`. 

---

## Environments

Conda environments isolate tool dependencies that don't co-exist (Python 2 vs 3):

| Env name              | Purpose                                              | Python |
|-----------------------|------------------------------------------------------|--------|
| `variant_benchmark`   | bwa-mem2, samtools, GATK4, bcftools, FreeBayes       | 3.10   |
| `picard`              | Picard MarkDuplicates / ValidateSamFile              | (own)  |
| `strelka`             | Strelka2 germline workflow                           | 2.7    |
| `happy`               | hap.py + RTG vcfeval (evaluation)                    | 2.7    |

YAML specs are in `envs/`.

---


## How to Reproduce

### Setup (one-time)

```bash
# Clone the repo
git clone https://github.com/Khodaee-Saeida/variant-calling-benchmark.git
cd variant-calling-benchmark

# Create conda envs (one-time, ~20 min)
for f in envs/*.yml; do conda env create -f "$f"; done

# Create your local config from the template
cp config/config.example.sh config/config.sh
$EDITOR config/config.sh    # fill in REF, PROJECT_DIR, RAW_AGILENT_BED
```

### Run the pipeline

Each script switches to the conda env it needs automatically — no manual `conda activate` between steps.

```bash
# Step 1 — alignment + MarkDuplicates
bash scripts/01_align_bwa.sh

# Step 2 — BAM QC (hard-fails if anything would corrupt downstream calling)
bash scripts/02_validate_bam.sh

# Step 3 — one-time: capture BED, GIAB truth, RTG SDF
bash scripts/03_prepare_inputs.sh

# Step 4 — variant calling (run all four, with timing summary)
bash scripts/run_all_callers.sh

# ...or individually:
#   bash scripts/04_call_gatk.sh
#   bash scripts/06_call_strelka2.sh
#   bash scripts/07_call_bcftools.sh
#   bash scripts/08_call_freebayes.sh

# Step 5 — normalize + PASS-filter all caller VCFs, build evaluable region
bash scripts/09_prepare_vcfs.sh

# Step 6 — hap.py benchmark + aggregated summary CSV
bash scripts/10_run_happy_benchmark.sh

# Step 7 — generate the comparison figure
python scripts/plot_benchmark.py results/all_callers_summary.csv docs/
```

## License

MIT
