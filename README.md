# Germline Variant Calling Benchmark on GIAB HG001 (Exome)

The pipeline for short-read germline variant calling and benchmarking against the Genome in a Bottle (GIAB) HG001 truth set, with a side-by-side comparison of four widely used callers: **GATK HaplotypeCaller**, **Strelka2**, **bcftools**, and **FreeBayes**.

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

![Benchmark comparison](docs/benchmark_comparison.png)

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

```
   FASTQ (HG001 Garvan exome)
            │
            ▼
   bwa-mem2 alignment to GRCh38_GIAB_noalt_masked
            │
            ▼
   Picard MarkDuplicates
            │
            ▼
   Picard ValidateSamFile  ─── BAM sanity check
            │
            ▼
   ┌────────┴────────────────────────────────┐
   │                                         │
   ▼              ▼              ▼           ▼
 GATK4         Strelka2       bcftools    FreeBayes
HaplotypeCaller germline      mpileup
 + GenotypeGVCFs                + call
 + VariantFiltration
   │              │              │           │
   └────────┬─────┴──────┬───────┴─────┬─────┘
            ▼            ▼             ▼
         bcftools norm (left-align, decompose multi-allelics)
                         │
                         ▼
              bcftools view -f PASS,.
                         │
                         ▼
            hap.py + RTG vcfeval (vs. GIAB v4.2.1)
                         │
                         ▼
                  per-caller summary.csv
                         │
                         ▼
                 all_callers_summary.csv
```

---

## Environments

Conda environments isolate tool dependencies that don't co-exist (Python 2 vs 3, JVM versions, etc.):

| Env name              | Purpose                                              | Python |
|-----------------------|------------------------------------------------------|--------|
| `variant_benchmark`   | bwa-mem2, samtools, GATK4, bcftools, FreeBayes       | 3.10   |
| `picard`              | Picard MarkDuplicates / ValidateSamFile              | (own)  |
| `strelka`             | Strelka2 germline workflow                           | 2.7    |
| `happy`               | hap.py + RTG vcfeval (evaluation)                    | 2.7    |

YAML specs are in `envs/`.

---

## Repository Structure

```
.
├── README.md                       <- this file
├── envs/                           <- conda env YAMLs
├── config/
│   └── config.sh                   <- centralized paths + parameters
├── scripts/
│   ├── 01_align_bwa.sh             <- bwa-mem2 alignment + sort
│   ├── 02_markdup_validate.sh      <- Picard MarkDuplicates + ValidateSamFile
│   ├── 03_prepare_inputs.sh        <- target BED + truth + reference + RTG SDF
│   ├── 04_call_gatk.sh             <- GATK HaplotypeCaller + filtering
│   ├── 06_call_strelka2.sh         <- Strelka2 germline workflow
│   ├── 07_call_bcftools.sh         <- bcftools mpileup + call
│   ├── 08_call_freebayes.sh        <- FreeBayes
│   ├── 09_run_happy_benchmark.sh   <- hap.py benchmarking + summary
│   ├── plot_benchmark.py           <- generates docs/benchmark_comparison.png
│   └── smoke_test_chr22.sh         <- chr22-only smoke test
├── results/
│   ├── all_callers_summary.csv     <- aggregated benchmark
│   └── per_caller/                 <- raw hap.py outputs
└── docs/
    └── benchmark_comparison.png    <- summary figure
```

---

## How to Reproduce

Clone, install conda envs, then run scripts in order:

```bash
git clone https://github.com/Khodaee-Saeida/variant-calling-benchmark.git
cd variant-calling-benchmark

# 1. Create envs (one-time, ~20 min)
for f in envs/*.yml; do conda env create -f "$f"; done

# 2. Edit config/config.sh with your paths
$EDITOR config/config.sh

# 3. Run the pipeline
conda activate variant_benchmark
bash scripts/01_align_bwa.sh
conda activate picard
bash scripts/02_markdup_validate.sh

bash scripts/03_prepare_inputs.sh        # target BED, truth, RTG SDF
bash scripts/04_call_gatk.sh             # ~1 hr exome
conda activate strelka
bash scripts/06_call_strelka2.sh
conda activate variant_benchmark
bash scripts/07_call_bcftools.sh
bash scripts/08_call_freebayes.sh

# 4. Benchmark
bash scripts/09_run_happy_benchmark.sh   # auto-handles env switching

# 5. Plot
python scripts/plot_benchmark.py
```

For a quick sanity check (no full-genome run), `scripts/smoke_test_chr22.sh` runs all four callers on chr22 only (~30 min total wall time).

---

## Notes on Methodology

**Why pre-normalize.** The four callers represent the same indel/MNP differently. `bcftools norm -f $REF -m -any` left-aligns indels and splits multi-allelics, giving consistent representation before benchmarking. `vcfeval` (used inside hap.py) handles a lot of this internally, but pre-normalizing reduces edge cases and makes downstream `bcftools` queries reproducible.

**Why restrict to the evaluable region.** The GIAB v4.2.1 truth set is highly accurate, but only inside its declared "high-confidence" regions. Variants called outside these regions can't be classified TP/FP/FN reliably. Intersecting with the AgilentV5 capture restricts the comparison to ~46.7 Mb where (a) we expect to call variants and (b) the truth is trustworthy.

**Why hap.py + vcfeval.** `hap.py` is the GA4GH-recommended benchmarking tool. Using `--engine=vcfeval` (RTG's haplotype-aware comparison) handles complex variant representations correctly, avoiding the false discrepancies that simple position-based comparison would produce.

**Why default GATK filters.** The GATK Best Practices hard-filter thresholds were originally tuned on WGS data. They are over-aggressive on exome capture data, particularly the `MQ < 40` and `ReadPosRankSum < -8` filters, because bait-edge reads frequently have marginal mapping quality. This benchmark uses the default Best Practices thresholds unchanged for reproducibility; tuning them for exome data closes most of the SNP recall gap.

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **bwa-mem2**: Vasimuddin et al., IPDPS 2019
- **GATK4**: Van der Auwera & O'Connor, *Genomics in the Cloud*, O'Reilly 2020
- **Strelka2**: Kim et al., Nat Methods 2018
- **bcftools**: Danecek et al., GigaScience 2021
- **FreeBayes**: Garrison & Marth, arXiv 2012
- **hap.py / vcfeval**: Krusche et al., Nat Biotech 2019
- **GIAB HG001**: Zook et al., Sci Data 2016

---

## License

MIT
