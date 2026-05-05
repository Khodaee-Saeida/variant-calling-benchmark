# Methodology

## Sample and data

We used **HG001 (NA12878)**, the most extensively characterised reference sample
in the NIST Genome-in-a-Bottle (GIAB) collection. Whole-exome reads were
captured with the **Agilent SureSelect Human All Exon V5** kit, producing
paired-end Illumina FASTQs across multiple lanes / barcodes.

## Reference

GRCh38 with alternate contigs removed and known false-duplication / decoy
regions masked, as recommended by GIAB for benchmarking. Using the un-masked
GRCh38 inflates false-positive counts in regions where reads multi-map across
duplicated contigs; using the masked reference makes per-caller comparisons
meaningful.

## Alignment and pre-processing

1. **Per-lane alignment with BWA-MEM2** using `-K 100000000 -M -Y`. The
   `-K` flag deterministically chunks input so the output BAM is byte-identical
   regardless of thread count, which matters for benchmarking reproducibility.
2. **Read groups** are constructed per FASTQ pair from the FASTQ header
   (flowcell, lane) and filename (barcode):
   - `ID = SM.FLOWCELL.LANE.BARCODE` (unique per lane)
   - `LB = BARCODE` (so duplicates are marked correctly per library)
   - `PU = FLOWCELL.LANE.BARCODE` (Platform Unit, used by BQSR)
3. **Lane BAMs are merged** with `samtools merge`.
4. **Duplicates marked** with **Picard MarkDuplicates** (`VALIDATION_STRINGENCY=LENIENT`,
   `CREATE_INDEX=true`).
5. *(Optional)* **GATK BQSR** using dbSNP and Mills indels as known-sites VCFs.

The dedup BAM (`${SM}.merged.dedup.bam`, or `…dedup.bqsr.bam` when BQSR is run)
is the **single shared input** to all five callers.

## Variant calling

Each caller is run on the **same BAM** and **restricted to the AgilentV5
capture intervals** (with a 100 bp pad where the caller supports it). This
ensures the comparison reflects genuine algorithmic differences rather than
differences in input, region selection, or filtering convention.

| Caller        | Filtering applied                                              |
|---------------|----------------------------------------------------------------|
| GATK HC       | Hard filters per Best-Practices (QD, FS, MQ, MQRankSum, …)     |
| DeepVariant   | Native quality scores; no extra filtering                      |
| Strelka2      | Native PASS / LowQual annotations                              |
| bcftools      | `QUAL<20 || INFO/DP<10` flagged as LowQual                     |
| FreeBayes     | `QUAL<20` flagged as LowQual                                   |

VCFs are normalised (`bcftools norm -f REF -m -any`) so multi-allelic sites
are decomposed and left-aligned identically across callers.

## Evaluation

Two independent benchmarking tools are used:

- **hap.py 0.3.15** — uses `vcfeval` as its match engine; reports Precision /
  Recall / F1 separately for SNP and INDEL.
- **RTG vcfeval 3.12.1** — independent ROC-aware match-pathfinding algorithm.

The **evaluation region** is the intersection of:

- the **AgilentV5 capture BED**, and
- the **GIAB v4.2.1 high-confidence BED** for HG001.

Restricting to this intersection eliminates two confounders:

1. variants outside the capture target (no reads, so trivially missed);
2. variants outside the high-confidence region (truth set unreliable).

All five callers are evaluated on the same evaluation region with the same
truth VCF, so any difference in F1 reflects caller behaviour.

When stratification BEDs from GIAB (low-complexity, segmental duplications,
homopolymer runs, etc.) are supplied via `STRAT_TSV` in the config, hap.py
reports per-stratum metrics — useful for identifying *where* each caller
struggles, not just how it ranks overall.

## Reproducibility

- All tool versions are pinned in `environment.yml`.
- All paths are centralized in `config/config.sh`; no magic constants in
  the scripts.
- Each script is idempotent: it skips work whose output already exists.
- BAMs are byte-deterministic (BWA-MEM2 `-K 100000000`).
- Container images for DeepVariant and Strelka2 are referenced by exact
  tag, not `latest`.

## Limitations

- Single sample (HG001). Extending to HG002–HG007 would let us report
  metric variance across samples.
- Single capture kit (AgilentV5). Other kits (IDT xGen, Twist Core, Agilent
  V7) have different probe density and would re-rank borderline callers.
- WES only. WGS would change the ranking, particularly for INDELs near
  homopolymers.
- Default thresholds. Caller-specific score thresholds were not tuned on
  a held-out set; the precision-recall trade-off can shift with tuned
  cutoffs.
