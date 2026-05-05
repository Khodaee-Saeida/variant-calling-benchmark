# Ready-to-paste blurbs for CV / cover letter

Pick whichever length fits the context. Replace the GitHub URL after you push.

---

## One-line (CV bullet)

> **WES Variant-Calling Benchmark on GIAB HG001** — Designed and implemented a
> reproducible end-to-end pipeline (BWA-MEM2 → MarkDuplicates → BQSR → 5
> callers: GATK HaplotypeCaller, DeepVariant, Strelka2, bcftools, FreeBayes;
> evaluated against GIAB v4.2.1 with hap.py and RTG vcfeval).
> [github.com/Khodaee-Saeida/variant-calling-benchmark](https://github.com/Khodaee-Saeida/variant-calling-benchmark)

---

## Two-line (CV project section)

> **WES Variant-Calling Benchmark on GIAB HG001** &nbsp;·&nbsp; *Bash, Python, GATK, DeepVariant, Strelka2, bcftools, FreeBayes, hap.py, RTG vcfeval, conda*
>
> Built a reproducible pipeline that aligns Agilent V5 exome reads with
> BWA-MEM2, calls variants with five widely used germline callers, and
> benchmarks each one against the GIAB v4.2.1 truth set inside the
> capture-target ∩ high-confidence region — cross-validated with two
> independent evaluation tools (hap.py and RTG vcfeval).
> [github.com/Khodaee-Saeida/variant-calling-benchmark](https://github.com/Khodaee-Saeida/variant-calling-benchmark)

---

## Paragraph (cover letter)

> To strengthen my hands-on experience with germline variant analysis, I
> built and published an end-to-end WES variant-calling benchmark on the
> NIST GIAB HG001 reference sample
> (github.com/Khodaee-Saeida/variant-calling-benchmark). The pipeline
> performs per-lane BWA-MEM2 alignment with correctly constructed read-group
> tags, merges and de-duplicates lane BAMs, optionally applies GATK BQSR,
> and runs five callers on the same input — GATK HaplotypeCaller,
> DeepVariant, Strelka2, bcftools, and FreeBayes. Each caller is evaluated
> inside the intersection of the AgilentV5 capture BED and the GIAB v4.2.1
> high-confidence BED, scored with both hap.py and RTG vcfeval, and the
> per-caller precision, recall, and F1 are aggregated into a single summary
> table and precision-recall plot. The repository documents methodology,
> data sources, and exact tool versions so the experiment can be reproduced
> on any GIAB sample, capture kit, or reference build.

---

## Skills demonstrated (use as keywords on LinkedIn / CV)

Bash scripting · Python · pandas · matplotlib · conda environment management ·
BWA-MEM2 · samtools · Picard · GATK · DeepVariant · Strelka2 · bcftools ·
FreeBayes · hap.py · RTG vcfeval · GIAB benchmarking · Singularity / Docker ·
Reproducible bioinformatics pipelines · WES analysis · GRCh38 · BED interval
arithmetic · Git / Gi