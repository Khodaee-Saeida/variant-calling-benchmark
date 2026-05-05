# Data sources and tool references

## Reference genome

- **GRCh38 (no-alt, masked)** — used as recommended by GIAB.
  GIAB FTP: `https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/references/GRCh38/`

## Truth set

- **GIAB v4.2.1 small-variant benchmark** for HG001 (NA12878):
  `https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/`
  - Truth VCF: `HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`
  - High-confidence BED: `HG001_GRCh38_1_22_v4.2.1_benchmark.bed`
- **GIAB stratification BEDs** (optional): `https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/`

## Capture intervals

- **Agilent SureSelect Human All Exon V5** — BED file from Agilent
  (require login); convert chromosome naming (`chr` prefix) and lift to
  GRCh38 if your copy is in GRCh37.

## Known-sites VCFs (BQSR)

Mirrored on GATK's resource bundle:
`https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0`

- `dbsnp_146.hg38.vcf.gz`
- `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz`

## Tool citations

| Tool         | Citation |
|--------------|----------|
| BWA-MEM2     | Vasimuddin et al., *IPDPS* 2019 |
| GATK         | Poplin et al., *bioRxiv* 2018; DePristo et al., *Nat Genet* 2011 |
| DeepVariant  | Poplin et al., *Nat Biotechnol* 2018 |
| Strelka2     | Kim et al., *Nat Methods* 2018 |
| bcftools / samtools | Danecek et al., *GigaScience* 2021 |
| FreeBayes    | Garrison & Marth, *arXiv* 2012 |
| hap.py       | Krusche et al., *Nat Biotechnol* 2019 (best-practices benchmarking) |
| RTG vcfeval  | Cleary et al., *bioRxiv* 2015 |
| GIAB         | Zook et al., *Sci Data* 2016; *Nat Biotechnol* 2019 |
