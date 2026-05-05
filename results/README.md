# results/

Files in this directory are **generated outputs** of the pipeline; they are
not committed to git (see `.gitignore`). Run the full pipeline, then
`scripts/11_aggregate_results.py`, and the following will appear here:

| File              | Description                                                    |
|-------------------|----------------------------------------------------------------|
| `summary.tsv`     | Tidy long-format table: caller × type × {hap.py, vcfeval} → TP, FP, FN, precision, recall, F1 |
| `pr_curve.png`    | Precision-recall scatter (SNPs and INDELs) across callers, from hap.py |

After your first successful run, commit a copy of `summary.tsv` and
`pr_curve.png` to the repo (override `.gitignore` for those two files) so
visitors to the GitHub page see real numbers without having to rerun the
pipeline.
