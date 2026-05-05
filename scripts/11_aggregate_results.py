#!/usr/bin/env python3
"""
11_aggregate_results.py

Walk the hap.py and vcfeval output trees, build a single tidy summary table
(`results/summary.tsv`), and produce a precision-recall plot
(`results/pr_curve.png`).

Reads:
    ${EVAL_DIR}/happy/<caller>/<caller>.summary.csv     (hap.py per-caller summary)
    ${EVAL_DIR}/vcfeval/<caller>/summary.txt            (RTG vcfeval per-caller summary)

Writes:
    results/summary.tsv
    results/pr_curve.png

Run:
    python scripts/11_aggregate_results.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# --- Resolve EVAL_DIR from the shell config (so this script and the bash ----
# pipeline agree on paths). Fall back to a relative default for offline review.
def resolve_eval_dir() -> Path:
    env = os.environ.get("EVAL_DIR")
    if env:
        return Path(env)
    return REPO_ROOT / "results" / "_eval"  # for local dry-runs


CALLERS = ["gatk", "deepvariant", "strelka2", "bcftools", "freebayes"]


def load_happy(eval_dir: Path) -> pd.DataFrame:
    """hap.py emits one CSV per caller; we keep the PASS rows for SNP and INDEL."""
    rows = []
    for caller in CALLERS:
        csv = eval_dir / "happy" / caller / f"{caller}.summary.csv"
        if not csv.exists():
            print(f"[WARN] missing hap.py output: {csv}", file=sys.stderr)
            continue
        df = pd.read_csv(csv)
        df = df[df["Filter"] == "PASS"]                 # exclude ALL row
        df = df[df["Type"].isin(["SNP", "INDEL"])]
        for _, r in df.iterrows():
            rows.append({
                "caller": caller,
                "tool": "hap.py",
                "type": r["Type"],
                "TP": int(r["TRUTH.TP"]),
                "FN": int(r["TRUTH.FN"]),
                "FP": int(r["QUERY.FP"]),
                "precision": float(r["METRIC.Precision"]),
                "recall":    float(r["METRIC.Recall"]),
                "f1":        float(r["METRIC.F1_Score"]),
            })
    return pd.DataFrame(rows)


def load_vcfeval(eval_dir: Path) -> pd.DataFrame:
    """vcfeval summary.txt is whitespace-aligned; the last row labelled
    'None' is the no-threshold (unfiltered) summary, the first row is per-type
    when available. We just grab the bottom 'None' row as the global summary."""
    rows = []
    for caller in CALLERS:
        f = eval_dir / "vcfeval" / caller / "summary.txt"
        if not f.exists():
            print(f"[WARN] missing vcfeval output: {f}", file=sys.stderr)
            continue
        txt = f.read_text().strip().splitlines()
        # find the data line (last non-header, non-divider line)
        data = [ln for ln in txt if ln and not ln.startswith("Threshold") and not ln.startswith("---")]
        if not data:
            continue
        cols = data[-1].split()
        # Layout: Threshold TP-base TP-call FP FN Precision Sensitivity F-measure
        try:
            tp_call = int(cols[2])
            fp      = int(cols[3])
            fn      = int(cols[4])
            prec    = float(cols[5])
            rec     = float(cols[6])
            f1      = float(cols[7])
        except (IndexError, ValueError):
            print(f"[WARN] could not parse {f}: {cols}", file=sys.stderr)
            continue
        rows.append({
            "caller": caller, "tool": "vcfeval", "type": "BOTH",
            "TP": tp_call, "FP": fp, "FN": fn,
            "precision": prec, "recall": rec, "f1": f1,
        })
    return pd.DataFrame(rows)


def plot_pr(df_happy: pd.DataFrame, out_png: Path) -> None:
    if df_happy.empty:
        print("[WARN] no hap.py results -> skipping plot", file=sys.stderr)
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, vtype in zip(axes, ["SNP", "INDEL"]):
        sub = df_happy[df_happy["type"] == vtype]
        ax.scatter(sub["recall"], sub["precision"], s=120)
        for _, r in sub.iterrows():
            ax.annotate(r["caller"], (r["recall"], r["precision"]),
                        textcoords="offset points", xytext=(6, 6), fontsize=10)
        ax.set_title(vtype)
        ax.set_xlabel("Recall")
        ax.set_ylabel("Precision")
        ax.grid(True, ls=":", alpha=0.4)
        ax.set_xlim(0.85, 1.005)
        ax.set_ylim(0.85, 1.005)
    fig.suptitle("Precision-Recall by Caller (hap.py, GIAB v4.2.1, AgilentV5 ∩ HC)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    print(f"[OK] wrote {out_png}")


def main() -> int:
    eval_dir = resolve_eval_dir()
    print(f"[INFO] EVAL_DIR = {eval_dir}")

    happy   = load_happy(eval_dir)
    vcfeval = load_vcfeval(eval_dir)
    summary = pd.concat([happy, vcfeval], ignore_index=True)

    if summary.empty:
        print("[ERROR] no evaluation results found yet; run 09 and 10 first.",
              file=sys.stderr)
        return 1

    out_tsv = RESULTS_DIR / "summary.tsv"
    summary = summary.sort_values(["tool", "type", "caller"]).reset_index(drop=True)
    summary.to_csv(out_tsv, sep="\t", index=False, float_format="%.4f")
    print(f"[OK] wrote {out_tsv} ({len(summary)} rows)")

    plot_pr(happy, RESULTS_DIR / "pr_curve.png")

    # Pretty print to stdout for quick inspection.
    with pd.option_context("display.max_rows", None, "display.width", 120):
        print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
