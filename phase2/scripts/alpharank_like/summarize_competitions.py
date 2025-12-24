#!/usr/bin/env python3
"""Summarize competition IoU results and generate tables/plots."""
import argparse
import logging
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PROJECT_ROOT = Path("/storage/kiran-stuff/aaRS/phase2")
RESULTS_DIR = PROJECT_ROOT / "results"
METRICS = RESULTS_DIR / "metrics_long.csv"
MANIFEST = RESULTS_DIR / "manifest_runs.csv"


def setup_logger(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def bootstrap_ci(data: np.ndarray, iters: int = 500, alpha: float = 0.05) -> Tuple[float, float]:
    if len(data) == 0:
        return np.nan, np.nan
    samples = []
    for _ in range(iters):
        resample = np.random.choice(data, size=len(data), replace=True)
        samples.append(np.mean(resample))
    lo = np.percentile(samples, 100 * (alpha / 2))
    hi = np.percentile(samples, 100 * (1 - alpha / 2))
    return float(lo), float(hi)


def summarize(metrics: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict] = []
    comp_rows = metrics[metrics["run_type"] == "competition"].copy()
    for (era, enzyme, lig1, lig2), grp in comp_rows.groupby(["era", "enzyme", "ligand_1", "ligand_2"]):
        rep = grp.drop_duplicates(subset=["replicate_id"])
        vals = rep["delta_iou"].dropna().values
        wins = rep["win"].dropna().values
        n = len(vals)
        win_rate = wins.mean() if len(wins) else np.nan
        mean_delta = vals.mean() if len(vals) else np.nan
        ci_lo, ci_hi = bootstrap_ci(vals, iters=1000) if len(vals) else (np.nan, np.nan)
        rows.append({
            "era": era,
            "enzyme": enzyme,
            "cognate": lig1,
            "competitor": lig2,
            "n": n,
            "win_rate": win_rate,
            "mean_delta": mean_delta,
            "delta_ci_lower": ci_lo,
            "delta_ci_upper": ci_hi,
        })
    return pd.DataFrame(rows)


def write_table(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["# TABLE 1 (Auto-generated)", "", "| Era | Enzyme | Cognate | Competitor | n | Win rate | Mean ΔIoU | 95% CI |", "|---|---|---|---|---|---|---|---|"]
    for _, r in df.iterrows():
        lines.append(f"| {r.era} | {r.enzyme} | {r.cognate} | {r.competitor} | {int(r.n)} | {r.win_rate:.2f} | {r.mean_delta:.2f} | [{r.delta_ci_lower:.2f}, {r.delta_ci_upper:.2f}] |")
    path.write_text("\n".join(lines))


def plot_bar(df: pd.DataFrame, path: Path) -> None:
    if df.empty:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 4))
    labels = [f"{row.enzyme}\n{row.cognate} vs {row.competitor}" for _, row in df.iterrows()]
    ax.bar(labels, df["win_rate"])
    ax.set_ylabel("Win rate (ΔIoU>0)")
    ax.set_ylim(0, 1)
    plt.xticks(rotation=45, ha="right")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def write_readme(manifest_path: Path, metrics_path: Path, summary_path: Path, table_path: Path, fig_path: Path) -> None:
    content = [
        "# AF3Score AlphaRank-like results",
        "",
        f"- Manifest: `{manifest_path}`",
        f"- Metrics: `{metrics_path}`",
        f"- Competition summary: `{summary_path}`",
        f"- Table: `{table_path}`",
        f"- Plot: `{fig_path}`",
        "",
        "Regenerate:",
        "```bash",
        "python scripts/alpharank_like/build_manifest.py",
        "python scripts/alpharank_like/extract_af3_metrics.py",
        "python scripts/alpharank_like/pocket_iou.py",
        "python scripts/alpharank_like/summarize_competitions.py",
        "```",
    ]
    (RESULTS_DIR / "README_results.md").write_text("\n".join(content))


def main() -> None:
    ap = argparse.ArgumentParser(description="Summarize competition IoU results.")
    ap.add_argument("-m", "--metrics", default=str(METRICS), help="metrics_long.csv with IoU columns")
    ap.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    args = ap.parse_args()
    setup_logger(args.verbose)
    metrics = pd.read_csv(args.metrics)
    summary = summarize(metrics)
    summary_path = RESULTS_DIR / "competition_summary.csv"
    summary.to_csv(summary_path, index=False)
    table_path = RESULTS_DIR / "TABLE1_AUTO.md"
    plot_path = RESULTS_DIR / "FIG_competition_barplot.png"
    write_table(summary, table_path)
    plot_bar(summary, plot_path)
    write_readme(MANIFEST, METRICS, summary_path, table_path, plot_path)
    logging.info("Wrote summary to %s", summary_path)


if __name__ == "__main__":
    main()
