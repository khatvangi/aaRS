"""
Final SWORD2 analysis for ProRS structures.

CLI:
python -m analysis.proRS_sword2_final_analysis \
    --sword2-summary sword2_results/proRS_sword2_summary.csv \
    --tex-out summary/proRS_sword2_table.tex \
    --fig-dir figures \
    --runs-dir sword2_results/runs
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import List, Optional
import json

import matplotlib.pyplot as plt
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(message)s")
logger = logging.getLogger(__name__)


def abbreviate_structure(name: str) -> str:
    return (
        name.replace("_model_chainA_A", "")
        .replace("deep_", "deep_")
        .replace("modern_", "modern_")
    )


def newest_run_dir(runs_dir: Path, prefix: str) -> Optional[Path]:
    matches = sorted([p for p in runs_dir.glob(f"{prefix}*") if p.is_dir()])
    if not matches:
        return None
    return matches[-1]


def extract_metrics_from_json(summary: dict) -> dict:
    """
    Extract n_domains, n_PUs, n_solutions from a SWORD2 summary JSON that
    contains "Optimal partition" / "Alternative partition N" blocks.
    """
    metrics = {"n_domains": None, "n_PUs": None, "n_solutions": None}

    partitions = {k: v for k, v in summary.items() if "partition" in k.lower()}
    metrics["n_solutions"] = len(partitions) if partitions else None

    if partitions:
        opt = partitions.get("Optimal partition")
        if opt is None:
            opt = partitions[sorted(partitions.keys())[0]]
        metrics["n_domains"] = opt.get("Nb. domains") or opt.get("n_domains")

        domains = opt.get("Domains", {}) if isinstance(opt, dict) else {}
        all_pus = []
        for dom in domains.values():
            pus = dom.get("PUs", {}) if isinstance(dom, dict) else {}
            all_pus.extend(pus.keys())
        metrics["n_PUs"] = len(all_pus)
    else:
        doms = summary.get("domains") or summary.get("Domains")
        pus = summary.get("PUs") or summary.get("Protein Units")
        sols = summary.get("solutions") or summary.get("Solutions")
        if doms and isinstance(doms, list):
            metrics["n_domains"] = len(doms)
        if pus and isinstance(pus, list):
            metrics["n_PUs"] = len(pus)
        if sols and isinstance(sols, list):
            metrics["n_solutions"] = len(sols)
    return metrics


def write_latex_table(df: pd.DataFrame, path: Path) -> None:
    cols = ["structure", "n_domains", "n_PUs", "n_solutions", "cat_PUs", "edit_PUs", "score"]
    latex_lines: List[str] = []
    latex_lines.append("\\begin{tabular}{lrrrrrr}")
    latex_lines.append("\\hline")
    latex_lines.append("Structure & Domains & PUs & Solutions & PU$_{\\text{cat}}$ & PU$_{\\text{edit}}$ & Score \\\\")
    latex_lines.append("\\hline")
    for _, row in df.iterrows():
        latex_lines.append(
            f"{row['structure']} & {row['n_domains']} & {row['n_PUs']} & {row['n_solutions']} & "
            f"{row['cat_PUs']} & {row['edit_PUs']} & {row['score']} \\\\"
        )
    latex_lines.append("\\hline")
    latex_lines.append("\\end{tabular}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(latex_lines))
    logger.info("Wrote LaTeX table: %s", path)


def barplot_grouped(df: pd.DataFrame, fig_dir: Path, out_name: str):
    fig, ax = plt.subplots(figsize=(7, 5))
    x = range(len(df))
    width = 0.35
    ax.bar([i - width / 2 for i in x], df["n_domains"], width=width, label="Domains", color="steelblue")
    ax.bar([i + width / 2 for i in x], df["n_PUs"], width=width, label="Protein Units", color="darkorange")
    ax.set_xticks(list(x))
    ax.set_xticklabels(df["structure"], rotation=20)
    ax.set_ylabel("Count")
    ax.set_title("Domains and Protein Units per structure")
    ax.legend()
    fig.tight_layout()
    fig_dir.mkdir(parents=True, exist_ok=True)
    path = fig_dir / out_name
    fig.savefig(path, dpi=300)
    plt.close(fig)
    logger.info("Wrote plot: %s", path)


def barplot_single(df: pd.DataFrame, col: str, title: str, ylabel: str, fig_dir: Path, out_name: str):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(df["structure"], df[col], color="slategray")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticklabels(df["structure"], rotation=20)
    fig.tight_layout()
    fig_dir.mkdir(parents=True, exist_ok=True)
    path = fig_dir / out_name
    fig.savefig(path, dpi=300)
    plt.close(fig)
    logger.info("Wrote plot: %s", path)


def main():
    parser = argparse.ArgumentParser(description="Final SWORD2 analysis for ProRS")
    parser.add_argument("--sword2-summary", type=Path, default=Path("sword2_results/proRS_sword2_summary.csv"))
    parser.add_argument("--tex-out", type=Path, default=Path("summary/proRS_sword2_table.tex"))
    parser.add_argument("--fig-dir", type=Path, default=Path("figures"))
    parser.add_argument("--runs-dir", type=Path, default=Path("sword2_results/runs"), help="SWORD2 runs directory")
    args = parser.parse_args()

    if not args.sword2_summary.exists():
        raise FileNotFoundError(f"SWORD2 summary CSV not found: {args.sword2_summary}")

    df = pd.read_csv(args.sword2_summary)
    for required in ["structure_id", "n_domains", "n_PUs", "n_solutions"]:
        if required not in df.columns:
            raise KeyError(f"Missing column '{required}' in {args.sword2_summary}")

    # Patch metrics from JSON (always refresh to ensure completeness)
    needs_update = False
    for idx, row in df.iterrows():
        prefix = row["structure_id"]
        run_dir = row.get("run_dir")
        json_path = None

        if isinstance(run_dir, str) and run_dir:
            cand = Path(run_dir) / "SWORD2_summary.json"
            if cand.exists():
                json_path = cand
        if json_path is None:
            latest = newest_run_dir(args.runs_dir, prefix)
            if latest:
                cand = latest / "SWORD2_summary.json"
                if cand.exists():
                    json_path = cand

        if not json_path:
            logger.warning("No SWORD2_summary.json found for %s", prefix)
            continue
        try:
            summary = json.loads(json_path.read_text())
            metrics = extract_metrics_from_json(summary)
            for key in ["n_domains", "n_PUs", "n_solutions"]:
                if metrics.get(key) is not None:
                    df.at[idx, key] = metrics[key]
                    needs_update = True
            logger.info("Patched metrics for %s from %s", prefix, json_path)
        except Exception as exc:
            logger.warning("Failed to parse %s: %s", json_path, exc)

    if needs_update:
        df.to_csv(args.sword2_summary, index=False)
        logger.info("Updated sword2 summary CSV: %s", args.sword2_summary)

    df["structure"] = df["structure_id"].apply(abbreviate_structure)
    df["cat_PUs"] = df.get("n_PUs_overlapping_catalytic", pd.Series([None] * len(df)))
    df["edit_PUs"] = df.get("n_PUs_overlapping_editing", pd.Series([None] * len(df)))
    df["score"] = df.get("best_solution_score", pd.Series([None] * len(df)))

    write_latex_table(df, args.tex_out)

    barplot_grouped(df, args.fig_dir, "proRS_sword2_domains_PUs.png")
    barplot_single(df, "n_solutions", "Number of alternative solutions", "Solutions", args.fig_dir, "proRS_sword2_solutions.png")

    if "n_PUs_overlapping_catalytic" in df.columns or "n_PUs_overlapping_editing" in df.columns:
        fig, ax = plt.subplots(figsize=(7, 5))
        x = range(len(df))
        width = 0.35
        ax.bar([i - width / 2 for i in x], df["cat_PUs"], width=width, label="Cat-PUs", color="seagreen")
        ax.bar([i + width / 2 for i in x], df["edit_PUs"], width=width, label="Edit-PUs", color="indianred")
        ax.set_xticks(list(x))
        ax.set_xticklabels(df["structure"], rotation=20)
        ax.set_ylabel("PUs overlapping region")
        ax.set_title("Protein Units overlapping catalytic / editing regions")
        ax.legend()
        fig.tight_layout()
        path = args.fig_dir / "proRS_sword2_PUs_overlap.png"
        fig.savefig(path, dpi=300)
        plt.close(fig)
        logger.info("Wrote plot: %s", path)


if __name__ == "__main__":
    main()
