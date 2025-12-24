"""
Final PoincaréMSA geometry analysis for ProRS.

Computes summary stats, significance tests, and produces cleaned figures.

CLI (defaults assume repo root = aaRS):
python -m analysis.proRS_final_geometry_analysis \
    --poincare-merged poincare_runs/proRS/merged_proRS_poincare_metadata.csv \
    --output-stats summary/proRS_poincare_stats.json \
    --output-text summary/proRS_poincare_stats.txt \
    --fig-dir figures
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# Configurable column names
SEQ_COL = "seq_id"
X_COL = "x"
Y_COL = "y"
KINGDOM_COL = "kingdom"
LUCA_COL = "is_LUCA"
EDIT_COL = "has_editing_domain"
PROMISCUITY_COL = "promiscuity_score"

logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(message)s")
logger = logging.getLogger(__name__)


def load_data(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Poincaré merged CSV not found: {path}")
    df = pd.read_csv(path)
    # Handle alternative column names if default not present
    col_map = {}
    if X_COL not in df.columns or Y_COL not in df.columns:
        if "pm1" in df.columns and "pm2" in df.columns:
            col_map[X_COL] = "pm1"
            col_map[Y_COL] = "pm2"
        else:
            raise KeyError(f"Missing coordinate columns in {path}. Columns: {list(df.columns)}")
    else:
        col_map[X_COL] = X_COL
        col_map[Y_COL] = Y_COL

    if SEQ_COL not in df.columns:
        if "proteins_id" in df.columns:
            col_map[SEQ_COL] = "proteins_id"
        else:
            col_map[SEQ_COL] = df.columns[0]  # fallback to first column
    else:
        col_map[SEQ_COL] = SEQ_COL

    df = df.rename(columns={v: k for k, v in col_map.items()})
    df["radius"] = np.sqrt(df[X_COL] ** 2 + df[Y_COL] ** 2)
    return df


def summarize_series(series: pd.Series) -> Dict[str, float]:
    series = series.dropna()
    return {
        "n": int(series.count()),
        "radius_mean": float(series.mean()) if len(series) else None,
        "radius_median": float(series.median()) if len(series) else None,
        "radius_std": float(series.std()) if len(series) else None,
        "radius_min": float(series.min()) if len(series) else None,
        "radius_max": float(series.max()) if len(series) else None,
        "radius_q25": float(series.quantile(0.25)) if len(series) else None,
        "radius_q75": float(series.quantile(0.75)) if len(series) else None,
    }


def mannwhitney_safe(a: pd.Series, b: pd.Series) -> Tuple[float, float]:
    a = a.dropna()
    b = b.dropna()
    if len(a) == 0 or len(b) == 0:
        return (np.nan, np.nan)
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return float(u), float(p)


def kruskal_safe(groups: Dict[str, pd.Series]) -> Tuple[float, float]:
    arrays = [g.dropna() for g in groups.values() if len(g.dropna()) > 0]
    if len(arrays) < 2:
        return (np.nan, np.nan)
    h, p = stats.kruskal(*arrays)
    return float(h), float(p)


def correlation_safe(x: pd.Series, y: pd.Series) -> Dict[str, float]:
    x = x.dropna()
    y = y.dropna()
    if len(x) == 0 or len(y) == 0 or len(x) != len(y):
        return {"available": False}
    pearson = stats.pearsonr(x, y)
    spearman = stats.spearmanr(x, y)
    return {
        "available": True,
        "pearson_r": float(pearson.statistic if hasattr(pearson, "statistic") else pearson[0]),
        "pearson_p": float(pearson.pvalue if hasattr(pearson, "pvalue") else pearson[1]),
        "spearman_rho": float(spearman.statistic if hasattr(spearman, "statistic") else spearman[0]),
        "spearman_p": float(spearman.pvalue if hasattr(spearman, "pvalue") else spearman[1]),
    }


def save_json(obj: Dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2))
    logger.info("Wrote stats JSON: %s", path)


def save_text(lines: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(lines)
    logger.info("Wrote stats TXT: %s", path)


def plot_scatter_by_category(
    df: pd.DataFrame,
    fig_dir: Path,
    color_col: str,
    title: str,
    out_name: str,
    special_mask: pd.Series | None = None,
    special_kwargs: Dict | None = None,
):
    fig, ax = plt.subplots(figsize=(6, 6))
    categories = sorted(df[color_col].dropna().unique())
    cmap = plt.get_cmap("tab10")
    for i, cat in enumerate(categories):
        mask = df[color_col] == cat
        ax.scatter(df.loc[mask, X_COL], df.loc[mask, Y_COL], s=18, color=cmap(i % 10), label=str(cat), alpha=0.9)
    if special_mask is not None and special_kwargs is not None:
        ax.scatter(
            df.loc[special_mask, X_COL],
            df.loc[special_mask, Y_COL],
            **special_kwargs,
        )
    ax.set_xlabel("Poincaré x")
    ax.set_ylabel("Poincaré y")
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig_dir.mkdir(parents=True, exist_ok=True)
    out_path = fig_dir / out_name
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    logger.info("Wrote plot: %s", out_path)


def plot_box_radius(df: pd.DataFrame, fig_dir: Path, out_name: str):
    kingdoms = list(sorted(df[KINGDOM_COL].dropna().unique()))
    data = [df.loc[df[KINGDOM_COL] == k, "radius"].dropna() for k in kingdoms]
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.boxplot(data, labels=kingdoms, showfliers=True)
    ax.set_ylabel("Radius")
    ax.set_title("Radius by kingdom")
    fig.tight_layout()
    out_path = fig_dir / out_name
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    logger.info("Wrote plot: %s", out_path)


def plot_radius_vs_promiscuity(df: pd.DataFrame, fig_dir: Path, out_name: str):
    if len(df) < 2:
        logger.warning("Not enough points for radius vs promiscuity plot; skipping.")
        return
    fig, ax = plt.subplots(figsize=(6, 5))
    cmap = plt.get_cmap("tab10")
    kingdoms = sorted(df[KINGDOM_COL].dropna().unique()) if KINGDOM_COL in df.columns else []
    if kingdoms:
        for i, k in enumerate(kingdoms):
            mask = df[KINGDOM_COL] == k
            ax.scatter(df.loc[mask, "radius"], df.loc[mask, PROMISCUITY_COL], s=18, color=cmap(i % 10), label=str(k))
        ax.legend(loc="best", fontsize=8)
    else:
        ax.scatter(df["radius"], df[PROMISCUITY_COL], s=18, color="steelblue")
    ax.set_xlabel("Radius")
    ax.set_ylabel("Promiscuity score")
    pearson = stats.pearsonr(df["radius"], df[PROMISCUITY_COL])
    spearman = stats.spearmanr(df["radius"], df[PROMISCUITY_COL])
    ax.set_title(f"Radius vs Promiscuity (pearson r={pearson.statistic:.2f}, spearman rho={spearman.statistic:.2f})")
    fig.tight_layout()
    out_path = fig_dir / out_name
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    logger.info("Wrote plot: %s", out_path)


def main():
    parser = argparse.ArgumentParser(description="Final Poincaré geometry analysis for ProRS.")
    parser.add_argument(
        "--poincare-merged",
        type=Path,
        default=Path("poincare_runs/proRS/merged_proRS_poincare_metadata.csv"),
    )
    parser.add_argument("--output-stats", type=Path, default=Path("summary/proRS_poincare_stats.json"))
    parser.add_argument("--output-text", type=Path, default=Path("summary/proRS_poincare_stats.txt"))
    parser.add_argument("--fig-dir", type=Path, default=Path("figures"))
    args = parser.parse_args()

    df = load_data(args.poincare_merged)

    stats_out: Dict = {}

    luca_mask = df[LUCA_COL].astype(bool) if LUCA_COL in df.columns else pd.Series([False] * len(df))
    luca_r = df.loc[luca_mask, "radius"]
    non_luca_r = df.loc[~luca_mask, "radius"]
    stats_out["LUCA"] = summarize_series(luca_r)
    stats_out["non_LUCA"] = summarize_series(non_luca_r)
    u, p = mannwhitney_safe(luca_r, non_luca_r)
    stats_out["luca_vs_non_luca"] = {"mannwhitney_U": u, "p_value": p}

    if KINGDOM_COL in df.columns:
        by_k = {}
        groups = {}
        for k in sorted(df[KINGDOM_COL].dropna().unique()):
            radii = df.loc[df[KINGDOM_COL] == k, "radius"]
            by_k[k] = summarize_series(radii)
            groups[k] = radii
        stats_out["by_kingdom"] = by_k
        h, p_kw = kruskal_safe(groups)
        stats_out["kingdom_kruskal"] = {"H": h, "p_value": p_kw}
    else:
        stats_out["by_kingdom"] = {}
        stats_out["kingdom_kruskal"] = {"H": np.nan, "p_value": np.nan}

    if EDIT_COL in df.columns:
        edit_mask = df[EDIT_COL].astype(bool)
        edit_r = df.loc[edit_mask, "radius"]
        no_edit_r = df.loc[~edit_mask, "radius"]
        stats_out["editing"] = {
            "has_editing": summarize_series(edit_r),
            "no_editing": summarize_series(no_edit_r),
        }
        u_edit, p_edit = mannwhitney_safe(edit_r, no_edit_r)
        stats_out["editing"]["mannwhitney_U"] = u_edit
        stats_out["editing"]["p_value"] = p_edit
    else:
        stats_out["editing"] = {"available": False}

    if PROMISCUITY_COL in df.columns:
        corr = correlation_safe(df["radius"], df[PROMISCUITY_COL])
        stats_out["radius_promiscuity_correlation"] = corr
    else:
        stats_out["radius_promiscuity_correlation"] = {"available": False}

    # LUCA neighborhood density (delta=0.1)
    if len(luca_r.dropna()) > 0:
        delta = 0.1
        threshold = float(luca_r.max()) + delta
        frac = float((df["radius"] <= threshold).mean())
        stats_out["luca_neighborhood"] = {"delta": delta, "fraction_within_delta": frac}

    save_json(stats_out, args.output_stats)

    lines = []
    lines.append(f"LUCA radius summary: {stats_out['LUCA']}")
    lines.append(f"Non-LUCA radius summary: {stats_out['non_LUCA']}")
    lines.append(f"Mann–Whitney LUCA vs non-LUCA: U={stats_out['luca_vs_non_luca']['mannwhitney_U']}, p={stats_out['luca_vs_non_luca']['p_value']}")
    if stats_out.get("by_kingdom"):
        lines.append(f"Kingdom summary: {stats_out['by_kingdom']}")
        lines.append(f"Kruskal–Wallis (kingdom): H={stats_out['kingdom_kruskal']['H']}, p={stats_out['kingdom_kruskal']['p_value']}")
    if isinstance(stats_out.get("editing"), dict) and stats_out["editing"].get("available", True):
        lines.append(f"Editing summary: {stats_out['editing']}")
    if stats_out["radius_promiscuity_correlation"].get("available"):
        lines.append(f"Radius vs promiscuity: {stats_out['radius_promiscuity_correlation']}")
    if stats_out.get("luca_neighborhood"):
        lines.append(f"LUCA neighborhood fraction (delta={stats_out['luca_neighborhood']['delta']}): {stats_out['luca_neighborhood']['fraction_within_delta']}")
    save_text("\n".join(lines), args.output_text)

    fig_dir = args.fig_dir

    # Plots
    if KINGDOM_COL in df.columns:
        plot_scatter_by_category(
            df,
            fig_dir,
            KINGDOM_COL,
            "Poincaré map by kingdom",
            "proRS_poincare_kingdom_clean.png",
            special_mask=luca_mask,
            special_kwargs={"s": 60, "color": "gold", "marker": "*", "edgecolors": "black", "label": "LUCA"},
        )
    else:
        plot_scatter_by_category(
            df.assign(dummy="all"),
            fig_dir,
            "dummy",
            "Poincaré map",
            "proRS_poincare_kingdom_clean.png",
            special_mask=luca_mask,
            special_kwargs={"s": 60, "color": "gold", "marker": "*", "edgecolors": "black", "label": "LUCA"},
        )

    # LUCA vs modern
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(df.loc[~luca_mask, X_COL], df.loc[~luca_mask, Y_COL], s=18, color="slategray", label="Modern")
    ax.scatter(df.loc[luca_mask, X_COL], df.loc[luca_mask, Y_COL], s=70, color="gold", marker="*", edgecolors="black", label="LUCA")
    ax.set_xlabel("Poincaré x")
    ax.set_ylabel("Poincaré y")
    ax.set_title("Poincaré: LUCA vs modern")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    out_path = fig_dir / "proRS_poincare_LUCA_vs_modern.png"
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    logger.info("Wrote plot: %s", out_path)

    # Editing overlay
    if EDIT_COL in df.columns:
        edit_mask = df[EDIT_COL].astype(bool)
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.scatter(df.loc[~edit_mask, X_COL], df.loc[~edit_mask, Y_COL], s=18, color="royalblue", label="No editing")
        ax.scatter(df.loc[edit_mask, X_COL], df.loc[edit_mask, Y_COL], s=18, color="darkorange", label="Has editing")
        ax.scatter(df.loc[luca_mask, X_COL], df.loc[luca_mask, Y_COL], s=70, color="gold", marker="*", edgecolors="black", label="LUCA")
        ax.set_xlabel("Poincaré x")
        ax.set_ylabel("Poincaré y")
        ax.set_title("Poincaré: editing overlay")
        ax.set_aspect("equal", adjustable="box")
        ax.legend(loc="best", fontsize=8)
        fig.tight_layout()
        out_path = fig_dir / "proRS_poincare_editing_overlay.png"
        fig.savefig(out_path, dpi=300)
        plt.close(fig)
        logger.info("Wrote plot: %s", out_path)

    # Radius by kingdom boxplot
    if KINGDOM_COL in df.columns:
        plot_box_radius(df, fig_dir, "proRS_radius_by_kingdom_boxplot.png")

    # Radius vs promiscuity
    if PROMISCUITY_COL in df.columns:
        plot_radius_vs_promiscuity(df.dropna(subset=["radius", PROMISCUITY_COL]), fig_dir, "proRS_radius_vs_promiscuity.png")


if __name__ == "__main__":
    main()
