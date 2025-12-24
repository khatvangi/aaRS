"""
Generate combined Poincaré + SWORD2 figure for the proRS_LUCA subset.

Outputs:
- figures/proRS_LUCA_SXX_poincare_sword2_combined.png
- figures/proRS_LUCA_SXX_poincare_sword2_combined.svg

Usage:
PYTHONPATH=. python -m analysis.proRS_LUCA_combined_figure \
  --poincare-merged poincare_runs/proRS_LUCA/merged_proRS_LUCA_poincare_metadata.csv \
  --sword2-summary sword2_results/proRS_sword2_summary.csv \
  --outdir figures
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_poincare(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "pm1" not in df.columns or "pm2" not in df.columns:
        raise ValueError("Poincaré merged CSV must contain pm1 and pm2 columns.")
    df["radius"] = np.sqrt(df["pm1"] ** 2 + df["pm2"] ** 2)
    return df


def load_sword2(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    needed = {"structure_id", "n_domains", "n_PUs", "n_solutions"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"SWORD2 summary missing columns: {missing}")
    return df


def make_colors(df: pd.DataFrame) -> List[str]:
    colors = []
    for _, row in df.iterrows():
        if row.get("is_LUCA", 0) == 1:
            colors.append("#d4a017")  # gold
        elif row.get("has_editing_domain", 0) == 1:
            colors.append("#d62728")  # red
        else:
            colors.append("#1f77b4")  # blue
    return colors


def plot_panel_a(ax, df: pd.DataFrame) -> None:
    ax.set_title("Poincaré disk: LUCA vs modern/extant", fontsize=12, fontweight="bold")
    circ = plt.Circle((0, 0), 1, edgecolor="lightgray", facecolor="none", lw=1.0)
    ax.add_patch(circ)

    colors = make_colors(df)
    markers = []
    for _, row in df.iterrows():
        markers.append("*" if row.get("is_LUCA", 0) == 1 else "o")

    for (pm1, pm2, c, m, sid) in zip(df["pm1"], df["pm2"], colors, markers, df["seq_id"]):
        size = 110 if m == "*" else 70
        edge = "black" if m == "*" else "white"
        ax.scatter(pm1, pm2, c=c, s=size, marker=m, edgecolors=edge, linewidths=0.8, alpha=0.9)
        if sid in {"LUCA_ProRS", "modern_prours_pro"}:
            ax.text(pm1 + 0.02, pm2, sid.replace("_", " "), fontsize=9, weight="bold")

    ax.set_xlabel("pm1")
    ax.set_ylabel("pm2")
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-1.05, 1.05)
    ax.set_aspect("equal", adjustable="box")
    ax.scatter([], [], c="#d4a017", marker="*", s=120, edgecolors="black", linewidths=0.8, label="LUCA")
    ax.scatter([], [], c="#d62728", marker="o", s=70, edgecolors="white", linewidths=0.8, label="Editing-lineage")
    ax.scatter([], [], c="#1f77b4", marker="o", s=70, edgecolors="white", linewidths=0.8, label="Non-editing extant")
    ax.legend(frameon=False, loc="upper right")


def _box_with_jitter(ax, data: Dict[str, List[float]], colors: Dict[str, str], ylim: Tuple[float, float], title: str) -> None:
    labels = list(data.keys())
    values = [data[k] for k in labels]
    bp = ax.boxplot(values, patch_artist=True, labels=labels, widths=0.5, showfliers=False)
    for patch, lab in zip(bp["boxes"], labels):
        patch.set_facecolor(colors.get(lab, "gray"))
        patch.set_alpha(0.7)
    for i, vals in enumerate(values, start=1):
        if len(vals) == 1:
            ax.scatter([i], vals, color=colors.get(labels[i - 1], "gray"), s=70, edgecolors="black", zorder=3)
        else:
            jitter = (np.random.rand(len(vals)) - 0.5) * 0.06
            ax.scatter(np.full(len(vals), i) + jitter, vals, color=colors.get(labels[i - 1], "gray"), alpha=0.6, s=40, zorder=3)
    ax.set_ylim(ylim)
    ax.set_ylabel("Hyperbolic radius")
    ax.set_title(title, fontsize=12, fontweight="bold")


def plot_panel_b(ax, df: pd.DataFrame, radius_ylim: Tuple[float, float]) -> None:
    lu = df[df.get("is_LUCA", 0) == 1]["radius"].tolist()
    non = df[df.get("is_LUCA", 0) != 1]["radius"].tolist()
    data = {"LUCA": lu, "Non-LUCA": non}
    colors = {"LUCA": "#d4a017", "Non-LUCA": "#666666"}
    _box_with_jitter(ax, data, colors, radius_ylim, "LUCA vs Non-LUCA radius\n(p \u2248 0.25)")


def plot_panel_c(ax, df: pd.DataFrame, radius_ylim: Tuple[float, float]) -> None:
    edit = df[df.get("has_editing_domain", 0) == 1]["radius"].tolist()
    non = df[df.get("has_editing_domain", 0) != 1]["radius"].tolist()
    data = {"Editing": edit, "Non-editing": non}
    colors = {"Editing": "#d62728", "Non-editing": "#1f77b4"}
    _box_with_jitter(ax, data, colors, radius_ylim, "Editing vs Non-editing radius\n(p \u2248 0.036)")


def plot_panel_d(ax1, ax2, df_s: pd.DataFrame) -> None:
    order = [
        "deep_editing_pro_model_chainA_A",
        "deep_domain_pro_model_chainA_A",
        "modern_prours_pro_model_chainA_A",
    ]
    labels = {
        "deep_editing_pro_model_chainA_A": "LUCA / editing",
        "deep_domain_pro_model_chainA_A": "deep domain",
        "modern_prours_pro_model_chainA_A": "modern",
    }
    subset = df_s.set_index("structure_id").loc[order].reset_index()
    x = np.arange(len(subset))
    width = 0.35
    ax1.bar(x - width / 2, subset["n_domains"], width=width, color="#4d4d4d", label="Domains")
    ax1.bar(x + width / 2, subset["n_PUs"], width=width, color="#b3b3b3", label="PUs")
    ax1.set_xticks(x)
    ax1.set_xticklabels([labels[i] for i in subset["structure_id"]], rotation=20)
    ax1.set_ylabel("Count")
    ax1.set_title("Domains & Protein Units", fontsize=11, fontweight="bold")
    ax1.legend(frameon=False)

    ax2.bar(x, subset["n_solutions"], color="#7b9acc")
    ax2.set_xticks(x)
    ax2.set_xticklabels([labels[i] for i in subset["structure_id"]], rotation=20)
    ax2.set_ylabel("# Solutions")
    ax2.set_title("Structural ambiguity (SWORD2)", fontsize=11, fontweight="bold")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build combined Poincaré+SWORD2 figure for proRS_LUCA run.")
    parser.add_argument("--poincare-merged", default="poincare_runs/proRS_LUCA/merged_proRS_LUCA_poincare_metadata.csv")
    parser.add_argument("--sword2-summary", default="sword2_results/proRS_sword2_summary.csv")
    parser.add_argument("--outdir", default="figures")
    args = parser.parse_args()

    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

    df_p = load_poincare(Path(args.poincare_merged))
    df_s = load_sword2(Path(args.sword2_summary))

    radius_min = df_p["radius"].min()
    radius_max = df_p["radius"].max()
    pad = 0.05
    radius_ylim = (max(0, radius_min - pad), radius_max + pad)

    fig = plt.figure(figsize=(12, 12))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05])
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0])
    gsD = gs[1, 1].subgridspec(1, 2, wspace=0.35)
    axD1 = fig.add_subplot(gsD[0, 0])
    axD2 = fig.add_subplot(gsD[0, 1])

    plot_panel_a(axA, df_p)
    plot_panel_b(axB, df_p, radius_ylim)
    plot_panel_c(axC, df_p, radius_ylim)
    plot_panel_d(axD1, axD2, df_s)

    fig.suptitle("Supplementary Figure SXX. Hyperbolic and Structural Evolution of ProRS", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    png_path = outdir / "proRS_LUCA_SXX_poincare_sword2_combined.png"
    svg_path = outdir / "proRS_LUCA_SXX_poincare_sword2_combined.svg"
    fig.savefig(png_path, dpi=300)
    fig.savefig(svg_path)
    print(f"Wrote {png_path}")
    print(f"Wrote {svg_path}")


if __name__ == "__main__":
    main()
