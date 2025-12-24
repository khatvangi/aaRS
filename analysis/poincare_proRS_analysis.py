"""
PoincaréMSA analysis for ProRS.

Tasks:
- locate the Poincaré CSV (first CSV*.csv under a run directory)
- optionally merge with metadata
- write merged CSV
- produce simple matplotlib scatter plots

CLI example:
python -m analysis.poincare_proRS_analysis \
    --poincare-dir aaRS/poincare_runs/proRS/poincare_knn5_gamma2 \
    --metadata-csv aaRS/metadata/proRS_metadata.csv \
    --id-col proteins_id --x-col pm1 --y-col pm2
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd

from .utils import ensure_dir, find_first_csv

logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(message)s")
logger = logging.getLogger(__name__)


def load_poincare_csv(
    run_dir: Path, id_col: str, x_col: str, y_col: str
) -> pd.DataFrame:
    csv_path = find_first_csv(run_dir, pattern="*.csv")
    logger.info("Loading Poincaré CSV: %s", csv_path)
    df = pd.read_csv(csv_path)
    for col in (id_col, x_col, y_col):
        if col not in df.columns:
            raise KeyError(f"Column '{col}' not found in {csv_path}. Columns: {list(df.columns)}")
    return df[[id_col, x_col, y_col]].copy()


def load_metadata(path: Optional[Path], id_col: str) -> Optional[pd.DataFrame]:
    if path is None:
        return None
    if not path.exists():
        raise FileNotFoundError(f"Metadata CSV not found: {path}")
    meta = pd.read_csv(path)
    if id_col not in meta.columns:
        raise KeyError(f"Metadata missing id column '{id_col}'. Columns: {list(meta.columns)}")
    return meta


def merge_metadata(coords: pd.DataFrame, meta: Optional[pd.DataFrame], id_col: str) -> pd.DataFrame:
    if meta is None:
        logger.info("No metadata provided; returning coordinates only.")
        return coords
    merged = coords.merge(meta, on=id_col, how="left")
    logger.info("Merged coordinates with metadata: %d rows", len(merged))
    return merged


def _plot_base(df: pd.DataFrame, x_col: str, y_col: str, color_col: Optional[str], title: str, out_path: Path):
    fig, ax = plt.subplots(figsize=(6, 6))
    if color_col and color_col in df.columns:
        categories = df[color_col].fillna("NA").unique()
        cmap = plt.get_cmap("tab20")
        for i, cat in enumerate(categories):
            subset = df[df[color_col] == cat]
            ax.scatter(subset[x_col], subset[y_col], label=str(cat), s=15, color=cmap(i % 20))
        ax.legend(loc="best", fontsize=8)
    else:
        ax.scatter(df[x_col], df[y_col], s=15, color="steelblue")
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    logger.info("Wrote plot: %s", out_path)


def plot_by_kingdom(df: pd.DataFrame, x_col: str, y_col: str, out_dir: Path, kingdom_col: str = "kingdom"):
    _plot_base(df, x_col, y_col, kingdom_col, "Poincaré: by kingdom", out_dir / "fig_proRS_poincare_kingdom.png")


def plot_luca_vs_modern(
    df: pd.DataFrame, x_col: str, y_col: str, out_dir: Path, luca_col: str = "is_LUCA"
):
    fig, ax = plt.subplots(figsize=(6, 6))
    luca_mask = df.get(luca_col, pd.Series([False] * len(df))).astype(bool)
    ax.scatter(df.loc[~luca_mask, x_col], df.loc[~luca_mask, y_col], s=20, color="slategray", label="Modern")
    ax.scatter(df.loc[luca_mask, x_col], df.loc[luca_mask, y_col], s=60, color="gold", marker="*", label="LUCA")
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title("Poincaré: LUCA vs modern")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "fig_proRS_poincare_LUCA.png", dpi=300)
    plt.close(fig)
    logger.info("Wrote plot: %s", out_dir / "fig_proRS_poincare_LUCA.png")


def plot_editing_overlay(
    df: pd.DataFrame, x_col: str, y_col: str, out_dir: Path, editing_col: str = "has_editing_domain"
):
    col_values = df.get(editing_col, pd.Series([False] * len(df))).astype(bool)
    colors = {True: "darkorange", False: "royalblue"}
    fig, ax = plt.subplots(figsize=(6, 6))
    for flag in [True, False]:
        mask = col_values == flag
        ax.scatter(df.loc[mask, x_col], df.loc[mask, y_col], s=20, color=colors[flag], label=str(flag))
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title("Poincaré: editing domain overlay")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "fig_proRS_poincare_editing.png", dpi=300)
    plt.close(fig)
    logger.info("Wrote plot: %s", out_dir / "fig_proRS_poincare_editing.png")


def main():
    parser = argparse.ArgumentParser(description="PoincaréMSA ProRS analysis")
    parser.add_argument("--poincare-dir", type=Path, required=True, help="Directory with Poincaré CSV outputs")
    parser.add_argument("--metadata-csv", type=Path, default=None, help="Optional metadata CSV to merge")
    parser.add_argument("--id-col", default="proteins_id", help="Identifier column name in Poincaré CSV and metadata")
    parser.add_argument("--x-col", default="pm1", help="X coordinate column")
    parser.add_argument("--y-col", default="pm2", help="Y coordinate column")
    parser.add_argument("--output-dir", type=Path, default=Path("aaRS/poincare_runs/proRS"), help="Output root dir")
    args = parser.parse_args()

    coords = load_poincare_csv(args.poincare_dir, args.id_col, args.x_col, args.y_col)
    meta = load_metadata(args.metadata_csv, args.id_col) if args.metadata_csv else None
    merged = merge_metadata(coords, meta, args.id_col)

    merged_path = args.output_dir / "merged_proRS_poincare_metadata.csv"
    ensure_dir(merged_path.parent)
    merged.to_csv(merged_path, index=False)
    logger.info("Wrote merged CSV: %s", merged_path)

    plot_by_kingdom(merged, args.x_col, args.y_col, args.output_dir)
    plot_luca_vs_modern(merged, args.x_col, args.y_col, args.output_dir)
    plot_editing_overlay(merged, args.x_col, args.y_col, args.output_dir)


if __name__ == "__main__":
    main()
