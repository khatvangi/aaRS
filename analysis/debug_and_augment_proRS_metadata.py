"""
Inspect and optionally augment ProRS Poincaré metadata.

Inspection:
  - Prints column names and head
  - Counts for is_LUCA and has_editing_domain if present
  - Saves inspection to summary/proRS_metadata_inspection.txt

Augmentation (when --config-json is provided):
  - Reads a JSON config describing id column, LUCA ids, editing ids, and promiscuity scores
  - Adds/updates is_LUCA, has_editing_domain, promiscuity_score columns
  - Writes augmented CSV to --metadata-out
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List

import pandas as pd


def inspect_metadata(df: pd.DataFrame) -> str:
    lines: List[str] = []
    lines.append("Columns:\n" + ", ".join(df.columns))
    lines.append("\nHead (first 10 rows):")
    lines.append(df.head(10).to_string(index=False))
    if "is_LUCA" in df.columns:
        lines.append(f"\nCount is_LUCA==1: {(df['is_LUCA']==1).sum()}")
    if "has_editing_domain" in df.columns:
        lines.append(f"Count has_editing_domain==1: {(df['has_editing_domain']==1).sum()}")
    lines.append(f"\nTotal rows: {len(df)}")
    return "\n".join(lines)


def apply_config(df: pd.DataFrame, config: Dict) -> pd.DataFrame:
    id_col = config.get("id_column")
    if not id_col:
        raise ValueError("Config JSON must contain 'id_column'.")
    if id_col not in df.columns:
        raise KeyError(f"id_column '{id_col}' not found in metadata. Columns: {list(df.columns)}")

    luca_ids = set(config.get("LUCA_ids", []))
    editing_ids = set(config.get("editing_ids", []))
    prom_dict = config.get("promiscuity", {})

    df = df.copy()
    if luca_ids:
        df["is_LUCA"] = df[id_col].isin(luca_ids).astype(int)
    if editing_ids:
        df["has_editing_domain"] = df[id_col].isin(editing_ids).astype(int)
    if prom_dict:
        existing = df.get("promiscuity_score")
        mapped = df[id_col].map(prom_dict)
        if existing is not None:
            df["promiscuity_score"] = mapped.combine_first(existing)
        else:
            df["promiscuity_score"] = mapped
    return df


def main():
    parser = argparse.ArgumentParser(description="Inspect and optionally augment ProRS metadata.")
    parser.add_argument(
        "--metadata-in",
        type=Path,
        default=Path("poincare_runs/proRS/merged_proRS_poincare_metadata.csv"),
        help="Input metadata CSV",
    )
    parser.add_argument(
        "--metadata-out",
        type=Path,
        default=Path("poincare_runs/proRS/merged_proRS_poincare_metadata_augmented.csv"),
        help="Output augmented metadata CSV",
    )
    parser.add_argument(
        "--config-json",
        type=Path,
        default=None,
        help="Annotation config JSON (if omitted, only inspection is performed)",
    )
    parser.add_argument(
        "--inspection-out",
        type=Path,
        default=Path("summary/proRS_metadata_inspection.txt"),
        help="Path to save inspection log",
    )
    args = parser.parse_args()

    if not args.metadata_in.exists():
        raise FileNotFoundError(f"Metadata input not found: {args.metadata_in}")
    df = pd.read_csv(args.metadata_in)

    inspection = inspect_metadata(df)
    args.inspection_out.parent.mkdir(parents=True, exist_ok=True)
    args.inspection_out.write_text(inspection)
    print(inspection)

    if args.config_json is None:
        return

    if not args.config_json.exists():
        raise FileNotFoundError(f"Config JSON not found: {args.config_json}")
    config = json.loads(args.config_json.read_text())
    df_aug = apply_config(df, config)
    args.metadata_out.parent.mkdir(parents=True, exist_ok=True)
    df_aug.to_csv(args.metadata_out, index=False)
    print(f"\nAugmented metadata written to {args.metadata_out}")


if __name__ == "__main__":
    main()
