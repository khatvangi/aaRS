"""
Merge the proRS_LUCA Poincaré output with manual annotations.

Example:
PYTHONPATH=. python -m analysis.build_proRS_LUCA_metadata \
  --poincare-dir poincare_runs/proRS_LUCA/poincare_knn5_gamma2 \
  --annotations poincare_runs/proRS_LUCA/proRS_LUCA_metadata_annotations.csv \
  --output poincare_runs/proRS_LUCA/merged_proRS_LUCA_poincare_metadata.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd


def _find_poincare_csv(pdir: Path) -> Path:
    csvs = sorted(pdir.glob("*.csv"))
    if not csvs:
        raise FileNotFoundError(f"No CSV files found in {pdir}")
    return csvs[0]


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge Poincaré CSV with manual annotations.")
    parser.add_argument("--poincare-dir", required=True, help="Directory containing Poincaré CSV")
    parser.add_argument("--annotations", required=True, help="Annotations CSV (proteins_id column)")
    parser.add_argument("--output", required=True, help="Output path for merged CSV")
    args = parser.parse_args()

    pdir = Path(args.poincare_dir)
    ann_path = Path(args.annotations)
    if not ann_path.exists():
        print(f"Annotation CSV not found: {ann_path}", file=sys.stderr)
        sys.exit(1)

    pcsv = _find_poincare_csv(pdir)
    df_p = pd.read_csv(pcsv)
    ann = pd.read_csv(ann_path)

    if "proteins_id" not in df_p.columns:
        print("Poincaré CSV missing 'proteins_id'", file=sys.stderr)
        sys.exit(1)
    if "proteins_id" not in ann.columns:
        print("Annotation CSV missing 'proteins_id'", file=sys.stderr)
        sys.exit(1)

    merged = df_p.merge(ann, on="proteins_id", how="left")
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_path, index=False)
    print(f"Wrote {out_path} with {len(merged)} rows (source {pcsv.name})")


if __name__ == "__main__":
    main()
