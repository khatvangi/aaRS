"""
Inspect augmented Poincaré metadata to list available IDs.

Reads: poincare_runs/proRS/merged_proRS_poincare_metadata_augmented.csv
Prints columns and unique IDs, and writes summary/proRS_metadata_ids.csv.
"""
from __future__ import annotations

from pathlib import Path
import pandas as pd


def main():
    path = Path("poincare_runs/proRS/merged_proRS_poincare_metadata_augmented.csv")
    if not path.exists():
        raise FileNotFoundError(f"Metadata file not found: {path}")
    df = pd.read_csv(path)

    cols = df.columns.tolist()
    print("Columns:", cols)

    id_cols = ["proteins_id"]
    for extra in ("seq_id", "name"):
        if extra in df.columns:
            id_cols.append(extra)

    uniq = df[id_cols].drop_duplicates().sort_values("proteins_id")
    print(uniq)

    out = Path("summary/proRS_metadata_ids.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    uniq.to_csv(out, index=False)
    print("Wrote", out)


if __name__ == "__main__":
    main()
