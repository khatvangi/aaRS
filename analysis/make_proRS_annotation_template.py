"""
Create a numeric-ID annotation template for ProRS.

Reads summary/proRS_metadata_ids.csv and writes
metadata/proRS_annotation_config_numeric.json.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


def main():
    ids_path = Path("summary/proRS_metadata_ids.csv")
    if not ids_path.exists():
        raise FileNotFoundError(f"Missing ID list: {ids_path}. Run analysis.inspect_proRS_metadata_ids first.")
    df = pd.read_csv(ids_path)
    ids = sorted(df["proteins_id"].unique().tolist())

    cfg = {
        "id_column": "proteins_id",
        "LUCA_ids": [],
        "editing_ids": [],
        "promiscuity": {str(i): None for i in ids},
    }

    out = Path("metadata/proRS_annotation_config_numeric.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(cfg, indent=2))
    print("Wrote", out)


if __name__ == "__main__":
    main()
