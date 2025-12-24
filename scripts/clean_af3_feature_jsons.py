#!/usr/bin/env python3
"""
Remove non-standard ligand keys from AF3 feature-cache JSONs and write cleaned copies.

Usage:
PYTHONPATH=. python -m scripts.clean_af3_feature_jsons
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import List

CLEAN_KEYS = ["ligand", "ligands", "small_molecule", "small_molecules"]
INPUT_GLOB = Path("phase2/af3_ligand_inputs").glob("proRS_*.json")


def clean_file(path: Path) -> Path:
    with path.open() as handle:
        data = json.load(handle)

    for key in CLEAN_KEYS:
        if key in data:
            data.pop(key)

    out_path = path.with_name(path.stem + "_clean.json")
    with out_path.open("w") as handle:
        json.dump(data, handle)
    return out_path


def main() -> None:
    input_paths = sorted(INPUT_GLOB)
    cleaned: List[Path] = []
    for path in input_paths:
        cleaned.append(clean_file(path))

    print(f"Cleaned {len(cleaned)} files:")
    for path in cleaned:
        print(f"  {path}")


if __name__ == "__main__":
    main()
