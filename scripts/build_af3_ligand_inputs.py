#!/usr/bin/env python3
"""
Build AF3 complex input JSONs by reusing cached apo features and adding ligands.

Run locally on boron:
PYTHONPATH=. python -m scripts.build_af3_ligand_inputs
"""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Dict, List

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
MODELS_YAML = REPO_ROOT / "metadata/af3_ligand_models.yaml"
OUT_DIR = REPO_ROOT / "phase2/af3_ligand_inputs"

LIGANDS: Dict[str, Dict[str, str]] = {
    "Pro": {"name": "proline", "smiles": "C1CC(NC1)C(=O)O"},
    "Thr": {"name": "threonine", "smiles": "CC(C(C(=O)O)N)O"},
}


def load_models_cfg(path: Path) -> Dict[str, Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Model map not found: {path}")
    with path.open() as handle:
        cfg = yaml.safe_load(handle)
    if not isinstance(cfg, dict):
        raise ValueError(f"Unexpected YAML structure in {path}")
    return cfg


def strip_existing_ligands(sequences: List[dict]) -> List[dict]:
    """Remove any existing ligand entries; keep protein/RNA blocks intact."""
    cleaned = []
    for entry in sequences:
        if "ligand" in entry:
            continue
        cleaned.append(entry)
    return cleaned


def add_ligand_block(data: dict, ligand_id: str, smiles: str) -> dict:
    """Return a copy of the AF3 input with a ligand appended."""
    out = copy.deepcopy(data)
    seqs = out.get("sequences", [])
    seqs = strip_existing_ligands(seqs)
    seqs.append({"ligand": {"id": ligand_id, "smiles": smiles}})
    out["sequences"] = seqs
    # Also add an explicit small_molecules block if the runner expects it.
    out["small_molecules"] = [{"id": ligand_id, "smiles": smiles}]
    return out


def build_for_model(label: str, model_cfg: Dict[str, str]) -> List[Path]:
    data_path = REPO_ROOT / model_cfg["data_json"]
    if not data_path.exists():
        raise FileNotFoundError(f"Missing apo data_json for {label}: {data_path}")
    with data_path.open() as handle:
        apo_data = json.load(handle)

    written: List[Path] = []
    for lig_key, lig_info in LIGANDS.items():
        ligand_id = lig_key
        complex_data = add_ligand_block(apo_data, ligand_id, lig_info["smiles"])
        out_path = OUT_DIR / f"{label}_{lig_key}.json"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w") as handle:
            json.dump(complex_data, handle, indent=2)
        written.append(out_path)
    return written


def main():
    models_cfg = load_models_cfg(MODELS_YAML)
    all_written: List[Path] = []
    for label, cfg in models_cfg.items():
        print(f"[build_af3_ligand_inputs] Using apo features for {label}: {cfg['data_json']}")
        all_written.extend(build_for_model(label, cfg))
    for path in all_written:
        print(f"Wrote {path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
