#!/usr/bin/env python3
"""
Build AF3 chain-list JSONs for domain constructs (protein + ligand).

Usage:
PYTHONPATH=. python -m scripts.build_domain_af3_chainlist_json
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List

LIGANDS: Dict[str, Dict[str, str]] = {
    "Pro": {
        "name": "proline",
        "sequence": "C1CC(NC1)C(=O)O",
    },
    "Thr": {
        "name": "threonine",
        "sequence": "CC(C(C(=O)O)N)O",
    },
}

DOMAIN_CONFIG = [
    {"label": "cat_LUCA_Pro", "fasta": "phase2/domain_fasta/cat_LUCA.fasta", "ligand": "Pro"},
    {"label": "cat_LUCA_Thr", "fasta": "phase2/domain_fasta/cat_LUCA.fasta", "ligand": "Thr"},
    {"label": "edit_LUCA_Pro", "fasta": "phase2/domain_fasta/edit_LUCA.fasta", "ligand": "Pro"},
    {"label": "edit_LUCA_Thr", "fasta": "phase2/domain_fasta/edit_LUCA.fasta", "ligand": "Thr"},
]

OUT_DIR = Path("phase2/domain_af3_inputs")


def read_single_fasta(path: Path) -> str:
    lines = path.read_text().splitlines()
    seq_lines = [ln.strip() for ln in lines if ln and not ln.startswith(">")]
    if not seq_lines:
        raise ValueError(f"No sequence found in FASTA: {path}")
    return "".join(seq_lines)


def build_chain_list(protein_seq: str, ligand_seq: str) -> List[Dict[str, object]]:
    return [
        {
            "chain_id": "A",
            "sequence": protein_seq,
            "is_protein": True,
            "is_na": False,
            "is_ligand": False,
        },
        {
            "chain_id": "L",
            "sequence": ligand_seq,
            "is_protein": False,
            "is_na": False,
            "is_ligand": True,
        },
    ]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []

    for entry in DOMAIN_CONFIG:
        label = entry["label"]
        fasta_path = Path(entry["fasta"])
        ligand_key = entry["ligand"]

        protein_seq = read_single_fasta(fasta_path)
        ligand_seq = LIGANDS[ligand_key]["sequence"]

        chain_list = build_chain_list(protein_seq, ligand_seq)

        out_path = OUT_DIR / f"{label}.json"
        with out_path.open("w") as handle:
            json.dump(chain_list, handle, indent=2)
        written.append(out_path)

    print(f"Wrote {len(written)} chain-list JSONs:")
    for path in written:
        print(f"  {path}")


if __name__ == "__main__":
    main()
