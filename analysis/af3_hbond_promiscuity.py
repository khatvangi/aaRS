#!/usr/bin/env python3
"""
Summarize ligand–pocket H-bonds and contacts for AF3 ProRS complexes.

Usage:
PYTHONPATH=. python -m analysis.af3_hbond_promiscuity \
    --config metadata/af3_promiscuity_complexes.json \
    --out-csv summary/af3_hbond_summary.csv \
    --out-fig figures/af3_hbond_bars.png
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Iterable, List

import matplotlib
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def load_structure(path: Path, structure_id: str):
    suffix = path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(structure_id, path)


def find_ligand_residues(structure) -> List:
    ligands = []
    for model in structure:
        for chain in model:
            for res in chain:
                hetflag = res.id[0].strip()
                if hetflag and res.get_resname().upper() not in {"HOH", "WAT"}:
                    ligands.append(res)
    return ligands


def collect_protein_atoms(structure, exclude: Iterable) -> List:
    exclude_set = set(exclude)
    atoms = []
    for model in structure:
        for chain in model:
            for res in chain:
                if res in exclude_set:
                    continue
                if res.id[0].strip():  # hetero
                    continue
                atoms.extend(res.get_atoms())
    return atoms


def compute_hbonds(structure, max_dist: float = 3.5, pocket_shell: float = 4.0):
    """
    Simple H-bond/contact proxy:
      - ligand atoms: any non-water hetero residues
      - pocket atoms: any protein atoms within pocket_shell of any ligand atom
      - H-bond: ligand N/O < max_dist to protein N/O
    """
    ligands = find_ligand_residues(structure)
    if not ligands:
        raise ValueError("No ligand residues found (non-water hetero atoms).")

    ligand_atoms = [a for res in ligands for a in res.get_atoms()]
    protein_atoms = collect_protein_atoms(structure, exclude=ligands)

    pocket_atoms = []
    for patom in protein_atoms:
        for latom in ligand_atoms:
            if patom - latom <= pocket_shell:
                pocket_atoms.append(patom)
                break
    pocket_atoms = list({a.get_full_id(): a for a in pocket_atoms}.values())

    hbonds = []
    for patom in pocket_atoms:
        if patom.element.upper() not in {"N", "O"}:
            continue
        for latom in ligand_atoms:
            if latom.element.upper() not in {"N", "O"}:
                continue
            dist = patom - latom
            if dist <= max_dist:
                hbonds.append((patom, latom, dist))

    return ligands, pocket_atoms, hbonds


def summarize_complex(label: str, pdb_path: Path):
    structure = load_structure(pdb_path, label)
    ligands, pocket_atoms, hbonds = compute_hbonds(structure)

    ligand_names = {res.get_resname().strip() for res in ligands}
    ligand_name = "/".join(sorted(ligand_names))

    per_res_counts = defaultdict(int)
    for patom, _, _ in hbonds:
        _, _, chain_id, res_id = patom.get_full_id()[0:4]
        resname = patom.get_parent().get_resname().strip()
        key = (chain_id, res_id[1], resname)
        per_res_counts[key] += 1

    rows = []
    for (chain_id, resseq, resname), count in per_res_counts.items():
        rows.append(
            {
                "label": label,
                "ligand_name": ligand_name,
                "chain": chain_id,
                "resseq": resseq,
                "resname": resname,
                "n_hbonds": count,
            }
        )

    summary = {
        "label": label,
        "ligand_name": ligand_name,
        "n_ligand_atoms": sum(1 for _ in (a for res in ligands for a in res.get_atoms())),
        "n_pocket_atoms": len(pocket_atoms),
        "n_hbonds": len(hbonds),
    }

    return rows, summary


def plot_summary(df_summary: pd.DataFrame, out_fig: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    x = np.arange(len(df_summary))
    colors = ["#D4A017", "#C06060", "#4A6FA5", "#4A6FA5"][: len(df_summary)]
    ax.bar(x, df_summary["n_hbonds"], color=colors)
    ax.set_xticks(x)
    ax.set_xticklabels(df_summary["label"], rotation=45, ha="right")
    ax.set_ylabel("Ligand–pocket H-bonds (AF3 proxy)")
    ax.set_title("Hydrogen-bond network size across ProRS complexes")
    fig.tight_layout()
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, dpi=300)


def main():
    ap = argparse.ArgumentParser(description="Summarize ligand–pocket H-bonds for AF3 complexes.")
    ap.add_argument("--config", required=True, help="JSON list of {label, pdb}")
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--out-fig", required=True)
    args = ap.parse_args()

    cfg = pd.read_json(args.config)
    all_rows: List[dict] = []
    summaries: List[dict] = []

    for _, row in cfg.iterrows():
        label = row["label"]
        pdb_path = Path(row["pdb"])
        if not pdb_path.exists():
            raise FileNotFoundError(f"PDB/CIF not found: {pdb_path}")
        print(f"[af3_hbond_promiscuity] Analyzing {label}: {pdb_path}")
        per_res, summary = summarize_complex(label, pdb_path)
        all_rows.extend(per_res)
        summaries.append(summary)

    df_contacts = pd.DataFrame(all_rows)
    df_summary = pd.DataFrame(summaries)

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df_contacts.to_csv(out_csv, index=False)

    plot_summary(df_summary, Path(args.out_fig))
    print(f"[af3_hbond_promiscuity] Wrote {out_csv} and {args.out_fig}")


if __name__ == "__main__":
    main()
