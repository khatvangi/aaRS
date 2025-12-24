#!/usr/bin/env python3
"""
Build AF3 ProRS complexes by transplanting the PYR ligand from a holo template.

Usage:
PYTHONPATH=. python -m analysis.make_af3_proRS_complexes \
    --template templates/1H4Q_proRS_holo.pdb \
    --ligand-resn PYR \
    --models-json metadata/af3_proRS_models.json \
    --out-dir phase2/af3_complexes
"""

from __future__ import annotations

import argparse
import json
from io import StringIO
from pathlib import Path
from typing import Dict, List

from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select, Superimposer


def load_structure(path: Path, structure_id: str):
    suffix = path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(structure_id, path)


def get_ca_atoms(structure) -> List:
    cas = []
    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0].strip():
                    continue
                if "CA" in res:
                    cas.append(res["CA"])
    if not cas:
        raise ValueError("No CA atoms found for alignment.")
    return cas


def get_ligand_residues(structure, ligand_resn: str):
    ligand_resn = ligand_resn.upper()
    ligands = []
    for res in structure.get_residues():
        if res.id[0].strip() and res.get_resname().upper() == ligand_resn:
            ligands.append(res)
    if not ligands:
        raise ValueError(f"Ligand {ligand_resn} not found in template.")
    return ligands


class ProteinSelect(Select):
    def accept_residue(self, residue):
        return not residue.id[0].strip()


class LigandSelect(Select):
    def __init__(self, ligand_resn: str):
        super().__init__()
        self.ligand_resn = ligand_resn.upper()

    def accept_residue(self, residue):
        return residue.id[0].strip() and residue.get_resname().upper() == self.ligand_resn


def superimpose_model(template_cas, model_cas, model_structure):
    n = min(len(template_cas), len(model_cas))
    if n == 0:
        raise ValueError("No overlapping CA atoms for alignment.")
    if len(template_cas) != len(model_cas):
        print(f"[make_af3_proRS_complexes] Truncating CA alignment to {n} residues.")
    sup = Superimposer()
    sup.set_atoms(template_cas[:n], model_cas[:n])
    sup.apply(list(model_structure.get_atoms()))
    return sup


def write_complex(
    model_structure, template_structure, ligand_resn: str, out_path: Path
) -> None:
    protein_buf = StringIO()
    ligand_buf = StringIO()

    protein_io = PDBIO()
    protein_io.set_structure(model_structure)
    protein_io.save(protein_buf, select=ProteinSelect())

    ligand_io = PDBIO()
    ligand_io.set_structure(template_structure)
    ligand_io.save(ligand_buf, select=LigandSelect(ligand_resn))

    protein_lines = protein_buf.getvalue().rstrip().splitlines()
    if protein_lines and protein_lines[-1].startswith("END"):
        protein_lines = protein_lines[:-1]
    ligand_lines = [ln for ln in ligand_buf.getvalue().splitlines() if ln.strip()]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as handle:
        for line in protein_lines:
            handle.write(f"{line}\n")
        for line in ligand_lines:
            if line.startswith("ATOM"):
                line = "HETATM" + line[6:]
            handle.write(f"{line}\n")
        handle.write("END\n")


def process_model(
    template_structure, ligand_resn: str, model_label: str, model_path: Path, out_dir: Path
) -> Path:
    print(f"[make_af3_proRS_complexes] Processing {model_label}: {model_path}")
    model_structure = load_structure(model_path, model_label)

    template_cas = get_ca_atoms(template_structure)
    model_cas = get_ca_atoms(model_structure)
    superimpose_model(template_cas, model_cas, model_structure)

    out_path = out_dir / f"{model_label}_pro_complex.pdb"
    write_complex(model_structure, template_structure, ligand_resn, out_path)
    print(f"[make_af3_proRS_complexes] Wrote {out_path}")
    return out_path


def parse_args():
    ap = argparse.ArgumentParser(description="Transplant template ligand into AF3 ProRS models.")
    ap.add_argument("--template", required=True, help="Template PDB with bound ligand.")
    ap.add_argument("--ligand-resn", required=True, help="Ligand residue name to transplant.")
    ap.add_argument("--models-json", required=True, help='JSON list of {"label", "pdb"} entries.')
    ap.add_argument("--out-dir", required=True, help="Output directory for complexes.")
    return ap.parse_args()


def main():
    args = parse_args()
    template_path = Path(args.template)
    models_json = Path(args.models_json)
    out_dir = Path(args.out_dir)

    if not template_path.exists():
        raise FileNotFoundError(f"Template not found: {template_path}")
    if not models_json.exists():
        raise FileNotFoundError(f"Model list not found: {models_json}")

    template_structure = load_structure(template_path, "template")
    ligand_residues = get_ligand_residues(template_structure, args.ligand_resn)
    print(
        f"[make_af3_proRS_complexes] Found {len(ligand_residues)} ligand residue(s) "
        f"named {args.ligand_resn} in template."
    )

    with models_json.open() as handle:
        models: List[Dict[str, str]] = json.load(handle)

    for entry in models:
        label = entry["label"]
        pdb_path = Path(entry["pdb"])
        if not pdb_path.exists():
            raise FileNotFoundError(f"Model file missing for {label}: {pdb_path}")
        process_model(template_structure, args.ligand_resn, label, pdb_path, out_dir)


if __name__ == "__main__":
    main()
