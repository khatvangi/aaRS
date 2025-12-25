#!/usr/bin/env python3
"""
make_balanced_jobs.py

Reads conditions.yaml and ligands.txt, generates:
  - inputs/json/*.json (AF3 input JSONs)
  - manifests/jobs.csv (job manifest)

Total jobs = 6 conditions x 12 ligands x 5 seeds = 360
"""

import yaml
import json
import csv
from pathlib import Path

BASE_DIR = Path(__file__).parent
CONFIGS_DIR = BASE_DIR / "configs"
INPUTS_DIR = BASE_DIR / "inputs" / "json"
MANIFESTS_DIR = BASE_DIR / "manifests"

SEEDS = [1, 2, 3, 4, 5]

# CCD codes for amino acid ligands
LIGAND_CCD = {
    "THR": "THR",
    "SER": "SER",
    "PRO": "PRO",
    "ALA": "ALA",
    "GLY": "GLY",
    "VAL": "VAL",
    "ILE": "ILE",
    "LEU": "LEU",
    "CYS": "CYS",
    "MET": "MET",
    "TYR": "TYR",
    "GLU": "GLU",
}


def load_conditions():
    with open(CONFIGS_DIR / "conditions.yaml") as f:
        data = yaml.safe_load(f)
    return data["conditions"]


def load_ligands():
    with open(CONFIGS_DIR / "ligands.txt") as f:
        return [line.strip() for line in f if line.strip()]


def read_fasta(fasta_path):
    """Read FASTA and return sequence (single chain assumed)."""
    seq_lines = []
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                seq_lines.append(line.strip())
    return "".join(seq_lines)


def make_af3_json(job_name, protein_seq, ligand_code, include_zn, seed):
    """
    Create AF3 input JSON structure.
    """
    sequences = [
        {
            "protein": {
                "id": "A",
                "sequence": protein_seq
            }
        },
        {
            "ligand": {
                "id": "B",
                "ccdCodes": [ligand_code]
            }
        }
    ]

    if include_zn:
        sequences.append({
            "ligand": {
                "id": "C",
                "ccdCodes": ["ZN"]
            }
        })

    return {
        "name": job_name,
        "modelSeeds": [seed],
        "sequences": sequences,
        "dialect": "alphafold3",
        "version": 1
    }


def main():
    conditions = load_conditions()
    ligands = load_ligands()

    # Check for placeholder paths
    for cond_name, cond_data in conditions.items():
        if "PLACEHOLDER" in cond_data["fasta"]:
            print(f"ERROR: {cond_name} still has placeholder FASTA path")
            print("Please update configs/conditions.yaml with real FASTA paths")
            return

    jobs = []

    for cond_name, cond_data in conditions.items():
        fasta_path = cond_data["fasta"]
        include_zn = cond_data.get("include_zn", False)

        try:
            protein_seq = read_fasta(fasta_path)
        except Exception as e:
            print(f"ERROR reading FASTA for {cond_name}: {e}")
            return

        for ligand in ligands:
            for seed in SEEDS:
                job_name = f"{cond_name}_{ligand}_seed{seed}"
                json_path = INPUTS_DIR / f"{job_name}.json"

                # Create JSON
                af3_input = make_af3_json(
                    job_name=job_name,
                    protein_seq=protein_seq,
                    ligand_code=LIGAND_CCD[ligand],
                    include_zn=include_zn,
                    seed=seed
                )

                with open(json_path, "w") as f:
                    json.dump(af3_input, f, indent=2)

                jobs.append({
                    "condition": cond_name,
                    "ligand": ligand,
                    "seed": seed,
                    "job_name": job_name,
                    "input_json": str(json_path),
                    "include_zn": include_zn
                })

    # Write manifest
    manifest_path = MANIFESTS_DIR / "jobs.csv"
    with open(manifest_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["condition", "ligand", "seed", "job_name", "input_json", "include_zn"])
        writer.writeheader()
        writer.writerows(jobs)

    print(f"Created {len(jobs)} jobs")
    print(f"JSONs written to: {INPUTS_DIR}")
    print(f"Manifest written to: {manifest_path}")

    # Verify count
    expected = 6 * 12 * 5
    if len(jobs) != expected:
        print(f"WARNING: Expected {expected} jobs, got {len(jobs)}")
    else:
        print(f"OK: Job count matches expected ({expected})")


if __name__ == "__main__":
    main()
