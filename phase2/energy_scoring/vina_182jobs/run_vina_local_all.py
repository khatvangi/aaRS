#!/usr/bin/env python3
"""
Run Vina with --local_only on ALL 182 jobs for consistency.
This performs local minimization to resolve minor clashes.
"""

import subprocess
import os
from pathlib import Path
import pandas as pd
import re

VINA = "/storage/kiran-stuff/conda-envs/vina/bin/vina"
OUTPUT_DIR = Path("/storage/kiran-stuff/aaRS/phase2/energy_scoring/vina_182jobs")

def get_ligand_center(pdb_path):
    """Get centroid of ligand"""
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
                except:
                    pass
    if not coords:
        return None, None, None
    cx = sum(c[0] for c in coords) / len(coords)
    cy = sum(c[1] for c in coords) / len(coords)
    cz = sum(c[2] for c in coords) / len(coords)
    return cx, cy, cz


def run_vina_local(work_dir):
    """Run Vina with local_only minimization"""
    receptor = work_dir / "receptor.pdbqt"
    ligand = work_dir / "ligand.pdbqt"

    if not receptor.exists() or not ligand.exists():
        return None, "Missing PDBQT files"

    cx, cy, cz = get_ligand_center(work_dir / "ligand.pdb")
    if cx is None:
        return None, "Could not get ligand center"

    cmd = [
        VINA,
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--center_x", f"{cx:.3f}",
        "--center_y", f"{cy:.3f}",
        "--center_z", f"{cz:.3f}",
        "--size_x", "25",
        "--size_y", "25",
        "--size_z", "25",
        "--out", str(work_dir / "ligand_min.pdbqt"),
        "--local_only",
        "--exhaustiveness", "1"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    # Parse energy
    energy = None
    for line in result.stdout.split('\n'):
        if 'Estimated Free Energy of Binding' in line:
            match = re.search(r':\s*([-\d.]+)', line)
            if match:
                energy = float(match.group(1))
            break

    # Save log
    with open(work_dir / "vina_local.log", "w") as f:
        f.write(result.stdout)
        if result.stderr:
            f.write("\n--- STDERR ---\n")
            f.write(result.stderr)

    return energy, None


def main():
    print("=" * 60)
    print("Vina LOCAL MINIMIZATION - All 182 Jobs")
    print("=" * 60)

    # Load original results to get job list
    df = pd.read_csv(OUTPUT_DIR / "vina_scores_182jobs_fixed.csv")

    results = []

    for i, row in df.iterrows():
        job_name = row['job_name']
        ligand = row['ligand']

        # Find work directory
        safe_name = f"{job_name}_{ligand}".replace("/", "_")
        work_dir = OUTPUT_DIR / safe_name

        if not work_dir.exists():
            print(f"[{i+1}/182] {job_name} + {ligand} ... SKIP (no dir)")
            results.append({
                'job_name': job_name,
                'ligand': ligand,
                'vina_local': None,
                'vina_original': row['vina_score'],
                'iptm': row.get('iptm'),
                'error': 'No directory'
            })
            continue

        print(f"[{i+1}/182] {job_name} + {ligand}", end=" ... ")

        energy, error = run_vina_local(work_dir)

        if energy is not None:
            delta = energy - row['vina_score'] if pd.notna(row['vina_score']) else None
            print(f"{energy:.2f} (was {row['vina_score']:.2f}, Δ={delta:+.2f})" if delta else f"{energy:.2f}")
        else:
            print(f"FAILED: {error}")

        results.append({
            'job_name': job_name,
            'ligand': ligand,
            'vina_local': energy,
            'vina_original': row['vina_score'],
            'iptm': row.get('iptm'),
            'ligand_chain': row.get('ligand_chain'),
            'error': error
        })

    # Save results
    results_df = pd.DataFrame(results)
    csv_path = OUTPUT_DIR / "vina_scores_LOCAL.csv"
    results_df.to_csv(csv_path, index=False)

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    n_success = results_df['vina_local'].notna().sum()
    print(f"Successful: {n_success}/182")

    if n_success > 0:
        local_scores = results_df['vina_local'].dropna()
        print(f"Score range: {local_scores.min():.2f} to {local_scores.max():.2f}")
        print(f"Mean: {local_scores.mean():.2f}")

        # Count positive scores
        n_positive = (local_scores > 0).sum()
        print(f"Positive scores (still clashing): {n_positive}")

    print(f"\nSaved: {csv_path}")


if __name__ == "__main__":
    main()
