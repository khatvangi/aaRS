#!/usr/bin/env python3
"""
run_vina.py

Run AutoDock Vina on the 4 GROMACS systems.
Prepares PDBQT files and calculates binding affinity.
"""

import subprocess
import os
from pathlib import Path
import json

BASE_DIR = Path("/storage/kiran-stuff/aaRS/phase2/energy_scoring/gromacs_4systems")
VINA = "/storage/kiran-stuff/conda-envs/vina/bin/vina"
OBABEL = "/usr/bin/obabel"

SYSTEMS = [
    "systemA_thrrs_THR",
    "systemB_thrrs_SER",
    "systemC_thrrs_ILE",
    "systemD_prors_edit"
]


def pdb_to_pdbqt(pdb_file, pdbqt_file, is_receptor=False):
    """Convert PDB to PDBQT using Open Babel"""
    cmd = [OBABEL, str(pdb_file), "-O", str(pdbqt_file)]
    if is_receptor:
        cmd.extend(["-xr"])  # receptor mode
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  Warning: {result.stderr}")
    return pdbqt_file.exists()


def get_ligand_center(pdb_file):
    """Calculate centroid of ligand from PDB"""
    x_coords, y_coords, z_coords = [], [], []

    with open(pdb_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                x_coords.append(x)
                y_coords.append(y)
                z_coords.append(z)

    if not x_coords:
        return None, None, None

    center_x = sum(x_coords) / len(x_coords)
    center_y = sum(y_coords) / len(y_coords)
    center_z = sum(z_coords) / len(z_coords)

    return center_x, center_y, center_z


def run_vina(system_dir, exhaustiveness=8):
    """Run Vina docking/scoring on a system"""

    receptor_pdb = system_dir / "protein.pdb"
    ligand_pdb = system_dir / "ligand.pdb"
    receptor_pdbqt = system_dir / "receptor.pdbqt"
    ligand_pdbqt = system_dir / "ligand.pdbqt"
    output_pdbqt = system_dir / "docked.pdbqt"

    # Convert to PDBQT
    print("  Converting receptor to PDBQT...")
    if not pdb_to_pdbqt(receptor_pdb, receptor_pdbqt, is_receptor=True):
        return None

    print("  Converting ligand to PDBQT...")
    if not pdb_to_pdbqt(ligand_pdb, ligand_pdbqt, is_receptor=False):
        return None

    # Get ligand center for search box
    cx, cy, cz = get_ligand_center(ligand_pdb)
    if cx is None:
        print("  ERROR: Could not calculate ligand center")
        return None

    print(f"  Ligand center: ({cx:.2f}, {cy:.2f}, {cz:.2f})")

    # Run Vina - score only mode (ligand already in place from AF3)
    cmd = [
        VINA,
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", str(cx),
        "--center_y", str(cy),
        "--center_z", str(cz),
        "--size_x", "20",
        "--size_y", "20",
        "--size_z", "20",
        "--score_only"  # Just score the current pose
    ]

    print(f"  Running Vina score_only...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Parse affinity from output
    affinity = None
    for line in result.stdout.split('\n'):
        if 'Affinity' in line:
            parts = line.split()
            for i, p in enumerate(parts):
                if p == 'Affinity:':
                    affinity = float(parts[i+1])
                    break

    # Also try local optimization
    cmd_local = [
        VINA,
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", str(cx),
        "--center_y", str(cy),
        "--center_z", str(cz),
        "--size_x", "20",
        "--size_y", "20",
        "--size_z", "20",
        "--out", str(output_pdbqt),
        "--local_only",  # Local optimization only
        "--exhaustiveness", str(exhaustiveness)
    ]

    print(f"  Running Vina local optimization...")
    result_local = subprocess.run(cmd_local, capture_output=True, text=True)

    # Parse best affinity from local optimization
    best_affinity = None
    for line in result_local.stdout.split('\n'):
        if line.strip().startswith('1'):
            parts = line.split()
            if len(parts) >= 2:
                try:
                    best_affinity = float(parts[1])
                except:
                    pass
                break

    return {
        "score_only_affinity": affinity,
        "local_opt_affinity": best_affinity,
        "center": (cx, cy, cz),
        "stdout": result.stdout,
        "stdout_local": result_local.stdout
    }


def main():
    print("="*60)
    print("AutoDock Vina Scoring for 4 Systems")
    print("="*60)

    results = {}

    for system in SYSTEMS:
        print(f"\n--- {system} ---")
        system_dir = BASE_DIR / system

        if not system_dir.exists():
            print(f"  ERROR: Directory not found")
            continue

        result = run_vina(system_dir)

        if result:
            results[system] = result
            print(f"  Score-only affinity: {result['score_only_affinity']} kcal/mol")
            print(f"  Local-opt affinity:  {result['local_opt_affinity']} kcal/mol")
        else:
            print(f"  FAILED")

    # Summary table
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"{'System':<25} {'Score-Only':<15} {'Local-Opt':<15}")
    print("-"*55)

    for system, result in results.items():
        score = result.get('score_only_affinity', 'N/A')
        local = result.get('local_opt_affinity', 'N/A')
        score_str = f"{score:.2f}" if isinstance(score, float) else str(score)
        local_str = f"{local:.2f}" if isinstance(local, float) else str(local)
        print(f"{system:<25} {score_str:<15} {local_str:<15}")

    # Save results
    results_file = BASE_DIR / "vina_results.json"
    with open(results_file, 'w') as f:
        # Convert tuples to lists for JSON
        json_results = {}
        for k, v in results.items():
            json_results[k] = {
                "score_only_affinity": v.get("score_only_affinity"),
                "local_opt_affinity": v.get("local_opt_affinity"),
                "center": list(v.get("center", []))
            }
        json.dump(json_results, f, indent=2)

    print(f"\nResults saved to: {results_file}")


if __name__ == "__main__":
    main()
