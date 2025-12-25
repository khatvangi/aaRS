#!/usr/bin/env python3
"""
fix_and_run_vina.py

Fix ligand separation (by chain, not residue name) and add hydrogens before Vina.
"""

import subprocess
import os
from pathlib import Path
import json

BASE_DIR = Path("/storage/kiran-stuff/aaRS/phase2/energy_scoring/gromacs_4systems")
VINA = "/storage/kiran-stuff/conda-envs/vina/bin/vina"
OBABEL = "/usr/bin/obabel"

SYSTEMS = {
    "systemA_thrrs_THR": {"protein_chain": "A", "ligand_chain": "B", "zn_chain": "C"},
    "systemB_thrrs_SER": {"protein_chain": "A", "ligand_chain": "B", "zn_chain": "C"},
    "systemC_thrrs_ILE": {"protein_chain": "A", "ligand_chain": "B", "zn_chain": "C"},
    "systemD_prors_edit": {"protein_chain": "A", "ligand_chain": "B", "zn_chain": None},
}


def separate_by_chain(pdb_file, output_dir, protein_chain, ligand_chain, zn_chain):
    """Separate PDB by chain ID"""

    protein_lines = []
    ligand_lines = []
    zn_lines = []

    with open(pdb_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21]

                if chain == protein_chain:
                    protein_lines.append(line)
                elif chain == ligand_chain:
                    ligand_lines.append(line)
                elif zn_chain and chain == zn_chain:
                    zn_lines.append(line)

    # Write protein
    protein_pdb = output_dir / "protein_chain.pdb"
    with open(protein_pdb, "w") as f:
        f.writelines(protein_lines)
        f.write("END\n")

    # Write ligand
    ligand_pdb = output_dir / "ligand_chain.pdb"
    with open(ligand_pdb, "w") as f:
        f.writelines(ligand_lines)
        f.write("END\n")

    # Write Zn if present
    if zn_lines:
        zn_pdb = output_dir / "zn_chain.pdb"
        with open(zn_pdb, "w") as f:
            f.writelines(zn_lines)
            f.write("END\n")

    return len(protein_lines), len(ligand_lines), len(zn_lines)


def add_hydrogens_and_convert(input_pdb, output_pdbqt, is_receptor=False):
    """Add hydrogens and convert to PDBQT"""

    # First add hydrogens
    temp_pdb = input_pdb.with_suffix('.h.pdb')

    cmd_h = [OBABEL, str(input_pdb), "-O", str(temp_pdb), "-h"]  # -h adds hydrogens
    result = subprocess.run(cmd_h, capture_output=True, text=True)

    if not temp_pdb.exists():
        print(f"    Warning: Failed to add H: {result.stderr}")
        temp_pdb = input_pdb  # Use original

    # Convert to PDBQT
    cmd_pdbqt = [OBABEL, str(temp_pdb), "-O", str(output_pdbqt)]
    if is_receptor:
        cmd_pdbqt.extend(["-xr"])

    result = subprocess.run(cmd_pdbqt, capture_output=True, text=True)

    return output_pdbqt.exists()


def get_ligand_center(pdb_file):
    """Calculate centroid of ligand"""
    coords = []

    with open(pdb_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append((x, y, z))

    if not coords:
        return None, None, None

    cx = sum(c[0] for c in coords) / len(coords)
    cy = sum(c[1] for c in coords) / len(coords)
    cz = sum(c[2] for c in coords) / len(coords)

    return cx, cy, cz


def run_vina_scoring(receptor_pdbqt, ligand_pdbqt, cx, cy, cz, output_dir):
    """Run Vina scoring"""

    # Score only (current pose)
    cmd_score = [
        VINA,
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", f"{cx:.3f}",
        "--center_y", f"{cy:.3f}",
        "--center_z", f"{cz:.3f}",
        "--size_x", "25",
        "--size_y", "25",
        "--size_z", "25",
        "--score_only"
    ]

    print(f"    Running: {' '.join(cmd_score[-8:])}")
    result = subprocess.run(cmd_score, capture_output=True, text=True)

    # Save output
    with open(output_dir / "vina_score.log", "w") as f:
        f.write(result.stdout)
        f.write("\n--- STDERR ---\n")
        f.write(result.stderr)

    # Parse affinity
    affinity = None
    for line in result.stdout.split('\n'):
        if 'Affinity:' in line:
            parts = line.split()
            for i, p in enumerate(parts):
                if p == 'Affinity:':
                    try:
                        affinity = float(parts[i+1])
                    except:
                        pass
                    break

    print(f"    Vina stdout:\n{result.stdout}")
    if result.stderr:
        print(f"    Vina stderr: {result.stderr[:200]}")

    return affinity


def process_system(system_name, config):
    """Process one system"""

    system_dir = BASE_DIR / system_name
    print(f"\n{'='*50}")
    print(f"Processing: {system_name}")
    print(f"{'='*50}")

    # 1. Separate by chain
    print("\n1. Separating by chain...")
    n_prot, n_lig, n_zn = separate_by_chain(
        system_dir / "complex.pdb",
        system_dir,
        config["protein_chain"],
        config["ligand_chain"],
        config["zn_chain"]
    )
    print(f"   Protein: {n_prot} atoms, Ligand: {n_lig} atoms, Zn: {n_zn} atoms")

    # 2. Add hydrogens and convert
    print("\n2. Adding hydrogens and converting to PDBQT...")

    receptor_pdbqt = system_dir / "receptor.pdbqt"
    ligand_pdbqt = system_dir / "ligand.pdbqt"

    print("   Converting receptor...")
    if not add_hydrogens_and_convert(system_dir / "protein_chain.pdb", receptor_pdbqt, is_receptor=True):
        print("   ERROR: Receptor conversion failed")
        return None

    print("   Converting ligand...")
    if not add_hydrogens_and_convert(system_dir / "ligand_chain.pdb", ligand_pdbqt, is_receptor=False):
        print("   ERROR: Ligand conversion failed")
        return None

    # Check H atoms in ligand
    h_count = subprocess.run(
        f"grep ' H ' {ligand_pdbqt} | wc -l",
        shell=True, capture_output=True, text=True
    )
    print(f"   Ligand H atoms: {h_count.stdout.strip()}")

    # 3. Get ligand center
    print("\n3. Calculating ligand center...")
    cx, cy, cz = get_ligand_center(system_dir / "ligand_chain.pdb")
    print(f"   Center: ({cx:.2f}, {cy:.2f}, {cz:.2f})")

    # 4. Run Vina
    print("\n4. Running Vina scoring...")
    affinity = run_vina_scoring(receptor_pdbqt, ligand_pdbqt, cx, cy, cz, system_dir)

    return {
        "affinity": affinity,
        "center": (cx, cy, cz),
        "n_protein_atoms": n_prot,
        "n_ligand_atoms": n_lig
    }


def main():
    print("="*60)
    print("AutoDock Vina - Fixed Pipeline")
    print("="*60)

    results = {}

    for system_name, config in SYSTEMS.items():
        result = process_system(system_name, config)
        if result:
            results[system_name] = result

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"{'System':<25} {'Affinity (kcal/mol)':<20}")
    print("-"*45)

    for system, result in results.items():
        aff = result.get('affinity')
        aff_str = f"{aff:.2f}" if aff else "N/A"
        print(f"{system:<25} {aff_str:<20}")

    # Save
    with open(BASE_DIR / "vina_results_fixed.json", "w") as f:
        json.dump({k: {**v, "center": list(v["center"])} for k, v in results.items()}, f, indent=2)

    print(f"\nResults saved to: {BASE_DIR / 'vina_results_fixed.json'}")


if __name__ == "__main__":
    main()
