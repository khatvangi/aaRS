#!/usr/bin/env python3
"""
run_vina_fixed.py

Run AutoDock Vina on all 182 AF3 jobs.
FIXED: Dynamically detect ligand chain (not hardcoded to chain B).
"""

import subprocess
import os
from pathlib import Path
import pandas as pd
import json
import re
from concurrent.futures import ProcessPoolExecutor, as_completed

VINA = "/storage/kiran-stuff/conda-envs/vina/bin/vina"
OBABEL = "/usr/bin/obabel"
OUTPUT_DIR = Path("/storage/kiran-stuff/aaRS/phase2/energy_scoring/vina_182jobs")
OUTPUT_DIR.mkdir(exist_ok=True)

DF_PATH = "/storage/kiran-stuff/aaRS/phase2/figures/data/comprehensive_ligand_analysis.csv"


def get_af3_jobs():
    """Get list of unique AF3 jobs from comprehensive analysis"""
    df = pd.read_csv(DF_PATH)
    df = df[~df['job_name'].str.contains('seed-', na=False)]
    df = df[df['job_name'].notna()]
    df['job_base'] = df['job_name'].str.replace(r'_\d{8}_\d{6}$', '', regex=True)
    df = df[df['AA_iptm'].notna()]
    df = df.drop_duplicates(subset=['job_base', 'ligand'], keep='first')
    df = df[~df['job_name'].str.contains('final_results', case=False, na=False)]

    jobs = []
    for _, row in df.iterrows():
        cif_path = row['cif_file']
        if pd.notna(cif_path) and Path(cif_path).exists():
            jobs.append({
                'job_name': row['job_base'],
                'ligand': row['ligand'],
                'cif_path': cif_path,
                'iptm': row['AA_iptm']
            })
    return jobs


def cif_to_pdb(cif_path, pdb_path):
    """Convert CIF to PDB using gemmi"""
    try:
        import gemmi
        structure = gemmi.read_structure(str(cif_path))
        structure.write_pdb(str(pdb_path))
        return True
    except Exception as e:
        print(f"CIF conversion error: {e}")
        return False


def detect_ligand_chain(pdb_path):
    """
    Detect ligand chain dynamically.
    Ligand is typically the smallest chain with 1 residue.
    Returns (protein_chains, ligand_chain)
    """
    chain_residue_counts = {}

    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21]
                resnum = line[22:26].strip()
                resname = line[17:20].strip()
                key = (chain, resnum, resname)
                if chain not in chain_residue_counts:
                    chain_residue_counts[chain] = set()
                chain_residue_counts[chain].add(key)

    # Find chains with only 1 residue (likely ligand)
    ligand_chains = []
    protein_chains = []

    for chain, residues in chain_residue_counts.items():
        if len(residues) == 1:
            ligand_chains.append(chain)
        else:
            protein_chains.append(chain)

    # If no single-residue chain found, check for chain L (common for ligands)
    if not ligand_chains:
        if 'L' in chain_residue_counts:
            ligand_chains = ['L']
            protein_chains = [c for c in chain_residue_counts if c != 'L']

    return protein_chains, ligand_chains


def separate_by_chain(pdb_path, work_dir, protein_chains, ligand_chain):
    """Separate protein from ligand by chain IDs"""
    protein_lines = []
    ligand_lines = []

    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21]
                if chain in protein_chains:
                    protein_lines.append(line)
                elif chain == ligand_chain:
                    ligand_lines.append(line)

    protein_pdb = work_dir / "protein.pdb"
    ligand_pdb = work_dir / "ligand.pdb"

    with open(protein_pdb, "w") as f:
        f.writelines(protein_lines)
        f.write("END\n")

    with open(ligand_pdb, "w") as f:
        f.writelines(ligand_lines)
        f.write("END\n")

    return len(protein_lines), len(ligand_lines)


def add_hydrogens_and_convert(input_pdb, output_pdbqt, is_receptor=False):
    """Add hydrogens using obabel and convert to PDBQT"""
    temp_h_pdb = input_pdb.with_suffix('.h.pdb')
    cmd_h = [OBABEL, str(input_pdb), "-O", str(temp_h_pdb), "-h", "--partialcharge", "gasteiger"]
    subprocess.run(cmd_h, capture_output=True, text=True)

    if not temp_h_pdb.exists():
        temp_h_pdb = input_pdb

    cmd_pdbqt = [OBABEL, str(temp_h_pdb), "-O", str(output_pdbqt)]
    if is_receptor:
        cmd_pdbqt.append("-xr")

    subprocess.run(cmd_pdbqt, capture_output=True, text=True)
    return output_pdbqt.exists()


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


def run_vina_score(receptor_pdbqt, ligand_pdbqt, cx, cy, cz):
    """Run Vina score_only"""
    cmd = [
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

    result = subprocess.run(cmd, capture_output=True, text=True)

    energy = None
    for line in result.stdout.split('\n'):
        if 'Estimated Free Energy of Binding' in line:
            match = re.search(r':\s*([-\d.]+)', line)
            if match:
                energy = float(match.group(1))
            break

    return energy, result.stdout


def process_job(job_info):
    """Process a single AF3 job"""
    job_name = job_info['job_name']
    ligand = job_info['ligand']
    cif_path = job_info['cif_path']

    safe_name = f"{job_name}_{ligand}".replace("/", "_")
    work_dir = OUTPUT_DIR / safe_name
    work_dir.mkdir(exist_ok=True)

    try:
        # 1. Convert CIF to PDB
        pdb_path = work_dir / "complex.pdb"
        if not cif_to_pdb(cif_path, pdb_path):
            return {**job_info, 'vina_score': None, 'error': 'CIF conversion failed', 'ligand_chain': None}

        # 2. Detect ligand chain dynamically
        protein_chains, ligand_chains = detect_ligand_chain(pdb_path)

        if not ligand_chains:
            return {**job_info, 'vina_score': None, 'error': 'No ligand chain detected', 'ligand_chain': None}

        ligand_chain = ligand_chains[0]  # Use first ligand chain

        # 3. Separate by chain
        n_prot, n_lig = separate_by_chain(pdb_path, work_dir, protein_chains, ligand_chain)

        if n_lig == 0:
            return {**job_info, 'vina_score': None, 'error': 'No ligand atoms found', 'ligand_chain': ligand_chain}

        # 4. Add hydrogens and convert
        receptor_pdbqt = work_dir / "receptor.pdbqt"
        ligand_pdbqt = work_dir / "ligand.pdbqt"

        if not add_hydrogens_and_convert(work_dir / "protein.pdb", receptor_pdbqt, is_receptor=True):
            return {**job_info, 'vina_score': None, 'error': 'Receptor PDBQT failed', 'ligand_chain': ligand_chain}

        if not add_hydrogens_and_convert(work_dir / "ligand.pdb", ligand_pdbqt, is_receptor=False):
            return {**job_info, 'vina_score': None, 'error': 'Ligand PDBQT failed', 'ligand_chain': ligand_chain}

        # 5. Get ligand center
        cx, cy, cz = get_ligand_center(work_dir / "ligand.pdb")
        if cx is None:
            return {**job_info, 'vina_score': None, 'error': 'Could not get ligand center', 'ligand_chain': ligand_chain}

        # 6. Run Vina
        energy, stdout = run_vina_score(receptor_pdbqt, ligand_pdbqt, cx, cy, cz)

        with open(work_dir / "vina.log", "w") as f:
            f.write(stdout)

        return {**job_info, 'vina_score': energy, 'error': None, 'ligand_chain': ligand_chain}

    except Exception as e:
        return {**job_info, 'vina_score': None, 'error': str(e), 'ligand_chain': None}


def main():
    print("=" * 60)
    print("AutoDock Vina - 182 AF3 Jobs (FIXED: Dynamic chain detection)")
    print("=" * 60)

    jobs = get_af3_jobs()
    print(f"\nFound {len(jobs)} jobs to process")

    results = []
    for i, job in enumerate(jobs):
        print(f"\n[{i+1}/{len(jobs)}] {job['job_name']} + {job['ligand']}", end=" ... ")
        result = process_job(job)

        if result['vina_score'] is not None:
            print(f"OK: {result['vina_score']:.2f} kcal/mol (chain {result['ligand_chain']})")
        else:
            print(f"FAILED: {result.get('error', 'Unknown')}")

        results.append(result)

    df = pd.DataFrame(results)
    csv_path = OUTPUT_DIR / "vina_scores_182jobs_fixed.csv"
    df.to_csv(csv_path, index=False)
    print(f"\n\nSaved results to: {csv_path}")

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    n_success = df['vina_score'].notna().sum()
    n_failed = df['vina_score'].isna().sum()
    print(f"Successful: {n_success}")
    print(f"Failed: {n_failed}")

    if n_success > 0:
        print(f"\nVina score range: {df['vina_score'].min():.2f} to {df['vina_score'].max():.2f} kcal/mol")
        print(f"Mean: {df['vina_score'].mean():.2f} kcal/mol")

    # Show chain distribution
    print("\nLigand chain distribution:")
    print(df['ligand_chain'].value_counts())


if __name__ == "__main__":
    main()
