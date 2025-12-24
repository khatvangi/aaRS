#!/usr/bin/env python3
"""
Pre-parameterize 20 amino acid ligands for MMPBSA.
Extract ligand PDB from one example CIF for each AA, then run antechamber + parmchk2.
"""
import os
import sys
import subprocess
import pandas as pd
import json

AA20 = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

def main():
    # Load manifest
    df = pd.read_csv('../mmgbsa_real/manifest_af3only.csv')

    # Get one example per standard AA
    examples = {}
    for aa in AA20:
        rows = df[df['lig_resname'] == aa]
        if len(rows) > 0:
            examples[aa] = rows.iloc[0]['file']
        else:
            print(f"WARNING: No example found for {aa}")

    print(f"Found examples for {len(examples)}/20 amino acids")

    # Create ligands directory
    os.makedirs('ligands', exist_ok=True)

    # Process each AA
    for aa, cif_path in examples.items():
        print(f"\n{'='*60}")
        print(f"Processing {aa}: {cif_path}")
        print('='*60)

        aa_dir = f"ligands/{aa}"
        os.makedirs(aa_dir, exist_ok=True)

        # Step 1: Extract ligand PDB using prep_from_cif.py
        prep_cmd = [
            sys.executable,
            '../prep_from_cif.py',
            f'../../{cif_path}',
            aa_dir
        ]
        result = subprocess.run(prep_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"ERROR extracting {aa}: {result.stderr}")
            continue

        # Check ligand.pdb exists
        ligand_pdb = f"{aa_dir}/ligand.pdb"
        if not os.path.exists(ligand_pdb):
            print(f"ERROR: {ligand_pdb} not created")
            continue

        # Step 2: Run antechamber (GAFF2, Gasteiger charges for speed)
        # Note: Using Gasteiger instead of AM1-BCC for faster parameterization
        #       AM1-BCC requires slow QM calculation; Gasteiger is empirical & fast
        print(f"  Running antechamber...")
        antechamber_exe = "/home/kiran/miniforge3/envs/amber_mmgbsa/bin/antechamber"
        antechamber_cmd = [
            antechamber_exe,
            "-i", ligand_pdb,
            "-fi", "pdb",
            "-o", f"{aa_dir}/lig.mol2",
            "-fo", "mol2",
            "-at", "gaff2",
            "-c", "gas",  # Gasteiger charges (fast)
            "-s", "2",
            "-rn", "LIG",
            "-nc", "0"  # zwitterion net charge = 0
        ]
        result = subprocess.run(antechamber_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"ERROR antechamber {aa}: {result.stderr[:500]}")
            continue

        # Step 3: Run parmchk2
        print(f"  Running parmchk2...")
        parmchk2_exe = "/home/kiran/miniforge3/envs/amber_mmgbsa/bin/parmchk2"
        parmchk_cmd = [
            parmchk2_exe,
            "-i", f"{aa_dir}/lig.mol2",
            "-f", "mol2",
            "-o", f"{aa_dir}/lig.frcmod"
        ]
        result = subprocess.run(parmchk_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"ERROR parmchk2 {aa}: {result.stderr[:500]}")
            continue

        print(f"  SUCCESS: {aa_dir}/lig.mol2 and lig.frcmod created")

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print('='*60)
    success = []
    for aa in AA20:
        if os.path.exists(f"ligands/{aa}/lig.mol2") and os.path.exists(f"ligands/{aa}/lig.frcmod"):
            success.append(aa)

    print(f"Successfully parameterized: {len(success)}/20")
    print(f"AAs: {', '.join(success)}")

    if len(success) < 20:
        missing = [aa for aa in AA20 if aa not in success]
        print(f"\nMissing: {', '.join(missing)}")

if __name__ == '__main__':
    main()
