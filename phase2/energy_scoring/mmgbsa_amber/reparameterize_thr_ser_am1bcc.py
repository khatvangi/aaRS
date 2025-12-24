#!/usr/bin/env python3
"""
Re-parameterize THR and SER with AM1-BCC charges for validation.
These two amino acids have hydroxyl groups that are sensitive to charge assignment.
"""
import os
import subprocess

LIGANDS_DIR = "ligands"

# Only THR and SER for AM1-BCC validation
amino_acids = ["THR", "SER"]

def parameterize_ligand(aa):
    """
    Parameterize a single amino acid ligand with AM1-BCC charges.
    """
    aa_dir = os.path.join(LIGANDS_DIR, aa)
    ligand_pdb = os.path.join(aa_dir, "ligand.pdb")

    if not os.path.exists(ligand_pdb):
        print(f"✗ {aa}: PDB not found at {ligand_pdb}")
        return False

    print(f"Parameterizing {aa} with AM1-BCC charges...")

    # AMBER tools paths
    antechamber_exe = "/home/kiran/miniforge3/envs/amber_mmgbsa/bin/antechamber"
    parmchk2_exe = "/home/kiran/miniforge3/envs/amber_mmgbsa/bin/parmchk2"

    # Antechamber: PDB -> Mol2 with GAFF2 + AM1-BCC charges
    # Note: sqm (AM1-BCC) can be slow, but more accurate than Gasteiger
    # Run antechamber in the ligands directory to allow sqm to create temp files
    # -j 5: fully judge bond types and add all hydrogens (PDB has no H)
    antechamber_cmd = [
        antechamber_exe,
        "-i", "ligand.pdb",
        "-fi", "pdb",
        "-o", "lig.mol2",
        "-fo", "mol2",
        "-at", "gaff2",
        "-c", "bcc",  # AM1-BCC charges (more accurate)
        "-j", "5",    # Judge bond types and add all H
        "-s", "2",
        "-rn", "LIG",
        "-nc", "0"  # zwitterion net charge = 0
    ]

    try:
        result = subprocess.run(
            antechamber_cmd,
            capture_output=True,
            text=True,
            timeout=600,  # 10 min timeout for AM1-BCC
            cwd=aa_dir  # Run in ligand directory
        )

        if result.returncode != 0:
            print(f"✗ {aa}: antechamber failed")
            print(f"  stderr: {result.stderr[:200]}")
            return False

    except subprocess.TimeoutExpired:
        print(f"✗ {aa}: antechamber timed out (AM1-BCC calculation took >10 min)")
        return False

    # Backup Gasteiger version first
    if os.path.exists(f"{aa_dir}/lig.mol2"):
        subprocess.run(["cp", "lig.mol2", "lig_gasteiger.mol2"], cwd=aa_dir)
    if os.path.exists(f"{aa_dir}/lig.frcmod"):
        subprocess.run(["cp", "lig.frcmod", "lig_gasteiger.frcmod"], cwd=aa_dir)

    # parmchk2: Generate missing force field parameters
    parmchk2_cmd = [
        parmchk2_exe,
        "-i", "lig.mol2",
        "-f", "mol2",
        "-o", "lig.frcmod"
    ]

    result = subprocess.run(parmchk2_cmd, capture_output=True, text=True, cwd=aa_dir)
    if result.returncode != 0:
        print(f"✗ {aa}: parmchk2 failed")
        return False

    # Check outputs exist
    if not os.path.exists(f"{aa_dir}/lig.mol2") or not os.path.exists(f"{aa_dir}/lig.frcmod"):
        print(f"✗ {aa}: Output files missing")
        return False

    print(f"✓ {aa}: AM1-BCC parameterization complete")
    return True

def main():
    print("="*60)
    print("Re-parameterizing THR and SER with AM1-BCC charges")
    print("="*60)

    success_count = 0
    for aa in amino_acids:
        if parameterize_ligand(aa):
            success_count += 1

    print("\n" + "="*60)
    print(f"Completed: {success_count}/{len(amino_acids)} amino acids")
    print("="*60)

    if success_count == len(amino_acids):
        print("\n✓ All amino acids re-parameterized successfully")
        print("  Gasteiger versions backed up as lig_gasteiger.mol2/frcmod")
    else:
        print(f"\n✗ {len(amino_acids) - success_count} failed")

if __name__ == '__main__':
    main()
