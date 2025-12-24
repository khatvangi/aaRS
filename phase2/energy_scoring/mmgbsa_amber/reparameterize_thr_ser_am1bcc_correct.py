#!/usr/bin/env python3
"""
Re-parameterize THR and SER with AM1-BCC charges using CORRECT net charge.

Previous error: Did not specify -nc 0, causing sqm to fail with "odd electrons" error.
Correct approach: Explicitly set -nc 0 for zwitterionic amino acids.
"""
import os
import subprocess

LIGANDS_DIR = "ligands"
amino_acids = ["THR", "SER"]

def parameterize_with_am1bcc(aa):
    """
    Parameterize with AM1-BCC charges, explicitly setting net charge = 0.
    """
    aa_dir = os.path.join(LIGANDS_DIR, aa)
    ligand_pdb = os.path.join(aa_dir, "ligand.pdb")

    if not os.path.exists(ligand_pdb):
        return f"{aa}: SKIP - no ligand.pdb"

    print(f"Parameterizing {aa} with AM1-BCC charges (-nc 0)...")

    # Backup Gasteiger version
    if os.path.exists(f"{aa_dir}/lig.mol2"):
        subprocess.run(["cp", "lig.mol2", "lig_gasteiger.mol2"], cwd=aa_dir, capture_output=True)
    if os.path.exists(f"{aa_dir}/lig.frcmod"):
        subprocess.run(["cp", "lig.frcmod", "lig_gasteiger.frcmod"], cwd=aa_dir, capture_output=True)

    # CRITICAL FIX: Add -nc 0 to specify net charge for zwitterion
    antechamber_exe = "/home/kiran/miniforge3/envs/amber_mmgbsa/bin/antechamber"
    antechamber_cmd = [
        antechamber_exe,
        "-i", "ligand.pdb",
        "-fi", "pdb",
        "-o", "lig_am1bcc.mol2",  # Separate file for AM1-BCC
        "-fo", "mol2",
        "-at", "gaff2",
        "-c", "bcc",      # AM1-BCC charges
        "-nc", "0",       # CRITICAL: net charge = 0 for zwitterion
        "-j", "5",        # Add all hydrogens
        "-s", "2",
        "-rn", "LIG"
    ]

    try:
        result = subprocess.run(
            antechamber_cmd,
            capture_output=True,
            text=True,
            timeout=600,  # AM1-BCC can take a few minutes
            cwd=aa_dir
        )

        if result.returncode != 0:
            print(f"  ✗ FAILED: {result.stderr[:200]}")
            return f"{aa}: FAIL antechamber"

    except subprocess.TimeoutExpired:
        print(f"  ✗ TIMEOUT (>10 min)")
        return f"{aa}: FAIL timeout"

    # parmchk2 for AM1-BCC version
    parmchk2_exe = "/home/kiran/miniforge3/envs/amber_mmgbsa/bin/parmchk2"
    parmchk2_cmd = [
        parmchk2_exe,
        "-i", "lig_am1bcc.mol2",
        "-f", "mol2",
        "-o", "lig_am1bcc.frcmod"
    ]

    result = subprocess.run(parmchk2_cmd, capture_output=True, text=True, cwd=aa_dir)
    if result.returncode != 0:
        print(f"  ✗ FAILED parmchk2")
        return f"{aa}: FAIL parmchk2"

    # Verify outputs
    if not os.path.exists(f"{aa_dir}/lig_am1bcc.mol2") or not os.path.exists(f"{aa_dir}/lig_am1bcc.frcmod"):
        print(f"  ✗ MISSING output files")
        return f"{aa}: FAIL missing output"

    # Count atoms
    with open(f"{aa_dir}/lig_am1bcc.mol2") as f:
        lines = f.readlines()
    for line in lines:
        if line.strip().startswith("@<TRIPOS>ATOM"):
            idx = lines.index(line)
            natoms = len([l for l in lines[idx+1:] if l.strip() and not l.startswith("@")])
            print(f"  ✓ SUCCESS ({natoms} atoms)")
            return f"{aa}: OK ({natoms} atoms)"

    return f"{aa}: OK"

def main():
    print("="*70)
    print("Re-parameterizing THR and SER with AM1-BCC charges")
    print("Using CORRECT net charge: -nc 0")
    print("="*70)
    print()

    results = []
    for aa in amino_acids:
        result = parameterize_with_am1bcc(aa)
        results.append(result)

    print()
    print("="*70)
    print("Results:")
    for res in results:
        print(f"  {res}")
    print("="*70)

    success = sum(1 for r in results if "OK" in r)
    if success == len(amino_acids):
        print()
        print("✓ AM1-BCC parameterization successful!")
        print("  Files created:")
        print("    - lig_am1bcc.mol2 (AM1-BCC charges)")
        print("    - lig_am1bcc.frcmod")
        print("    - lig_gasteiger.mol2 (backup)")
        print("    - lig_gasteiger.frcmod (backup)")
        print()
        print("Next: Run validation on ~10 ThrRS+Zn models")
    else:
        print()
        print(f"✗ {len(amino_acids) - success} failed")

if __name__ == '__main__':
    main()
