#!/usr/bin/env python3
"""
Correctly parameterize THR and SER with AM1-BCC charges.

Fix: Pre-add hydrogens with obabel before running AM1-BCC to avoid odd electron count.
"""
import os
import subprocess

LIGANDS_DIR = "ligands"
amino_acids = ["THR", "SER"]

def parameterize_am1bcc(aa):
    """
    Parameterize with AM1-BCC using obabel to add H first.
    """
    aa_dir = os.path.join(LIGANDS_DIR, aa)
    ligand_pdb = os.path.join(aa_dir, "ligand.pdb")

    if not os.path.exists(ligand_pdb):
        return f"{aa}: SKIP - no ligand.pdb"

    print(f"Parameterizing {aa} with AM1-BCC...")

    # Step 1: Add hydrogens with obabel
    print(f"  1. Adding hydrogens with obabel...")
    obabel_cmd = [
        "/usr/bin/obabel",
        "ligand.pdb",
        "-O", "ligand_H.pdb",
        "-h",           # Add hydrogens
        "-p", "7.0"     # pH 7
    ]

    result = subprocess.run(obabel_cmd, capture_output=True, text=True, cwd=aa_dir)
    if result.returncode != 0:
        print(f"  ✗ FAILED obabel: {result.stderr[:100]}")
        return f"{aa}: FAIL obabel"

    if not os.path.exists(f"{aa_dir}/ligand_H.pdb"):
        print(f"  ✗ FAILED to create ligand_H.pdb")
        return f"{aa}: FAIL no H pdb"

    # Count atoms
    with open(f"{aa_dir}/ligand_H.pdb") as f:
        natoms_H = len([l for l in f if l.startswith("HETATM") or l.startswith("ATOM")])
    print(f"  2. Hydrogens added: {natoms_H} total atoms")

    # Step 2: Run antechamber with AM1-BCC on hydrogenated structure
    print(f"  3. Running antechamber with AM1-BCC...")
    antechamber_exe = "/home/kiran/miniforge3/envs/amber_mmgbsa/bin/antechamber"
    antechamber_cmd = [
        antechamber_exe,
        "-i", "ligand_H.pdb",
        "-fi", "pdb",
        "-o", "lig_am1bcc.mol2",
        "-fo", "mol2",
        "-at", "gaff2",
        "-c", "bcc",      # AM1-BCC
        "-nc", "0",       # Net charge 0
        "-s", "2",
        "-rn", "LIG"
    ]

    try:
        result = subprocess.run(
            antechamber_cmd,
            capture_output=True,
            text=True,
            timeout=600,
            cwd=aa_dir
        )

        if result.returncode != 0:
            print(f"  ✗ FAILED antechamber")
            # Check sqm.out
            if os.path.exists(f"{aa_dir}/sqm.out"):
                with open(f"{aa_dir}/sqm.out") as f:
                    sqm_lines = f.readlines()
                    for line in sqm_lines[-10:]:
                        if "ERROR" in line or "QMMM:" in line:
                            print(f"     {line.strip()}")
            return f"{aa}: FAIL antechamber"

    except subprocess.TimeoutExpired:
        print(f"  ✗ TIMEOUT")
        return f"{aa}: FAIL timeout"

    # Step 3: parmchk2
    print(f"  4. Running parmchk2...")
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
        print(f"  ✗ Missing output files")
        return f"{aa}: FAIL missing output"

    print(f"  ✓ SUCCESS!")
    return f"{aa}: OK (AM1-BCC charges)"

def main():
    print("="*70)
    print("Parameterizing THR and SER with AM1-BCC charges")
    print("Method: obabel (add H) → antechamber (AM1-BCC, -nc 0)")
    print("="*70)
    print()

    results = []
    for aa in amino_acids:
        result = parameterize_am1bcc(aa)
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
        print()
        print("Files created:")
        print("  - ligand_H.pdb (hydrogenated input)")
        print("  - lig_am1bcc.mol2 (AM1-BCC charges)")
        print("  - lig_am1bcc.frcmod")
        print()
        print("Existing Gasteiger files:")
        print("  - lig.mol2 (Gasteiger charges)")
        print("  - lig.frcmod")
        print()
        print("Next: Run validation comparison on ~10 ThrRS+Zn models")
    else:
        print()
        print(f"✗ {len(amino_acids) - success} failed")

if __name__ == '__main__':
    main()
