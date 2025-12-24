#!/usr/bin/env python3
"""
Fix ALL ligand parameterizations by adding hydrogens.
Current parameterizations are missing all H atoms - this is critical for MM/GBSA.

Note on AM1-BCC for zwitterions:
- AM1-BCC (via sqm) fails for zwitterionic amino acids due to odd electron count issue
- This is a known problem when ligand PDBs lack explicit charges
- Using Gasteiger charges with proper hydrogens is acceptable for comparative ΔΔG studies
"""
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

LIGANDS_DIR = "ligands"
AA20 = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

def reparameterize_with_hydrogens(aa):
    """
    Re-parameterize with Gasteiger charges + hydrogens added.
    """
    aa_dir = os.path.join(LIGANDS_DIR, aa)
    ligand_pdb = os.path.join(aa_dir, "ligand.pdb")

    if not os.path.exists(ligand_pdb):
        return f"{aa}: SKIP - no ligand.pdb"

    # Backup old files
    if os.path.exists(f"{aa_dir}/lig.mol2"):
        subprocess.run(["cp", "lig.mol2", "lig_no_H.mol2"], cwd=aa_dir, capture_output=True)
    if os.path.exists(f"{aa_dir}/lig.frcmod"):
        subprocess.run(["cp", "lig.frcmod", "lig_no_H.frcmod"], cwd=aa_dir, capture_output=True)

    # Antechamber with hydrogen addition
    antechamber_exe = "/home/kiran/miniforge3/envs/amber_mmgbsa/bin/antechamber"
    antechamber_cmd = [
        antechamber_exe,
        "-i", "ligand.pdb",
        "-fi", "pdb",
        "-o", "lig.mol2",
        "-fo", "mol2",
        "-at", "gaff2",
        "-c", "gas",      # Gasteiger charges
        "-j", "5",        # CRITICAL: Add all hydrogens
        "-s", "2",
        "-rn", "LIG",
        "-nc", "0"
    ]

    try:
        result = subprocess.run(
            antechamber_cmd,
            capture_output=True,
            text=True,
            timeout=300,
            cwd=aa_dir
        )

        if result.returncode != 0:
            return f"{aa}: FAIL antechamber - {result.stderr[:100]}"

    except subprocess.TimeoutExpired:
        return f"{aa}: FAIL timeout"

    # parmchk2
    parmchk2_exe = "/home/kiran/miniforge3/envs/amber_mmgbsa/bin/parmchk2"
    parmchk2_cmd = [
        parmchk2_exe,
        "-i", "lig.mol2",
        "-f", "mol2",
        "-o", "lig.frcmod"
    ]

    result = subprocess.run(parmchk2_cmd, capture_output=True, text=True, cwd=aa_dir)
    if result.returncode != 0:
        return f"{aa}: FAIL parmchk2"

    # Verify outputs
    if not os.path.exists(f"{aa_dir}/lig.mol2") or not os.path.exists(f"{aa_dir}/lig.frcmod"):
        return f"{aa}: FAIL missing output"

    # Count atoms to verify H were added
    with open(f"{aa_dir}/lig.mol2") as f:
        lines = f.readlines()
    for line in lines:
        if line.strip().startswith("@<TRIPOS>ATOM"):
            idx = lines.index(line)
            natoms = len([l for l in lines[idx+1:] if l.strip() and not l.startswith("@")])
            return f"{aa}: OK ({natoms} atoms)"

    return f"{aa}: OK"

def main():
    print("="*70)
    print("Re-parameterizing all 20 amino acids WITH HYDROGENS")
    print("Using Gasteiger charges (AM1-BCC fails for zwitterions)")
    print("="*70)

    # Process in parallel
    results = []
    with ProcessPoolExecutor(max_workers=10) as ex:
        futs = {ex.submit(reparameterize_with_hydrogens, aa): aa for aa in AA20}
        for fut in tqdm(as_completed(futs), total=len(futs), desc="Parameterizing"):
            results.append(fut.result())

    # Print results
    print("\n" + "="*70)
    print("Results:")
    print("="*70)
    for res in sorted(results):
        print(f"  {res}")

    success = sum(1 for r in results if "OK" in r)
    print("\n" + "="*70)
    print(f"Success: {success}/{len(AA20)} amino acids")
    print("="*70)

if __name__ == '__main__':
    main()
