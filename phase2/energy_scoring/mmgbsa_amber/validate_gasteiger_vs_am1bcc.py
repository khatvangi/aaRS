#!/usr/bin/env python3
"""
Validation: Compare Gasteiger vs AM1-BCC charges for THR and SER.

Run MM/GBSA on selected ThrRS+Zn jobs with both charge sets.
Compare ΔΔG rankings to ensure Gasteiger is trustworthy.

If rankings differ significantly, we must abort main batch and redo with AM1-BCC.
"""
import os
import sys
import json
import subprocess
import pandas as pd
import numpy as np

# Validation jobs: ThrRS+Zn with THR or SER
VALIDATION_JOBS = [
    "anc_thrrs_cat_zn_THR",
    "anc_thrrs_cat_zn_SER",
    "modern_thrrs_ecoli_zn_THR",
    "modern_thrrs_ecoli_zn_SER",
    "modern_thrrs_ecoli_THR_zinc",
    # Skip COMPETITION jobs for now - they have 2 ligands
]

JOBS_DIR = "jobs"
LIGANDS_DIR = "ligands"
VALIDATION_DIR = "validation_am1bcc"

def run_mmpbsa_with_params(job_name, ligand_mol2, ligand_frcmod, charge_method):
    """
    Run full MM/GBSA pipeline for one job with specified ligand parameters.
    Returns ΔG_bind or None if failed.
    """
    job_dir = os.path.join(JOBS_DIR, job_name)

    if not os.path.exists(job_dir):
        print(f"  ✗ {job_name}: job directory not found")
        return None

    # Check if minimization already done
    if not os.path.exists(os.path.join(job_dir, "min.rst")):
        print(f"  ✗ {job_name}: minimization not complete yet")
        return None

    # Create validation subdirectory
    val_dir = os.path.join(job_dir, f"validate_{charge_method}")
    os.makedirs(val_dir, exist_ok=True)

    # Copy ligand parameters
    subprocess.run(["cp", ligand_mol2, os.path.join(val_dir, "LIG.mol2")], check=True)
    subprocess.run(["cp", ligand_frcmod, os.path.join(val_dir, "LIG.frcmod")], check=True)
    subprocess.run(["cp", os.path.join(job_dir, "zn.mol2"), os.path.join(val_dir, "zn.mol2")], check=True)
    subprocess.run(["cp", os.path.join(job_dir, "zn.frcmod"), os.path.join(val_dir, "zn.frcmod")], check=True)
    subprocess.run(["cp", os.path.join(job_dir, "receptor.pdb"), val_dir], check=True)
    subprocess.run(["cp", os.path.join(job_dir, "ligand.pdb"), val_dir], check=True)
    subprocess.run(["cp", os.path.join(job_dir, "complex.pdb"), val_dir], check=True)

    # Create tleap input
    tleap_in = """source leaprc.protein.ff14SB
source leaprc.gaff2

loadAmberParams LIG.frcmod
LIG = loadMol2 LIG.mol2

loadAmberParams zn.frcmod
ZN = loadMol2 zn.mol2

REC = loadPdb receptor.pdb
L   = loadPdb ligand.pdb
COM = loadPdb complex.pdb

saveAmberParm REC receptor.prmtop receptor.inpcrd
saveAmberParm L   ligand.prmtop   ligand.inpcrd
saveAmberParm COM complex.prmtop  complex.inpcrd
quit
"""
    with open(os.path.join(val_dir, "tleap.in"), "w") as f:
        f.write(tleap_in)

    # Run tleap
    tleap_exe = "/home/kiran/miniforge3/envs/amber_mmgbsa/bin/tleap"
    result = subprocess.run([tleap_exe, "-f", "tleap.in"],
                          capture_output=True, text=True, cwd=val_dir)
    if result.returncode != 0 or not os.path.exists(os.path.join(val_dir, "complex.prmtop")):
        print(f"  ✗ {job_name} {charge_method}: tleap failed")
        return None

    # Use existing minimized structure (min.rst from parent job)
    # Create trajectory from it
    cpptraj_in = """trajin ../min.rst
trajout oneframe.nc netcdf
"""
    with open(os.path.join(val_dir, "cpptraj.in"), "w") as f:
        f.write(cpptraj_in)

    cpptraj_exe = "/home/kiran/miniforge3/envs/amber_mmgbsa/bin/cpptraj"
    result = subprocess.run([cpptraj_exe, "-p", "complex.prmtop", "-i", "cpptraj.in"],
                          capture_output=True, text=True, cwd=val_dir)
    if result.returncode != 0:
        print(f"  ✗ {job_name} {charge_method}: cpptraj failed")
        return None

    # Run MM/GBSA
    mmpbsa_in = """&general
  startframe=1, endframe=1, interval=1,
  verbose=1,
/
&gb
  igb=5, saltcon=0.150
/
"""
    with open(os.path.join(val_dir, "mmpbsa.in"), "w") as f:
        f.write(mmpbsa_in)

    mmpbsa_exe = "/home/kiran/miniforge3/envs/amber_mmgbsa/bin/MMPBSA.py"
    env = os.environ.copy()
    env['AMBERHOME'] = '/home/kiran/miniforge3/envs/amber_mmgbsa'
    result = subprocess.run([
        mmpbsa_exe, "-O",
        "-i", "mmpbsa.in",
        "-cp", "complex.prmtop",
        "-rp", "receptor.prmtop",
        "-lp", "ligand.prmtop",
        "-y", "oneframe.nc",
        "-o", "FINAL_RESULTS_MMPBSA.dat"
    ], capture_output=True, text=True, cwd=val_dir, env=env)

    if result.returncode != 0:
        print(f"  ✗ {job_name} {charge_method}: MMPBSA.py failed")
        return None

    # Parse results
    results_file = os.path.join(val_dir, "FINAL_RESULTS_MMPBSA.dat")
    if not os.path.exists(results_file):
        print(f"  ✗ {job_name} {charge_method}: no results file")
        return None

    with open(results_file) as f:
        lines = f.readlines()

    delta_total = None
    for line in lines:
        if "DELTA TOTAL" in line:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    delta_total = float(parts[2])
                except ValueError:
                    pass
            break

    if delta_total is None:
        print(f"  ✗ {job_name} {charge_method}: failed to parse ΔG")
        return None

    print(f"  ✓ {job_name} {charge_method}: ΔG = {delta_total:.2f} kcal/mol")
    return delta_total

def main():
    print("="*70)
    print("Validation: Gasteiger vs AM1-BCC charges")
    print("="*70)
    print()

    # Check if AM1-BCC params exist
    thr_am1bcc_mol2 = os.path.join(LIGANDS_DIR, "THR", "lig_am1bcc.mol2")
    thr_am1bcc_frcmod = os.path.join(LIGANDS_DIR, "THR", "lig_am1bcc.frcmod")
    ser_am1bcc_mol2 = os.path.join(LIGANDS_DIR, "SER", "lig_am1bcc.mol2")
    ser_am1bcc_frcmod = os.path.join(LIGANDS_DIR, "SER", "lig_am1bcc.frcmod")

    for f in [thr_am1bcc_mol2, thr_am1bcc_frcmod, ser_am1bcc_mol2, ser_am1bcc_frcmod]:
        if not os.path.exists(f):
            print(f"ERROR: Missing AM1-BCC parameter file: {f}")
            print("Run parameterize_am1bcc_correct.py first!")
            sys.exit(1)

    # Gasteiger params
    thr_gas_mol2 = os.path.join(LIGANDS_DIR, "THR", "lig.mol2")
    thr_gas_frcmod = os.path.join(LIGANDS_DIR, "THR", "lig.frcmod")
    ser_gas_mol2 = os.path.join(LIGANDS_DIR, "SER", "lig.mol2")
    ser_gas_frcmod = os.path.join(LIGANDS_DIR, "SER", "lig.frcmod")

    # Run validation
    results = []

    for job_name in VALIDATION_JOBS:
        # Determine ligand type
        meta_file = os.path.join(JOBS_DIR, job_name, "meta.json")
        if not os.path.exists(meta_file):
            print(f"Skipping {job_name}: no metadata")
            continue

        with open(meta_file) as f:
            meta = json.load(f)

        lig_resname = meta['ligand_original_resname']

        if lig_resname not in ["THR", "SER"]:
            print(f"Skipping {job_name}: ligand is {lig_resname}, not THR/SER")
            continue

        print(f"\n{job_name} ({lig_resname}):")

        # Gasteiger
        if lig_resname == "THR":
            dG_gas = run_mmpbsa_with_params(job_name, thr_gas_mol2, thr_gas_frcmod, "gasteiger")
        else:  # SER
            dG_gas = run_mmpbsa_with_params(job_name, ser_gas_mol2, ser_gas_frcmod, "gasteiger")

        # AM1-BCC
        if lig_resname == "THR":
            dG_am1bcc = run_mmpbsa_with_params(job_name, thr_am1bcc_mol2, thr_am1bcc_frcmod, "am1bcc")
        else:  # SER
            dG_am1bcc = run_mmpbsa_with_params(job_name, ser_am1bcc_mol2, ser_am1bcc_frcmod, "am1bcc")

        if dG_gas is not None and dG_am1bcc is not None:
            diff = dG_am1bcc - dG_gas
            results.append({
                'job_name': job_name,
                'ligand': lig_resname,
                'dG_gasteiger': dG_gas,
                'dG_am1bcc': dG_am1bcc,
                'difference': diff
            })

    # Summary
    print()
    print("="*70)
    print("Validation Results:")
    print("="*70)

    if len(results) == 0:
        print("No results - jobs may not be ready yet. Wait for main batch to complete.")
        sys.exit(0)

    df = pd.DataFrame(results)
    print(df.to_string(index=False))

    # Calculate ΔΔG
    print()
    print("="*70)
    print("ΔΔG Analysis (relative to THR):")
    print("="*70)

    # Group by ligand and average
    avg_gas = df.groupby('ligand')['dG_gasteiger'].mean()
    avg_am1bcc = df.groupby('ligand')['dG_am1bcc'].mean()

    if 'THR' in avg_gas.index and 'SER' in avg_gas.index:
        ddG_gas = avg_gas['SER'] - avg_gas['THR']
        ddG_am1bcc = avg_am1bcc['SER'] - avg_am1bcc['THR']

        print(f"\nGasteiger:")
        print(f"  THR (cognate):    {avg_gas['THR']:.2f} kcal/mol")
        print(f"  SER (similar):    {avg_gas['SER']:.2f} kcal/mol")
        print(f"  ΔΔG (SER vs THR): {ddG_gas:.2f} kcal/mol")

        print(f"\nAM1-BCC:")
        print(f"  THR (cognate):    {avg_am1bcc['THR']:.2f} kcal/mol")
        print(f"  SER (similar):    {avg_am1bcc['SER']:.2f} kcal/mol")
        print(f"  ΔΔG (SER vs THR): {ddG_am1bcc:.2f} kcal/mol")

        print()
        print("="*70)
        print("VERDICT:")
        print("="*70)

        # Check if sign matches
        if (ddG_gas > 0 and ddG_am1bcc > 0) or (ddG_gas < 0 and ddG_am1bcc < 0):
            print("✓ SIGN MATCHES - Gasteiger and AM1-BCC agree on ranking")
            print(f"  Both predict SER is {'worse' if ddG_gas > 0 else 'better'} than THR")

            # Check magnitude
            diff_pct = abs(ddG_am1bcc - ddG_gas) / abs(ddG_am1bcc) * 100 if ddG_am1bcc != 0 else 0
            print(f"  ΔΔG difference: {abs(ddG_am1bcc - ddG_gas):.2f} kcal/mol ({diff_pct:.1f}%)")

            if diff_pct < 30:
                print()
                print("✓ GASTEIGER IS ACCEPTABLE")
                print("  Proceed with current batch using Gasteiger charges.")
            else:
                print()
                print("⚠ WARNING: Large ΔΔG difference")
                print("  Consider re-running with AM1-BCC for higher accuracy.")
        else:
            print("✗ SIGN MISMATCH - Gasteiger and AM1-BCC DISAGREE on ranking!")
            print("  Gasteiger predicts SER is " + ('worse' if ddG_gas > 0 else 'better') + " than THR")
            print("  AM1-BCC predicts SER is " + ('worse' if ddG_am1bcc > 0 else 'better') + " than THR")
            print()
            print("✗ GASTEIGER IS UNRELIABLE")
            print("  ABORT current batch and re-run with AM1-BCC charges!")
    else:
        print("Insufficient data for ΔΔG comparison")

    # Save results
    df.to_csv("validation_gasteiger_vs_am1bcc.csv", index=False)
    print()
    print(f"Saved: validation_gasteiger_vs_am1bcc.csv")
    print("="*70)

if __name__ == '__main__':
    main()
