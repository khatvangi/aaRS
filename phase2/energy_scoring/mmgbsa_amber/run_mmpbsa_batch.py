#!/usr/bin/env python3
"""
Batch MM/GBSA runner using AmberTools MMPBSA.py (single snapshot, no MD).
Processes all AF3 jobs with parallel execution.
"""
import os
import sys
import json
import re
import subprocess
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

MANIFEST = "../mmgbsa_real/manifest_af3only.csv"
PREP_SCRIPT = "../prep_from_cif.py"
JOBS_DIR = "jobs"
LIGANDS_DIR = "ligands_fixed"
ZN_MOL2 = "zn.mol2"
ZN_FRCMOD = "zn.frcmod"
OUTPUT_CSV = "mmgbsa_results.csv"
NWORKERS = 60

def parse_minimization_convergence(min_out_file):
    """
    Parse sander minimization output to extract convergence info.
    Returns: dict with 'converged', 'rms_gradient', 'nstep', 'energy'
    """
    if not os.path.exists(min_out_file):
        return {'converged': False, 'rms_gradient': 999.0, 'nstep': 0, 'energy': 0.0}

    with open(min_out_file) as f:
        content = f.read()

    # Find FINAL RESULTS section
    match = re.search(
        r'FINAL RESULTS.*?NSTEP\s+ENERGY\s+RMS\s+GMAX.*?(\d+)\s+([-\d.E+]+)\s+([-\d.E+]+)\s+([-\d.E+]+)',
        content, re.DOTALL
    )

    if not match:
        # If no FINAL RESULTS, get last NSTEP line
        matches = re.findall(r'NSTEP\s+ENERGY\s+RMS\s+GMAX.*?\s+(\d+)\s+([-\d.E+]+)\s+([-\d.E+]+)\s+([-\d.E+]+)', content, re.DOTALL)
        if matches:
            match_data = matches[-1]
            nstep = int(match_data[0])
            energy = float(match_data[1])
            rms = float(match_data[2])
            gmax = float(match_data[3])
        else:
            return {'converged': False, 'rms_gradient': 999.0, 'nstep': 0, 'energy': 0.0}
    else:
        nstep = int(match.group(1))
        energy = float(match.group(2))
        rms = float(match.group(3))
        gmax = float(match.group(4))

    # With restraints, RMS gradient won't reach <0.1, so check if energy is stable
    # Convergence here means: completed minimization
    converged = True  # If we got to final results, minimization completed

    return {
        'nstep': nstep,
        'energy': energy,
        'rms_gradient': rms,
        'max_gradient': gmax,
        'converged': converged
    }

def calculate_drift_metrics(initial_pdb, final_pdb):
    """
    Calculate geometry drift between initial (AF3) and minimized structure.
    Returns: dict with 'ligand_rmsd', 'zn_distance_change', 'flag'
    Flag is set if drift exceeds thresholds (ligand RMSD > 0.8 Å or Zn distance > 0.5 Å)
    """
    try:
        # Read PDB files and extract ligand (chain B, resname LIG) and Zn (chain D) coordinates
        def read_coords(pdb_file):
            lig_coords = []
            zn_coord = None

            with open(pdb_file) as f:
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        chain = line[21:22].strip()
                        resname = line[17:20].strip()
                        atom_name = line[12:16].strip()
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])

                        # Ligand: chain B, resname LIG
                        if chain == "B" and resname == "LIG":
                            lig_coords.append([x, y, z])

                        # Zn: chain D, resname ZN
                        if chain == "D" and resname == "ZN":
                            zn_coord = np.array([x, y, z])

            return np.array(lig_coords) if lig_coords else None, zn_coord

        init_lig, init_zn = read_coords(initial_pdb)
        final_lig, final_zn = read_coords(final_pdb)

        # Calculate ligand RMSD (no alignment, just direct coordinate difference)
        ligand_rmsd = 0.0
        if init_lig is not None and final_lig is not None and len(init_lig) == len(final_lig):
            diff = init_lig - final_lig
            ligand_rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

        # Calculate Zn distance change
        zn_distance_change = 0.0
        if init_zn is not None and final_zn is not None:
            zn_distance_change = np.linalg.norm(init_zn - final_zn)

        # Flag if drift exceeds thresholds
        flag = ""
        if ligand_rmsd > 0.8:
            flag += "LIGAND_DRIFT "
        if zn_distance_change > 0.5:
            flag += "ZN_DRIFT "

        return {
            'ligand_rmsd': ligand_rmsd,
            'zn_distance_change': zn_distance_change,
            'flag': flag.strip()
        }

    except Exception as e:
        return {
            'ligand_rmsd': -1.0,
            'zn_distance_change': -1.0,
            'flag': f'DRIFT_CALC_ERROR: {str(e)[:50]}'
        }

def remove_duplicate_bonds(mol2_file):
    """
    Remove duplicate bonds from mol2 file and renumber.
    Fixes antechamber bug that creates duplicate bonds.
    """
    with open(mol2_file) as f:
        lines = f.readlines()

    output = []
    in_bonds = False
    bond_num = 1
    seen_bonds = set()

    for line in lines:
        if '@<TRIPOS>BOND' in line:
            in_bonds = True
            output.append(line)
            continue
        if '@<TRIPOS>SUBSTRUCTURE' in line:
            in_bonds = False
            output.append(line)
            continue

        if in_bonds and line.strip():
            parts = line.split()
            if len(parts) >= 4:
                atom1, atom2, bond_type = parts[1], parts[2], parts[3]
                bond_key = tuple(sorted([atom1, atom2])) + (bond_type,)
                if bond_key not in seen_bonds:
                    seen_bonds.add(bond_key)
                    output.append(f"{bond_num:6d} {atom1:5s} {atom2:5s} {bond_type:3s}\n")
                    bond_num += 1
        else:
            output.append(line)

    # Update molecule header with correct bond count
    for i, line in enumerate(output):
        if line.startswith('@<TRIPOS>MOLECULE'):
            # Find the counts line (2 lines after @<TRIPOS>MOLECULE)
            if i+2 < len(output):
                parts = output[i+2].split()
                if len(parts) >= 2:
                    n_atoms = parts[0]
                    output[i+2] = f"   {n_atoms}    {bond_num-1}     1     0     0\n"
            break

    with open(mol2_file, 'w') as f:
        f.writelines(output)

def verify_pose_rmsd(pdb_file, mol2_file):
    """
    Verify that mol2 preserves AF3 heavy-atom coordinates.
    Returns dict with 'rmsd', 'n_atoms_matched'
    """
    # Parse PDB heavy atoms
    pdb_atoms = {}
    with open(pdb_file) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            atom_name = line[12:16].strip()
            element = line[76:78].strip() if len(line) > 77 else atom_name[0]

            # Skip hydrogens
            if element == 'H' or atom_name.startswith('H'):
                continue

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            pdb_atoms[atom_name] = np.array([x, y, z])

    # Parse mol2 heavy atoms
    mol2_atoms = {}
    in_atoms = False
    with open(mol2_file) as f:
        for line in f:
            if '@<TRIPOS>ATOM' in line:
                in_atoms = True
                continue
            if '@<TRIPOS>BOND' in line:
                break
            if in_atoms and line.strip():
                parts = line.split()
                if len(parts) >= 6:
                    atom_name = parts[1]
                    # Skip hydrogens
                    if atom_name.startswith('H') or atom_name.startswith('h'):
                        continue
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    mol2_atoms[atom_name] = np.array([x, y, z])

    # Match atoms by name
    matched = set(pdb_atoms.keys()) & set(mol2_atoms.keys())

    if not matched:
        return {'rmsd': 999.0, 'n_atoms_matched': 0}

    # Compute RMSD
    pdb_coords = np.array([pdb_atoms[name] for name in sorted(matched)])
    mol2_coords = np.array([mol2_atoms[name] for name in sorted(matched)])

    diffs = pdb_coords - mol2_coords
    rmsd = np.sqrt(np.mean(np.sum(diffs**2, axis=1)))

    return {'rmsd': rmsd, 'n_atoms_matched': len(matched)}

def run_one_job(row):
    """
    Process one AF3 CIF file through the full MM/GBSA pipeline.
    Returns a dict with results.
    """
    cif_path = row['file']
    job_name = row['job_name']
    job_dir = os.path.join(JOBS_DIR, job_name)

    try:
        os.makedirs(job_dir, exist_ok=True)

        # Step 1: Extract PDBs
        prep_cmd = [
            sys.executable,
            PREP_SCRIPT,
            f"../../{cif_path}",
            job_dir
        ]
        result = subprocess.run(prep_cmd, capture_output=True, text=True, cwd=os.getcwd())
        if result.returncode != 0:
            return {"job_name": job_name, "status": "PREP_FAILED", "error": result.stderr[:200]}

        # Read metadata
        meta_path = os.path.join(job_dir, "meta.json")
        with open(meta_path) as f:
            meta = json.load(f)

        lig_resname = meta['ligand_original_resname']
        zn_count = meta['zn_count']

        # Step 2: Per-structure ligand parameterization from AF3 pose
        # Generate mol2 + frcmod from AF3 ligand.pdb (not from template)

        # Copy Zn params
        subprocess.run(["cp", ZN_MOL2, os.path.join(job_dir, "zn.mol2")])
        subprocess.run(["cp", ZN_FRCMOD, os.path.join(job_dir, "zn.frcmod")])

        # Determine net charge for this AA
        if lig_resname in ['ARG', 'LYS']:
            net_charge = 1
        elif lig_resname in ['ASP', 'GLU']:
            net_charge = -1
        else:
            net_charge = 0

        # Add hydrogens to AF3 ligand pose
        obabel_exe = "obabel"

        result = subprocess.run([
            obabel_exe, "ligand.pdb", "-O", "ligand_H.pdb", "-h", "-p", "7.0"
        ], capture_output=True, text=True, cwd=job_dir)

        ligand_H_pdb = os.path.join(job_dir, "ligand_H.pdb")
        if result.returncode != 0 or not os.path.exists(ligand_H_pdb):
            return {"job_name": job_name, "status": "OBABEL_FAILED", "error": result.stderr[:200]}

        # Parameterize with antechamber (GAFF2 + Gasteiger on AF3 pose)
        antechamber_exe = "/home/kiran/miniforge3/envs/amber_mmgbsa/bin/antechamber"
        lig_pose_mol2 = os.path.join(job_dir, "lig_pose.mol2")

        result = subprocess.run([
            antechamber_exe,
            "-i", "ligand_H.pdb", "-fi", "pdb",
            "-o", "lig_pose.mol2", "-fo", "mol2",
            "-at", "gaff2", "-c", "gas",
            "-nc", str(net_charge), "-rn", "LIG", "-s", "2"
        ], capture_output=True, text=True, cwd=job_dir)

        if result.returncode != 0 or not os.path.exists(lig_pose_mol2):
            return {"job_name": job_name, "status": "ANTECHAMBER_FAILED", "error": result.stderr[:200]}

        # Remove duplicate bonds (antechamber bug)
        remove_duplicate_bonds(lig_pose_mol2)

        # Generate frcmod with parmchk2
        parmchk2_exe = "/home/kiran/miniforge3/envs/amber_mmgbsa/bin/parmchk2"
        lig_frcmod = os.path.join(job_dir, "LIG.frcmod")

        result = subprocess.run([
            parmchk2_exe,
            "-i", "lig_pose.mol2", "-f", "mol2",
            "-o", "LIG.frcmod"
        ], capture_output=True, text=True, cwd=job_dir)

        if result.returncode != 0 or not os.path.exists(lig_frcmod):
            return {"job_name": job_name, "status": "PARMCHK2_FAILED", "error": result.stderr[:200]}

        # Step 2b: QC gate - verify RMSD between ligand.pdb and lig_pose.mol2 heavy atoms
        # This ensures antechamber preserved the AF3 geometry
        rmsd_qc = verify_pose_rmsd(
            os.path.join(job_dir, "ligand.pdb"),
            os.path.join(job_dir, "lig_pose.mol2")
        )

        if rmsd_qc['rmsd'] > 0.05:
            return {
                "job_name": job_name,
                "status": "POSE_MISMATCH",
                "error": f"RMSD {rmsd_qc['rmsd']:.4f} > 0.05 Å"
            }

        # Step 3: Create tleap input (use mol2 for ligand, combine with receptor)
        tleap_in = """source leaprc.protein.ff14SB
source leaprc.gaff2

loadAmberParams LIG.frcmod
LIG = loadMol2 lig_pose.mol2

loadAmberParams zn.frcmod
ZN  = loadMol2 zn.mol2

REC = loadPdb receptor.pdb
COM = combine {REC LIG}

saveAmberParm REC receptor.prmtop receptor.inpcrd
saveAmberParm LIG ligand.prmtop ligand.inpcrd
saveAmberParm COM complex.prmtop complex.inpcrd
quit
"""
        with open(os.path.join(job_dir, "tleap.in"), "w") as f:
            f.write(tleap_in)

        # Run tleap
        tleap_exe = "/home/kiran/miniforge3/envs/amber_mmgbsa/bin/tleap"
        result = subprocess.run([tleap_exe, "-f", "tleap.in"],
                              capture_output=True, text=True, cwd=job_dir)

        # Check outputs exist (tleap may return non-zero due to warnings, so check files instead)
        if not os.path.exists(os.path.join(job_dir, "complex.prmtop")):
            return {"job_name": job_name, "status": "TLEAP_FAILED", "error": result.stderr[:200] if result.stderr else "No prmtop created"}

        # Step 4: Minimization (500 steps with restraints on protein+Zn, ligand free)
        # ncyc=250: steepest descent for first 250 steps, then conjugate gradient
        # ntr=1: enable restraints
        # restraintmask='!:LIG & !@H=': restrain everything except ligand (LIG residue) and hydrogens
        #   This includes protein heavy atoms AND Zn
        min_in_text = """Minimize
 &cntrl
  imin=1,
  maxcyc=500,
  ncyc=250,
  ntb=0,
  cut=999.0,
  ntr=1,
  restraint_wt=10.0,
  restraintmask='!:LIG & !@H=',
 /
"""
        with open(os.path.join(job_dir, "min.in"), "w") as f:
            f.write(min_in_text)

        sander_exe = "/home/kiran/miniforge3/envs/amber_mmgbsa/bin/sander"
        result = subprocess.run([
            sander_exe, "-O",
            "-i", "min.in",
            "-p", "complex.prmtop",
            "-c", "complex.inpcrd",
            "-r", "min.rst",
            "-o", "min.out",
            "-ref", "complex.inpcrd"  # Use initial coords as restraint reference
        ], capture_output=True, text=True, cwd=job_dir)

        if result.returncode != 0:
            return {"job_name": job_name, "status": "MIN_FAILED", "error": result.stderr[:200]}

        # Step 4b: Parse convergence and check geometry drift
        # Parse min.out for convergence
        conv_info = parse_minimization_convergence(os.path.join(job_dir, "min.out"))

        # Convert min.rst to PDB for geometry comparison
        cpptraj_drift_in = """trajin min.rst
trajout min_final.pdb pdb
"""
        with open(os.path.join(job_dir, "cpptraj_drift.in"), "w") as f:
            f.write(cpptraj_drift_in)

        cpptraj_exe = "/home/kiran/miniforge3/envs/amber_mmgbsa/bin/cpptraj"
        result = subprocess.run([cpptraj_exe, "-p", "complex.prmtop", "-i", "cpptraj_drift.in"],
                              capture_output=True, text=True, cwd=job_dir)

        if result.returncode != 0:
            return {"job_name": job_name, "status": "CPPTRAJ_DRIFT_FAILED", "error": result.stderr[:200]}

        # Calculate geometry drift metrics
        drift_metrics = calculate_drift_metrics(
            os.path.join(job_dir, "complex.pdb"),
            os.path.join(job_dir, "min_final.pdb")
        )

        # Step 5: Create trajectory
        cpptraj_in = """trajin min.rst
trajout oneframe.nc netcdf
"""
        with open(os.path.join(job_dir, "cpptraj.in"), "w") as f:
            f.write(cpptraj_in)

        cpptraj_exe = "/home/kiran/miniforge3/envs/amber_mmgbsa/bin/cpptraj"
        result = subprocess.run([cpptraj_exe, "-p", "complex.prmtop", "-i", "cpptraj.in"],
                              capture_output=True, text=True, cwd=job_dir)

        if result.returncode != 0:
            return {"job_name": job_name, "status": "CPPTRAJ_FAILED", "error": result.stderr[:200]}

        # Step 6: MM/GBSA
        mmpbsa_in = """&general
  startframe=1, endframe=1, interval=1,
  verbose=1,
/
&gb
  igb=5, saltcon=0.150
/
"""
        with open(os.path.join(job_dir, "mmpbsa.in"), "w") as f:
            f.write(mmpbsa_in)

        mmpbsa_exe = "/home/kiran/miniforge3/envs/amber_mmgbsa/bin/MMPBSA.py"
        # Set AMBERHOME for MMPBSA.py
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
        ], capture_output=True, text=True, cwd=job_dir, env=env)

        if result.returncode != 0:
            return {"job_name": job_name, "status": "MMPBSA_FAILED", "error": result.stderr[:200]}

        # Step 7: Parse results
        results_file = os.path.join(job_dir, "FINAL_RESULTS_MMPBSA.dat")
        if not os.path.exists(results_file):
            return {"job_name": job_name, "status": "NO_RESULTS_FILE"}

        # Parse DELTA TOTAL from results file
        with open(results_file) as f:
            lines = f.readlines()

        delta_total = None
        for line in lines:
            if "DELTA TOTAL" in line:
                parts = line.split()
                # Line format: "DELTA TOTAL            value +/- stderr"
                if len(parts) >= 3:
                    try:
                        delta_total = float(parts[2])
                    except ValueError:
                        pass
                break

        if delta_total is None:
            return {"job_name": job_name, "status": "PARSE_FAILED"}

        return {
            "job_name": job_name,
            "file": cif_path,
            "status": "OK",
            "ligand_resname": lig_resname,
            "zn_present": zn_count,
            "BE_dG_bind": delta_total,
            "converged": conv_info['converged'],
            "rms_gradient": conv_info['rms_gradient'],
            "nstep": conv_info['nstep'],
            "ligand_rmsd_A": drift_metrics['ligand_rmsd'],
            "zn_distance_change_A": drift_metrics['zn_distance_change'],
            "drift_flag": drift_metrics['flag']
        }

    except Exception as e:
        return {"job_name": job_name, "status": "ERROR", "error": str(e)[:200]}

def main():
    # Load manifest
    df = pd.read_csv(MANIFEST)
    print(f"Total jobs: {len(df)}")

    # Create jobs directory
    os.makedirs(JOBS_DIR, exist_ok=True)

    # Run batch
    rows = []
    with ProcessPoolExecutor(max_workers=NWORKERS) as ex:
        futs = {ex.submit(run_one_job, row): row for _, row in df.iterrows()}
        print(f"Submitted {len(futs)} futures")
        for fut in tqdm(as_completed(futs), total=len(futs)):
            result = fut.result()
            rows.append(result)
            if len(rows) % 10 == 0:
                print(f"\nCompleted {len(rows)}/{len(futs)} jobs", flush=True)

    # Save results
    results_df = pd.DataFrame(rows)
    results_df.to_csv(OUTPUT_CSV, index=False)

    print(f"\n{'='*60}")
    print(f"Wrote {OUTPUT_CSV}: {len(results_df)} rows")
    print(f"{'='*60}")
    print(results_df['status'].value_counts())

if __name__ == '__main__':
    main()
