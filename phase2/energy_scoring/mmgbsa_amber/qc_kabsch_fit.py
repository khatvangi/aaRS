#!/usr/bin/env python3
"""
QC Step 2: Verify Kabsch fit integrity.
For all 171 successful jobs, compute:
- RMSD between AF3 heavy atoms (ligand.pdb) and lig_pose.mol2 heavy atoms
- Max per-atom deviation

Pass criteria:
- RMSD median < 0.05 Å
- Max deviation < 0.2 Å
"""
import pandas as pd
import numpy as np
import os
import sys

def parse_pdb_heavy_atoms(pdb_file):
    """Parse PDB and return {atom_name: (x, y, z)} for heavy atoms."""
    atoms = {}
    with open(pdb_file) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue

            atom_name = line[12:16].strip()
            element = line[76:78].strip() if len(line) > 77 else atom_name[0]

            # Skip hydrogens
            if element == 'H' or atom_name.startswith('H') or atom_name.startswith('h'):
                continue

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            atoms[atom_name] = (x, y, z)

    return atoms

def parse_mol2_heavy_atoms(mol2_file):
    """Parse mol2 and return {atom_name: (x, y, z)} for heavy atoms."""
    atoms = {}
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

                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4])

                    atoms[atom_name] = (x, y, z)

    return atoms

def compute_rmsd_and_max_dev(pdb_atoms, mol2_atoms):
    """
    Compute RMSD and max deviation between matched atoms.
    Returns: (rmsd, max_dev, n_matched)
    """
    # Match by atom name
    matched_names = set(pdb_atoms.keys()) & set(mol2_atoms.keys())

    if not matched_names:
        return np.nan, np.nan, 0

    pdb_coords = np.array([pdb_atoms[name] for name in sorted(matched_names)])
    mol2_coords = np.array([mol2_atoms[name] for name in sorted(matched_names)])

    # Compute per-atom distances
    diffs = pdb_coords - mol2_coords
    distances = np.sqrt((diffs**2).sum(axis=1))

    rmsd = np.sqrt((distances**2).mean())
    max_dev = distances.max()

    return rmsd, max_dev, len(matched_names)

def analyze_job_fit(job_dir, job_name):
    """Analyze Kabsch fit quality for one job."""
    pdb_file = os.path.join(job_dir, 'ligand.pdb')
    mol2_file = os.path.join(job_dir, 'lig_pose.mol2')

    if not os.path.exists(pdb_file):
        return {
            'job_name': job_name,
            'rmsd': np.nan,
            'max_deviation': np.nan,
            'n_atoms_matched': 0,
            'error': 'ligand.pdb missing'
        }

    if not os.path.exists(mol2_file):
        return {
            'job_name': job_name,
            'rmsd': np.nan,
            'max_deviation': np.nan,
            'n_atoms_matched': 0,
            'error': 'lig_pose.mol2 missing'
        }

    # Parse atoms
    pdb_atoms = parse_pdb_heavy_atoms(pdb_file)
    mol2_atoms = parse_mol2_heavy_atoms(mol2_file)

    # Compute RMSD
    rmsd, max_dev, n_matched = compute_rmsd_and_max_dev(pdb_atoms, mol2_atoms)

    return {
        'job_name': job_name,
        'rmsd': rmsd,
        'max_deviation': max_dev,
        'n_atoms_matched': n_matched
    }

def main():
    # Load results
    df = pd.read_csv('mmgbsa_results.csv')
    ok_df = df[df['status'] == 'OK'].copy()

    print(f"Total OK jobs: {len(ok_df)}")
    print(f"\nAnalyzing Kabsch fit quality for all {len(ok_df)} jobs...")

    # Analyze each job
    results = []
    for idx, row in ok_df.iterrows():
        job_name = row['job_name']
        job_dir = os.path.join('jobs', job_name)

        result = analyze_job_fit(job_dir, job_name)
        result['ligand_resname'] = row['ligand_resname']
        result['BE_dG_bind'] = row['BE_dG_bind']

        results.append(result)

        if (len(results) % 50 == 0):
            print(f"  Processed {len(results)}/{len(ok_df)}...")

    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    results_df.to_csv('qc/kabsch_fit_qc.csv', index=False)

    print(f"\n✓ Wrote qc/kabsch_fit_qc.csv")

    # Summary statistics
    valid = results_df[~results_df['rmsd'].isna()]

    print(f"\n{'='*70}")
    print("Summary Statistics:")
    print(f"{'='*70}")
    print(f"Valid fits: {len(valid)}/{len(results_df)}")

    if len(valid) > 0:
        print(f"\nRMSD (Å):")
        print(f"  Median: {valid['rmsd'].median():.6f}")
        print(f"  Mean:   {valid['rmsd'].mean():.6f}")
        print(f"  Std:    {valid['rmsd'].std():.6f}")
        print(f"  Min:    {valid['rmsd'].min():.6f}")
        print(f"  Max:    {valid['rmsd'].max():.6f}")

        print(f"\nMax deviation (Å):")
        print(f"  Median: {valid['max_deviation'].median():.6f}")
        print(f"  Mean:   {valid['max_deviation'].mean():.6f}")
        print(f"  Max:    {valid['max_deviation'].max():.6f}")

        print(f"\nAtoms matched:")
        print(f"  Mean: {valid['n_atoms_matched'].mean():.1f}")
        print(f"  Min:  {valid['n_atoms_matched'].min():.0f}")
        print(f"  Max:  {valid['n_atoms_matched'].max():.0f}")

        # Pass/Fail criteria
        median_rmsd = valid['rmsd'].median()
        max_max_dev = valid['max_deviation'].max()

        print(f"\n{'='*70}")
        print("Pass/Fail Criteria:")
        print(f"{'='*70}")

        if median_rmsd < 0.05:
            print(f"✓ PASS: Median RMSD = {median_rmsd:.6f} Å < 0.05 Å")
        else:
            print(f"✗ FAIL: Median RMSD = {median_rmsd:.6f} Å ≥ 0.05 Å")

        if max_max_dev < 0.2:
            print(f"✓ PASS: Max deviation = {max_max_dev:.6f} Å < 0.2 Å")
        else:
            print(f"✗ FAIL: Max deviation = {max_max_dev:.6f} Å ≥ 0.2 Å")

        # Flag outliers
        outliers = valid[valid['rmsd'] > 0.5]
        if len(outliers) > 0:
            print(f"\n⚠ WARNING: {len(outliers)} jobs with RMSD > 0.5 Å:")
            print(outliers[['job_name', 'ligand_resname', 'rmsd', 'max_deviation']].head(10))

    # Create histogram
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(12, 4))

        # RMSD histogram
        axes[0].hist(valid['rmsd'], bins=50, edgecolor='black')
        axes[0].axvline(0.05, color='red', linestyle='--', label='Target median (0.05 Å)')
        axes[0].axvline(valid['rmsd'].median(), color='green', linestyle='-', label=f'Actual median ({valid["rmsd"].median():.4f} Å)')
        axes[0].set_xlabel('RMSD (Å)')
        axes[0].set_ylabel('Count')
        axes[0].set_title('Kabsch Fit RMSD Distribution')
        axes[0].legend()
        axes[0].set_xlim(left=0)

        # Max deviation histogram
        axes[1].hist(valid['max_deviation'], bins=50, edgecolor='black')
        axes[1].axvline(0.2, color='red', linestyle='--', label='Target max (0.2 Å)')
        axes[1].set_xlabel('Max Deviation (Å)')
        axes[1].set_ylabel('Count')
        axes[1].set_title('Max Per-Atom Deviation Distribution')
        axes[1].legend()
        axes[1].set_xlim(left=0)

        plt.tight_layout()
        plt.savefig('qc/kabsch_fit_distribution.png', dpi=150)
        print(f"\n✓ Saved qc/kabsch_fit_distribution.png")

    except ImportError:
        print("\n(matplotlib not available - skipping plots)")

if __name__ == '__main__':
    main()
