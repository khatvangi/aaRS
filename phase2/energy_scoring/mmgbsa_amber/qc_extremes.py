#!/usr/bin/env python3
"""
QC Step 1: Analyze extreme binding energies for pathological contacts.
For top 20 most positive and bottom 20 most negative BE:
- Extract minimum ligand-protein heavy-atom distance
- Extract minimum ligand-Zn distance (if Zn present)
- Count ligand atoms within 2.0 Å of any protein heavy atom (clash count)
"""
import pandas as pd
import numpy as np
import os
import subprocess
import sys

def parse_pdb_coords(pdb_file, chain=None, resname=None, exclude_H=True):
    """
    Parse PDB and return coordinates for specified chain/resname.
    Returns: list of (atom_name, x, y, z, element)
    """
    coords = []
    with open(pdb_file) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue

            atom_chain = line[21:22].strip()
            atom_resname = line[17:20].strip()
            atom_name = line[12:16].strip()
            element = line[76:78].strip() if len(line) > 77 else atom_name[0]

            # Skip hydrogens if requested
            if exclude_H and element == 'H':
                continue

            # Filter by chain/resname if specified
            if chain and atom_chain != chain:
                continue
            if resname and atom_resname != resname:
                continue

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            coords.append((atom_name, x, y, z, element))

    return coords

def min_distance(coords1, coords2):
    """Calculate minimum distance between two sets of coordinates."""
    if not coords1 or not coords2:
        return np.nan

    pos1 = np.array([(x, y, z) for _, x, y, z, _ in coords1])
    pos2 = np.array([(x, y, z) for _, x, y, z, _ in coords2])

    # Pairwise distances
    dists = np.sqrt(((pos1[:, np.newaxis, :] - pos2[np.newaxis, :, :])**2).sum(axis=2))

    return dists.min()

def count_clashes(coords1, coords2, threshold=2.0):
    """Count atoms in coords1 within threshold of any atom in coords2."""
    if not coords1 or not coords2:
        return 0

    pos1 = np.array([(x, y, z) for _, x, y, z, _ in coords1])
    pos2 = np.array([(x, y, z) for _, x, y, z, _ in coords2])

    # For each atom in coords1, check if min distance to coords2 < threshold
    dists = np.sqrt(((pos1[:, np.newaxis, :] - pos2[np.newaxis, :, :])**2).sum(axis=2))
    min_dists = dists.min(axis=1)

    return (min_dists < threshold).sum()

def convert_rst_to_pdb(job_dir, prmtop="complex.prmtop", rst="min.rst", output="min_final.pdb"):
    """Convert Amber restart to PDB using cpptraj."""
    cpptraj_in = f"""trajin {rst}
trajout {output} pdb
"""
    cpptraj_in_file = os.path.join(job_dir, "cpptraj_qc.in")
    with open(cpptraj_in_file, 'w') as f:
        f.write(cpptraj_in)

    cpptraj_exe = "/home/kiran/miniforge3/envs/amber_mmgbsa/bin/cpptraj"
    result = subprocess.run(
        [cpptraj_exe, "-p", prmtop, "-i", "cpptraj_qc.in"],
        capture_output=True, text=True, cwd=job_dir
    )

    return result.returncode == 0

def analyze_job_contacts(job_dir, job_name, meta):
    """Analyze contacts for one job."""
    # Convert rst to PDB if needed
    min_pdb = os.path.join(job_dir, "min_final.pdb")
    if not os.path.exists(min_pdb):
        success = convert_rst_to_pdb(job_dir)
        if not success:
            return {
                'job_name': job_name,
                'min_lig_prot_dist': np.nan,
                'min_lig_zn_dist': np.nan,
                'clash_count_2A': np.nan,
                'error': 'cpptraj_failed'
            }

    # Parse minimized structure
    # Ligand: chain B or resname LIG
    lig_coords = parse_pdb_coords(min_pdb, resname='LIG', exclude_H=True)

    # Protein: chain A (exclude LIG and ZN)
    all_coords = parse_pdb_coords(min_pdb, exclude_H=True)
    prot_coords = [(n, x, y, z, e) for n, x, y, z, e in all_coords
                   if e not in ['H'] and n not in ['ZN']]

    # Zn: element Zn or resname ZN
    zn_coords = [(n, x, y, z, e) for n, x, y, z, e in all_coords if e == 'Zn' or n == 'ZN']

    # Remove ligand atoms from protein coords
    lig_atom_names = set([n for n, _, _, _, _ in lig_coords])
    prot_coords = [(n, x, y, z, e) for n, x, y, z, e in prot_coords if n not in lig_atom_names]

    # Calculate metrics
    min_lig_prot = min_distance(lig_coords, prot_coords)
    min_lig_zn = min_distance(lig_coords, zn_coords) if zn_coords else np.nan
    clash_count = count_clashes(lig_coords, prot_coords, threshold=2.0)

    return {
        'job_name': job_name,
        'min_lig_prot_dist': min_lig_prot,
        'min_lig_zn_dist': min_lig_zn,
        'clash_count_2A': clash_count
    }

def main():
    # Load results
    df = pd.read_csv('mmgbsa_results.csv')
    ok_df = df[df['status'] == 'OK'].copy()

    print(f"Total OK jobs: {len(ok_df)}")

    # Get top 20 most positive and bottom 20 most negative
    top20_pos = ok_df.nlargest(20, 'BE_dG_bind')
    top20_neg = ok_df.nsmallest(20, 'BE_dG_bind')

    extremes = pd.concat([top20_pos, top20_neg])

    print(f"\nAnalyzing {len(extremes)} extreme cases:")
    print(f"  Top 20 positive: {top20_pos['BE_dG_bind'].min():.2f} to {top20_pos['BE_dG_bind'].max():.2f}")
    print(f"  Top 20 negative: {top20_neg['BE_dG_bind'].min():.2f} to {top20_neg['BE_dG_bind'].max():.2f}")

    # Analyze each
    results = []
    for idx, row in extremes.iterrows():
        job_name = row['job_name']
        job_dir = os.path.join('jobs', job_name)

        print(f"  Analyzing {job_name}...")

        result = analyze_job_contacts(job_dir, job_name, row)
        result['BE_dG_bind'] = row['BE_dG_bind']
        result['ligand_resname'] = row['ligand_resname']
        result['zn_present'] = row.get('zn_present', 0)

        results.append(result)

    # Save results
    results_df = pd.DataFrame(results)
    results_df = results_df.merge(extremes[['job_name', 'BE_dG_bind']], on='job_name', suffixes=('', '_check'))
    results_df.drop('BE_dG_bind_check', axis=1, errors='ignore')
    results_df.to_csv('qc/extremes_diagnostics.csv', index=False)

    print(f"\n✓ Wrote qc/extremes_diagnostics.csv")

    # Summary statistics
    print(f"\n{'='*70}")
    print("Summary Statistics:")
    print(f"{'='*70}")
    print(f"Min ligand-protein distance (Å):")
    print(f"  Mean: {results_df['min_lig_prot_dist'].mean():.3f}")
    print(f"  Min:  {results_df['min_lig_prot_dist'].min():.3f}")
    print(f"  Max:  {results_df['min_lig_prot_dist'].max():.3f}")

    print(f"\nMin ligand-Zn distance (Å, where Zn present):")
    zn_subset = results_df[~results_df['min_lig_zn_dist'].isna()]
    if len(zn_subset) > 0:
        print(f"  Mean: {zn_subset['min_lig_zn_dist'].mean():.3f}")
        print(f"  Min:  {zn_subset['min_lig_zn_dist'].min():.3f}")
        print(f"  Max:  {zn_subset['min_lig_zn_dist'].max():.3f}")

    print(f"\nClash count (<2.0 Å):")
    print(f"  Mean: {results_df['clash_count_2A'].mean():.1f}")
    print(f"  Max:  {results_df['clash_count_2A'].max():.0f}")

    # Flag pathological cases
    pathological = results_df[
        (results_df['min_lig_prot_dist'] < 1.6) |
        (results_df['min_lig_zn_dist'] < 1.6)
    ]

    if len(pathological) > 0:
        print(f"\n⚠ WARNING: {len(pathological)} cases with distances < 1.6 Å:")
        print(pathological[['job_name', 'BE_dG_bind', 'min_lig_prot_dist', 'min_lig_zn_dist']])
    else:
        print(f"\n✓ No pathological contacts (< 1.6 Å) detected")

    # Create histograms
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(12, 4))

        # Min distance histogram
        axes[0].hist(results_df['min_lig_prot_dist'].dropna(), bins=20, edgecolor='black')
        axes[0].axvline(1.6, color='red', linestyle='--', label='Pathology threshold (1.6 Å)')
        axes[0].set_xlabel('Min Ligand-Protein Distance (Å)')
        axes[0].set_ylabel('Count')
        axes[0].set_title('Minimum Heavy-Atom Distances')
        axes[0].legend()

        # Clash count histogram
        axes[1].hist(results_df['clash_count_2A'].dropna(), bins=20, edgecolor='black')
        axes[1].set_xlabel('Number of Clashes (< 2.0 Å)')
        axes[1].set_ylabel('Count')
        axes[1].set_title('Clash Counts (Extreme Energy Cases)')

        plt.tight_layout()
        plt.savefig('qc/extremes_histograms.png', dpi=150)
        print(f"\n✓ Saved qc/extremes_histograms.png")

    except ImportError:
        print("\n(matplotlib not available - skipping plots)")

if __name__ == '__main__':
    main()
