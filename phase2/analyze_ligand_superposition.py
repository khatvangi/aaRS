#!/usr/bin/env python3
"""
Ligand superposition analysis: Do all ligands bind in the same spot?
Or do cognate ligands cluster separately from non-cognate?

Approach:
1. For each enzyme variant (job_name base)
2. Load all structures with different ligands
3. Align protein backbones (CA atoms)
4. Extract ligand heavy-atom coordinates
5. Compute pairwise ligand RMSD
6. Test: Do cognate ligands cluster separately?
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

def parse_cif_atoms(cif_file, chain_id):
    """Extract atoms from a specific chain in CIF file."""
    atoms = []

    with open(cif_file) as f:
        lines = f.readlines()

    in_atom_site = False
    for line in lines:
        if line.startswith('_atom_site.'):
            in_atom_site = True
            continue
        if line.startswith('#') and in_atom_site:
            in_atom_site = False
            continue

        if in_atom_site and (line.startswith('ATOM') or line.startswith('HETATM')):
            parts = line.split()
            if len(parts) < 13:
                continue

            atom_chain = parts[6]
            if atom_chain != chain_id:
                continue

            atom_name = parts[3]
            resname = parts[5]
            resnum_str = parts[8]
            x = float(parts[10])
            y = float(parts[11])
            z = float(parts[12])
            element = parts[2] if len(parts) > 2 else atom_name[0]

            # Skip hydrogens
            if element == 'H' or atom_name.startswith('H'):
                continue

            # Handle missing resnum
            try:
                resnum = int(resnum_str)
            except ValueError:
                resnum = 1  # default for ligands

            atoms.append({
                'atom_name': atom_name,
                'resname': resname,
                'resnum': resnum,
                'element': element,
                'coords': np.array([x, y, z])
            })

    return atoms

def get_ca_atoms(atoms):
    """Extract CA atoms for alignment."""
    ca_atoms = [a for a in atoms if a['atom_name'] == 'CA']
    return ca_atoms

def kabsch_align(coords_mobile, coords_ref):
    """
    Kabsch algorithm for optimal rotation.
    Returns rotation matrix and translation vector.
    """
    # Center coordinates
    centroid_mobile = np.mean(coords_mobile, axis=0)
    centroid_ref = np.mean(coords_ref, axis=0)

    mobile_centered = coords_mobile - centroid_mobile
    ref_centered = coords_ref - centroid_ref

    # Compute rotation matrix
    H = mobile_centered.T @ ref_centered
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T

    # Ensure right-handed coordinate system
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    # Translation vector
    t = centroid_ref - R @ centroid_mobile

    return R, t

def apply_transform(coords, R, t):
    """Apply rotation and translation to coordinates."""
    return (R @ coords.T).T + t

def compute_ligand_rmsd(lig1_atoms, lig2_atoms):
    """
    Compute RMSD between two ligands (already aligned).
    Match atoms by name.
    """
    # Match atoms by name
    atoms1_dict = {a['atom_name']: a['coords'] for a in lig1_atoms}
    atoms2_dict = {a['atom_name']: a['coords'] for a in lig2_atoms}

    common_atoms = set(atoms1_dict.keys()) & set(atoms2_dict.keys())

    if len(common_atoms) < 3:
        return np.nan

    coords1 = np.array([atoms1_dict[name] for name in sorted(common_atoms)])
    coords2 = np.array([atoms2_dict[name] for name in sorted(common_atoms)])

    diff = coords1 - coords2
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

    return rmsd

def analyze_enzyme_variant(structures_df, base_name, protein_chain='A'):
    """
    Analyze ligand superposition for one enzyme variant.
    structures_df: rows with same enzyme variant but different ligands
    """
    if len(structures_df) < 2:
        return None

    print(f"\nAnalyzing {base_name}: {len(structures_df)} ligands")

    # Load reference structure (first one)
    ref_row = structures_df.iloc[0]
    ref_file = ref_row['file']
    ref_protein_chain = ref_row['protein_chain'] if 'protein_chain' in ref_row else protein_chain
    ref_ligand_chain = ref_row['ligand_chain'] if 'ligand_chain' in ref_row else 'B'

    print(f"  Reference: {ref_row['ligand']} from {os.path.basename(ref_file)}")

    # Load reference protein CA atoms
    ref_protein = parse_cif_atoms(ref_file, ref_protein_chain)
    ref_ca = get_ca_atoms(ref_protein)

    if len(ref_ca) < 50:
        print(f"  ERROR: Too few CA atoms ({len(ref_ca)})")
        return None

    ref_ca_coords = np.array([a['coords'] for a in ref_ca])

    # Load reference ligand
    ref_ligand = parse_cif_atoms(ref_file, ref_ligand_chain)

    # Store aligned ligands
    aligned_ligands = {}
    aligned_ligands[ref_row['ligand']] = ref_ligand

    # Align all other structures
    for idx, row in structures_df.iloc[1:].iterrows():
        ligand_name = row['ligand']
        cif_file = row['file']
        prot_chain = row['protein_chain'] if 'protein_chain' in row else protein_chain
        lig_chain = row['ligand_chain'] if 'ligand_chain' in row else 'B'

        # Load protein CA atoms
        protein_atoms = parse_cif_atoms(cif_file, prot_chain)
        ca_atoms = get_ca_atoms(protein_atoms)

        if len(ca_atoms) != len(ref_ca):
            print(f"  WARNING: {ligand_name} has {len(ca_atoms)} CA atoms (ref: {len(ref_ca)})")
            continue

        ca_coords = np.array([a['coords'] for a in ca_atoms])

        # Compute alignment
        R, t = kabsch_align(ca_coords, ref_ca_coords)

        # Load and align ligand
        ligand_atoms = parse_cif_atoms(cif_file, lig_chain)

        # Apply transformation to ligand
        aligned_lig = []
        for atom in ligand_atoms:
            new_coords = apply_transform(atom['coords'].reshape(1, -1), R, t)[0]
            aligned_lig.append({
                **atom,
                'coords': new_coords
            })

        aligned_ligands[ligand_name] = aligned_lig

    print(f"  Successfully aligned {len(aligned_ligands)} ligands")

    # Compute pairwise RMSD matrix
    ligand_names = list(aligned_ligands.keys())
    n = len(ligand_names)
    rmsd_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i+1, n):
            rmsd = compute_ligand_rmsd(
                aligned_ligands[ligand_names[i]],
                aligned_ligands[ligand_names[j]]
            )
            rmsd_matrix[i, j] = rmsd
            rmsd_matrix[j, i] = rmsd

    return {
        'base_name': base_name,
        'ligand_names': ligand_names,
        'rmsd_matrix': rmsd_matrix,
        'aligned_ligands': aligned_ligands,
        'n_ligands': len(ligand_names)
    }

def identify_cognate(base_name):
    """Identify cognate ligand for each enzyme."""
    base_lower = base_name.lower()

    if 'prors' in base_lower or 'prours' in base_lower:
        return 'PRO'
    elif 'thrrs' in base_lower:
        return 'THR'
    else:
        return None

def analyze_specificity(result):
    """
    Test if cognate ligands cluster separately from non-cognate.
    """
    base_name = result['base_name']
    cognate = identify_cognate(base_name)

    if cognate is None or cognate not in result['ligand_names']:
        return None

    ligand_names = result['ligand_names']
    rmsd_matrix = result['rmsd_matrix']

    cognate_idx = ligand_names.index(cognate)

    # RMSD from cognate to cognate (should be 0)
    # RMSD from cognate to non-cognate
    cognate_to_noncognate = []

    for i, lig in enumerate(ligand_names):
        if i != cognate_idx:
            cognate_to_noncognate.append(rmsd_matrix[cognate_idx, i])

    # RMSD among non-cognate ligands
    noncognate_to_noncognate = []
    for i in range(len(ligand_names)):
        if i == cognate_idx:
            continue
        for j in range(i+1, len(ligand_names)):
            if j == cognate_idx:
                continue
            noncognate_to_noncognate.append(rmsd_matrix[i, j])

    mean_cognate_to_noncognate = np.mean(cognate_to_noncognate) if cognate_to_noncognate else np.nan
    mean_noncognate_to_noncognate = np.mean(noncognate_to_noncognate) if noncognate_to_noncognate else np.nan

    return {
        'base_name': base_name,
        'cognate': cognate,
        'n_ligands': len(ligand_names),
        'mean_cognate_to_noncognate_rmsd': mean_cognate_to_noncognate,
        'mean_noncognate_to_noncognate_rmsd': mean_noncognate_to_noncognate,
        'cognate_more_distinct': mean_cognate_to_noncognate > mean_noncognate_to_noncognate
    }

def main():
    # Load geometry metrics to get structure list
    df = pd.read_csv('geometry_metrics.csv')
    print(f"Loaded {len(df)} structures")

    # Load manifest for chain info
    manifest = pd.read_csv('manifest.csv')

    # Merge to get chain assignments
    df = df.merge(
        manifest[['file', 'protein_chain', 'ligand_chain']].drop_duplicates(),
        on='file',
        how='left'
    )

    # Filter to main enzymes and sufficient data
    df_main = df[df['enzyme'].isin(['ProRS', 'ThrRS'])].copy()

    # Group by base structure (enzyme variant)
    # Extract base name (remove ligand-specific suffix)
    def get_base_name(job_name):
        """Extract base structure name (enzyme variant without ligand)."""
        # Remove trailing _XXX (amino acid codes)
        # Common patterns: modern_prors_ALA, anc_thrrs_cat_THR
        parts = job_name.split('_')

        # Check if last part is 3-letter AA code
        aa_codes = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                   'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                   'THR', 'TRP', 'TYR', 'VAL']

        if parts[-1] in aa_codes:
            return '_'.join(parts[:-1])
        else:
            return job_name

    df_main['base_name'] = df_main['job_name'].apply(get_base_name)

    # Count ligands per base structure
    base_counts = df_main.groupby('base_name').size()

    # Filter to bases with multiple ligands
    bases_multi = base_counts[base_counts >= 3].index.tolist()

    print(f"\nFound {len(bases_multi)} enzyme variants with ≥3 ligands")

    # Analyze each base structure
    results = []
    specificity_results = []

    for base_name in bases_multi:
        structures = df_main[df_main['base_name'] == base_name].copy()

        # Analyze superposition
        result = analyze_enzyme_variant(structures, base_name)

        if result:
            results.append(result)

            # Analyze specificity
            spec = analyze_specificity(result)
            if spec:
                specificity_results.append(spec)

    # Save results
    print("\n" + "="*70)
    print("LIGAND SUPERPOSITION ANALYSIS SUMMARY")
    print("="*70)

    # Overall statistics
    all_rmsds = []
    for r in results:
        matrix = r['rmsd_matrix']
        # Get upper triangle (excluding diagonal)
        rmsds = matrix[np.triu_indices_from(matrix, k=1)]
        all_rmsds.extend(rmsds[~np.isnan(rmsds)])

    print(f"\nOverall ligand RMSD statistics (all pairs):")
    print(f"  Mean:   {np.mean(all_rmsds):.2f} Å")
    print(f"  Median: {np.median(all_rmsds):.2f} Å")
    print(f"  Std:    {np.std(all_rmsds):.2f} Å")
    print(f"  Min:    {np.min(all_rmsds):.2f} Å")
    print(f"  Max:    {np.max(all_rmsds):.2f} Å")

    # Specificity analysis
    if specificity_results:
        print("\n" + "="*70)
        print("COGNATE vs NON-COGNATE CLUSTERING")
        print("="*70)

        spec_df = pd.DataFrame(specificity_results)
        print(spec_df.to_string(index=False))

        # Summary
        print(f"\nStructures where cognate is MORE distinct from non-cognate:")
        print(f"  {spec_df['cognate_more_distinct'].sum()}/{len(spec_df)} ({100*spec_df['cognate_more_distinct'].sum()/len(spec_df):.1f}%)")

        spec_df.to_csv('ligand_specificity_clustering.csv', index=False)
        print(f"\nSaved ligand_specificity_clustering.csv")

    # Create visualization: RMSD distribution
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Plot 1: Overall RMSD distribution
    ax = axes[0]
    ax.hist(all_rmsds, bins=50, alpha=0.7, edgecolor='black')
    ax.axvline(np.median(all_rmsds), color='red', linestyle='--', linewidth=2, label=f'Median = {np.median(all_rmsds):.2f} Å')
    ax.set_xlabel('Ligand-ligand RMSD (Å)', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title('Distribution of Ligand RMSD\n(all pairwise comparisons)', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)

    # Plot 2: Cognate vs non-cognate
    if specificity_results:
        ax = axes[1]
        spec_df = pd.DataFrame(specificity_results)

        cognate_vals = spec_df['mean_cognate_to_noncognate_rmsd'].dropna()
        noncognate_vals = spec_df['mean_noncognate_to_noncognate_rmsd'].dropna()

        positions = [1, 2]
        bp = ax.boxplot([cognate_vals, noncognate_vals], positions=positions,
                        widths=0.6, patch_artist=True)

        bp['boxes'][0].set_facecolor('lightcoral')
        bp['boxes'][1].set_facecolor('lightblue')

        ax.set_xticks(positions)
        ax.set_xticklabels(['Cognate to\nNon-cognate', 'Non-cognate to\nNon-cognate'])
        ax.set_ylabel('Mean RMSD (Å)', fontsize=12)
        ax.set_title('Cognate Ligand Clustering', fontsize=14, fontweight='bold')
        ax.grid(alpha=0.3, axis='y')

        # Add significance test
        from scipy import stats
        if len(cognate_vals) > 1 and len(noncognate_vals) > 1:
            t_stat, p_val = stats.ttest_ind(cognate_vals, noncognate_vals)
            ax.text(0.5, 0.95, f'p = {p_val:.3f}', transform=ax.transAxes,
                   ha='center', va='top', fontsize=10)

    plt.tight_layout()
    plt.savefig('ligand_superposition_analysis.png', dpi=300, bbox_inches='tight')
    print("\nSaved ligand_superposition_analysis.png")
    plt.close()

    # Create heatmap for a few selected structures
    print("\n" + "="*70)
    print("Creating RMSD heatmaps for selected structures...")
    print("="*70)

    # Select top structures with most ligands
    top_structures = sorted(results, key=lambda x: x['n_ligands'], reverse=True)[:4]

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()

    for idx, result in enumerate(top_structures):
        ax = axes[idx]

        matrix = result['rmsd_matrix']
        ligands = result['ligand_names']

        # Plot heatmap
        im = ax.imshow(matrix, cmap='viridis', aspect='auto')
        ax.set_xticks(range(len(ligands)))
        ax.set_yticks(range(len(ligands)))
        ax.set_xticklabels(ligands, rotation=45, ha='right', fontsize=8)
        ax.set_yticklabels(ligands, fontsize=8)

        # Add colorbar
        plt.colorbar(im, ax=ax, label='RMSD (Å)')

        # Highlight cognate
        cognate = identify_cognate(result['base_name'])
        if cognate and cognate in ligands:
            cognate_idx = ligands.index(cognate)
            # Draw box around cognate row/column
            ax.axhline(cognate_idx - 0.5, color='red', linewidth=2)
            ax.axhline(cognate_idx + 0.5, color='red', linewidth=2)
            ax.axvline(cognate_idx - 0.5, color='red', linewidth=2)
            ax.axvline(cognate_idx + 0.5, color='red', linewidth=2)

        ax.set_title(f"{result['base_name']}\n({result['n_ligands']} ligands)",
                    fontsize=11, fontweight='bold')

    plt.tight_layout()
    plt.savefig('ligand_rmsd_heatmaps.png', dpi=300, bbox_inches='tight')
    print("Saved ligand_rmsd_heatmaps.png")
    plt.close()

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print("\nGenerated files:")
    print("  - ligand_specificity_clustering.csv")
    print("  - ligand_superposition_analysis.png")
    print("  - ligand_rmsd_heatmaps.png")

if __name__ == '__main__':
    main()
