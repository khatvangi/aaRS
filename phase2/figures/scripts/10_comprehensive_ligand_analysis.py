#!/usr/bin/env python3
"""
Comprehensive Ligand Interaction Analysis
Analyzes ALL amino acid ligands in ALL AF3 predictions to understand:
- Ligand-Zn coordination distances
- Ligand position in active site
- Protein-ligand contacts
- H-bond networks

This reveals the molecular mechanism behind ipTM scores.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.vectors import Vector

# Distance thresholds
ZN_COORD_CUTOFF = 3.0  # Å - coordination bond
CONTACT_CUTOFF = 4.0   # Å - van der Waals contact
HBOND_CUTOFF = 3.5     # Å - hydrogen bond

def find_all_cif_files():
    """Find all CIF files from AF3 predictions."""
    base_dir = Path('/storage/kiran-stuff/aaRS/phase2')
    cif_files = []

    # Search in outputs directory
    for cif_path in base_dir.glob('outputs/**/*_model.cif'):
        cif_files.append(cif_path)

    # Search in af3_modern_matrix directory
    for cif_path in base_dir.glob('af3_modern_matrix/**/*_model.cif'):
        cif_files.append(cif_path)

    # Search in other potential locations
    for cif_path in base_dir.glob('**/*_model.cif'):
        if cif_path not in cif_files:
            cif_files.append(cif_path)

    return sorted(list(set(cif_files)))


def identify_chains(structure):
    """Identify protein, ligand, and Zn chains."""
    chains = {}

    for model in structure:
        for chain in model:
            chain_id = chain.id
            residues = list(chain.get_residues())

            if len(residues) == 0:
                continue

            # Check first residue
            first_res = residues[0]
            resname = first_res.resname.strip()

            # Classify chain type
            if len(residues) == 1 and len(resname) == 3 and resname.isupper():
                # Single residue with 3-letter code = amino acid ligand
                chains['ligand'] = chain_id
                chains['ligand_name'] = resname
            elif resname == 'ZN':
                chains['zn'] = chain_id
            elif len(residues) > 10:
                # Long chain = protein
                if 'protein' not in chains:
                    chains['protein'] = chain_id

    return chains


def calculate_zn_distances(structure, chains):
    """Calculate distances from ligand to Zn ion."""
    if 'ligand' not in chains or 'zn' not in chains:
        return {}

    results = {}

    for model in structure:
        ligand_chain = model[chains['ligand']]
        zn_chain = model[chains['zn']]

        # Get Zn coordinates
        zn_atoms = [atom for atom in zn_chain.get_atoms()]
        if len(zn_atoms) == 0:
            return {}
        zn_coord = zn_atoms[0].coord

        # Calculate distances from each ligand atom to Zn
        for residue in ligand_chain:
            for atom in residue:
                atom_coord = atom.coord
                distance = np.linalg.norm(atom_coord - zn_coord)

                atom_name = atom.get_name()
                results[f'zn_to_{atom_name}'] = float(distance)

        # Find minimum distance
        if results:
            results['zn_min_distance'] = min(results.values())

            # Identify coordinating atoms (< 3.0 Å)
            coordinating = [k for k, v in results.items()
                          if k.startswith('zn_to_') and v < ZN_COORD_CUTOFF]
            results['coordinating_atoms'] = [k.replace('zn_to_', '')
                                            for k in coordinating]

    return results


def calculate_active_site_distance(structure, chains):
    """Calculate distance from ligand to active site center."""
    if 'ligand' not in chains or 'protein' not in chains:
        return {}

    results = {}

    for model in structure:
        protein_chain = model[chains['protein']]
        ligand_chain = model[chains['ligand']]

        # Calculate ligand centroid
        ligand_coords = []
        for residue in ligand_chain:
            for atom in residue:
                ligand_coords.append(atom.coord)

        if len(ligand_coords) == 0:
            return {}

        ligand_centroid = np.mean(ligand_coords, axis=0)

        # Find nearest protein residues to define active site
        protein_atoms = [atom for atom in protein_chain.get_atoms()]
        if len(protein_atoms) == 0:
            return {}

        # Find 10 closest protein atoms to ligand
        distances = []
        for atom in protein_atoms:
            dist = np.linalg.norm(atom.coord - ligand_centroid)
            distances.append(dist)

        distances = sorted(distances)
        results['nearest_protein_distance'] = float(distances[0]) if distances else None
        results['avg_10_nearest_distance'] = float(np.mean(distances[:10])) if len(distances) >= 10 else None
        results['ligand_centroid_depth'] = float(np.mean(distances[:10])) if len(distances) >= 10 else None

    return results


def count_contacts(structure, chains):
    """Count protein atoms within contact distance of ligand."""
    if 'ligand' not in chains or 'protein' not in chains:
        return {}

    results = {}

    for model in structure:
        protein_chain = model[chains['protein']]
        ligand_chain = model[chains['ligand']]

        # Get all atoms
        protein_atoms = [atom for atom in protein_chain.get_atoms()]
        ligand_atoms = [atom for atom in ligand_chain.get_atoms()]

        if len(protein_atoms) == 0 or len(ligand_atoms) == 0:
            return {}

        # Use NeighborSearch for efficiency
        ns = NeighborSearch(protein_atoms)

        # Count contacts within 4Å
        contact_atoms = set()
        contact_residues = set()

        for ligand_atom in ligand_atoms:
            nearby = ns.search(ligand_atom.coord, CONTACT_CUTOFF, level='A')
            contact_atoms.update(nearby)

            # Get residues
            for atom in nearby:
                contact_residues.add(atom.get_parent())

        results['contact_atoms_4A'] = len(contact_atoms)
        results['contact_residues_4A'] = len(contact_residues)

        # Also count closer contacts (3Å)
        contact_atoms_3A = set()
        for ligand_atom in ligand_atoms:
            nearby = ns.search(ligand_atom.coord, 3.0, level='A')
            contact_atoms_3A.update(nearby)

        results['contact_atoms_3A'] = len(contact_atoms_3A)

        # Identify contacting residues
        residue_names = []
        for residue in contact_residues:
            res_id = f"{residue.resname}{residue.id[1]}"
            residue_names.append(res_id)

        results['contacting_residues'] = sorted(residue_names)[:20]  # Top 20

    return results


def find_hbonds(structure, chains):
    """Find hydrogen bonds between ligand and protein."""
    if 'ligand' not in chains or 'protein' not in chains:
        return {}

    results = {}
    hbonds = []

    for model in structure:
        protein_chain = model[chains['protein']]
        ligand_chain = model[chains['ligand']]

        # Get potential H-bond donors/acceptors
        protein_atoms = [atom for atom in protein_chain.get_atoms()
                        if atom.element in ['O', 'N']]
        ligand_atoms = [atom for atom in ligand_chain.get_atoms()
                       if atom.element in ['O', 'N']]

        if len(protein_atoms) == 0 or len(ligand_atoms) == 0:
            return {}

        # Find H-bonds (simple distance-based)
        for ligand_atom in ligand_atoms:
            for protein_atom in protein_atoms:
                distance = ligand_atom - protein_atom

                if distance < HBOND_CUTOFF:
                    parent_res = protein_atom.get_parent()
                    hbond_info = {
                        'ligand_atom': ligand_atom.get_name(),
                        'protein_atom': protein_atom.get_name(),
                        'protein_residue': f"{parent_res.resname}{parent_res.id[1]}",
                        'distance': float(distance)
                    }
                    hbonds.append(hbond_info)

        results['hbond_count'] = len(hbonds)
        results['hbonds'] = hbonds

        if hbonds:
            results['avg_hbond_distance'] = float(np.mean([h['distance'] for h in hbonds]))
            results['min_hbond_distance'] = float(min([h['distance'] for h in hbonds]))

    return results


def analyze_cif_file(cif_path):
    """Comprehensive analysis of a single CIF file."""
    parser = MMCIFParser(QUIET=True)

    try:
        structure = parser.get_structure('structure', str(cif_path))
    except Exception as e:
        return {'error': str(e)}

    # Identify chains
    chains = identify_chains(structure)

    if 'ligand' not in chains:
        return {'error': 'No ligand found'}

    results = {
        'cif_file': str(cif_path),
        'job_name': cif_path.parent.name,
        'ligand': chains.get('ligand_name', 'UNK'),
        'has_zn': 'zn' in chains,
    }

    # Zn coordination analysis
    if 'zn' in chains:
        zn_data = calculate_zn_distances(structure, chains)
        results.update(zn_data)

    # Active site position
    position_data = calculate_active_site_distance(structure, chains)
    results.update(position_data)

    # Contact analysis
    contact_data = count_contacts(structure, chains)
    results.update(contact_data)

    # H-bond analysis
    hbond_data = find_hbonds(structure, chains)
    results.update(hbond_data)

    return results


def main():
    print("="*80)
    print("COMPREHENSIVE LIGAND INTERACTION ANALYSIS")
    print("="*80)

    # Find all CIF files
    print("\n1. Finding CIF files...")
    cif_files = find_all_cif_files()
    print(f"   Found {len(cif_files)} CIF files")

    # Analyze each file
    print("\n2. Analyzing structures...")
    all_results = []

    for i, cif_path in enumerate(cif_files, 1):
        if i % 10 == 0:
            print(f"   Processed {i}/{len(cif_files)}...")

        result = analyze_cif_file(cif_path)
        all_results.append(result)

    # Convert to DataFrame
    df = pd.DataFrame(all_results)

    # Load AF3 scores to merge
    print("\n3. Merging with AF3 scores...")
    af3_df = pd.read_csv('AF3_RESULTS_CORRECTED.csv')

    # Merge on job_name
    merged_df = df.merge(
        af3_df[['job_name', 'ptm', 'AA_iptm', 'protein_len', 'ligands']],
        on='job_name',
        how='left'
    )

    # Save results
    output_file = 'figures/data/comprehensive_ligand_analysis.csv'
    merged_df.to_csv(output_file, index=False)
    print(f"   Saved: {output_file}")

    # Save detailed JSON (includes hbond lists)
    json_file = 'figures/data/comprehensive_ligand_analysis.json'
    with open(json_file, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"   Saved: {json_file}")

    # Generate summary statistics
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)

    # Filter for Zn-containing structures
    zn_df = merged_df[merged_df['has_zn'] == True].copy()

    if len(zn_df) > 0:
        print(f"\n{len(zn_df)} structures with Zn²⁺:")

        # Group by ligand
        ligand_summary = zn_df.groupby('ligand').agg({
            'zn_min_distance': 'mean',
            'contact_atoms_4A': 'mean',
            'hbond_count': 'mean',
            'AA_iptm': 'mean'
        }).round(3)

        ligand_summary = ligand_summary.sort_values('AA_iptm', ascending=False)
        print("\nLigand Summary (Zn-containing structures):")
        print(ligand_summary.to_string())

        # Highlight key ligands
        print("\n" + "="*80)
        print("KEY COMPARISONS")
        print("="*80)

        key_ligands = ['THR', 'SER', 'ILE', 'VAL', 'PRO', 'ALA']
        for ligand in key_ligands:
            ligand_data = zn_df[zn_df['ligand'] == ligand]
            if len(ligand_data) > 0:
                row = ligand_data.iloc[0]
                print(f"\n{ligand}:")
                print(f"  ipTM: {row.get('AA_iptm', 'N/A'):.3f}")
                print(f"  Zn min distance: {row.get('zn_min_distance', 'N/A'):.2f} Å")
                print(f"  Coordinating atoms: {row.get('coordinating_atoms', [])}")
                print(f"  Contact atoms (4Å): {row.get('contact_atoms_4A', 0):.0f}")
                print(f"  H-bonds: {row.get('hbond_count', 0):.0f}")

        # Identify coordination vs non-coordination
        print("\n" + "="*80)
        print("COORDINATION ANALYSIS")
        print("="*80)

        zn_df['is_coordinating'] = zn_df['zn_min_distance'] < ZN_COORD_CUTOFF

        coordinating = zn_df[zn_df['is_coordinating'] == True]
        non_coordinating = zn_df[zn_df['is_coordinating'] == False]

        print(f"\nCoordinating ligands (< {ZN_COORD_CUTOFF} Å): {len(coordinating)}")
        if len(coordinating) > 0:
            print(coordinating[['ligand', 'zn_min_distance', 'AA_iptm', 'contact_atoms_4A']].to_string(index=False))

        print(f"\nNon-coordinating ligands (≥ {ZN_COORD_CUTOFF} Å): {len(non_coordinating)}")
        if len(non_coordinating) > 0:
            print(non_coordinating[['ligand', 'zn_min_distance', 'AA_iptm', 'contact_atoms_4A']].to_string(index=False))

    # Summary for all structures
    print("\n" + "="*80)
    print("ALL STRUCTURES SUMMARY")
    print("="*80)

    print(f"\nTotal structures analyzed: {len(merged_df)}")
    print(f"With Zn: {len(merged_df[merged_df['has_zn'] == True])}")
    print(f"Without Zn: {len(merged_df[merged_df['has_zn'] == False])}")
    print(f"Errors: {len(merged_df[merged_df.get('error').notna()])}")

    print("\n" + "="*80)
    print("COMPLETE!")
    print("="*80)


if __name__ == '__main__':
    main()
