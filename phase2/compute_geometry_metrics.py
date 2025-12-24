#!/usr/bin/env python3
"""
Compute pure geometry/contact metrics for AF3 structures.
No force fields, no MM/GBSA - just physical geometry that AF3 predicts.
"""
import os
import sys
import json
import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Element data for atom type classification
POLAR_ATOMS = {'N', 'O', 'S'}
NONPOLAR_ATOMS = {'C'}
HBOND_DONORS = {'N', 'O'}  # simplified
HBOND_ACCEPTORS = {'N', 'O', 'S'}

def parse_cif_atoms(cif_file):
    """Parse CIF file and extract atom coordinates by chain."""
    atoms_by_chain = {}

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

        if in_atom_site and line.startswith('ATOM') or line.startswith('HETATM'):
            parts = line.split()
            if len(parts) < 13:
                continue

            atom_type = parts[0]
            atom_id = parts[1]
            atom_name = parts[3]
            resname = parts[5]
            chain_id = parts[6]
            resnum = parts[8]
            x = float(parts[10])
            y = float(parts[11])
            z = float(parts[12])
            element = parts[2] if len(parts) > 2 else atom_name[0]

            # Skip hydrogens
            if element == 'H' or atom_name.startswith('H'):
                continue

            if chain_id not in atoms_by_chain:
                atoms_by_chain[chain_id] = []

            atoms_by_chain[chain_id].append({
                'atom_name': atom_name,
                'resname': resname,
                'resnum': resnum,
                'element': element,
                'coords': np.array([x, y, z])
            })

    return atoms_by_chain

def compute_distance(coord1, coord2):
    """Euclidean distance between two coordinates."""
    return np.linalg.norm(coord1 - coord2)

def is_hbond(donor_atom, acceptor_atom, donor_neighbors):
    """
    Check if two atoms form a hydrogen bond.
    Simplified: D-H...A distance < 3.5 Å
    """
    dist = compute_distance(donor_atom['coords'], acceptor_atom['coords'])

    # Simple distance criterion (proper would include angle)
    if dist < 3.5:
        # Check if donor element can donate
        if donor_atom['element'] in HBOND_DONORS and acceptor_atom['element'] in HBOND_ACCEPTORS:
            return True
    return False

def compute_contacts(protein_atoms, ligand_atoms, cutoff=4.0):
    """Count contacts within cutoff distance."""
    contacts = 0
    for p_atom in protein_atoms:
        for l_atom in ligand_atoms:
            dist = compute_distance(p_atom['coords'], l_atom['coords'])
            if dist <= cutoff:
                contacts += 1
    return contacts

def compute_clashes(protein_atoms, ligand_atoms, clash_cutoff=2.2):
    """Count steric clashes (contacts < clash_cutoff)."""
    clashes = 0
    for p_atom in protein_atoms:
        for l_atom in ligand_atoms:
            dist = compute_distance(p_atom['coords'], l_atom['coords'])
            if dist < clash_cutoff:
                clashes += 1
    return clashes

def compute_polar_close_contacts(protein_atoms, ligand_atoms):
    """
    Count polar close contacts between protein and ligand.

    NOTE: This is NOT true H-bonding (no angle check, no explicit hydrogens).
    Counts N/O atom pairs within 3.5 Å as a proxy for polar interactions.

    Distance cutoff: 3.5 Å
    """
    polar_contacts = 0

    # Protein polar atoms to ligand polar atoms
    for p_atom in protein_atoms:
        if p_atom['element'] in HBOND_DONORS:
            for l_atom in ligand_atoms:
                if l_atom['element'] in HBOND_ACCEPTORS:
                    dist = compute_distance(p_atom['coords'], l_atom['coords'])
                    if dist < 3.5:
                        polar_contacts += 1

    # Ligand polar atoms to protein polar atoms
    for l_atom in ligand_atoms:
        if l_atom['element'] in HBOND_DONORS:
            for p_atom in protein_atoms:
                if p_atom['element'] in HBOND_ACCEPTORS:
                    dist = compute_distance(l_atom['coords'], p_atom['coords'])
                    if dist < 3.5:
                        polar_contacts += 1

    return polar_contacts

def compute_contact_types(protein_atoms, ligand_atoms, cutoff=4.0):
    """Count polar vs nonpolar contacts."""
    polar_contacts = 0
    nonpolar_contacts = 0

    for p_atom in protein_atoms:
        for l_atom in ligand_atoms:
            dist = compute_distance(p_atom['coords'], l_atom['coords'])
            if dist <= cutoff:
                # Classify contact
                p_polar = p_atom['element'] in POLAR_ATOMS
                l_polar = l_atom['element'] in POLAR_ATOMS

                if p_polar or l_polar:
                    polar_contacts += 1
                else:
                    nonpolar_contacts += 1

    return polar_contacts, nonpolar_contacts

def compute_zn_metrics(protein_atoms, ligand_atoms, zn_atoms):
    """
    Compute Zn coordination metrics.
    Returns: min_zn_ligand_dist, zn_coordination_count, zn_ligand_contacts
    """
    if not zn_atoms:
        return None, 0, 0

    min_dist = 999.9
    coord_count = 0  # atoms within 2.4 Å of Zn
    ligand_contacts = 0  # ligand atoms near Zn

    for zn_atom in zn_atoms:
        zn_coord = zn_atom['coords']

        # Min distance to ligand
        for l_atom in ligand_atoms:
            dist = compute_distance(zn_coord, l_atom['coords'])
            if dist < min_dist:
                min_dist = dist
            if dist < 3.5:  # ligand engaged with Zn
                ligand_contacts += 1

        # Coordination count (protein + ligand atoms within 2.4 Å)
        for p_atom in protein_atoms:
            if p_atom['element'] in {'N', 'O', 'S'}:  # coordinating atoms
                dist = compute_distance(zn_coord, p_atom['coords'])
                if dist < 2.4:
                    coord_count += 1

        for l_atom in ligand_atoms:
            if l_atom['element'] in {'N', 'O', 'S'}:
                dist = compute_distance(zn_coord, l_atom['coords'])
                if dist < 2.4:
                    coord_count += 1

    return min_dist if min_dist < 999 else None, coord_count, ligand_contacts

def classify_structure(job_name):
    """
    Classify structure by enzyme, epoch, domain.
    Returns: enzyme, epoch, domain
    """
    name_lower = job_name.lower()

    # Enzyme
    if 'prors' in name_lower or 'prours' in name_lower:
        enzyme = 'ProRS'
    elif 'thrrs' in name_lower:
        enzyme = 'ThrRS'
    else:
        enzyme = 'Other'

    # Epoch
    if 'anc' in name_lower:
        epoch = 'Ancestral'
    elif 'modern' in name_lower:
        epoch = 'Modern'
    elif 'deep' in name_lower:
        epoch = 'Deep'
    else:
        epoch = 'Unknown'

    # Domain
    if 'cat' in name_lower and 'catalytic' not in name_lower:
        domain = 'Catalytic'
    elif 'edit' in name_lower:
        domain = 'Editing'
    elif 'full' in name_lower or 'ecoli' in name_lower or 'human' in name_lower:
        domain = 'Full-length'
    else:
        domain = 'Domain'

    return enzyme, epoch, domain

def process_structure(cif_file, protein_chain, ligand_chain, zn_chain=None):
    """Process one CIF structure and compute all geometry metrics."""
    try:
        atoms_by_chain = parse_cif_atoms(cif_file)

        if protein_chain not in atoms_by_chain or ligand_chain not in atoms_by_chain:
            return None

        protein_atoms = atoms_by_chain[protein_chain]
        ligand_atoms = atoms_by_chain[ligand_chain]

        # Get Zn atoms if present
        zn_atoms = []
        if zn_chain and zn_chain in atoms_by_chain:
            zn_atoms = [a for a in atoms_by_chain[zn_chain] if a['element'] == 'ZN']

        # Compute all metrics
        contacts_4A = compute_contacts(protein_atoms, ligand_atoms, cutoff=4.0)
        clashes = compute_clashes(protein_atoms, ligand_atoms, clash_cutoff=2.2)
        polar_close_contacts_3p5A = compute_polar_close_contacts(protein_atoms, ligand_atoms)
        polar_contacts, nonpolar_contacts = compute_contact_types(protein_atoms, ligand_atoms)

        # Zn metrics
        zn_min_dist, zn_coord_count, zn_lig_contacts = compute_zn_metrics(
            protein_atoms, ligand_atoms, zn_atoms
        )

        return {
            'contacts_4A': contacts_4A,
            'clashes_2.2A': clashes,
            'polar_close_contacts_3p5A': polar_close_contacts_3p5A,
            'polar_contacts': polar_contacts,
            'nonpolar_contacts': nonpolar_contacts,
            'total_contacts': polar_contacts + nonpolar_contacts,
            'polar_fraction': polar_contacts / (polar_contacts + nonpolar_contacts) if (polar_contacts + nonpolar_contacts) > 0 else 0,
            'zn_min_dist': zn_min_dist,
            'zn_coordination_count': zn_coord_count,
            'zn_ligand_contacts': zn_lig_contacts,
            'ligand_heavy_atoms': len(ligand_atoms),
            'protein_heavy_atoms': len(protein_atoms)
        }

    except Exception as e:
        print(f"Error processing {cif_file}: {e}")
        return None

def main():
    # Load manifest
    manifest = pd.read_csv('manifest.csv')
    print(f"Loaded manifest: {len(manifest)} structures")

    # Load AF3 metrics
    af3_metrics = pd.read_csv('MASTER_AF3_RESULTS.csv')
    print(f"Loaded AF3 metrics: {len(af3_metrics)} structures")

    # Filter to top-ranked models only
    manifest_top = manifest[~manifest['file'].str.contains('sample-', na=False)]
    print(f"Top-ranked models: {len(manifest_top)}")

    results = []

    for idx, row in manifest_top.iterrows():
        cif_file = row['file']
        job_name = os.path.basename(cif_file).replace('_model.cif', '')

        print(f"Processing {idx+1}/{len(manifest_top)}: {job_name}")

        # Parse chains
        protein_chain = row['protein_chain']
        ligand_chain = row['ligand_chain']
        zn_present = row['zn_present']

        # Determine Zn chain
        zn_chain = None
        if zn_present > 0:
            # Usually chain C or D for Zn
            chains_str = row['chains']
            if 'C:1:1' in chains_str or 'D:1:1' in chains_str:
                # Find the Zn chain
                for ch in ['C', 'D', 'E']:
                    if f'{ch}:1:1' in chains_str:
                        zn_chain = ch
                        break

        # Compute geometry metrics
        metrics = process_structure(cif_file, protein_chain, ligand_chain, zn_chain)

        if metrics is None:
            continue

        # Classify structure
        enzyme, epoch, domain = classify_structure(job_name)

        # Combine all data
        result = {
            'job_name': job_name,
            'file': cif_file,
            'ligand': row['lig_resname'],
            'zn_present': zn_present > 0,
            'enzyme': enzyme,
            'epoch': epoch,
            'domain': domain,
            **metrics
        }

        results.append(result)

    # Create DataFrame
    df = pd.DataFrame(results)

    # Merge with AF3 metrics
    df_merged = df.merge(
        af3_metrics[['job_name', 'iptm', 'ptm', 'pocket_iptm', 'has_clash', 'ranking_score']],
        on='job_name',
        how='left'
    )

    # Reorder columns
    column_order = [
        'job_name', 'file', 'enzyme', 'epoch', 'domain', 'ligand', 'zn_present',
        'iptm', 'ptm', 'pocket_iptm', 'ranking_score',
        'contacts_4A', 'clashes_2.2A', 'polar_close_contacts_3p5A',
        'polar_contacts', 'nonpolar_contacts', 'polar_fraction',
        'zn_min_dist', 'zn_coordination_count', 'zn_ligand_contacts',
        'ligand_heavy_atoms', 'protein_heavy_atoms', 'has_clash'
    ]

    df_merged = df_merged[column_order]

    # Save
    output_file = 'geometry_metrics.csv'
    df_merged.to_csv(output_file, index=False)
    print(f"\nSaved {len(df_merged)} structures to {output_file}")

    # Print summary
    print("\n" + "="*70)
    print("SUMMARY STATISTICS")
    print("="*70)
    print(f"\nBy enzyme:")
    print(df_merged.groupby('enzyme').size())
    print(f"\nBy epoch:")
    print(df_merged.groupby('epoch').size())
    print(f"\nBy domain:")
    print(df_merged.groupby('domain').size())
    print(f"\nWith Zn: {df_merged['zn_present'].sum()}/{len(df_merged)}")

    print("\n" + "="*70)
    print("GEOMETRY METRICS SUMMARY")
    print("="*70)
    print(df_merged[['contacts_4A', 'clashes_2.2A', 'polar_close_contacts_3p5A', 'polar_fraction', 'pocket_iptm']].describe())

if __name__ == '__main__':
    main()
