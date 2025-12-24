#!/usr/bin/env python3
"""
Recompute geometry metrics from CIF files using gemmi.
No energies, no MD - pure geometry only.
"""
import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

try:
    import gemmi
    HAS_GEMMI = True
except ImportError:
    HAS_GEMMI = False
    print("WARNING: gemmi not available, falling back to manual CIF parsing")

def parse_cif_gemmi(cif_path, protein_chain, ligand_chain, zn_chain=None):
    """Parse CIF using gemmi library."""
    structure = gemmi.read_structure(str(cif_path))
    model = structure[0]

    protein_atoms = []
    ligand_atoms = []
    zn_atoms = []

    for chain in model:
        if chain.name == protein_chain:
            for residue in chain:
                for atom in residue:
                    if atom.element.name != 'H':
                        protein_atoms.append({
                            'element': atom.element.name,
                            'coords': np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                        })
        elif chain.name == ligand_chain:
            for residue in chain:
                for atom in residue:
                    if atom.element.name != 'H':
                        ligand_atoms.append({
                            'element': atom.element.name,
                            'coords': np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                        })
        elif zn_chain and chain.name == zn_chain:
            for residue in chain:
                for atom in residue:
                    if atom.element.name == 'Zn':
                        zn_atoms.append({
                            'element': 'Zn',
                            'coords': np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                        })

    return protein_atoms, ligand_atoms, zn_atoms

def parse_cif_manual(cif_path, protein_chain, ligand_chain, zn_chain=None):
    """Manual CIF parsing fallback."""
    atoms_by_chain = {}

    with open(cif_path) as f:
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

            element = parts[2]
            chain_id = parts[6]
            x, y, z = float(parts[10]), float(parts[11]), float(parts[12])

            if element == 'H':
                continue

            if chain_id not in atoms_by_chain:
                atoms_by_chain[chain_id] = []

            atoms_by_chain[chain_id].append({
                'element': element,
                'coords': np.array([x, y, z])
            })

    protein_atoms = atoms_by_chain.get(protein_chain, [])
    ligand_atoms = atoms_by_chain.get(ligand_chain, [])
    zn_atoms = []
    if zn_chain and zn_chain in atoms_by_chain:
        zn_atoms = [a for a in atoms_by_chain[zn_chain] if a['element'] == 'Zn']

    return protein_atoms, ligand_atoms, zn_atoms

def compute_distance(coord1, coord2):
    """Euclidean distance."""
    return np.linalg.norm(coord1 - coord2)

def compute_contacts(protein_atoms, ligand_atoms, cutoff):
    """Count contacts within cutoff."""
    contacts = 0
    for p in protein_atoms:
        for l in ligand_atoms:
            if compute_distance(p['coords'], l['coords']) <= cutoff:
                contacts += 1
    return contacts

def compute_clashes(protein_atoms, ligand_atoms, cutoff):
    """Count clashes below cutoff."""
    clashes = 0
    for p in protein_atoms:
        for l in ligand_atoms:
            if compute_distance(p['coords'], l['coords']) < cutoff:
                clashes += 1
    return clashes

def compute_polar_close_contacts(protein_atoms, ligand_atoms, cutoff):
    """
    Count N/O pairs within cutoff (no angle check).
    This is NOT true H-bonds.
    """
    polar_atoms = {'N', 'O', 'S'}
    count = 0

    for p in protein_atoms:
        if p['element'] in polar_atoms:
            for l in ligand_atoms:
                if l['element'] in polar_atoms:
                    if compute_distance(p['coords'], l['coords']) < cutoff:
                        count += 1

    return count

def compute_polar_contacts(protein_atoms, ligand_atoms, cutoff):
    """Count polar contacts (at least one polar atom) within cutoff."""
    polar_atoms = {'N', 'O', 'S'}
    count = 0

    for p in protein_atoms:
        for l in ligand_atoms:
            if compute_distance(p['coords'], l['coords']) <= cutoff:
                if p['element'] in polar_atoms or l['element'] in polar_atoms:
                    count += 1

    return count

def compute_zn_metrics(ligand_atoms, zn_atoms, zn_engaged_cutoff=3.0):
    """
    Compute Zn-ligand distance to heteroatoms (O/N/S) only.
    Returns: min_dist, engaged flag.
    """
    if not zn_atoms or not ligand_atoms:
        return None, False

    hetero_atoms = {'O', 'N', 'S'}
    min_dist = 999.9

    for zn in zn_atoms:
        for lig in ligand_atoms:
            if lig['element'] in hetero_atoms:
                dist = compute_distance(zn['coords'], lig['coords'])
                if dist < min_dist:
                    min_dist = dist

    engaged = (min_dist <= zn_engaged_cutoff) if min_dist < 999 else False
    min_dist = min_dist if min_dist < 999 else None

    return min_dist, engaged

def recompute_metrics_for_row(row, params_dict=None,
                               contacts_cutoff=4.0,
                               clashes_cutoff=2.2,
                               polar_close_cutoff=3.5,
                               zn_engaged_cutoff=3.0):
    """
    Recompute all metrics for one structure.

    Args:
        row: dict with 'file' and optionally 'zn_present', 'protein_chain', 'ligand_chain', 'chains'
        params_dict: optional dict with keys: contact_cutoff, polar_cutoff, polar_close_cutoff,
                     clash_cutoff, zn_engaged_cutoff
        (keyword args are used if params_dict is None)
    """
    # Extract params from dict if provided
    if params_dict is not None:
        contacts_cutoff = params_dict.get('contact_cutoff', contacts_cutoff)
        # polar_cutoff is used for polar contacts (same as contact cutoff typically)
        polar_cutoff = params_dict.get('polar_cutoff', contacts_cutoff)
        polar_close_cutoff = params_dict.get('polar_close_cutoff', polar_close_cutoff)
        clashes_cutoff = params_dict.get('clash_cutoff', clashes_cutoff)
        zn_engaged_cutoff = params_dict.get('zn_engaged_cutoff', zn_engaged_cutoff)
    else:
        polar_cutoff = contacts_cutoff

    cif_path = row['file']
    if not Path(cif_path).exists():
        print(f"WARNING: {cif_path} not found")
        return None

    # Determine chains
    protein_chain = row.get('protein_chain', 'A')
    ligand_chain = row.get('ligand_chain', 'B')

    # Zn chain (if present)
    zn_chain = None
    if row.get('zn_present', False):
        # Try to find Zn chain from original data
        if 'chains' in row:
            chains_str = str(row['chains'])
            for ch in ['C', 'D', 'E']:
                if f'{ch}:1:1' in chains_str:
                    zn_chain = ch
                    break

    # Parse CIF
    try:
        if HAS_GEMMI:
            protein_atoms, ligand_atoms, zn_atoms = parse_cif_gemmi(
                cif_path, protein_chain, ligand_chain, zn_chain
            )
        else:
            protein_atoms, ligand_atoms, zn_atoms = parse_cif_manual(
                cif_path, protein_chain, ligand_chain, zn_chain
            )
    except Exception as e:
        print(f"ERROR parsing {cif_path}: {e}")
        return None

    if not ligand_atoms:
        print(f"WARNING: No ligand atoms found in {cif_path}")
        return None

    # Compute metrics
    contacts = compute_contacts(protein_atoms, ligand_atoms, contacts_cutoff)
    clashes = compute_clashes(protein_atoms, ligand_atoms, clashes_cutoff)
    polar_close = compute_polar_close_contacts(protein_atoms, ligand_atoms, polar_close_cutoff)
    polar_contacts = compute_polar_contacts(protein_atoms, ligand_atoms, polar_cutoff)

    n_ligand_atoms = len(ligand_atoms)

    # Normalized
    contacts_per_atom = contacts / max(1, n_ligand_atoms)
    polar_contacts_per_atom = polar_contacts / max(1, n_ligand_atoms)
    clash_rate = clashes / max(1, n_ligand_atoms)
    polar_close_per_atom = polar_close / max(1, n_ligand_atoms)

    # Zn metrics
    zn_min_dist_hetero, zn_engaged = compute_zn_metrics(
        ligand_atoms, zn_atoms, zn_engaged_cutoff
    )

    return {
        f'contacts_{contacts_cutoff}A': contacts,
        f'clashes_{clashes_cutoff}A': clashes,
        f'polar_close_{polar_close_cutoff}A': polar_close,
        'polar_contacts': polar_contacts,
        'ligand_heavy_atoms': n_ligand_atoms,
        'contacts_per_atom': contacts_per_atom,
        'polar_contacts_per_atom': polar_contacts_per_atom,
        'clash_rate': clash_rate,
        'polar_close_contacts_per_atom': polar_close_per_atom,
        'zn_min_dist_hetero': zn_min_dist_hetero,
        'zn_engaged': zn_engaged,
        'zn_floating': (row.get('zn_present', False) and not zn_engaged)
    }

def main(input_csv='geometry_metrics_clean.csv',
         output_csv='geometry_metrics_recomputed.csv',
         contacts_cutoff=4.0,
         clashes_cutoff=2.2,
         polar_close_cutoff=3.5,
         zn_engaged_cutoff=3.0):
    """Main recomputation function."""

    print(f"Loading {input_csv}...")
    df = pd.read_csv(input_csv)
    print(f"Loaded {len(df)} structures")

    # Load manifest if available for chain info
    if Path('manifest.csv').exists():
        manifest = pd.read_csv('manifest.csv')
        df = df.merge(
            manifest[['file', 'protein_chain', 'ligand_chain', 'chains']],
            on='file',
            how='left',
            suffixes=('', '_manifest')
        )

    results = []

    for idx, row in df.iterrows():
        if (idx + 1) % 20 == 0:
            print(f"Processing {idx+1}/{len(df)}...")

        metrics = recompute_metrics_for_row(
            row,
            contacts_cutoff=contacts_cutoff,
            clashes_cutoff=clashes_cutoff,
            polar_close_cutoff=polar_close_cutoff,
            zn_engaged_cutoff=zn_engaged_cutoff
        )

        if metrics is None:
            continue

        # Merge with original row
        result = row.to_dict()
        result.update(metrics)
        results.append(result)

    df_out = pd.DataFrame(results)
    df_out.to_csv(output_csv, index=False)
    print(f"\nSaved {len(df_out)} structures to {output_csv}")

    # Summary stats
    print("\n" + "="*70)
    print("RECOMPUTED METRICS SUMMARY")
    print("="*70)
    print(f"\nSettings:")
    print(f"  contacts_cutoff: {contacts_cutoff} Å")
    print(f"  clashes_cutoff: {clashes_cutoff} Å")
    print(f"  polar_close_cutoff: {polar_close_cutoff} Å")
    print(f"  zn_engaged_cutoff: {zn_engaged_cutoff} Å")

    print(f"\nZn classification:")
    print(f"  Zn engaged: {df_out['zn_engaged'].sum()}")
    print(f"  Zn floating: {df_out['zn_floating'].sum()}")

    print("\nMetrics:")
    print(df_out[['contacts_per_atom', 'polar_contacts_per_atom',
                  'clash_rate', 'polar_close_contacts_per_atom']].describe())

    return df_out

if __name__ == '__main__':
    main()
