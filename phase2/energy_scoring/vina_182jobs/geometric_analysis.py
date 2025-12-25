#!/usr/bin/env python3
"""
Geometric Analysis of Vina Minimized Poses

Validates energy findings with structural measurements:
1. Ligand burial / contact surface
2. Zn coordination distances
3. H-bond distances
4. Pocket compactness
"""

import numpy as np
from pathlib import Path
import pandas as pd
from collections import defaultdict

OUTPUT_DIR = Path("/storage/kiran-stuff/aaRS/phase2/energy_scoring/vina_182jobs")

def parse_pdbqt(filepath):
    """Parse PDBQT file, return list of (atom_name, x, y, z, atom_type)"""
    atoms = []
    if not filepath.exists():
        return atoms
    with open(filepath) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    atom_name = line[12:16].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    atom_type = line[77:79].strip() if len(line) > 77 else ""
                    atoms.append((atom_name, x, y, z, atom_type))
                except:
                    pass
    return atoms


def parse_pdb(filepath):
    """Parse PDB file"""
    atoms = []
    if not filepath.exists():
        return atoms
    with open(filepath) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    chain = line[21]
                    res_num = int(line[22:26])
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    element = line[76:78].strip() if len(line) > 76 else atom_name[0]
                    atoms.append({
                        'name': atom_name, 'res': res_name, 'chain': chain,
                        'resnum': res_num, 'x': x, 'y': y, 'z': z, 'element': element
                    })
                except:
                    pass
    return atoms


def distance(p1, p2):
    """Euclidean distance"""
    return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)


def get_ligand_center(atoms):
    """Get centroid of atoms"""
    if not atoms:
        return None
    coords = [(a[1], a[2], a[3]) if isinstance(a, tuple) else (a['x'], a['y'], a['z']) for a in atoms]
    return tuple(np.mean(coords, axis=0))


def count_contacts(ligand_atoms, protein_atoms, cutoff=4.0):
    """Count protein atoms within cutoff of any ligand atom"""
    contacts = 0
    contact_residues = set()
    for la in ligand_atoms:
        lcoord = (la[1], la[2], la[3]) if isinstance(la, tuple) else (la['x'], la['y'], la['z'])
        for pa in protein_atoms:
            pcoord = (pa['x'], pa['y'], pa['z'])
            if distance(lcoord, pcoord) < cutoff:
                contacts += 1
                contact_residues.add((pa['res'], pa['resnum']))
    return contacts, len(contact_residues)


def find_zn_distance(ligand_atoms, protein_atoms):
    """Find distance from ligand N (amino) to nearest Zn"""
    # Find Zn in protein
    zn_atoms = [a for a in protein_atoms if a['element'] == 'ZN' or a['name'] == 'ZN']
    if not zn_atoms:
        return None, None

    # Find ligand N (amino group)
    lig_n = [a for a in ligand_atoms if a[0] == 'N' or (isinstance(a, dict) and a['name'] == 'N')]
    if not lig_n:
        return None, None

    # Get closest Zn to ligand N
    min_dist = float('inf')
    for zn in zn_atoms:
        zn_coord = (zn['x'], zn['y'], zn['z'])
        for n in lig_n:
            n_coord = (n[1], n[2], n[3]) if isinstance(n, tuple) else (n['x'], n['y'], n['z'])
            d = distance(zn_coord, n_coord)
            if d < min_dist:
                min_dist = d

    return min_dist if min_dist < 10 else None, len(zn_atoms)


def find_hbond_distances(ligand_atoms, protein_atoms):
    """Find potential H-bond distances (N-O, O-O within 3.5A)"""
    hbond_candidates = []

    # Ligand polar atoms (N, O)
    lig_polar = [(a[0], a[1], a[2], a[3]) for a in ligand_atoms
                 if isinstance(a, tuple) and a[0].startswith(('N', 'O'))]

    # Protein polar atoms
    prot_polar = [(a['name'], a['x'], a['y'], a['z'], a['res'], a['resnum'])
                  for a in protein_atoms if a['name'].startswith(('N', 'O'))]

    for la in lig_polar:
        for pa in prot_polar:
            d = distance((la[1], la[2], la[3]), (pa[1], pa[2], pa[3]))
            if d < 3.5:
                hbond_candidates.append({
                    'lig_atom': la[0],
                    'prot_atom': pa[0],
                    'prot_res': f"{pa[4]}{pa[5]}",
                    'distance': d
                })

    return hbond_candidates


def analyze_job(job_name, ligand):
    """Analyze a single job's geometry"""
    safe_name = f"{job_name}_{ligand}".replace("/", "_")
    work_dir = OUTPUT_DIR / safe_name

    if not work_dir.exists():
        return None

    # Load minimized ligand
    lig_min_path = work_dir / "ligand_min.pdbqt"
    if not lig_min_path.exists():
        lig_min_path = work_dir / "ligand.pdbqt"

    ligand_atoms = parse_pdbqt(lig_min_path)
    protein_atoms = parse_pdb(work_dir / "protein.pdb")

    if not ligand_atoms or not protein_atoms:
        return None

    # Measurements
    lig_center = get_ligand_center(ligand_atoms)
    n_contacts, n_contact_res = count_contacts(ligand_atoms, protein_atoms, cutoff=4.0)
    zn_dist, n_zn = find_zn_distance(ligand_atoms, protein_atoms)
    hbonds = find_hbond_distances(ligand_atoms, protein_atoms)

    return {
        'job_name': job_name,
        'ligand': ligand,
        'n_lig_atoms': len(ligand_atoms),
        'n_contacts_4A': n_contacts,
        'n_contact_residues': n_contact_res,
        'has_zn': n_zn > 0 if n_zn else False,
        'zn_N_distance': zn_dist,
        'n_hbonds_3.5A': len(hbonds),
        'best_hbond': min([h['distance'] for h in hbonds]) if hbonds else None
    }


def main():
    print("=" * 70)
    print("GEOMETRIC ANALYSIS OF VINA MINIMIZED POSES")
    print("=" * 70)

    # Load Vina scores
    df = pd.read_csv(OUTPUT_DIR / "vina_scores_LOCAL.csv")

    results = []
    for _, row in df.iterrows():
        result = analyze_job(row['job_name'], row['ligand'])
        if result:
            result['vina_local'] = row['vina_local']
            result['iptm'] = row.get('iptm')
            results.append(result)

    results_df = pd.DataFrame(results)

    # === ANALYSIS 1: Zn Coordination ===
    print("\n" + "=" * 70)
    print("1. ZINC COORDINATION ANALYSIS")
    print("=" * 70)

    zn_jobs = results_df[results_df['has_zn'] == True]
    no_zn_jobs = results_df[results_df['has_zn'] == False]

    print(f"\nJobs with Zn: {len(zn_jobs)}")
    print(f"Jobs without Zn: {len(no_zn_jobs)}")

    if len(zn_jobs) > 0:
        print(f"\nZn-N distance (ligand amino to Zn):")
        zn_dists = zn_jobs['zn_N_distance'].dropna()
        print(f"  Mean: {zn_dists.mean():.2f} Å")
        print(f"  Range: {zn_dists.min():.2f} - {zn_dists.max():.2f} Å")
        print(f"  Coordination (<2.5Å): {(zn_dists < 2.5).sum()}")
        print(f"  Close (<3.0Å): {(zn_dists < 3.0).sum()}")

    # === ANALYSIS 2: ThrRS Zn Paradox ===
    print("\n" + "=" * 70)
    print("2. ThrRS ZINC PARADOX - GEOMETRIC VALIDATION")
    print("=" * 70)

    # Compare modern_thrrs (no Zn) vs modern_thrrs_ecoli_zn (with Zn)
    thrrs_no_zn = results_df[results_df['job_name'] == 'modern_thrrs']
    thrrs_zn = results_df[results_df['job_name'] == 'modern_thrrs_ecoli_zn']

    if len(thrrs_no_zn) > 0 and len(thrrs_zn) > 0:
        print("\nmodern_thrrs (NO Zn):")
        print(f"  Mean contacts (4Å): {thrrs_no_zn['n_contacts_4A'].mean():.1f}")
        print(f"  Mean H-bonds (3.5Å): {thrrs_no_zn['n_hbonds_3.5A'].mean():.1f}")
        print(f"  Mean Vina: {thrrs_no_zn['vina_local'].mean():.2f} kcal/mol")

        print("\nmodern_thrrs_ecoli_zn (WITH Zn):")
        print(f"  Mean contacts (4Å): {thrrs_zn['n_contacts_4A'].mean():.1f}")
        print(f"  Mean H-bonds (3.5Å): {thrrs_zn['n_hbonds_3.5A'].mean():.1f}")
        print(f"  Mean Zn-N distance: {thrrs_zn['zn_N_distance'].mean():.2f} Å")
        print(f"  Mean Vina: {thrrs_zn['vina_local'].mean():.2f} kcal/mol")

        # Cognate comparison
        thr_no_zn = thrrs_no_zn[thrrs_no_zn['ligand'] == 'THR']
        thr_zn = thrrs_zn[thrrs_zn['ligand'] == 'THR']

        if len(thr_no_zn) > 0 and len(thr_zn) > 0:
            print("\n--- COGNATE (THR) ---")
            print(f"  No Zn:   {thr_no_zn['n_contacts_4A'].values[0]} contacts, {thr_no_zn['n_hbonds_3.5A'].values[0]} H-bonds, Vina={thr_no_zn['vina_local'].values[0]:.2f}")
            print(f"  With Zn: {thr_zn['n_contacts_4A'].values[0]} contacts, {thr_zn['n_hbonds_3.5A'].values[0]} H-bonds, Zn-N={thr_zn['zn_N_distance'].values[0]:.2f}Å, Vina={thr_zn['vina_local'].values[0]:.2f}")

    # === ANALYSIS 3: ProRS Evolution ===
    print("\n" + "=" * 70)
    print("3. ProRS EVOLUTIONARY TRAJECTORY - GEOMETRIC VALIDATION")
    print("=" * 70)

    prors_constructs = ['anc_prors_cat', 'anc_prors_edit', 'modern_prors']

    for construct in prors_constructs:
        subset = results_df[results_df['job_name'] == construct]
        if len(subset) == 0:
            continue

        pro_data = subset[subset['ligand'] == 'PRO']

        print(f"\n{construct}:")
        print(f"  All ligands - Mean contacts: {subset['n_contacts_4A'].mean():.1f}, H-bonds: {subset['n_hbonds_3.5A'].mean():.1f}")
        if len(pro_data) > 0:
            print(f"  COGNATE (PRO) - Contacts: {pro_data['n_contacts_4A'].values[0]}, H-bonds: {pro_data['n_hbonds_3.5A'].values[0]}, Vina: {pro_data['vina_local'].values[0]:.2f}")

    # === ANALYSIS 4: Contact vs Energy Correlation ===
    print("\n" + "=" * 70)
    print("4. CONTACTS vs VINA ENERGY CORRELATION")
    print("=" * 70)

    from scipy import stats

    valid = results_df[results_df['vina_local'].notna() & results_df['n_contacts_4A'].notna()]
    r_contacts, p_contacts = stats.pearsonr(valid['n_contacts_4A'], valid['vina_local'])
    r_hbonds, p_hbonds = stats.pearsonr(valid['n_hbonds_3.5A'], valid['vina_local'])

    print(f"Contacts (4Å) vs Vina: r = {r_contacts:.3f}, p = {p_contacts:.2e}")
    print(f"H-bonds (3.5Å) vs Vina: r = {r_hbonds:.3f}, p = {p_hbonds:.2e}")

    # Save results
    results_df.to_csv(OUTPUT_DIR / "geometric_analysis.csv", index=False)
    print(f"\nSaved: {OUTPUT_DIR / 'geometric_analysis.csv'}")

    return results_df


if __name__ == "__main__":
    df = main()
