#!/usr/bin/env python3
"""
Geometric Analysis v2 - Reads Zn from complex.pdb
Validates the Zinc Paradox and ProRS Evolution geometrically
"""

import numpy as np
from pathlib import Path
import pandas as pd
from scipy import stats

OUTPUT_DIR = Path("/storage/kiran-stuff/aaRS/phase2/energy_scoring/vina_182jobs")


def parse_pdb_atoms(filepath):
    """Parse PDB, return list of atom dicts"""
    atoms = []
    if not filepath.exists():
        return atoms
    with open(filepath) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    atoms.append({
                        'name': line[12:16].strip(),
                        'res': line[17:20].strip(),
                        'chain': line[21],
                        'resnum': int(line[22:26]),
                        'x': float(line[30:38]),
                        'y': float(line[38:46]),
                        'z': float(line[46:54]),
                        'element': line[76:78].strip() if len(line) > 76 else line[12:14].strip()
                    })
                except:
                    pass
    return atoms


def distance(a1, a2):
    """3D distance between two atom dicts"""
    return np.sqrt((a1['x']-a2['x'])**2 + (a1['y']-a2['y'])**2 + (a1['z']-a2['z'])**2)


def analyze_job(job_name, ligand):
    """Comprehensive geometric analysis"""
    safe_name = f"{job_name}_{ligand}".replace("/", "_")
    work_dir = OUTPUT_DIR / safe_name

    if not work_dir.exists():
        return None

    # Load files
    complex_atoms = parse_pdb_atoms(work_dir / "complex.pdb")
    ligand_atoms = parse_pdb_atoms(work_dir / "ligand.pdb")

    # Use minimized ligand if available
    lig_min_path = work_dir / "ligand_min.pdbqt"
    if lig_min_path.exists():
        with open(lig_min_path) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        # Update ligand positions from minimized file
                        name = line[12:16].strip()
                        for la in ligand_atoms:
                            if la['name'] == name:
                                la['x'] = float(line[30:38])
                                la['y'] = float(line[38:46])
                                la['z'] = float(line[46:54])
                                break
                    except:
                        pass

    protein_atoms = parse_pdb_atoms(work_dir / "protein.pdb")

    if not ligand_atoms or not protein_atoms:
        return None

    # Find Zn in complex
    zn_atoms = [a for a in complex_atoms if a['res'] == 'ZN' or a['element'] == 'ZN']
    has_zn = len(zn_atoms) > 0

    # Find ligand N (amino group)
    lig_n = [a for a in ligand_atoms if a['name'] == 'N']
    lig_o = [a for a in ligand_atoms if a['name'].startswith('O')]

    # Zn-N distance
    zn_n_dist = None
    if zn_atoms and lig_n:
        zn_n_dist = min([distance(zn, n) for zn in zn_atoms for n in lig_n])

    # Contacts within 4Å
    contacts_4A = 0
    contact_residues = set()
    for la in ligand_atoms:
        for pa in protein_atoms:
            if distance(la, pa) < 4.0:
                contacts_4A += 1
                contact_residues.add((pa['res'], pa['resnum']))

    # H-bond candidates (N-O, O-O within 3.5Å)
    hbonds = []
    lig_polar = [a for a in ligand_atoms if a['name'].startswith(('N', 'O'))]
    prot_polar = [a for a in protein_atoms if a['name'].startswith(('N', 'O'))]

    for la in lig_polar:
        for pa in prot_polar:
            d = distance(la, pa)
            if d < 3.5:
                hbonds.append({
                    'lig_atom': la['name'],
                    'prot_atom': pa['name'],
                    'prot_res': f"{pa['res']}{pa['resnum']}",
                    'distance': d
                })

    # Ligand burial (fraction of atoms with >10 contacts)
    burial_score = 0
    for la in ligand_atoms:
        nearby = sum(1 for pa in protein_atoms if distance(la, pa) < 5.0)
        if nearby > 10:
            burial_score += 1
    burial_frac = burial_score / len(ligand_atoms) if ligand_atoms else 0

    return {
        'job_name': job_name,
        'ligand': ligand,
        'has_zn': has_zn,
        'zn_n_dist': zn_n_dist,
        'contacts_4A': contacts_4A,
        'contact_residues': len(contact_residues),
        'n_hbonds': len(hbonds),
        'best_hbond': min([h['distance'] for h in hbonds]) if hbonds else None,
        'burial_frac': burial_frac,
        'hbond_details': hbonds[:5] if hbonds else []
    }


def main():
    print("=" * 70)
    print("GEOMETRIC ANALYSIS v2 - Vina Minimized Poses")
    print("=" * 70)

    df = pd.read_csv(OUTPUT_DIR / "vina_scores_LOCAL.csv")

    results = []
    for _, row in df.iterrows():
        result = analyze_job(row['job_name'], row['ligand'])
        if result:
            result['vina_local'] = row['vina_local']
            result['iptm'] = row.get('iptm')
            results.append(result)

    rdf = pd.DataFrame(results)

    # ================================================================
    # 1. ZINC COORDINATION
    # ================================================================
    print("\n" + "=" * 70)
    print("1. ZINC COORDINATION ANALYSIS")
    print("=" * 70)

    zn_jobs = rdf[rdf['has_zn'] == True]
    no_zn_jobs = rdf[rdf['has_zn'] == False]

    print(f"\nJobs with Zn: {len(zn_jobs)}")
    print(f"Jobs without Zn: {len(no_zn_jobs)}")

    if len(zn_jobs) > 0:
        zn_dists = zn_jobs['zn_n_dist'].dropna()
        print(f"\nZn-N(amino) distances:")
        print(f"  Mean: {zn_dists.mean():.2f} Å")
        print(f"  Range: {zn_dists.min():.2f} - {zn_dists.max():.2f} Å")
        print(f"  Coordinated (<2.5Å): {(zn_dists < 2.5).sum()}")

    # ================================================================
    # 2. ThrRS ZINC PARADOX
    # ================================================================
    print("\n" + "=" * 70)
    print("2. ThrRS ZINC PARADOX - GEOMETRIC VALIDATION")
    print("=" * 70)

    thrrs_zn = rdf[rdf['job_name'] == 'modern_thrrs_ecoli_zn']
    thrrs_no_zn = rdf[rdf['job_name'] == 'modern_thrrs']

    if len(thrrs_zn) > 0 and len(thrrs_no_zn) > 0:
        print("\n--- ALL LIGANDS ---")
        print(f"{'Metric':<25} {'No Zn':>12} {'With Zn':>12} {'Δ':>10}")
        print("-" * 60)

        metrics = ['contacts_4A', 'n_hbonds', 'burial_frac', 'vina_local']
        labels = ['Contacts (4Å)', 'H-bonds (3.5Å)', 'Burial fraction', 'Vina (kcal/mol)']

        for m, l in zip(metrics, labels):
            no_zn_val = thrrs_no_zn[m].mean()
            zn_val = thrrs_zn[m].mean()
            delta = zn_val - no_zn_val
            print(f"{l:<25} {no_zn_val:>12.2f} {zn_val:>12.2f} {delta:>+10.2f}")

        # Cognate only
        print("\n--- COGNATE (THR) ONLY ---")
        thr_no_zn = thrrs_no_zn[thrrs_no_zn['ligand'] == 'THR'].iloc[0] if len(thrrs_no_zn[thrrs_no_zn['ligand'] == 'THR']) > 0 else None
        thr_zn = thrrs_zn[thrrs_zn['ligand'] == 'THR'].iloc[0] if len(thrrs_zn[thrrs_zn['ligand'] == 'THR']) > 0 else None

        if thr_no_zn is not None and thr_zn is not None:
            print(f"{'Metric':<25} {'No Zn':>12} {'With Zn':>12}")
            print("-" * 50)
            print(f"{'Contacts (4Å)':<25} {thr_no_zn['contacts_4A']:>12} {thr_zn['contacts_4A']:>12}")
            print(f"{'H-bonds (3.5Å)':<25} {thr_no_zn['n_hbonds']:>12} {thr_zn['n_hbonds']:>12}")
            print(f"{'Zn-N distance (Å)':<25} {'N/A':>12} {thr_zn['zn_n_dist']:>12.2f}")
            print(f"{'Vina (kcal/mol)':<25} {thr_no_zn['vina_local']:>12.2f} {thr_zn['vina_local']:>12.2f}")

        # Key insight
        print("\n*** KEY FINDING ***")
        print(f"Zn provides {thrrs_zn['contacts_4A'].mean() - thrrs_no_zn['contacts_4A'].mean():.0f} more contacts on average")
        print(f"Zn improves Vina by {thrrs_no_zn['vina_local'].mean() - thrrs_zn['vina_local'].mean():.2f} kcal/mol (more negative = better)")

    # ================================================================
    # 3. ProRS EVOLUTION
    # ================================================================
    print("\n" + "=" * 70)
    print("3. ProRS EVOLUTIONARY TRAJECTORY")
    print("=" * 70)

    constructs = ['anc_prors_cat', 'anc_prors_edit', 'modern_prors']
    labels = ['Ancestral Catalytic', 'Ancestral + Editing', 'Modern ProRS']

    print(f"\n{'Construct':<22} {'Contacts':>10} {'H-bonds':>10} {'Burial':>10} {'Vina(PRO)':>12}")
    print("-" * 70)

    for c, l in zip(constructs, labels):
        subset = rdf[rdf['job_name'] == c]
        if len(subset) == 0:
            continue

        pro = subset[subset['ligand'] == 'PRO']
        if len(pro) == 0:
            continue

        pro = pro.iloc[0]
        print(f"{l:<22} {pro['contacts_4A']:>10} {pro['n_hbonds']:>10} {pro['burial_frac']:>10.2f} {pro['vina_local']:>12.2f}")

    # ================================================================
    # 4. CORRELATIONS
    # ================================================================
    print("\n" + "=" * 70)
    print("4. STRUCTURE-ENERGY CORRELATIONS")
    print("=" * 70)

    valid = rdf[rdf['vina_local'].notna()]

    r_cont, p_cont = stats.pearsonr(valid['contacts_4A'], valid['vina_local'])
    r_hb, p_hb = stats.pearsonr(valid['n_hbonds'], valid['vina_local'])
    r_burial, p_burial = stats.pearsonr(valid['burial_frac'], valid['vina_local'])

    print(f"\nCorrelation with Vina score:")
    print(f"  Contacts (4Å):  r = {r_cont:.3f}, p = {p_cont:.2e}")
    print(f"  H-bonds:        r = {r_hb:.3f}, p = {p_hb:.2e}")
    print(f"  Burial:         r = {r_burial:.3f}, p = {p_burial:.2e}")

    # Save
    rdf.to_csv(OUTPUT_DIR / "geometric_analysis_v2.csv", index=False)
    print(f"\nSaved: {OUTPUT_DIR / 'geometric_analysis_v2.csv'}")


if __name__ == "__main__":
    main()
