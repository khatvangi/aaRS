#!/usr/bin/env python3
"""
Geometric Validation - Direct measurements for paper
Validates: 1) Zinc Paradox, 2) ProRS Evolution, 3) Contacts vs Energy
"""

import numpy as np
from pathlib import Path
import pandas as pd

OUTPUT_DIR = Path("/storage/kiran-stuff/aaRS/phase2/energy_scoring/vina_182jobs")


def parse_pdb(filepath):
    """Parse PDB file"""
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


def dist(a1, a2):
    return np.sqrt((a1['x']-a2['x'])**2 + (a1['y']-a2['y'])**2 + (a1['z']-a2['z'])**2)


def analyze_structure(job_dir):
    """Analyze a single structure"""
    complex_atoms = parse_pdb(job_dir / "complex.pdb")
    ligand_atoms = parse_pdb(job_dir / "ligand.pdb")
    protein_atoms = parse_pdb(job_dir / "protein.pdb")

    if not ligand_atoms:
        return None

    # Find Zn
    zn = [a for a in complex_atoms if a['res'] == 'ZN' or a['element'] == 'ZN']

    # Ligand N (amino group)
    lig_n = [a for a in ligand_atoms if a['name'] == 'N']
    lig_o = [a for a in ligand_atoms if a['name'].startswith('O')]

    # Zn-N distance
    zn_n_dist = None
    if zn and lig_n:
        zn_n_dist = min([dist(z, n) for z in zn for n in lig_n])

    # Contacts
    contacts = sum(1 for la in ligand_atoms for pa in protein_atoms if dist(la, pa) < 4.0)

    # H-bonds
    lig_polar = [a for a in ligand_atoms if a['name'].startswith(('N', 'O'))]
    prot_polar = [a for a in protein_atoms if a['name'].startswith(('N', 'O'))]
    hbonds = sum(1 for la in lig_polar for pa in prot_polar if dist(la, pa) < 3.5)

    return {
        'has_zn': len(zn) > 0,
        'zn_n_dist': zn_n_dist,
        'contacts': contacts,
        'hbonds': hbonds
    }


def main():
    df = pd.read_csv(OUTPUT_DIR / "vina_scores_LOCAL.csv")

    print("=" * 80)
    print("GEOMETRIC VALIDATION FOR PAPER")
    print("=" * 80)

    # ================================================================
    # 1. ZINC PARADOX - ThrRS comparison
    # ================================================================
    print("\n" + "=" * 80)
    print("1. ZINC PARADOX: modern_thrrs (no Zn) vs modern_thrrs_ecoli_zn (with Zn)")
    print("=" * 80)

    # Identify the two construct groups
    no_zn_jobs = df[df['job_name'].str.startswith('modern_thrrs_') &
                    ~df['job_name'].str.contains('ecoli_zn') &
                    ~df['job_name'].str.contains('ecoli_THR')]
    zn_jobs = df[df['job_name'].str.startswith('modern_thrrs_ecoli_zn_')]

    print(f"\nNo-Zn jobs: {len(no_zn_jobs)}")
    print(f"With-Zn jobs: {len(zn_jobs)}")

    # Analyze each
    results_no_zn = []
    results_zn = []

    for _, row in no_zn_jobs.iterrows():
        lig = row['ligand']
        job_dir = OUTPUT_DIR / f"modern_thrrs_{lig}_{lig}"
        data = analyze_structure(job_dir)
        if data:
            data['ligand'] = lig
            data['vina'] = row['vina_local']
            results_no_zn.append(data)

    for _, row in zn_jobs.iterrows():
        lig = row['ligand']
        job_dir = OUTPUT_DIR / f"modern_thrrs_ecoli_zn_{lig}_{lig}"
        data = analyze_structure(job_dir)
        if data:
            data['ligand'] = lig
            data['vina'] = row['vina_local']
            results_zn.append(data)

    df_no_zn = pd.DataFrame(results_no_zn)
    df_zn = pd.DataFrame(results_zn)

    print(f"\n{'Metric':<25} {'No Zn':>15} {'With Zn':>15} {'Δ':>12}")
    print("-" * 70)

    for metric, label in [('contacts', 'Contacts (4Å)'), ('hbonds', 'H-bonds (3.5Å)'), ('vina', 'Vina (kcal/mol)')]:
        no_zn_val = df_no_zn[metric].mean()
        zn_val = df_zn[metric].mean()
        delta = zn_val - no_zn_val
        print(f"{label:<25} {no_zn_val:>15.2f} {zn_val:>15.2f} {delta:>+12.2f}")

    # Cognate comparison
    print("\n--- COGNATE (THR) COMPARISON ---")
    thr_no_zn = df_no_zn[df_no_zn['ligand'] == 'THR'].iloc[0] if len(df_no_zn[df_no_zn['ligand'] == 'THR']) > 0 else None
    thr_zn = df_zn[df_zn['ligand'] == 'THR'].iloc[0] if len(df_zn[df_zn['ligand'] == 'THR']) > 0 else None

    if thr_no_zn is not None and thr_zn is not None:
        print(f"  No Zn:   Contacts={thr_no_zn['contacts']}, H-bonds={thr_no_zn['hbonds']}, Vina={thr_no_zn['vina']:.2f}")
        print(f"  With Zn: Contacts={thr_zn['contacts']}, H-bonds={thr_zn['hbonds']}, Zn-N={thr_zn['zn_n_dist']:.2f}Å, Vina={thr_zn['vina']:.2f}")

    # Zn-N distances for all ligands with Zn
    print("\n--- Zn-N DISTANCES (coordination) ---")
    zn_dists = df_zn[['ligand', 'zn_n_dist']].dropna()
    coordinated = zn_dists[zn_dists['zn_n_dist'] < 2.5]
    print(f"  Coordinated (<2.5Å): {len(coordinated)}/{len(zn_dists)}")
    print(f"  Mean Zn-N distance: {zn_dists['zn_n_dist'].mean():.2f} Å")

    # ================================================================
    # 2. ProRS EVOLUTION
    # ================================================================
    print("\n" + "=" * 80)
    print("2. ProRS EVOLUTIONARY TRAJECTORY")
    print("=" * 80)

    constructs = [
        ('anc_prors_cat', 'Ancestral Catalytic'),
        ('anc_prors_edit', 'Ancestral + Editing'),
        ('modern_prors', 'Modern ProRS')
    ]

    print(f"\n{'Construct':<25} {'PRO Contacts':>12} {'PRO H-bonds':>12} {'PRO Vina':>12}")
    print("-" * 65)

    evolution_data = []
    for prefix, label in constructs:
        pro_row = df[df['job_name'] == f"{prefix}_PRO"]
        if len(pro_row) == 0:
            continue

        job_dir = OUTPUT_DIR / f"{prefix}_PRO_PRO"
        data = analyze_structure(job_dir)
        if data:
            vina = pro_row['vina_local'].values[0]
            print(f"{label:<25} {data['contacts']:>12} {data['hbonds']:>12} {vina:>12.2f}")
            evolution_data.append({
                'construct': label,
                'contacts': data['contacts'],
                'hbonds': data['hbonds'],
                'vina': vina
            })

    if len(evolution_data) >= 3:
        print("\n*** EVOLUTIONARY PATTERN ***")
        print(f"  Anc Cat → Anc+Edit: Contacts {evolution_data[0]['contacts']} → {evolution_data[1]['contacts']} (Δ={evolution_data[1]['contacts']-evolution_data[0]['contacts']:+d})")
        print(f"  Anc+Edit → Modern:  Contacts {evolution_data[1]['contacts']} → {evolution_data[2]['contacts']} (Δ={evolution_data[2]['contacts']-evolution_data[1]['contacts']:+d})")
        print(f"  Vina: {evolution_data[0]['vina']:.2f} → {evolution_data[1]['vina']:.2f} → {evolution_data[2]['vina']:.2f}")

    # ================================================================
    # 3. SUMMARY TABLE FOR PAPER
    # ================================================================
    print("\n" + "=" * 80)
    print("3. SUMMARY TABLE FOR PAPER")
    print("=" * 80)

    # Get selectivity for each construct
    constructs = [
        ('modern_prors', 'Modern ProRS', 'PRO'),
        ('anc_prors_cat', 'Anc ProRS (cat)', 'PRO'),
        ('anc_prors_edit', 'Anc ProRS (edit)', 'PRO'),
        ('modern_thrrs_ecoli_zn', 'Modern ThrRS +Zn', 'THR'),
        ('modern_thrrs', 'Modern ThrRS -Zn', 'THR'),
        ('anc_thrrs_cat_zn', 'Anc ThrRS +Zn', 'THR'),
        ('anc_thrrs_cat', 'Anc ThrRS -Zn', 'THR'),
    ]

    print(f"\n{'Construct':<22} {'Cognate':>8} {'Vina(cog)':>10} {'Contacts':>10} {'H-bonds':>8} {'Zn-N(Å)':>8}")
    print("-" * 75)

    for prefix, label, cognate in constructs:
        # Get cognate row
        cog_row = df[df['job_name'] == f"{prefix}_{cognate}"]
        if len(cog_row) == 0:
            continue

        job_dir = OUTPUT_DIR / f"{prefix}_{cognate}_{cognate}"
        data = analyze_structure(job_dir)
        if data:
            vina = cog_row['vina_local'].values[0]
            zn_str = f"{data['zn_n_dist']:.2f}" if data['zn_n_dist'] else "N/A"
            print(f"{label:<22} {cognate:>8} {vina:>10.2f} {data['contacts']:>10} {data['hbonds']:>8} {zn_str:>8}")

    print("\n" + "=" * 80)
    print("VALIDATION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
