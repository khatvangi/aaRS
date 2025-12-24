#!/usr/bin/env python3
"""
CORRECTED: Extract protein-ligand ipTM from [0][1], not [0][2]
For Zn runs: [0][1] = AA binding, [0][2] = Zn binding
"""
import json
import glob
import pandas as pd
import os

BASE_DIR = "/storage/kiran-stuff/aaRS/phase2"

def extract_all():
    conf_files = glob.glob(f"{BASE_DIR}/**/*summary_confidences.json", recursive=True)
    print(f"Found {len(conf_files)} files")

    records = []
    seen = set()

    for fp in conf_files:
        try:
            with open(fp) as f:
                data = json.load(f)

            cp = data.get('chain_pair_iptm', [])
            if not cp or len(cp) < 2:
                continue

            # CORRECTED: Always use [0][1] for protein-ligand
            aa_iptm = cp[0][1] if len(cp[0]) > 1 else None
            zn_iptm = cp[0][2] if len(cp[0]) > 2 else None

            # Find data.json for ligand info
            data_files = glob.glob(os.path.join(os.path.dirname(fp), "*_data.json"))
            if not data_files:
                data_files = glob.glob(os.path.join(os.path.dirname(os.path.dirname(fp)), "*_data.json"))

            ligands = []
            protein_len = 0
            has_rna = False
            job_name = ""

            if data_files:
                with open(data_files[0]) as f:
                    d = json.load(f)
                job_name = d.get('name', '')
                for seq in d.get('sequences', []):
                    if 'ligand' in seq:
                        ligands.extend(seq['ligand'].get('ccdCodes', []))
                    if 'rna' in seq:
                        has_rna = True
                    if 'protein' in seq:
                        protein_len = len(seq['protein'].get('sequence', ''))

            # Unique key
            key = f"{job_name}_{','.join(ligands)}_{protein_len}"
            if key in seen:
                continue
            seen.add(key)

            records.append({
                'job_name': job_name,
                'ligands': ','.join(ligands),
                'AA_iptm': aa_iptm,  # CORRECTED protein-AA binding
                'Zn_iptm': zn_iptm,  # Protein-Zn (if present)
                'global_iptm': data.get('iptm'),
                'ptm': data.get('ptm'),
                'protein_len': protein_len,
                'has_rna': has_rna,
                'chain_pair_iptm': str(cp),
            })
        except Exception as e:
            pass

    df = pd.DataFrame(records)
    df = df.sort_values('AA_iptm', ascending=False)
    df.to_csv(f"{BASE_DIR}/AF3_RESULTS_CORRECTED.csv", index=False)

    print(f"\nExtracted {len(df)} unique results")
    print("\n" + "="*80)
    print("TOP 30 BY AMINO ACID BINDING (AA_iptm = chain_pair_iptm[0][1])")
    print("="*80)
    print(df[['job_name', 'ligands', 'AA_iptm', 'Zn_iptm', 'global_iptm', 'protein_len']].head(30).to_string())

    # Summary by ligand (excluding Zn)
    df['primary_ligand'] = df['ligands'].str.split(',').str[0]
    print("\n" + "="*80)
    print("AVERAGE AA_iptm BY PRIMARY LIGAND")
    print("="*80)
    stats = df.groupby('primary_ligand')['AA_iptm'].agg(['mean', 'std', 'count']).round(3)
    print(stats.sort_values('mean', ascending=False).to_string())

    return df

if __name__ == '__main__':
    df = extract_all()
