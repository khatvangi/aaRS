#!/usr/bin/env python3
"""
Extract ALL AF3 results from nested output directories.
Trust the actual JSON content, not filenames.
"""
import json
import os
import glob
import pandas as pd
from pathlib import Path

BASE_DIR = "/storage/kiran-stuff/aaRS/phase2"

def extract_ligand_from_data_json(data_path):
    """Extract actual ligand from _data.json file."""
    try:
        with open(data_path) as f:
            data = json.load(f)

        ligands = []
        has_rna = False
        protein_len = 0

        for seq in data.get('sequences', []):
            if 'ligand' in seq:
                codes = seq['ligand'].get('ccdCodes', [])
                ligands.extend(codes)
            if 'rna' in seq:
                has_rna = True
            if 'protein' in seq:
                protein_len = len(seq['protein'].get('sequence', ''))

        return {
            'ligands': ligands,
            'has_rna': has_rna,
            'protein_length': protein_len,
            'job_name': data.get('name', '')
        }
    except Exception as e:
        return None

def extract_scores(conf_path):
    """Extract scores from summary_confidences.json."""
    try:
        with open(conf_path) as f:
            data = json.load(f)

        scores = {
            'iptm': data.get('iptm'),
            'ptm': data.get('ptm'),
            'ranking_score': data.get('ranking_score'),
            'fraction_disordered': data.get('fraction_disordered'),
            'has_clash': data.get('has_clash'),
        }

        # Chain pair iptm - critical for binding
        cp = data.get('chain_pair_iptm', [])
        if cp:
            scores['chain_pair_iptm'] = str(cp)
            # Extract protein-ligand interface score
            if len(cp) >= 2:
                if len(cp) == 2:  # [protein, ligand]
                    scores['pocket_iptm'] = cp[0][1]
                elif len(cp) >= 3:  # [protein, rna, ligand] or [protein, ligand, zinc]
                    scores['pocket_iptm'] = cp[0][2] if len(cp[0]) > 2 else cp[0][1]

        # pLDDT from chain_ptm
        chain_ptm = data.get('chain_ptm', [])
        if chain_ptm:
            scores['protein_ptm'] = chain_ptm[0]

        return scores
    except Exception as e:
        return None

def find_data_json(conf_path):
    """Find corresponding _data.json for a confidence file."""
    # Try same directory
    conf_dir = os.path.dirname(conf_path)

    # Pattern 1: In same directory
    data_files = glob.glob(os.path.join(conf_dir, "*_data.json"))
    if data_files:
        return data_files[0]

    # Pattern 2: Parent directory
    parent = os.path.dirname(conf_dir)
    data_files = glob.glob(os.path.join(parent, "*_data.json"))
    if data_files:
        return data_files[0]

    # Pattern 3: Two levels up
    grandparent = os.path.dirname(parent)
    data_files = glob.glob(os.path.join(grandparent, "*_data.json"))
    if data_files:
        return data_files[0]

    return None

def main():
    print("Scanning for AF3 outputs...")

    # Find all confidence files
    conf_files = glob.glob(f"{BASE_DIR}/**/*summary_confidences.json", recursive=True)
    print(f"Found {len(conf_files)} confidence files")

    records = []
    seen_jobs = set()  # Track unique jobs to avoid duplicates

    for i, conf_path in enumerate(conf_files):
        if i % 100 == 0:
            print(f"Processing {i}/{len(conf_files)}...")

        # Extract scores
        scores = extract_scores(conf_path)
        if not scores:
            continue

        # Find and extract ligand info
        data_path = find_data_json(conf_path)
        input_info = extract_ligand_from_data_json(data_path) if data_path else None

        # Build record
        record = {
            'filepath': conf_path,
            'rel_path': os.path.relpath(conf_path, BASE_DIR),
            **scores
        }

        if input_info:
            record['job_name'] = input_info['job_name']
            record['ligands'] = ','.join(input_info['ligands'])
            record['has_rna'] = input_info['has_rna']
            record['protein_length'] = input_info['protein_length']

            # Create unique key
            job_key = f"{input_info['job_name']}_{','.join(input_info['ligands'])}_{input_info['protein_length']}"
        else:
            # Fallback: extract from path
            record['job_name'] = os.path.basename(os.path.dirname(os.path.dirname(conf_path)))
            job_key = conf_path

        # Skip exact duplicates (keep best ranking_score)
        if job_key in seen_jobs:
            continue
        seen_jobs.add(job_key)

        records.append(record)

    print(f"\nExtracted {len(records)} unique results")

    # Create DataFrame
    df = pd.DataFrame(records)

    # Sort by pocket_iptm
    if 'pocket_iptm' in df.columns:
        df = df.sort_values('pocket_iptm', ascending=False)

    # Save full results
    df.to_csv(f"{BASE_DIR}/MASTER_AF3_RESULTS.csv", index=False)
    print(f"Saved to {BASE_DIR}/MASTER_AF3_RESULTS.csv")

    # Print summary
    print("\n" + "="*80)
    print("TOP 30 BY POCKET ipTM (protein-ligand binding)")
    print("="*80)
    cols = ['job_name', 'ligands', 'pocket_iptm', 'iptm', 'ptm', 'protein_length', 'has_rna']
    cols = [c for c in cols if c in df.columns]
    print(df[cols].head(30).to_string())

    # Group by ligand
    if 'ligands' in df.columns:
        print("\n" + "="*80)
        print("AVERAGE pocket_iptm BY LIGAND")
        print("="*80)
        ligand_stats = df.groupby('ligands').agg({
            'pocket_iptm': ['mean', 'std', 'count'],
            'iptm': 'mean'
        }).round(3)
        print(ligand_stats.sort_values(('pocket_iptm', 'mean'), ascending=False).to_string())

    return df

if __name__ == '__main__':
    df = main()
