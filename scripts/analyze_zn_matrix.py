#!/usr/bin/env python3
"""
Zinc-Aware AlphaFold3 Analyzer
Specifically designed to differentiate between Amino Acid binding (Chain B)
and Zinc binding (Chain C) in multi-ligand AF3 runs.
"""

import json
import pandas as pd
from pathlib import Path
import sys

def parse_summary_confidences(json_path):
    """Extract key metrics, specifically separating AA vs Zinc ipTM."""
    try:
        with open(json_path) as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error reading {json_path}: {e}")
        return {}
    
    metrics = {
        'iptm': data.get('iptm', 0.0),
        'ptm': data.get('ptm', 0.0),
        'ranking_score': data.get('ranking_score', 0.0),
    }
    
    # Extract Matrix: Protein(0) vs Others
    if 'chain_pair_iptm' in data:
        matrix = data['chain_pair_iptm']
        n_chains = len(matrix)
        
        # Logic for Zinc Runs: Protein(0) -> AA(1) -> Zinc(2)
        if n_chains >= 3:
            metrics['aa_iptm'] = matrix[0][1]   # Protein vs Amino Acid
            metrics['zinc_iptm'] = matrix[0][2] # Protein vs Zinc
        
        # Logic for Standard Runs: Protein(0) -> AA(1)
        elif n_chains == 2:
            metrics['aa_iptm'] = matrix[0][1]
            metrics['zinc_iptm'] = 0.0
            
        else:
            metrics['aa_iptm'] = 0.0
            metrics['zinc_iptm'] = 0.0
            
    return metrics

def analyze_model(model_dir):
    """Analyze a single AF3 model output directory"""
    model_dir = Path(model_dir)
    model_name = model_dir.name
    
    results = {'model': model_name}
    
    # Scan for summary file (robust search)
    summary_files = list(model_dir.glob('**/*_summary_confidences.json'))
    
    if not summary_files:
        return None

    # Use the first valid summary found
    summary_json = summary_files[0]
    
    # Parse
    metrics = parse_summary_confidences(summary_json)
    results.update(metrics)
    
    return results

def main():
    base_dir = Path.cwd()
    print("Zinc-Aware AlphaFold3 Analyzer")
    print("=" * 100)
    print(f"Scanning directory: {base_dir}")
    
    subdirs = [d for d in base_dir.iterdir() if d.is_dir()]
    results = []
    
    print(f"Found {len(subdirs)} folders. extracting data...")
    
    for d in sorted(subdirs):
        res = analyze_model(d)
        if res:
            results.append(res)
            print(f"  + Processed: {d.name}")
            
    if not results:
        print("No results found.")
        return

    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Sort by Amino Acid ipTM (The biological specificity score)
    if 'aa_iptm' in df.columns:
        df = df.sort_values(by='aa_iptm', ascending=False)
    
    # Save CSV
    output_csv = base_dir / 'zinc_specificity_results.csv'
    df.to_csv(output_csv, index=False)
    
    # Print Table
    print("\n" + "="*115)
    print(f"{'Model Name':<40} | {'AA ipTM (Specificity)':<22} | {'Zinc ipTM':<11} | {'pTM':<6}")
    print("-" * 115)
    
    for _, row in df.iterrows():
        name = row['model']
        if len(name) > 38: name = name[:35] + "..."
        
        aa_score = row.get('aa_iptm', 0.0)
        zn_score = row.get('zinc_iptm', 0.0)
        ptm = row.get('ptm', 0.0)
        
        print(f"{name:<40} | {aa_score:.4f}                 | {zn_score:.4f}      | {ptm:.4f}")
        
    print("="*115)
    print(f"Saved to: {output_csv}")
    print("Column 'AA ipTM' is the Protein-Ligand score.")
    print("Column 'Zinc ipTM' is the Protein-Zinc score (Control).")

if __name__ == '__main__':
    main()
