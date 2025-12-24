#!/usr/bin/env python3
"""
Universal AlphaFold3 Analyzer
Scans current directory for ANY AF3 run and extracts key metrics.
Specifically targets the Protein-Ligand Interface (chain_pair_iptm).
"""

import json
import pandas as pd
from pathlib import Path
import sys
import glob

def parse_summary_confidences(json_path):
    """Extract key confidence metrics from summary_confidences.json"""
    try:
        with open(json_path) as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error reading {json_path}: {e}")
        return {}
    
    metrics = {
        'iptm': data.get('iptm', None),
        'ptm': data.get('ptm', None),
        'ranking_score': data.get('ranking_score', None),
        'fraction_disordered': data.get('fraction_disordered', None),
        'has_clash': data.get('has_clash', None),
    }
    
    # Extract detailed matrices
    if 'chain_ptm' in data:
        metrics['chain_ptm'] = str(data['chain_ptm'])
    if 'chain_iptm' in data:
        metrics['chain_iptm'] = str(data['chain_iptm'])
        
    # CRITICAL: Extract Protein-Ligand Interaction specifically
    # Assuming standard AF3 input order: Protein (0) -> tRNA (1) -> Ligand (2)
    # If no tRNA, Ligand might be (1). We try to grab the last interaction.
    if 'chain_pair_iptm' in data:
        matrix = data['chain_pair_iptm']
        metrics['chain_pair_iptm'] = str(matrix)
        
        try:
            # Attempt to grab Protein(0) vs Last Chain (Ligand)
            # This is the "Pocket Score"
            n_chains = len(matrix)
            if n_chains >= 2:
                # Interaction between first chain (Protein) and last chain (Ligand)
                metrics['pocket_iptm'] = matrix[0][-1]
            else:
                metrics['pocket_iptm'] = 0.0
        except:
            metrics['pocket_iptm'] = 0.0
    
    return metrics

def parse_confidences(json_path):
    """Extract detailed per-residue confidence scores (pLDDT)"""
    try:
        with open(json_path) as f:
            data = json.load(f)
        
        # Try different possible field names
        plddt = None
        if 'confidenceScore' in data:
            plddt = data['confidenceScore']
        elif 'plddt' in data:
            plddt = data['plddt']
        elif 'atom_plddts' in data:
            plddt = data['atom_plddts']
        
        if plddt and len(plddt) > 0:
            return {
                'mean_plddt': sum(plddt) / len(plddt),
                'min_plddt': min(plddt),
                'max_plddt': max(plddt),
            }
    except:
        pass
    
    return {'mean_plddt': None, 'min_plddt': None, 'max_plddt': None}

def analyze_model(model_dir):
    """Analyze a single AF3 model output directory"""
    model_dir = Path(model_dir)
    model_name = model_dir.name
    
    results = {
        'model': model_name,
    }
    
    # We look for the "best" seed (usually seed-1_sample-0 or similar)
    # But since folder structures vary, we scan for the summary file directly
    summary_files = list(model_dir.glob('**/seed-1_sample-0/*_summary_confidences.json'))
    if not summary_files:
        # Try finding ANY summary file if standard structure doesn't exist
        summary_files = list(model_dir.glob('**/*_summary_confidences.json'))
        
    if not summary_files:
        # print(f"Skipping {model_name}: No summary JSON found.")
        return None

    # Use the first valid summary file found (usually seed 1 sample 0)
    summary_json = summary_files[0]
    # Deduce the confidences file path based on summary path
    conf_json = summary_json.parent / summary_json.name.replace('_summary_confidences.json', '_confidences.json')
    
    # Parse Summary
    summary_metrics = parse_summary_confidences(summary_json)
    results.update(summary_metrics)
    
    # Parse pLDDT
    if conf_json.exists():
        conf_metrics = parse_confidences(conf_json)
        results.update(conf_metrics)
    
    return results

def main():
    # Detect current directory
    base_dir = Path.cwd()
    
    print("AlphaFold3 Universal Analyzer")
    print("=" * 80)
    print(f"Scanning directory: {base_dir}")
    
    # Find all subdirectories
    subdirs = [d for d in base_dir.iterdir() if d.is_dir()]
    
    results = []
    print(f"Found {len(subdirs)} folders. Scanning for AF3 data...")
    
    for d in sorted(subdirs):
        res = analyze_model(d)
        if res:
            results.append(res)
            print(f"  + Processed: {d.name}")
            
    if not results:
        print("No valid AF3 output folders found containing summary_confidences.json")
        return

    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Sort by Pocket Score (if available) or Global ipTM
    sort_col = 'pocket_iptm' if 'pocket_iptm' in df.columns else 'iptm'
    df = df.sort_values(by=sort_col, ascending=False)
    
    # Save CSV
    output_csv = base_dir / 'comprehensive_af3_results.csv'
    df.to_csv(output_csv, index=False)
    
    # Print Summary Table to Console
    print("\n" + "="*100)
    print(f"{'Model Name':<40} | {'Pocket ipTM':<11} | {'Global ipTM':<11} | {'pTM':<6}")
    print("-" * 100)
    
    for _, row in df.iterrows():
        name = row['model']
        # Truncate name if too long
        if len(name) > 38: name = name[:35] + "..."
        
        pocket = row.get('pocket_iptm', 0.0)
        iptm = row.get('iptm', 0.0)
        ptm = row.get('ptm', 0.0)
        
        # Handle None types
        pocket = 0.0 if pocket is None else pocket
        iptm = 0.0 if iptm is None else iptm
        ptm = 0.0 if ptm is None else ptm
        
        print(f"{name:<40} | {pocket:.4f}      | {iptm:.4f}      | {ptm:.4f}")
        
    print("="*100)
    print(f"Saved full details to: {output_csv}")
    print(f"Note: 'Pocket ipTM' is chain_pair_iptm[0][-1] (Protein vs Ligand)")

if __name__ == '__main__':
    main()
