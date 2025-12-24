#!/usr/bin/env python3
"""
COMPREHENSIVE AF3 DATA AUDIT
============================
Run this on Boron to verify ALL AF3 outputs against manuscript Table 1 claims.

Usage:
    python audit_all_af3_data.py /storage/kiran-stuff/aaRS

Output:
    - complete_af3_audit.csv: All metrics from all runs
    - table1_verification.csv: Cross-reference with manuscript claims
    - audit_report.txt: Summary of findings
"""

import json
import os
import glob
import sys
from pathlib import Path
from datetime import datetime
import csv

# Table 1 claimed values (from manuscript v3)
TABLE1_CLAIMS = {
    'LUCA_ProRS_Pro_Cat': {'ipTM': 0.75, 'context': 'Catalytic domain'},
    'LUCA_ProRS_Thr_Cat': {'ipTM': 0.62, 'context': 'Catalytic domain', 'cross_reactivity': '83%'},
    'LUCA_ProRS_Pro_Full': {'ipTM': 0.75, 'context': 'Full-length'},
    'LUCA_ProRS_Thr_Full': {'ipTM': 0.67, 'context': 'Full-length', 'cross_reactivity': '90%'},
    'LUCA_ThrRS_Thr_Cat': {'ipTM': 0.89, 'context': 'Catalytic domain'},
    'LUCA_ThrRS_Pro_Cat': {'ipTM': 0.88, 'context': 'Catalytic domain', 'cross_reactivity': '99%'},
    'Euk_ProRS_Pro_Full': {'ipTM': 0.83, 'context': 'Full-length'},
    'Euk_ProRS_Thr_Full': {'ipTM': 0.70, 'context': 'Full-length', 'cross_reactivity': '85%'},
    'Ecoli_ProRS_Pro_Cat': {'ipTM': 0.80, 'context': 'Catalytic domain'},
    'Ecoli_ProRS_Thr_Cat': {'ipTM': 0.78, 'context': 'Catalytic domain', 'cross_reactivity': '97.5%'},
    'Human_ProRS_Pro_Full': {'ipTM': 0.80, 'context': 'Full-length'},
    'Human_ProRS_Thr_Full': {'ipTM': 0.78, 'context': 'Full-length', 'cross_reactivity': '97.5%'},
    'Human_ThrRS_Thr_Cat': {'ipTM': 0.84, 'context': 'Catalytic domain'},
    'Human_ThrRS_Pro_Cat': {'ipTM': 0.57, 'context': 'Catalytic domain', 'cross_reactivity': '68%'},
    'LUCA_ProRS_Edit_Pro': {'ipTM': 0.14, 'context': 'Editing domain'},
    'LUCA_ProRS_Edit_Thr': {'ipTM': 0.45, 'context': 'Editing domain', 'cross_reactivity': '321%'},
    'LUCA_ProRS_Trp_Control': {'ipTM': 0.13, 'context': 'Control'},
    'LUCA_ProRS_Phe_Control': {'ipTM': 0.15, 'context': 'Control'},
}


def find_all_af3_outputs(base_dir):
    """Find all AF3 summary_confidences.json files"""
    patterns = [
        '**/summary_confidences.json',
        '**/*summary_confidences.json',
        '**/*_summary_confidences.json',
    ]
    
    all_files = set()
    for pattern in patterns:
        files = glob.glob(os.path.join(base_dir, pattern), recursive=True)
        all_files.update(files)
    
    return sorted(list(all_files))


def extract_all_metrics(json_path):
    """Extract EVERY metric from a summary_confidences.json file"""
    try:
        with open(json_path) as f:
            data = json.load(f)
    except Exception as e:
        return {'file': json_path, 'error': str(e)}
    
    metrics = {
        'file': json_path,
        'parent_dir': os.path.basename(os.path.dirname(os.path.dirname(json_path))),
        'run_name': os.path.basename(os.path.dirname(json_path)),
        'iptm': data.get('iptm'),
        'ptm': data.get('ptm'),
        'ranking_score': data.get('ranking_score'),
        'fraction_disordered': data.get('fraction_disordered'),
        'has_clash': data.get('has_clash'),
    }
    
    # Chain-level metrics
    chain_iptm = data.get('chain_iptm', [])
    chain_ptm = data.get('chain_ptm', [])
    
    for i, v in enumerate(chain_iptm):
        metrics[f'chain_iptm_{i}'] = v
    
    for i, v in enumerate(chain_ptm):
        metrics[f'chain_ptm_{i}'] = v
    
    # Chain pair matrix (THIS IS THE KEY - gives the 0.75 values!)
    chain_pair_iptm = data.get('chain_pair_iptm', [])
    for i, row in enumerate(chain_pair_iptm):
        for j, v in enumerate(row):
            metrics[f'chain_pair_iptm_{i}_{j}'] = v
    
    chain_pair_pae_min = data.get('chain_pair_pae_min', [])
    for i, row in enumerate(chain_pair_pae_min):
        for j, v in enumerate(row):
            metrics[f'chain_pair_pae_min_{i}_{j}'] = v
    
    return metrics


def analyze_directory_structure(files):
    """Understand the directory structure of all runs"""
    structure = {}
    for f in files:
        parts = f.split('/')
        # Find the key directories
        for i, part in enumerate(parts):
            if 'phase' in part.lower() or 'af3' in part.lower() or 'output' in part.lower():
                key = '/'.join(parts[i:i+2])
                if key not in structure:
                    structure[key] = []
                structure[key].append(f)
                break
    return structure


def identify_run_type(filepath, metrics):
    """Try to identify what type of run this was based on path and metrics"""
    path_lower = filepath.lower()
    
    run_info = {
        'enzyme': 'unknown',
        'substrate': 'unknown', 
        'context': 'unknown',
        'age': 'unknown'
    }
    
    # Enzyme
    if 'prors' in path_lower or 'pro_rs' in path_lower or '_pro' in path_lower:
        run_info['enzyme'] = 'ProRS'
    elif 'thrrs' in path_lower or 'thr_rs' in path_lower:
        run_info['enzyme'] = 'ThrRS'
    
    # Substrate
    if '_pro.' in path_lower or '_pro_' in path_lower or 'proline' in path_lower:
        run_info['substrate'] = 'Pro'
    elif '_thr.' in path_lower or '_thr_' in path_lower or 'threonine' in path_lower:
        run_info['substrate'] = 'Thr'
    elif 'trp' in path_lower or 'tryptophan' in path_lower:
        run_info['substrate'] = 'Trp'
    elif 'phe' in path_lower or 'phenylalanine' in path_lower:
        run_info['substrate'] = 'Phe'
    
    # Context
    if 'full' in path_lower or 'fulllength' in path_lower:
        run_info['context'] = 'Full-length'
    elif 'domain' in path_lower or 'cat' in path_lower or 'catalytic' in path_lower:
        run_info['context'] = 'Catalytic'
    elif 'edit' in path_lower:
        run_info['context'] = 'Editing'
    
    # Age
    if 'luca' in path_lower or 'deep' in path_lower or 'ancestral' in path_lower:
        run_info['age'] = 'LUCA'
    elif 'shallow' in path_lower or 'euk' in path_lower:
        run_info['age'] = 'Eukaryotic'
    elif 'ecoli' in path_lower or 'modern' in path_lower:
        run_info['age'] = 'Modern'
    elif 'human' in path_lower:
        run_info['age'] = 'Human'
    
    return run_info


def find_matching_claims(metrics, run_info):
    """Try to match this run to a Table 1 claim"""
    matches = []
    
    # Get the likely "binding" metric (chain_pair_iptm[0][2] if 3 chains)
    binding_metric = metrics.get('chain_pair_iptm_0_2')
    overall_iptm = metrics.get('iptm')
    
    for claim_name, claim_data in TABLE1_CLAIMS.items():
        claimed_iptm = claim_data['ipTM']
        
        # Check if this run's metrics are close to the claim
        if binding_metric is not None:
            if abs(binding_metric - claimed_iptm) < 0.05:
                matches.append({
                    'claim': claim_name,
                    'claimed_value': claimed_iptm,
                    'found_value': binding_metric,
                    'metric': 'chain_pair_iptm_0_2',
                    'delta': abs(binding_metric - claimed_iptm)
                })
        
        if overall_iptm is not None:
            if abs(overall_iptm - claimed_iptm) < 0.05:
                matches.append({
                    'claim': claim_name,
                    'claimed_value': claimed_iptm,
                    'found_value': overall_iptm,
                    'metric': 'overall_iptm',
                    'delta': abs(overall_iptm - claimed_iptm)
                })
    
    return matches


def main(base_dir):
    print("=" * 80)
    print("COMPREHENSIVE AF3 DATA AUDIT")
    print(f"Base directory: {base_dir}")
    print(f"Timestamp: {datetime.now().isoformat()}")
    print("=" * 80)
    
    # Find all files
    print("\n[1/5] Finding all AF3 output files...")
    files = find_all_af3_outputs(base_dir)
    print(f"Found {len(files)} summary_confidences.json files")
    
    # Analyze structure
    print("\n[2/5] Analyzing directory structure...")
    structure = analyze_directory_structure(files)
    for key, file_list in structure.items():
        print(f"  {key}: {len(file_list)} files")
    
    # Extract all metrics
    print("\n[3/5] Extracting all metrics from all files...")
    all_metrics = []
    for f in files:
        metrics = extract_all_metrics(f)
        run_info = identify_run_type(f, metrics)
        metrics.update(run_info)
        all_metrics.append(metrics)
    
    # Save raw audit data
    if all_metrics:
        # Determine all columns
        all_cols = set()
        for m in all_metrics:
            all_cols.update(m.keys())
        all_cols = sorted(list(all_cols))
        
        with open('complete_af3_audit.csv', 'w', newline='') as csvf:
            writer = csv.DictWriter(csvf, fieldnames=all_cols)
            writer.writeheader()
            writer.writerows(all_metrics)
        print(f"  Saved: complete_af3_audit.csv")
    
    # Match to Table 1 claims
    print("\n[4/5] Cross-referencing with Table 1 claims...")
    verification_results = []
    
    for metrics in all_metrics:
        run_info = {k: metrics.get(k) for k in ['enzyme', 'substrate', 'context', 'age']}
        matches = find_matching_claims(metrics, run_info)
        
        for match in matches:
            verification_results.append({
                'file': metrics['file'],
                'run_info': str(run_info),
                **match
            })
    
    if verification_results:
        with open('table1_verification.csv', 'w', newline='') as csvf:
            writer = csv.DictWriter(csvf, fieldnames=verification_results[0].keys())
            writer.writeheader()
            writer.writerows(verification_results)
        print(f"  Saved: table1_verification.csv")
    
    # Generate report
    print("\n[5/5] Generating audit report...")
    
    report_lines = [
        "=" * 80,
        "AF3 DATA AUDIT REPORT",
        f"Generated: {datetime.now().isoformat()}",
        f"Base directory: {base_dir}",
        "=" * 80,
        "",
        "SUMMARY",
        "-" * 40,
        f"Total AF3 output files found: {len(files)}",
        f"Unique run directories: {len(structure)}",
        "",
        "DIRECTORY STRUCTURE",
        "-" * 40,
    ]
    
    for key, file_list in structure.items():
        report_lines.append(f"  {key}: {len(file_list)} samples")
    
    report_lines.extend([
        "",
        "CRITICAL METRICS SUMMARY",
        "-" * 40,
    ])
    
    # Summarize key findings
    low_iptm_runs = [m for m in all_metrics if m.get('iptm') and m['iptm'] < 0.5]
    high_iptm_runs = [m for m in all_metrics if m.get('iptm') and m['iptm'] >= 0.5]
    
    report_lines.extend([
        f"Runs with ipTM < 0.5 (FAILED): {len(low_iptm_runs)}",
        f"Runs with ipTM >= 0.5 (OK): {len(high_iptm_runs)}",
        "",
    ])
    
    # Check chain_pair_iptm values
    chain_pair_runs = [m for m in all_metrics if m.get('chain_pair_iptm_0_2')]
    if chain_pair_runs:
        report_lines.extend([
            f"Runs with chain_pair_iptm data: {len(chain_pair_runs)}",
            "",
            "RUNS MATCHING TABLE 1 VALUES (using chain_pair_iptm[0][2]):",
            "-" * 40,
        ])
        
        for m in chain_pair_runs:
            val = m.get('chain_pair_iptm_0_2')
            if val and 0.5 <= val <= 0.95:
                report_lines.append(f"  {m['file']}")
                report_lines.append(f"    chain_pair_iptm[0][2] = {val:.3f}")
                report_lines.append(f"    overall iptm = {m.get('iptm', 'N/A')}")
    
    report_lines.extend([
        "",
        "POTENTIAL ISSUES",
        "-" * 40,
    ])
    
    # Flag issues
    if low_iptm_runs:
        report_lines.append("⚠️  Multiple runs have ipTM < 0.5 (indicates failed docking)")
        for m in low_iptm_runs[:5]:  # Show first 5
            report_lines.append(f"    - {m['file']}: ipTM = {m.get('iptm', 'N/A')}")
    
    # Check for full-length LUCA runs
    luca_full = [m for m in all_metrics 
                 if 'deep' in m.get('file', '').lower() 
                 and 'full' in m.get('file', '').lower()]
    if luca_full:
        report_lines.extend([
            "",
            "LUCA FULL-LENGTH RUNS (critical for manuscript):",
        ])
        for m in luca_full:
            report_lines.append(f"  {m['file']}")
            report_lines.append(f"    ipTM = {m.get('iptm', 'N/A')}")
            report_lines.append(f"    chain_pair_iptm[0][2] = {m.get('chain_pair_iptm_0_2', 'N/A')}")
    
    report_lines.extend([
        "",
        "=" * 80,
        "END OF AUDIT REPORT",
        "=" * 80,
    ])
    
    report_text = '\n'.join(report_lines)
    
    with open('audit_report.txt', 'w') as f:
        f.write(report_text)
    print(f"  Saved: audit_report.txt")
    
    # Print report to console
    print("\n" + report_text)
    
    return all_metrics, verification_results


if __name__ == '__main__':
    if len(sys.argv) > 1:
        base_dir = sys.argv[1]
    else:
        base_dir = '/storage/kiran-stuff/aaRS'
    
    main(base_dir)
