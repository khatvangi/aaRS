#!/usr/bin/env python3
"""Generate verification CSV showing all data sources."""

import json
import os
from pathlib import Path
import pandas as pd

base_path = Path('/storage/kiran-stuff/aaRS/phase2')

models = {
    'deep_domain_pro': {'enzyme': 'LUCA ProRS', 'ligand': 'PRO', 'type': 'domain'},
    'deep_domain_thr': {'enzyme': 'LUCA ProRS', 'ligand': 'THR', 'type': 'domain'},
    'deep_catalytic_trp': {'enzyme': 'LUCA ProRS', 'ligand': 'TRP', 'type': 'domain'},
    'deep_catalytic_phe': {'enzyme': 'LUCA ProRS', 'ligand': 'PHE', 'type': 'domain'},
    'deep_editing_pro': {'enzyme': 'LUCA ProRS Editing', 'ligand': 'PRO', 'type': 'domain'},
    'deep_editing_thr': {'enzyme': 'LUCA ProRS Editing', 'ligand': 'THR', 'type': 'domain'},
    'deep_thrrs_pro': {'enzyme': 'LUCA ThrRS', 'ligand': 'PRO', 'type': 'domain'},
    'deep_thrrs_thr': {'enzyme': 'LUCA ThrRS', 'ligand': 'THR', 'type': 'domain'},
    'shallow_domain_pro': {'enzyme': 'Eukaryotic ProRS', 'ligand': 'PRO', 'type': 'domain'},
    'shallow_domain_thr': {'enzyme': 'Eukaryotic ProRS', 'ligand': 'THR', 'type': 'domain'},
    'modern_prours_pro': {'enzyme': 'Modern ProRS', 'ligand': 'PRO', 'type': 'domain'},
    'modern_prours_thr': {'enzyme': 'Modern ProRS', 'ligand': 'THR', 'type': 'domain'},
    'modern_thrrs_pro': {'enzyme': 'Modern ThrRS', 'ligand': 'PRO', 'type': 'domain'},
    'modern_thrrs_thr': {'enzyme': 'Modern ThrRS', 'ligand': 'THR', 'type': 'domain'},
    'fulllength_deep_pro': {'enzyme': 'LUCA ProRS', 'ligand': 'PRO', 'type': 'fulllength'},
    'fulllength_deep_thr': {'enzyme': 'LUCA ProRS', 'ligand': 'THR', 'type': 'fulllength'},
    'fulllength_shallow_pro': {'enzyme': 'Eukaryotic ProRS', 'ligand': 'PRO', 'type': 'fulllength'},
    'fulllength_shallow_thr': {'enzyme': 'Eukaryotic ProRS', 'ligand': 'THR', 'type': 'fulllength'},
}

results = []

for model_name, info in models.items():
    paths_to_try = [
        base_path / 'outputs' / model_name / 'seed-1_sample-0' / f'{model_name}_seed-1_sample-0_summary_confidences.json',
        base_path / 'outputs' / model_name / model_name / 'seed-1_sample-0' / f'{model_name}_seed-1_sample-0_summary_confidences.json',
        base_path / 'af3_output_full' / model_name / 'seed-1_sample-0' / f'{model_name}_seed-1_sample-0_summary_confidences.json',
    ]

    found = False
    for json_path in paths_to_try:
        if json_path.exists():
            with open(json_path) as f:
                data = json.load(f)

            results.append({
                'model_name': model_name,
                'enzyme': info['enzyme'],
                'ligand': info['ligand'],
                'model_type': info['type'],
                'json_file_path': str(json_path),
                'iptm': data.get('iptm', 'N/A'),
                'ptm': data.get('ptm', 'N/A'),
                'chain_iptm': str(data.get('chain_iptm', 'N/A')),
                'chain_ptm': str(data.get('chain_ptm', 'N/A')),
            })
            found = True
            break

    if not found:
        results.append({
            'model_name': model_name,
            'enzyme': info['enzyme'],
            'ligand': info['ligand'],
            'model_type': info['type'],
            'json_file_path': 'NOT FOUND',
            'iptm': 'N/A',
            'ptm': 'N/A',
            'chain_iptm': 'N/A',
            'chain_ptm': 'N/A',
        })

df = pd.DataFrame(results)
output_path = '/storage/kiran-stuff/aaRS/figures/data_source_verification.csv'
df.to_csv(output_path, index=False)
print(f"Saved: {output_path}")
print("\n" + "="*80)
print("DATA SOURCE VERIFICATION")
print("="*80)
print(df.to_string(index=False))
