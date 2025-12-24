#!/usr/bin/env python3
"""Generate validation tests (all can run locally - no OOM)."""

import json, re
from pathlib import Path

def clean_fasta(path, is_rna=False):
    with open(path) as f:
        seq = ''.join(l.strip() for l in f if not l.startswith('>'))
    seq = re.sub(r'[-.\*X]', '', seq)
    if is_rna:
        seq = seq.replace('T', 'U')
    return seq

# Ancestral sequences
deep_full = clean_fasta('../phase1b/results/Anc-ProThrRS-LUCA.fasta')
deep_thrrs = clean_fasta('../phase1b/results/Anc-ThrRS-LUCA.fasta')
rna = clean_fasta('../phase1/results/Anc-tRNA-ProThr.fasta', True)

# Modern sequences
modern_prours = clean_fasta('inputs/Modern_Human_ProRS.fasta')
modern_thrrs = clean_fasta('inputs/Modern_Human_ThrRS.fasta')

print(f"Loaded sequences:")
print(f"  Deep ProRS: {len(deep_full)} aa")
print(f"  Deep ThrRS: {len(deep_thrrs)} aa")
print(f"  Modern ProRS: {len(modern_prours)} aa")
print(f"  Modern ThrRS: {len(modern_thrrs)} aa")

# Extract domains
deep_editing = deep_full[1400:1700]  # 300 aa
deep_thrrs_cat = deep_thrrs[600:878]  # 278 aa
modern_pro_cat = modern_prours[200:700] if len(modern_prours) > 700 else modern_prours[:500]
modern_thr_cat = modern_thrrs[200:700] if len(modern_thrrs) > 700 else modern_thrrs[:500]
deep_cat = deep_full[200:700]

print("\n" + "="*70)
print("VALIDATION SUITE GENERATION")
print("="*70)

tests = [
    # Modern controls
    {'name': 'modern_prours_pro', 'protein': modern_pro_cat, 'ligand': 'PRO', 
     'priority': 'P1', 'desc': 'Modern ProRS + PRO (cognate)'},
    {'name': 'modern_prours_thr', 'protein': modern_pro_cat, 'ligand': 'THR', 
     'priority': 'P1', 'desc': 'Modern ProRS + THR (specificity test)'},
    {'name': 'modern_thrrs_thr', 'protein': modern_thr_cat, 'ligand': 'THR', 
     'priority': 'P1', 'desc': 'Modern ThrRS + THR (cognate)'},
    {'name': 'modern_thrrs_pro', 'protein': modern_thr_cat, 'ligand': 'PRO', 
     'priority': 'P1', 'desc': 'Modern ThrRS + PRO (specificity test)'},
    
    # Editing domain
    {'name': 'deep_editing_pro', 'protein': deep_editing, 'ligand': 'PRO', 
     'priority': 'P2', 'desc': 'LUCA editing domain + PRO'},
    {'name': 'deep_editing_thr', 'protein': deep_editing, 'ligand': 'THR', 
     'priority': 'P2', 'desc': 'LUCA editing domain + THR'},
    
    # Reverse test
    {'name': 'deep_thrrs_thr', 'protein': deep_thrrs_cat, 'ligand': 'THR', 
     'priority': 'P3', 'desc': 'LUCA ThrRS + THR (cognate)'},
    {'name': 'deep_thrrs_pro', 'protein': deep_thrrs_cat, 'ligand': 'PRO', 
     'priority': 'P3', 'desc': 'LUCA ThrRS + PRO (reverse promiscuity)'},
    
    # Negative controls
    {'name': 'deep_cat_trp', 'protein': deep_cat, 'ligand': 'TRP', 
     'priority': 'P4', 'desc': 'Negative control (Trp)'},
    {'name': 'deep_cat_phe', 'protein': deep_cat, 'ligand': 'PHE', 
     'priority': 'P4', 'desc': 'Negative control (Phe)'},
]

output_dir = Path('inputs/af3_jsons_validation')
output_dir.mkdir(exist_ok=True)

print(f"\nGenerating {len(tests)} validation JSONs...\n")

for t in tests:
    data = {
        "name": t['name'],
        "dialect": "alphafold3",
        "version": 1,
        "sequences": [
            {"protein": {"id": ["A"], "sequence": t['protein']}},
            {"rna": {"id": ["B"], "sequence": rna}},
            {"ligand": {"id": ["C"], "ccdCodes": [t['ligand']]}}
        ],
        "modelSeeds": [1]
    }
    
    path = output_dir / f"{t['name']}.json"
    with open(path, 'w') as f:
        json.dump(data, f, indent=2)
    
    print(f"✓ [{t['priority']}] {t['name']}: {len(t['protein'])} aa - {t['desc']}")

print("\n" + "="*70)
print(f"✓ {len(tests)} validation tests ready")
print("  Estimated time: ~4 hours (all local, no OOM)")
print("="*70)
