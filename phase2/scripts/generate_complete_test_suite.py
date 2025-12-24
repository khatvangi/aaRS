#!/usr/bin/env python3
"""
Generate complete test suite for publication-quality validation.
"""

import json, re
from pathlib import Path

def clean_fasta(path, is_rna=False):
    with open(path) as f:
        seq = ''.join(l.strip() for l in f if not l.startswith('>'))
    seq = re.sub(r'[-.\*X]', '', seq)
    if is_rna:
        seq = seq.replace('T', 'U')
    return seq

# Load ancestral sequences
shallow_full = clean_fasta('../phase1/results/Anc-ProThrRS.fasta')
deep_full = clean_fasta('../phase1b/results/Anc-ProThrRS-LUCA.fasta')
deep_thrrs = clean_fasta('../phase1b/results/Anc-ThrRS-LUCA.fasta')
rna = clean_fasta('../phase1/results/Anc-tRNA-ProThr.fasta', True)

# Extract domains
shallow_catalytic = shallow_full[200:700]
deep_catalytic = deep_full[200:700]
deep_editing = deep_full[1400:1700]
deep_thrrs_catalytic = deep_thrrs[600:878]  # Class II catalytic from Pfam

print("="*70)
print("COMPLETE TEST SUITE GENERATION")
print("="*70)

# Modern sequences (you'll need to download these)
modern_prours = "PASTE_HUMAN_PRORS_SEQUENCE_HERE"  # UniProt P07814
modern_thrrs = "PASTE_HUMAN_THRRS_SEQUENCE_HERE"  # UniProt P26639

tests = [
    # P2: Modern controls
    {'name': 'modern_prours_pro', 'protein': modern_prours[200:700] if len(modern_prours) > 700 else None, 
     'ligand': 'PRO', 'desc': 'Modern control (cognate)'},
    {'name': 'modern_prours_thr', 'protein': modern_prours[200:700] if len(modern_prours) > 700 else None, 
     'ligand': 'THR', 'desc': 'Modern control (should NOT bind)'},
    
    # P3: Editing domain
    {'name': 'deep_editing_pro', 'protein': deep_editing, 'ligand': 'PRO', 
     'desc': 'Editing domain + PRO'},
    {'name': 'deep_editing_thr', 'protein': deep_editing, 'ligand': 'THR', 
     'desc': 'Editing domain + THR (test discrimination)'},
    
    # P4: Reverse test
    {'name': 'deep_thrrs_pro', 'protein': deep_thrrs_catalytic, 'ligand': 'PRO', 
     'desc': 'ThrRS + PRO (reverse promiscuity test)'},
    {'name': 'deep_thrrs_thr', 'protein': deep_thrrs_catalytic, 'ligand': 'THR', 
     'desc': 'ThrRS + THR (cognate control)'},
    
    # P5: Negative controls
    {'name': 'deep_catalytic_trp', 'protein': deep_catalytic, 'ligand': 'TRP', 
     'desc': 'Negative control (unrelated AA)'},
    {'name': 'deep_catalytic_phe', 'protein': deep_catalytic, 'ligand': 'PHE', 
     'desc': 'Negative control (unrelated AA)'},
]

output_dir = Path('inputs/af3_jsons_validation')
output_dir.mkdir(exist_ok=True)

print("\nGenerating validation test JSONs...\n")

for t in tests:
    if t['protein'] is None:
        print(f"⏭️  {t['name']}: SKIPPED (need modern sequence)")
        continue
    
    if len(t['protein']) < 50:
        print(f"⏭️  {t['name']}: SKIPPED (sequence too short)")
        continue
    
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
    
    print(f"✓ {t['name']}: {len(t['protein'])} aa - {t['desc']}")

print("\n" + "="*70)
print("NEXT STEPS")
print("="*70)
print("""
1. Download modern sequences:
   - Human ProRS: https://www.uniprot.org/uniprotkb/P07814/entry
   - Human ThrRS: https://www.uniprot.org/uniprotkb/P26639/entry
   
2. Paste sequences into script (replace PASTE_HERE)

3. Rerun to generate modern controls

4. Launch validation suite:
   bash run_validation_suite.sh
""")
