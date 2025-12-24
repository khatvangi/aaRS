#!/usr/bin/env python3
"""
Extract catalytic domains (aa 200-700) and generate AF3 JSONs.
This region contains the amino acid binding site for promiscuity testing.
"""

import json
from pathlib import Path
import re

def read_fasta_clean(fasta_path, is_rna=False):
    """Read and clean FASTA."""
    with open(fasta_path) as f:
        lines = f.readlines()
    seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
    seq = re.sub(r'[-.\*X]', '', seq)
    if is_rna:
        seq = seq.replace('T', 'U')
    return seq

# Load sequences
print("Loading sequences...")
shallow_full = read_fasta_clean('../phase1/results/Anc-ProThrRS.fasta')
deep_full = read_fasta_clean('../phase1b/results/Anc-ProThrRS-LUCA.fasta')
ancestral_rna = read_fasta_clean('../phase1/results/Anc-tRNA-ProThr.fasta', is_rna=True)

# Extract catalytic domains (200-700)
START = 200
END = 700

shallow_domain = shallow_full[START:END]
deep_domain = deep_full[START:END]

print(f"✓ Shallow full-length: {len(shallow_full)} aa")
print(f"✓ Deep full-length:    {len(deep_full)} aa")
print(f"✓ Extracted region:    aa {START}-{END}")
print(f"✓ Shallow domain:      {len(shallow_domain)} aa")
print(f"✓ Deep domain:         {len(deep_domain)} aa")
print(f"✓ tRNA:                {len(ancestral_rna)} nt")

# Ligands (CCD codes only)
ligands = {'PRO': 'PRO', 'THR': 'THR'}

# Test cases
test_cases = [
    {'name': 'shallow_domain_pro', 'protein': shallow_domain, 'ligand': 'PRO', 'desc': 'Shallow catalytic domain + PRO'},
    {'name': 'shallow_domain_thr', 'protein': shallow_domain, 'ligand': 'THR', 'desc': 'Shallow catalytic domain + THR (promiscuity test)'},
    {'name': 'deep_domain_pro', 'protein': deep_domain, 'ligand': 'PRO', 'desc': 'Deep catalytic domain + PRO'},
    {'name': 'deep_domain_thr', 'protein': deep_domain, 'ligand': 'THR', 'desc': 'Deep catalytic domain + THR (promiscuity test)'}
]

# Create domain JSONs directory
output_dir = Path('inputs/af3_jsons_domain')
output_dir.mkdir(parents=True, exist_ok=True)

print("\n" + "="*60)
print("GENERATING CATALYTIC DOMAIN JSONS")
print("="*60)

for case in test_cases:
    af3_input = {
        "name": case['name'],
        "dialect": "alphafold3",
        "version": 1,
        "sequences": [
            {"protein": {"id": ["A"], "sequence": case['protein']}},
            {"rna": {"id": ["B"], "sequence": ancestral_rna}},
            {"ligand": {"id": ["C"], "ccdCodes": [ligands[case['ligand']]]}}
        ],
        "modelSeeds": [1]
    }
    
    output_path = output_dir / f"{case['name']}.json"
    with open(output_path, 'w') as f:
        json.dump(af3_input, f, indent=2)
    
    print(f"\n✓ {case['name']}.json")
    print(f"  {case['desc']}")
    print(f"  Size: {len(case['protein'])} aa")

print("\n" + "="*60)
print("✓ ALL DOMAIN JSONS GENERATED!")
print("="*60)
print(f"\nEstimated GPU memory: ~5-8 GB per job (vs 23+ GB for full-length)")
print(f"Estimated time per job: 20-30 minutes (vs 45+ min)")
print(f"Total time for 4 jobs: ~1.5-2 hours")
