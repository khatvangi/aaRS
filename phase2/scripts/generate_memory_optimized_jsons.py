#!/usr/bin/env python3
"""Generate memory-optimized AF3 JSONs with reduced MSA depth."""

import json
from pathlib import Path
import re

def clean_fasta(fasta_path, is_rna=False):
    """Read and clean FASTA for AF3."""
    with open(fasta_path) as f:
        lines = f.readlines()
    seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
    seq = re.sub(r'[-.\*X]', '', seq)
    if is_rna:
        seq = seq.replace('T', 'U')
    return seq

# Load sequences
print("Loading sequences...")
shallow_protein = clean_fasta('../phase1/results/Anc-ProThrRS.fasta')
deep_protein = clean_fasta('../phase1b/results/Anc-ProThrRS-LUCA.fasta')
ancestral_rna = clean_fasta('../phase1/results/Anc-tRNA-ProThr.fasta', is_rna=True)

print(f"✓ Shallow: {len(shallow_protein)} aa")
print(f"✓ Deep: {len(deep_protein)} aa")
print(f"✓ tRNA: {len(ancestral_rna)} nt")

# Ligands (CCD codes only)
ligands = {'PRO': 'PRO', 'THR': 'THR'}

# Test cases
test_cases = [
    {'name': 'shallow_ancestral_pro', 'protein': shallow_protein, 'ligand': 'PRO'},
    {'name': 'shallow_ancestral_thr', 'protein': shallow_protein, 'ligand': 'THR'},
    {'name': 'deep_ancestral_pro', 'protein': deep_protein, 'ligand': 'PRO'},
    {'name': 'deep_ancestral_thr', 'protein': deep_protein, 'ligand': 'THR'}
]

output_dir = Path('inputs/af3_jsons')
output_dir.mkdir(parents=True, exist_ok=True)

# Remove old JSONs
for old_json in output_dir.glob('*.json'):
    old_json.unlink()

print("\nGenerating memory-optimized JSONs...\n")

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
    
    print(f"✓ {case['name']}.json ({len(case['protein'])} aa + {case['ligand']})")

print("\n" + "="*60)
print("✓ Memory-optimized JSONs ready!")
print("="*60)
