#!/usr/bin/env python3
"""
Regenerate AF3 JSONs with clean sequences.
- Remove gaps from protein sequences
- Convert DNA (T) to RNA (U) notation
"""

import json
from pathlib import Path
import re

def read_fasta_clean(fasta_path, is_rna=False):
    """Read FASTA and clean for AF3."""
    with open(fasta_path) as f:
        lines = f.readlines()
    
    # Get sequence (skip header)
    seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
    
    # Remove ALL gap characters
    clean_seq = re.sub(r'[-.\*X]', '', seq)
    
    # Convert T→U for RNA
    if is_rna:
        clean_seq = clean_seq.replace('T', 'U')
    
    return clean_seq

# Load and clean sequences
print("Loading and cleaning ancestral sequences...")
shallow_protein = read_fasta_clean('../phase1/results/Anc-ProThrRS.fasta', is_rna=False)
deep_protein = read_fasta_clean('../phase1b/results/Anc-ProThrRS-LUCA.fasta', is_rna=False)
ancestral_rna = read_fasta_clean('../phase1/results/Anc-tRNA-ProThr.fasta', is_rna=True)

print(f"✓ Shallow ancestor: {len(shallow_protein)} aa (gaps removed)")
print(f"✓ Deep (LUCA) ancestor: {len(deep_protein)} aa (gaps removed)")
print(f"✓ Ancestral tRNA: {len(ancestral_rna)} nt (T→U converted)")

# Verify no invalid characters remain
def verify_sequence(seq, seq_type):
    if seq_type == 'protein':
        valid_chars = set('ACDEFGHIKLMNPQRSTVWY')
    else:  # RNA
        valid_chars = set('ACGU')
    
    invalid = set(seq) - valid_chars
    if invalid:
        print(f"  ⚠ WARNING: Found invalid characters: {invalid}")
        return False
    return True

print("\nValidating sequences...")
assert verify_sequence(shallow_protein, 'protein'), "Shallow protein has invalid chars"
assert verify_sequence(deep_protein, 'protein'), "Deep protein has invalid chars"
assert verify_sequence(ancestral_rna, 'rna'), "RNA has invalid chars"
print("✓ All sequences validated!")

# Define ligands
ligands = {
    'PRO': {
        'smiles': 'C1CC(NC1)C(=O)O',
        'ccd': 'PRO'
    },
    'THR': {
        'smiles': 'CC(C(C(=O)O)N)O',
        'ccd': 'THR'
    }
}

# Define all test cases
test_cases = [
    {
        'name': 'shallow_ancestral_pro',
        'protein': shallow_protein,
        'ligand': 'PRO',
        'description': 'Shallow ancestor + Proline (cognate)'
    },
    {
        'name': 'shallow_ancestral_thr',
        'protein': shallow_protein,
        'ligand': 'THR',
        'description': 'Shallow ancestor + Threonine (promiscuity test)'
    },
    {
        'name': 'deep_ancestral_pro',
        'protein': deep_protein,
        'ligand': 'PRO',
        'description': 'Deep (LUCA) ancestor + Proline (cognate)'
    },
    {
        'name': 'deep_ancestral_thr',
        'protein': deep_protein,
        'ligand': 'THR',
        'description': 'Deep (LUCA) ancestor + Threonine (promiscuity test)'
    }
]

output_dir = Path('inputs/af3_jsons')
output_dir.mkdir(parents=True, exist_ok=True)

print(f"\nGenerating {len(test_cases)} clean AF3 JSONs...\n")

for case in test_cases:
    ligand_info = ligands[case['ligand']]
    
    af3_input = {
        "name": case['name'],
        "dialect": "alphafold3",
        "version": 1,
        "sequences": [
            {
                "protein": {
                    "id": ["A"],
                    "sequence": case['protein']
                }
            },
            {
                "rna": {
                    "id": ["B"],
                    "sequence": ancestral_rna
                }
            },
            {
                "ligand": {
                    "id": ["C"],
                    "smiles": ligand_info['smiles'],
                    "ccdCodes": [ligand_info['ccd']]
                }
            }
        ],
        "modelSeeds": [1]
    }
    
    output_path = output_dir / f"{case['name']}.json"
    with open(output_path, 'w') as f:
        json.dump(af3_input, f, indent=2)
    
    print(f"✓ {output_path.name}")
    print(f"  {case['description']}")
    print(f"  Protein: {len(case['protein'])} aa")

print(f"\n{'='*60}")
print("✓ All clean JSONs generated!")
print(f"{'='*60}")
print("\nReady to launch:")
print("  nohup bash run_dual_test.sh > logs/phase2_final.log 2>&1 &")
