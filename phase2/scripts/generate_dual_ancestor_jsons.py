#!/usr/bin/env python3
"""
Generate AF3 JSONs for both shallow and deep ancestral reconstructions.
Tests the promiscuity hypothesis at two evolutionary timepoints.
"""

import json
from pathlib import Path

def read_fasta(fasta_path):
    """Read single sequence from FASTA file."""
    with open(fasta_path) as f:
        lines = f.readlines()
    seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return seq

# Load sequences
print("Loading ancestral sequences...")
shallow_protein = read_fasta('../phase1/results/Anc-ProThrRS.fasta')
deep_protein = read_fasta('../phase1b/results/Anc-ProThrRS-LUCA.fasta')
ancestral_rna = read_fasta('../phase1/results/Anc-tRNA-ProThr.fasta')

print(f"✓ Shallow ancestor: {len(shallow_protein)} aa")
print(f"✓ Deep (LUCA) ancestor: {len(deep_protein)} aa")
print(f"✓ Ancestral tRNA: {len(ancestral_rna)} nt")

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

# Define all test cases (2 ancestors × 2 ligands = 4 tests)
test_cases = [
    # Shallow ancestor tests
    {
        'name': 'shallow_ancestral_pro',
        'protein': shallow_protein,
        'ligand': 'PRO',
        'description': 'Shallow (eukaryotic) ancestor + Proline (cognate)'
    },
    {
        'name': 'shallow_ancestral_thr',
        'protein': shallow_protein,
        'ligand': 'THR',
        'description': 'Shallow (eukaryotic) ancestor + Threonine (promiscuity test)'
    },
    # Deep (LUCA) ancestor tests
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

print(f"\nGenerating {len(test_cases)} AF3 input JSONs...\n")

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

print(f"\n{'='*60}")
print("✓ All dual-ancestor JSONs generated!")
print(f"{'='*60}")
print("\nReady to launch Phase 2 AF3 modeling:")
print("  nohup bash run_dual_test.sh > logs/phase2_dual.log 2>&1 &")
