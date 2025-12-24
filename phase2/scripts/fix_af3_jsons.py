#!/usr/bin/env python3
"""
Fix AF3 JSON inputs to include required dialect/version fields
and use correct sequence schema format.
"""

import json
from pathlib import Path

# Read the ancestral sequences
def read_fasta(fasta_path):
    """Read single sequence from FASTA file."""
    with open(fasta_path) as f:
        lines = f.readlines()
    seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return seq

# Load sequences
protein_seq = read_fasta('../phase1/results/Anc-ProThrRS.fasta')
rna_seq = read_fasta('../phase1/results/Anc-tRNA-ProThr.fasta')

# Define ligands (using SMILES + CCD codes)
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

# Generate corrected JSONs
test_cases = [
    {
        'name': 'ancestral_pro',
        'ligand': 'PRO',
        'description': 'Ancestral aaRS + tRNA + Proline (cognate)'
    },
    {
        'name': 'ancestral_thr',
        'ligand': 'THR',
        'description': 'Ancestral aaRS + tRNA + Threonine (test promiscuity)'
    }
]

output_dir = Path('inputs/af3_jsons')
output_dir.mkdir(parents=True, exist_ok=True)

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
                    "sequence": protein_seq
                }
            },
            {
                "rna": {
                    "id": ["B"],
                    "sequence": rna_seq
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
    
    print(f"✓ Generated: {output_path}")
    print(f"  Description: {case['description']}")

print("\n✓ All JSONs regenerated with correct schema!")
