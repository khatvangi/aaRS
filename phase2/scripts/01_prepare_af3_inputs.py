#!/usr/bin/env python3
"""Generate AlphaFold3 input JSON files for all test cases"""

import json
from pathlib import Path
from Bio import SeqIO

def create_af3_input(name, protein_seq, rna_seq, ligand_smiles, ligand_name):
    """Create AF3 JSON input format"""
    return {
        "name": name,
        "modelSeeds": [1],
        "sequences": [
            {
                "proteinChain": {
                    "sequence": protein_seq,
                    "count": 1
                }
            },
            {
                "rnaChain": {
                    "sequence": rna_seq,
                    "count": 1
                }
            },
            {
                "ligand": {
                    "smiles": ligand_smiles,
                    "ccdCodes": [ligand_name]
                }
            }
        ]
    }

# Amino acid SMILES
AA_SMILES = {
    "Pro": "C1CC(NC1)C(=O)O",  # Proline
    "Thr": "CC(C(C(=O)O)N)O"   # Threonine
}

# Load ancestral sequences
anc_aars = str(SeqIO.read("inputs/Anc-ProThrRS.fasta", "fasta").seq)
anc_trna = str(SeqIO.read("inputs/Anc-tRNA-ProThr.fasta", "fasta").seq).replace('T', 'U')

# TODO: You need to add modern sequences here
# For now, using ancestral as placeholder
modern_prors = anc_aars  # Replace with real modern ProRS
modern_trna_pro = anc_trna  # Replace with real modern tRNA-Pro

test_cases = [
    # Ancestral (should be promiscuous)
    ("ancestral_pro", anc_aars, anc_trna, AA_SMILES["Pro"], "PRO"),
    ("ancestral_thr", anc_aars, anc_trna, AA_SMILES["Thr"], "THR"),
    
    # Modern (should be specific)
    ("modern_pro_cognate", modern_prors, modern_trna_pro, AA_SMILES["Pro"], "PRO"),
    ("modern_thr_noncognate", modern_prors, modern_trna_pro, AA_SMILES["Thr"], "THR"),
]

Path("inputs/af3_jsons").mkdir(exist_ok=True)

for name, prot, rna, smiles, ccd in test_cases:
    input_data = create_af3_input(name, prot, rna, smiles, ccd)
    
    with open(f"inputs/af3_jsons/{name}.json", 'w') as f:
        json.dump(input_data, f, indent=2)
    
    print(f"✓ Created {name}.json")

print(f"\n✓ Generated {len(test_cases)} AF3 input files")
print("⚠ WARNING: Modern sequences are placeholders - update with real sequences!")
