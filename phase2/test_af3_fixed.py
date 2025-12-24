#!/usr/bin/env python3
"""Test AF3 with correct JSON format"""

import json
from pathlib import Path
from Bio import SeqIO

# Load ancestral sequences
anc_aars = str(SeqIO.read("../phase1/results/Anc-ProThrRS.fasta", "fasta").seq)
anc_trna = str(SeqIO.read("../phase1/results/Anc-tRNA-ProThr.fasta", "fasta").seq).replace('T', 'U')

# Truncate for faster test
test_aars = anc_aars[:200]

# CORRECT AF3 JSON FORMAT (v1.0 dialect)
test_input = {
    "name": "test_aars_trna",
    "dialect": "alphafold3",
    "version": 1,
    "sequences": [
        {
            "protein": {
                "id": ["A"],
                "sequence": test_aars
            }
        },
        {
            "rna": {
                "id": ["B"],
                "sequence": anc_trna
            }
        }
    ]
}

Path("test_input").mkdir(exist_ok=True)
with open("test_input/test_fixed.json", 'w') as f:
    json.dump(test_input, f, indent=2)

print(f"✓ Created test_input/test_fixed.json")
print(f"  aaRS: {len(test_aars)} aa")
print(f"  tRNA: {len(anc_trna)} nt")
