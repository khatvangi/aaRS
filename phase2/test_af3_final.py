#!/usr/bin/env python3
"""AF3 test with correct JSON format (matching your sample)"""

import json
from pathlib import Path
from Bio import SeqIO

anc_aars = str(SeqIO.read("../phase1/results/Anc-ProThrRS.fasta", "fasta").seq)
anc_trna = str(SeqIO.read("../phase1/results/Anc-tRNA-ProThr.fasta", "fasta").seq).replace('T', 'U')

test_aars = anc_aars[:200]

# CORRECT FORMAT (matching your working example)
test_input = {
    "name": "test_aars_trna_complex",
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
    ],
    "modelSeeds": [1],
    "dialect": "alphafold3",
    "version": 1
}

Path("test_input").mkdir(exist_ok=True)
with open("test_input/test_final.json", 'w') as f:
    json.dump(test_input, f, indent=2)

print("✓ Created test_input/test_final.json")
print(f"  aaRS: {len(test_aars)} aa")
print(f"  tRNA: {len(anc_trna)} nt")
print("\nRun with: cd test_input && af3 --json_path=test_final.json")
