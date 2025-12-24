#!/usr/bin/env python3
"""Test AF3 wrapper with a minimal aaRS-tRNA complex"""

import json
from pathlib import Path
from Bio import SeqIO

# Load ancestral sequences
anc_aars = str(SeqIO.read("../phase1/results/Anc-ProThrRS.fasta", "fasta").seq)
anc_trna = str(SeqIO.read("../phase1/results/Anc-tRNA-ProThr.fasta", "fasta").seq).replace('T', 'U')

# Truncate for faster test (first 200 aa of aaRS)
test_aars = anc_aars[:200]

print(f"Test sequences:")
print(f"  aaRS: {len(test_aars)} aa")
print(f"  tRNA: {len(anc_trna)} nt")

# Create minimal AF3 JSON (no ligand for first test)
test_input = {
    "name": "test_aars_trna",
    "modelSeeds": [1],
    "sequences": [
        {
            "proteinChain": {
                "sequence": test_aars,
                "count": 1
            }
        },
        {
            "rnaChain": {
                "sequence": anc_trna,
                "count": 1
            }
        }
    ]
}

Path("test_input").mkdir(exist_ok=True)
with open("test_input/test.json", 'w') as f:
    json.dump(test_input, f, indent=2)

print("✓ Created test_input/test.json")
print("\nRun with:")
print("  cd test_input")
print("  af3 --json_path=test.json")
