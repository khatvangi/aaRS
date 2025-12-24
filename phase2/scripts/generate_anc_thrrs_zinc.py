import json
import os

# CONFIGURATION
# We need to extract the sequence from an existing file first
source_file = "anc_thrrs_cat_VAL.json" 
output_prefix = "anc_thrrs_cat_zn"

amino_acids = [
    "ALA", "ARG", "ASN", "ASP", "CYS", 
    "GLN", "GLU", "GLY", "HIS", "ILE", 
    "LEU", "LYS", "MET", "PHE", "PRO", 
    "SER", "THR", "TRP", "TYR", "VAL"
]

# 1. Get the Sequence
if not os.path.exists(source_file):
    print(f"❌ Error: Could not find {source_file} to extract sequence.")
    exit(1)

with open(source_file, 'r') as f:
    data = json.load(f)
    # Extract protein sequence (Chain A)
    protein_seq = data['sequences'][0]['protein']['sequence']
    print(f"✅ Extracted Ancestral Sequence ({len(protein_seq)} residues)")

# 2. Generate Zinc-Included Inputs
print(f"Generating 20 inputs with Zinc...")

for aa in amino_acids:
    run_name = f"{output_prefix}_{aa}"
    filename = f"{run_name}.json"
    
    new_data = {
        "name": run_name,
        "dialect": "alphafold3",
        "version": 1,
        "modelSeeds": [1],
        "sequences": [
            {
                "protein": {
                    "id": ["A"],
                    "sequence": protein_seq
                }
            },
            {
                "ligand": {
                    "id": ["B"],
                    "ccdCodes": [aa]
                }
            },
            {
                "ligand": {
                    "id": ["C"],
                    "ccdCodes": ["ZN"]
                }
            }
        ]
    }
    
    with open(filename, 'w') as f:
        json.dump(new_data, f, indent=2)

print("✅ Done! created 20 files starting with 'anc_thrrs_cat_zn_'")
