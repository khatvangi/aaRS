import json
import os

# ---------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------
output_dir = "af3_inputs_ancient"
os.makedirs(output_dir, exist_ok=True)

# 1. Sources for Full Length (JSONs that we know work)
json_sources = {
    "anc_prors_full": "af3_jsons_fulllength/fulllength_deep_pro.json",
    "anc_thrrs_full": "af3_jsons_fulllength/fulllength_deep_thr.json"
}

# 2. Sources for Domains (FASTAs are safer)
# Based on your tree structure
fasta_sources = {
    "anc_prors_cat":  "domain_fasta/cat_LUCA.fasta",
    "anc_prors_edit": "domain_fasta/edit_LUCA.fasta",
    "anc_thrrs_cat":  "domain_fasta/cat_LUCA_Thr.fasta" # Assuming this is ThrRS Cat
}

amino_acids = [
    "ALA", "ARG", "ASN", "ASP", "CYS", 
    "GLN", "GLU", "GLY", "HIS", "ILE", 
    "LEU", "LYS", "MET", "PHE", "PRO", 
    "SER", "THR", "TRP", "TYR", "VAL"
]

sequences = {}

# ---------------------------------------------------------
# EXTRACT SEQUENCES
# ---------------------------------------------------------
print("Extracting sequences...")

# Read JSONs
for name, path in json_sources.items():
    if os.path.exists(path):
        try:
            with open(path, 'r') as f:
                data = json.load(f)
                # Handle list wrapper if present
                if isinstance(data, list): data = data[0]
                seq = data['sequences'][0]['protein']['sequence']
                sequences[name] = seq
                print(f"  ✅ Loaded {name} from JSON")
        except Exception as e:
            print(f"  ❌ Failed {name} (JSON): {e}")
    else:
        print(f"  ⚠️ Missing file: {path}")

# Read FASTAs
def read_fasta(path):
    seq = ""
    with open(path, 'r') as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()
    return seq

for name, path in fasta_sources.items():
    if os.path.exists(path):
        try:
            seq = read_fasta(path)
            if len(seq) > 0:
                sequences[name] = seq
                print(f"  ✅ Loaded {name} from FASTA")
            else:
                print(f"  ❌ Empty FASTA: {path}")
        except Exception as e:
            print(f"  ❌ Failed {name} (FASTA): {e}")
    else:
        # Fallback: Try finding it in the broken JSON by raw string search
        # (Only if FASTA is missing)
        print(f"  ⚠️ FASTA missing: {path}. Skipping.")

# ---------------------------------------------------------
# GENERATE INPUTS
# ---------------------------------------------------------
print("-" * 40)
print(f"Generating inputs in: {output_dir}")
count = 0

for name, seq in sequences.items():
    for aa in amino_acids:
        run_name = f"{name}_{aa}"
        
        data = {
            "name": run_name,
            "dialect": "alphafold3",
            "version": 1,
            "modelSeeds": [1],
            "sequences": [
                { "protein": { "id": ["A"], "sequence": seq } },
                { "ligand": { "id": ["B"], "ccdCodes": [aa] } }
            ]
        }
        
        with open(os.path.join(output_dir, f"{run_name}.json"), 'w') as f:
            json.dump(data, f, indent=2)
        count += 1

print(f"Done! Generated {count} files.")
