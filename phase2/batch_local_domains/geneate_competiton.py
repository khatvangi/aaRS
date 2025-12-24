import json
import os

def create_competition(source_json, protein_name, lig1, lig2):
    if not os.path.exists(source_json):
        print(f"❌ Missing source: {source_json}")
        return

    with open(source_json, 'r') as f:
        data = json.load(f)
        seq = data['sequences'][0]['protein']['sequence']

    # Create the 4-chain setup: Protein, Lig1, Lig2, Zinc
    comp_data = {
        "name": f"COMPETITION_{protein_name}_{lig1}_vs_{lig2}",
        "dialect": "alphafold3",
        "version": 1,
        "modelSeeds": [1],
        "sequences": [
            {"protein": {"id": ["A"], "sequence": seq}},       # Chain A
            {"ligand": {"id": ["B"], "ccdCodes": [lig1]}},     # Chain B (Combatant 1)
            {"ligand": {"id": ["C"], "ccdCodes": [lig2]}},     # Chain C (Combatant 2)
            {"ligand": {"id": ["D"], "ccdCodes": ["ZN"]}}      # Chain D (The Prize)
        ]
    }

    filename = f"COMPETITION_{protein_name}_{lig1}_vs_{lig2}.json"
    with open(filename, 'w') as f:
        json.dump(comp_data, f, indent=2)
    print(f"✅ Created Arena: {filename}")

# --- DEFINING THE MATCHES ---

# 1. Ancestral ThrRS: Threonine (Cognate) vs Isoleucine (Hydrophobic Imposter)
create_competition(
    source_json="anc_thrrs_cat_zn_ILE.json", # Source for sequence
    protein_name="anc_thrrs",
    lig1="THR",
    lig2="ILE"
)

# 2. Modern ThrRS: Threonine (Cognate) vs Isoleucine (Hydrophobic Imposter)
# We need to find your modern sequence file. 
# Assuming you have one nearby or we can use the sequence from the Ancestor 
# (Wait, Modern seq is different. I will assume 'modern_thrrs_ecoli_zn_THR.json' exists in your path or similar).
# If you don't have the modern source file here, copy it in.
create_competition(
    source_json="../af3_modern_matrix/thrrs_ecoli_zn_jobs/modern_thrrs_ecoli_zn_THR.json", 
    protein_name="modern_thrrs",
    lig1="THR",
    lig2="ILE"
)
