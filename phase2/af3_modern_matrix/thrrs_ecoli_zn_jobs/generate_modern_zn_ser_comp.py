import json
import os

# CONFIGURATION
# We use the Modern ThrRS sequence file you likely used for the previous run
source_file = "modern_thrrs_ecoli_zn_THR.json" 

def create_competition(source_json, protein_name, lig1, lig2):
    if not os.path.exists(source_json):
        # Fallback: Try to find any modern file if specific one is missing
        print(f"⚠️ Warning: {source_json} not found. Checking current directory...")
        possible_files = [f for f in os.listdir('.') if 'modern' in f and 'json' in f]
        if possible_files:
            source_json = possible_files[0]
            print(f"✅ Found alternative: {source_json}")
        else:
            print("❌ Error: No source JSON found for Modern ThrRS sequence.")
            return

    with open(source_json, 'r') as f:
        data = json.load(f)
        seq = data['sequences'][0]['protein']['sequence']

    # The Arena: Protein + Thr + Ser + Zinc
    comp_data = {
        "name": f"COMPETITION_{protein_name}_{lig1}_vs_{lig2}",
        "dialect": "alphafold3",
        "version": 1,
        "modelSeeds": [1],
        "sequences": [
            {"protein": {"id": ["A"], "sequence": seq}},       # Chain A
            {"ligand": {"id": ["B"], "ccdCodes": [lig1]}},     # Chain B (Thr)
            {"ligand": {"id": ["C"], "ccdCodes": [lig2]}},     # Chain C (Ser)
            {"ligand": {"id": ["D"], "ccdCodes": ["ZN"]}}      # Chain D (Zinc)
        ]
    }

    filename = f"COMPETITION_{protein_name}_{lig1}_vs_{lig2}.json"
    with open(filename, 'w') as f:
        json.dump(comp_data, f, indent=2)
    print(f"✅ Created Zinc Trap Arena: {filename}")

# Run the generator
create_competition(
    source_json=source_file,
    protein_name="modern_thrrs",
    lig1="THR",
    lig2="SER"
)
