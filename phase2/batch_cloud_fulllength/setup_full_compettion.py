import json
import os

# 1. READ THE SOURCE SEQUENCE
source_file = "anc_prors_full_PRO.json"

if not os.path.exists(source_file):
    print(f"❌ Error: Could not find {source_file}")
    exit()

with open(source_file, 'r') as f:
    data = json.load(f)
    # Extract the full length sequence from Chain A
    full_sequence = data['sequences'][0]['protein']['sequence']
    print(f"✅ Loaded Sequence ({len(full_sequence)} residues)")

# 2. DEFINITION FUNCTION
def create_arena(filename, sequence, lig1, lig2, zinc=False):
    # Base structure
    comp_data = {
        "name": filename.replace(".json", ""),
        "dialect": "alphafold3",
        "version": 1,
        "modelSeeds": [1],
        "sequences": [
            # Chain A: The Massive Ancestral Enzyme
            {"protein": {"id": ["A"], "sequence": sequence}},
            # Chain B: Combatant 1
            {"ligand": {"id": ["B"], "ccdCodes": [lig1]}},
            # Chain C: Combatant 2
            {"ligand": {"id": ["C"], "ccdCodes": [lig2]}}
        ]
    }
    
    # Add Chain D: Zinc (Only if requested, ProRS usually doesn't need it, ThrRS does)
    if zinc:
        comp_data["sequences"].append({"ligand": {"id": ["D"], "ccdCodes": ["ZN"]}})

    with open(filename, 'w') as f:
        json.dump(comp_data, f, indent=2)
    print(f"✅ Created Arena: {filename}")

# 3. GENERATE THE GLADIATOR MATCHES

# Match A: The "Bucket Check" (Pro vs Glu)
# Does the active site reject Glu? Or does it bind in the Editing Domain?
create_arena(
    filename="COMPETITION_FULL_anc_prors_PRO_vs_GLU.json",
    sequence=full_sequence,
    lig1="PRO",
    lig2="GLU",
    zinc=False # Class IIa ProRS generally Mg2+ dependent, Zinc usually not structural here
)

# Match B: The "LUCA Check" (Pro vs Thr)
# Is this actually a generalist ProThrRS?
create_arena(
    filename="COMPETITION_FULL_anc_prors_PRO_vs_THR.json",
    sequence=full_sequence,
    lig1="PRO",
    lig2="THR",
    zinc=False 
)
