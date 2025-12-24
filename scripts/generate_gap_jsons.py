import json
import urllib.request
import time

# --- Configuration ---

# Ligand Definitions (SMILES + CCD)
# Note: For standard amino acids as ligands, ccdCodes are preferred.
LIGANDS = {
    "PRO": {"smiles": "C1CC(NC1)C(=O)O", "ccd": "PRO"},
    "THR": {"smiles": "CC(C(C(=O)O)N)O", "ccd": "THR"},
    "TRP": {"smiles": "C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N", "ccd": "TRP"},
    "PHE": {"smiles": "C1=CC=C(C=C1)CC(C(=O)O)N", "ccd": "PHE"}
}

# Job Definitions
JOBS = [
    # --- E. coli ProRS (P16659) ---
    {
        "filename": "modern_ecoli_cat_trp.json",
        "uniprot": "P16659",
        "ligand_code": "TRP",
        "truncate_to": 470
    },
    {
        "filename": "modern_ecoli_cat_phe.json",
        "uniprot": "P16659",
        "ligand_code": "PHE",
        "truncate_to": 470
    },
    {
        "filename": "modern_ecoli_full_pro.json",
        "uniprot": "P16659",
        "ligand_code": "PRO",
        "truncate_to": None
    },
    {
        "filename": "modern_ecoli_full_thr.json",
        "uniprot": "P16659",
        "ligand_code": "THR",
        "truncate_to": None
    },
    
    # --- Human ThrRS (P26639) ---
    {
        "filename": "modern_human_full_thr.json",
        "uniprot": "P26639",
        "ligand_code": "THR",
        "truncate_to": None
    },
    {
        "filename": "modern_human_full_pro.json",
        "uniprot": "P26639",
        "ligand_code": "PRO",
        "truncate_to": None
    }
]

def get_sequence(uniprot_id):
    """Fetches sequence from UniProt."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        with urllib.request.urlopen(url) as response:
            fasta_data = response.read().decode('utf-8')
            lines = fasta_data.split('\n')
            # Join lines, skipping the header (first line starting with >)
            sequence = "".join(line.strip() for line in lines if line and not line.startswith('>'))
            print(f"Fetched {uniprot_id}: {len(sequence)} residues.")
            return sequence
    except Exception as e:
        print(f"Error fetching {uniprot_id}: {e}")
        return None

def create_json(job, full_sequence):
    """Generates the AlphaFold 3 input JSON with required dialect/version fields."""
    
    # Handle Truncation
    if job["truncate_to"]:
        sequence = full_sequence[:job["truncate_to"]]
        print(f"  > Truncated to {len(sequence)} residues")
    else:
        sequence = full_sequence
        print(f"  > Using full length ({len(sequence)} residues)")

    ligand_def = LIGANDS[job["ligand_code"]]
    
    # The JSON structure strictly required by run_alphafold.py
    content = {
        "name": job["filename"].replace(".json", ""),
        "dialect": "alphafold3",
        "version": 1,
        "modelSeeds": [1],
        "sequences": [
            {
                "protein": {
                    "id": ["A"],
                    "sequence": sequence
                }
            },
            {
                "ligand": {
                    "id": ["B"],
                    "smiles": ligand_def["smiles"],
                    "ccdCodes": [ligand_def["ccd"]]
                }
            }
        ]
    }
    
    with open(job["filename"], 'w') as f:
        json.dump(content, f, indent=2)
    print(f"  > Saved: {job['filename']}")

# --- Main Execution ---
cache = {}

print("Starting generation of fixed JSONs...")

for job in JOBS:
    uid = job["uniprot"]
    if uid not in cache:
        cache[uid] = get_sequence(uid)
        time.sleep(0.5) # Be nice to UniProt API
    
    if cache[uid]:
        print(f"\nProcessing: {job['filename']}")
        create_json(job, cache[uid])
    else:
        print(f"Skipping {job['filename']} due to download error.")

print("\nDone.")

