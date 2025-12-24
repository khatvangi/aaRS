import json
import os

# ---------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------
base_dir = "af3_modern_matrix"
os.makedirs(base_dir, exist_ok=True)

# Sequences (Confirmed Modern E. coli)
seq_pro = "MRTSQYLLSTLKETPADAEVISHQLMLRAGMIRKLASGLYTWLPTGVRVLKKVENIVREEMNNAGAIEVSMPVVQPADLWQESGRWEQYGPELLRFVDRGERPFVLGPTHEEVITDLIRNELSSYKQLPLNFYQIQTKFRDEVRPRFGVMRSREFLMKDAYSFHTSQESLQETYDAMYAAYSKIFSRMGLDFRAVQADTGSIGGSASHEFQVLAQSGEDDVVFSDTSDYAANIELAEAIAPKEPRAAATQEMTLVDTPNAKTIAELVEQFNLPIEKTVKTLLVKAVEGSSFPQVALLVRGDHELNEVKAEKLPQVASPLTFATEEEIRAVVKAGPGSLGPVNMPIPVVIDRTVAAMSDFAAGANIDGKHYFGINWDRDVATPEVADIRNVVAGDPSPDGQGRLLIKRGIEVGHIFQLGTKYSEALKASVQGEDGRNQILTMGCYGIGVTRVVAAAIEQNYDERGIVWPDAIAPFQVAILPMNMHKSFRVQELAEKLYSELRAQGIEVLLDDRKERPGVMFADMELIGIPHTIVLGDRNLDNDDIEYKYRRNGEKQLIKTGDIVEYLVKQIKG"
seq_thr = "MPVITLPDGSQRHYDHAVSPMDVALDIGPGLAKACIIGRNQYLVDDSDRGFHICRAIDGVTRNDTEIMPFINDPKSEAHIAKCVAAWPEVADRIGANAAALPAVLELPEDVVEDIKRRGVLPRRMPVIKREHSEEIGRDLASEIGLVAQKYLGEGYASENLELHIVRGPMEIRLRQNLEEALKIAQAKGDEVESLKSRVAAELQALKEALVEQEVKEQVELKAELDAKLAQLEAKREQLRDEALAKARAELEKLLYGSKKPGEPDELVAEIADRREKAAALEAALAQSEEQFSLDAQAQKLPLDKLWAQANRAPSDFADQVIRRKAERGLYPPLMVTARGDLAASITPLRWRLPELFLKQGMPIGAYKGRKVFVKEFMQNPAYGRGIMIGIHQFNQDVVHRIVPGEGYFLFGPTPYQIMPLQVALGALQGLNPLNMATVFVGVIKPTNYNPHNGSEIPIKVSASEPEGDVSIVGFLVLPQEIISNIMNKEQELAHSIEVLKAGDNHHGEIAANIGGFKDIFGDFFRFTEQHAQLREIGIQPFYEDRGYIDLPDLVKEGKALIVEKGQIEQASGTSHNFVLELQNASENLVGSVGARVDRIGAGLPGGIIEDTLPLLKNAQATVLDGKKVTVKPGNTVVAKRDADGYARLQVAGTVDGRLQVKP"

# The 20 Standard Amino Acids
amino_acids = [
    "ALA", "ARG", "ASN", "ASP", "CYS", 
    "GLN", "GLU", "GLY", "HIS", "ILE", 
    "LEU", "LYS", "MET", "PHE", "PRO", 
    "SER", "THR", "TRP", "TYR", "VAL"
]

# Define the runs
enzymes = [
    ("modern_prors", seq_pro),
    ("modern_thrrs", seq_thr)
]

print(f"Setting up run matrix in {base_dir}...")

for enz_name, seq in enzymes:
    # 1. Define a 'Template' run (we use ALA as the default template)
    template_run_name = f"{enz_name}_ALA"
    
    for aa in amino_acids:
        run_name = f"{enz_name}_{aa}"
        
        # Create output subdirectory for this run
        # AF3 expects the JSON in the input_dir, but we are prepping the structure 
        # so we can inject MSAs into the output_dir
        run_output_folder = os.path.join(base_dir, run_name)
        os.makedirs(run_output_folder, exist_ok=True)
        
        # Create the JSON input
        json_data = {
            "name": run_name,
            "dialect": "alphafold3",
            "version": 1,
            "modelSeeds": [1],
            "sequences": [
                { "protein": { "id": ["A"], "sequence": seq } },
                { "ligand": { "id": ["B"], "ccdCodes": [aa] } }
            ]
        }
        
        # Save JSON in the base dir (standard AF3 input location)
        # OR inside the folder if your runner supports it. 
        # We will save it in the folder to keep things self-contained.
        with open(os.path.join(run_output_folder, f"{run_name}.json"), 'w') as f:
            json.dump(json_data, f, indent=2)

print("\nMatrix Setup Complete.")
print("Run the following two 'Template' jobs first to generate the MSAs:")
print(f"  1. modern_prors_ALA")
print(f"  2. modern_thrrs_ALA")
