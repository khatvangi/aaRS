import json
import os

# ---------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------
output_dir = "af3_inputs_modern_allAAs"
os.makedirs(output_dir, exist_ok=True)

# Confirmed Sequences from previous steps
seq_modern_pro = "MRTSQYLLSTLKETPADAEVISHQLMLRAGMIRKLASGLYTWLPTGVRVLKKVENIVREEMNNAGAIEVSMPVVQPADLWQESGRWEQYGPELLRFVDRGERPFVLGPTHEEVITDLIRNELSSYKQLPLNFYQIQTKFRDEVRPRFGVMRSREFLMKDAYSFHTSQESLQETYDAMYAAYSKIFSRMGLDFRAVQADTGSIGGSASHEFQVLAQSGEDDVVFSDTSDYAANIELAEAIAPKEPRAAATQEMTLVDTPNAKTIAELVEQFNLPIEKTVKTLLVKAVEGSSFPQVALLVRGDHELNEVKAEKLPQVASPLTFATEEEIRAVVKAGPGSLGPVNMPIPVVIDRTVAAMSDFAAGANIDGKHYFGINWDRDVATPEVADIRNVVAGDPSPDGQGRLLIKRGIEVGHIFQLGTKYSEALKASVQGEDGRNQILTMGCYGIGVTRVVAAAIEQNYDERGIVWPDAIAPFQVAILPMNMHKSFRVQELAEKLYSELRAQGIEVLLDDRKERPGVMFADMELIGIPHTIVLGDRNLDNDDIEYKYRRNGEKQLIKTGDIVEYLVKQIKG"
seq_modern_thr = "MPVITLPDGSQRHYDHAVSPMDVALDIGPGLAKACIIGRNQYLVDDSDRGFHICRAIDGVTRNDTEIMPFINDPKSEAHIAKCVAAWPEVADRIGANAAALPAVLELPEDVVEDIKRRGVLPRRMPVIKREHSEEIGRDLASEIGLVAQKYLGEGYASENLELHIVRGPMEIRLRQNLEEALKIAQAKGDEVESLKSRVAAELQALKEALVEQEVKEQVELKAELDAKLAQLEAKREQLRDEALAKARAELEKLLYGSKKPGEPDELVAEIADRREKAAALEAALAQSEEQFSLDAQAQKLPLDKLWAQANRAPSDFADQVIRRKAERGLYPPLMVTARGDLAASITPLRWRLPELFLKQGMPIGAYKGRKVFVKEFMQNPAYGRGIMIGIHQFNQDVVHRIVPGEGYFLFGPTPYQIMPLQVALGALQGLNPLNMATVFVGVIKPTNYNPHNGSEIPIKVSASEPEGDVSIVGFLVLPQEIISNIMNKEQELAHSIEVLKAGDNHHGEIAANIGGFKDIFGDFFRFTEQHAQLREIGIQPFYEDRGYIDLPDLVKEGKALIVEKGQIEQASGTSHNFVLELQNASENLVGSVGARVDRIGAGLPGGIIEDTLPLLKNAQATVLDGKKVTVKPGNTVVAKRDADGYARLQVAGTVDGRLQVKP"

# The 20 Standard Amino Acids
amino_acids = [
    "ALA", "ARG", "ASN", "ASP", "CYS", 
    "GLN", "GLU", "GLY", "HIS", "ILE", 
    "LEU", "LYS", "MET", "PHE", "PRO", 
    "SER", "THR", "TRP", "TYR", "VAL"
]

enzymes = [
    ("modern_ecoli_full_prors", seq_modern_pro),
    ("modern_ecoli_full_thrrs", seq_modern_thr)
]

# ---------------------------------------------------------
# GENERATOR
# ---------------------------------------------------------
print(f"Generating {len(enzymes) * len(amino_acids)} input files in {output_dir}...")

for enz_name, enz_seq in enzymes:
    for aa in amino_acids:
        run_name = f"{enz_name}_{aa}"
        
        data = {
            "name": run_name,
            "dialect": "alphafold3",
            "version": 1,
            "modelSeeds": [1],
            "sequences": [
                {
                    "protein": {
                        "id": ["A"],
                        "sequence": enz_seq
                    }
                },
                {
                    "ligand": {
                        "id": ["B"],
                        "ccdCodes": [aa]
                    }
                }
            ]
        }
        
        filename = os.path.join(output_dir, f"{run_name}.json")
        with open(filename, 'w') as f:
            json.dump(data, f, indent=2)

print("Done! Modern inputs ready.")
