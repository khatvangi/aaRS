#!/usr/bin/env python3
"""
Generate AF3 JSON input files for E. coli ThrRS with Zn and all 20 amino acids.
Uses the correct E. coli sequence from PDB 1EVL.
"""

import json
import os

# E. coli ThrRS sequence from PDB 1EVL (with proper Zn-binding site)
ECOLI_THRRS_SEQUENCE = """RDHRKIGKQLDLYHMQEEAPGMVFWHNDGWTIFRELEVFVRSKLKEYQYQEVKGPFMMDRVLWEKTGHWDNYKDAMFTTSSENREYCIKPMNCPGHVQIFNQGLKSYRDLPLRMAEFGSCHRNEPSGSLHGLMRVRGFTQDDAHIFCTEEQIRDEVNGCIRLVYDMYSTFGFEKIVVKLSTRPEKRIGSDEMWDRAEADLAVALEENNIPFEYQLGEGAFYGPKIEFTLYDCLDRAWQCGTVQLDFSLPSRLSASYVGEDNERKVPVMIHRAILGSMERFIGILTEEFAGFFPTWLAPVQVVIMNITDSQSEYVNELTQKLSNAGIRVKADLRNEKIGFKIREHTLRRVPYMLVCGDKEVESGKVAVRTRRGKDLGSMDVNEVIEKLQQEIRSRSLKQLEE""".replace("\n", "")

# All 20 standard amino acids (3-letter codes for CCD)
AMINO_ACIDS = [
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL"
]

def generate_af3_json(amino_acid, output_dir):
    """Generate AF3 JSON for ThrRS + Zn + amino acid ligand."""
    
    job_name = f"modern_thrrs_ecoli_zn_{amino_acid}"
    
    af3_input = {
        "name": job_name,
        "dialect": "alphafold3",
        "version": 1,
        "modelSeeds": [2],
        "sequences": [
            {
                "protein": {
                    "id": ["A"],
                    "sequence": ECOLI_THRRS_SEQUENCE
                }
            },
            {
                "ligand": {
                    "id": ["B"],
                    "ccdCodes": [amino_acid]
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
    
    # Create output directory if needed
    os.makedirs(output_dir, exist_ok=True)
    
    # Write JSON file
    output_path = os.path.join(output_dir, f"{job_name}.json")
    with open(output_path, 'w') as f:
        json.dump(af3_input, f, indent=2)
    
    print(f"Generated: {output_path}")
    return output_path

def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Generate AF3 JSONs for E. coli ThrRS + Zn with all amino acids"
    )
    parser.add_argument(
        "-o", "--output-dir",
        default="./thrrs_ecoli_zn_jobs",
        help="Output directory for JSON files (default: ./thrrs_ecoli_zn_jobs)"
    )
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("Generating AF3 inputs for E. coli ThrRS + Zn")
    print(f"Sequence: PDB 1EVL ({len(ECOLI_THRRS_SEQUENCE)} residues)")
    print(f"Output directory: {args.output_dir}")
    print("=" * 60)
    
    generated_files = []
    for aa in AMINO_ACIDS:
        path = generate_af3_json(aa, args.output_dir)
        generated_files.append(path)
    
    print("=" * 60)
    print(f"Generated {len(generated_files)} JSON files")
    print("\nTo run with AF3:")
    print(f"  for f in {args.output_dir}/*.json; do")
    print(f"    alphafold3 --json_path $f --output_dir af3_output/")
    print(f"  done")
    print("=" * 60)

if __name__ == "__main__":
    main()
