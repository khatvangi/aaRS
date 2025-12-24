#!/usr/bin/env python3
"""
Complete AF3 JSON generation pipeline:
1. Read source FASTAs
2. Remove all gaps and invalid characters
3. Convert DNA→RNA notation
4. Validate sequences
5. Generate AF3 JSONs
"""

import json
from pathlib import Path
import re

def clean_fasta(fasta_path, is_rna=False):
    """Read FASTA and thoroughly clean for AF3."""
    print(f"  Reading: {fasta_path}")
    
    with open(fasta_path) as f:
        lines = f.readlines()
    
    # Extract sequence (skip header lines)
    seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
    
    print(f"    Original length: {len(seq)}")
    
    # Count gaps before removal
    gap_count = seq.count('-')
    if gap_count > 0:
        print(f"    Found {gap_count} gaps - removing...")
    
    # Remove ALL non-standard characters
    if is_rna:
        # For RNA: keep only ACGTU, then convert T→U
        seq = re.sub(r'[^ACGTU]', '', seq.upper())
        seq = seq.replace('T', 'U')
        valid_chars = set('ACGU')
    else:
        # For protein: keep only 20 standard amino acids
        seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', seq.upper())
        valid_chars = set('ACDEFGHIKLMNPQRSTVWY')
    
    print(f"    Clean length: {len(seq)}")
    
    # Final validation
    invalid = set(seq) - valid_chars
    if invalid:
        raise ValueError(f"Invalid characters remaining: {invalid}")
    
    print(f"    ✓ Validated")
    return seq

# Load and clean all sequences
print("="*60)
print("CLEANING SOURCE SEQUENCES")
print("="*60)

print("\n1. Shallow ancestor protein:")
shallow_protein = clean_fasta('../phase1/results/Anc-ProThrRS.fasta', is_rna=False)

print("\n2. Deep (LUCA) ancestor protein:")
deep_protein = clean_fasta('../phase1b/results/Anc-ProThrRS-LUCA.fasta', is_rna=False)

print("\n3. Ancestral tRNA:")
ancestral_rna = clean_fasta('../phase1/results/Anc-tRNA-ProThr.fasta', is_rna=True)

print("\n" + "="*60)
print("FINAL SEQUENCE SUMMARY")
print("="*60)
print(f"Shallow protein:  {len(shallow_protein)} aa")
print(f"Deep protein:     {len(deep_protein)} aa")
print(f"Ancestral tRNA:   {len(ancestral_rna)} nt")
print(f"Pairwise identity: {sum(a==b for a,b in zip(shallow_protein, deep_protein))/min(len(shallow_protein), len(deep_protein))*100:.1f}%")

# Define ligands
ligands = {
    'PRO': {'smiles': 'C1CC(NC1)C(=O)O', 'ccd': 'PRO'},
    'THR': {'smiles': 'CC(C(C(=O)O)N)O', 'ccd': 'THR'}
}

# Define all test cases
test_cases = [
    {'name': 'shallow_ancestral_pro', 'protein': shallow_protein, 'ligand': 'PRO'},
    {'name': 'shallow_ancestral_thr', 'protein': shallow_protein, 'ligand': 'THR'},
    {'name': 'deep_ancestral_pro', 'protein': deep_protein, 'ligand': 'PRO'},
    {'name': 'deep_ancestral_thr', 'protein': deep_protein, 'ligand': 'THR'}
]

# Generate JSONs
output_dir = Path('inputs/af3_jsons')
output_dir.mkdir(parents=True, exist_ok=True)

# Clean directory first
for old_json in output_dir.glob('*.json'):
    old_json.unlink()
    print(f"  Removed old: {old_json.name}")

print("\n" + "="*60)
print("GENERATING AF3 INPUT JSONS")
print("="*60)

for case in test_cases:
    ligand_info = ligands[case['ligand']]
    
    af3_input = {
        "name": case['name'],
        "dialect": "alphafold3",
        "version": 1,
        "sequences": [
            {"protein": {"id": ["A"], "sequence": case['protein']}},
            {"rna": {"id": ["B"], "sequence": ancestral_rna}},
            {"ligand": {"id": ["C"], "smiles": ligand_info['smiles'], "ccdCodes": [ligand_info['ccd']]}}
        ],
        "modelSeeds": [1]
    }
    
    output_path = output_dir / f"{case['name']}.json"
    with open(output_path, 'w') as f:
        json.dump(af3_input, f, indent=2)
    
    print(f"\n✓ {output_path.name}")
    print(f"  Protein: {len(case['protein'])} aa")
    print(f"  Ligand:  {case['ligand']}")

print("\n" + "="*60)
print("✓ ALL JSONS GENERATED - READY TO LAUNCH!")
print("="*60)
print("\nLaunch command:")
print("  nohup bash run_dual_test.sh > logs/phase2_clean.log 2>&1 &")
