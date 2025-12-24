#!/usr/bin/env python3
"""
Test script to inspect CIF files and identify ligand residues
"""

import pymol
from pymol import cmd
import os

# Structure file paths
BASE_DIR = "/storage/kiran-stuff/aaRS/phase2/outputs"

STRUCTURES = {
    'luca_cat_pro': f"{BASE_DIR}/deep_domain_pro/deep_domain_pro_model.cif",
    'luca_cat_thr': f"{BASE_DIR}/deep_domain_thr/deep_domain_thr_model.cif",
}

def inspect_structure(name, path):
    """Load structure and print detailed information"""
    print(f"\n{'='*60}")
    print(f"INSPECTING: {name}")
    print(f"File: {path}")
    print('='*60)

    if not os.path.exists(path):
        print(f"✗ ERROR: File not found!")
        return

    cmd.reinitialize()
    cmd.load(path, name)

    # Get basic info
    print(f"\nTotal atoms: {cmd.count_atoms(name)}")

    # List all chains
    chains = []
    cmd.iterate(name, "chains.append(chain)", space={'chains': chains})
    unique_chains = sorted(set(chains))
    print(f"Chains: {unique_chains}")

    # List all residue types
    resnames = []
    cmd.iterate(name, "resnames.append(resn)", space={'resnames': resnames})
    unique_resnames = sorted(set(resnames))
    print(f"\nResidue types found ({len(unique_resnames)}): {unique_resnames}")

    # Check for small molecules (potential ligands)
    print("\nSearching for ligands...")

    # Method 1: Look for PRO/THR as ligands (not in polymer)
    pro_count = cmd.count_atoms(f"{name} and resn PRO")
    thr_count = cmd.count_atoms(f"{name} and resn THR")
    print(f"  PRO residues: {pro_count}")
    print(f"  THR residues: {thr_count}")

    # Method 2: Look for hetatm
    hetatm_count = cmd.count_atoms(f"{name} and hetatm")
    print(f"  HETATM atoms: {hetatm_count}")

    if hetatm_count > 0:
        hetatm_res = []
        cmd.iterate(f"{name} and hetatm", "hetatm_res.append((resi, resn))", space={'hetatm_res': hetatm_res})
        unique_hetatm = sorted(set(hetatm_res))
        print(f"  HETATM residues: {unique_hetatm}")

    # Method 3: Look for non-standard amino acids
    standard_aa = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                   'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

    non_standard = [r for r in unique_resnames if r not in standard_aa and r not in ['HOH', 'WAT']]
    if non_standard:
        print(f"  Non-standard residues: {non_standard}")

    # Check last few residues in each chain (ligands often at the end)
    for chain in unique_chains:
        chain_res = []
        cmd.iterate(f"{name} and chain {chain}", "chain_res.append((resi, resn))", space={'chain_res': chain_res})
        if chain_res:
            print(f"\n  Chain {chain} - Last 5 residues: {chain_res[-5:]}")

def main():
    print("\n" + "#"*60)
    print("# LIGAND DETECTION TEST")
    print("#"*60)

    for name, path in STRUCTURES.items():
        inspect_structure(name, path)

    print("\n" + "#"*60)
    print("# INSPECTION COMPLETE")
    print("#"*60)

if __name__ == '__main__':
    main()
