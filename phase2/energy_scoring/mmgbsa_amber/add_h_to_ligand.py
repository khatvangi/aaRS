#!/usr/bin/env python3
"""
Add hydrogens to ligand PDB using mol2 template for atom names and types.
Matches heavy atoms by name, then adds hydrogens with correct names from mol2.
"""
import sys
import re

def parse_mol2_atoms(mol2_file):
    """Parse mol2 file and return dict of {atom_name: (x, y, z, type)}"""
    atoms = {}
    in_atoms = False
    with open(mol2_file) as f:
        for line in f:
            if "@<TRIPOS>ATOM" in line:
                in_atoms = True
                continue
            if "@<TRIPOS>BOND" in line:
                break
            if in_atoms and line.strip():
                parts = line.split()
                if len(parts) >= 6:
                    atom_name = parts[1]
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    atom_type = parts[5]
                    atoms[atom_name] = (x, y, z, atom_type)
    return atoms

def parse_pdb_atoms(pdb_file):
    """Parse PDB file and return list of (name, x, y, z) for HETATM/ATOM"""
    atoms = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append((name, x, y, z))
    return atoms

def write_pdb_with_h(output_file, mol2_atoms, pdb_heavy_atoms, resname="LIG", chain="B"):
    """Write PDB with heavy atoms from PDB and H from mol2"""
    # Create mapping of heavy atom names from PDB to mol2 coordinates
    # We'll use PDB coordinates for heavy atoms, mol2 coordinates for H

    with open(output_file, 'w') as f:
        serial = 1

        # First, write heavy atoms from PDB
        for pdb_name, x, y, z in pdb_heavy_atoms:
            if pdb_name in mol2_atoms:
                atom_type = mol2_atoms[pdb_name][3]
                # Guess element from atom name
                element = pdb_name[0] if pdb_name[0].isalpha() else "C"

                f.write(f"HETATM{serial:5d}  {pdb_name:4s} {resname:3s} {chain:1s}   1    ")
                f.write(f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 50.00          {element:>2s}\n")
                serial += 1

        # Then, write H atoms from mol2
        for atom_name, (x, y, z, atom_type) in mol2_atoms.items():
            if atom_name.startswith('H') or atom_name.startswith('h'):
                f.write(f"HETATM{serial:5d}  {atom_name:4s} {resname:3s} {chain:1s}   1    ")
                f.write(f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 50.00           H\n")
                serial += 1

        f.write("END\n")

def main():
    if len(sys.argv) < 4:
        print("Usage: python add_h_to_ligand.py ligand.pdb template.mol2 output.pdb")
        sys.exit(1)

    pdb_file = sys.argv[1]
    mol2_file = sys.argv[2]
    output_file = sys.argv[3]

    mol2_atoms = parse_mol2_atoms(mol2_file)
    pdb_atoms = parse_pdb_atoms(pdb_file)

    print(f"Mol2 has {len(mol2_atoms)} atoms")
    print(f"PDB has {len(pdb_atoms)} atoms")

    write_pdb_with_h(output_file, mol2_atoms, pdb_atoms)
    print(f"Wrote {output_file}")

if __name__ == '__main__':
    main()
