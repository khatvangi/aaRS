#!/usr/bin/env python3
"""
Rigid-fit a template mol2 (with H) onto an AF3 ligand PDB (heavy atoms only).
Uses Kabsch algorithm on heavy-atom coordinates, applies same transform to all atoms.

Usage: python fit_ligand_mol2_to_af3.py template.mol2 af3_ligand.pdb output_pose.mol2
"""
import sys
import numpy as np

def parse_mol2_atoms(mol2_file):
    """
    Parse mol2 file and return:
    - atoms: list of dicts with {name, x, y, z, type, charge, ...}
    - bonds: list of bond lines
    - header: molecule info lines
    """
    atoms = []
    bonds = []
    header_lines = []

    in_atoms = False
    in_bonds = False
    in_molecule = False

    with open(mol2_file) as f:
        for line in f:
            if "@<TRIPOS>MOLECULE" in line:
                in_molecule = True
                header_lines.append(line)
                continue

            if "@<TRIPOS>ATOM" in line:
                in_molecule = False
                in_atoms = True
                in_bonds = False
                continue

            if "@<TRIPOS>BOND" in line:
                in_atoms = False
                in_bonds = True
                continue

            if "@<TRIPOS>SUBSTRUCTURE" in line:
                in_bonds = False
                continue

            if in_molecule:
                header_lines.append(line)

            if in_atoms and line.strip():
                parts = line.split()
                if len(parts) >= 9:
                    atoms.append({
                        'serial': int(parts[0]),
                        'name': parts[1],
                        'x': float(parts[2]),
                        'y': float(parts[3]),
                        'z': float(parts[4]),
                        'type': parts[5],
                        'resid': parts[6],
                        'resname': parts[7],
                        'charge': float(parts[8]),
                        'original_line': line
                    })

            if in_bonds and line.strip():
                bonds.append(line)

    return atoms, bonds, header_lines

def parse_pdb_atoms(pdb_file):
    """
    Parse PDB file and return dict of {atom_name: (x, y, z)}
    """
    atoms = {}
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms[name] = (x, y, z)
    return atoms

def kabsch_alignment(P, Q):
    """
    Kabsch algorithm: find rotation matrix R and translation t such that
    Q_transformed = R @ P + t minimizes RMSD.

    P: (N, 3) array - source coordinates (mol2 heavy atoms)
    Q: (N, 3) array - target coordinates (AF3 PDB heavy atoms)

    Returns: (R, t) where R is 3x3 rotation matrix, t is translation vector
    """
    # Center both sets
    P_center = P.mean(axis=0)
    Q_center = Q.mean(axis=0)

    P_centered = P - P_center
    Q_centered = Q - Q_center

    # Compute covariance matrix
    H = P_centered.T @ Q_centered

    # SVD
    U, S, Vt = np.linalg.svd(H)

    # Rotation matrix
    R = Vt.T @ U.T

    # Ensure proper rotation (det(R) = 1)
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    # Translation
    t = Q_center - R @ P_center

    return R, t

def write_mol2(output_file, atoms, bonds, header_lines):
    """
    Write mol2 file with updated coordinates
    """
    with open(output_file, 'w') as f:
        # Header
        for line in header_lines:
            f.write(line)

        # Atoms
        f.write("@<TRIPOS>ATOM\n")
        for atom in atoms:
            f.write(f"{atom['serial']:7d} {atom['name']:<8s} "
                   f"{atom['x']:9.4f} {atom['y']:9.4f} {atom['z']:9.4f} "
                   f"{atom['type']:<6s} {atom['resid']:>3s} {atom['resname']:<8s} "
                   f"{atom['charge']:9.6f}\n")

        # Bonds
        f.write("@<TRIPOS>BOND\n")
        for bond_line in bonds:
            f.write(bond_line)

        # Substructure (minimal)
        f.write("@<TRIPOS>SUBSTRUCTURE\n")
        f.write("     1 LIG         1 TEMP              0 ****  ****    0 ROOT\n")

def main():
    if len(sys.argv) != 4:
        print("Usage: python fit_ligand_mol2_to_af3.py template.mol2 af3_ligand.pdb output_pose.mol2")
        sys.exit(1)

    template_mol2 = sys.argv[1]
    af3_pdb = sys.argv[2]
    output_mol2 = sys.argv[3]

    # Parse inputs
    mol2_atoms, bonds, header = parse_mol2_atoms(template_mol2)
    pdb_atoms = parse_pdb_atoms(af3_pdb)

    print(f"Template mol2: {len(mol2_atoms)} atoms")
    print(f"AF3 PDB: {len(pdb_atoms)} atoms")

    # Match heavy atoms by name
    matched_pairs = []
    mol2_heavy_names = []

    for atom in mol2_atoms:
        name = atom['name']
        # Heavy atoms: not H or h
        if not (name.startswith('H') or name.startswith('h')):
            if name in pdb_atoms:
                matched_pairs.append((atom, pdb_atoms[name]))
                mol2_heavy_names.append(name)
            else:
                print(f"Warning: Heavy atom {name} in mol2 not found in PDB")

    print(f"Matched {len(matched_pairs)} heavy atoms: {', '.join(mol2_heavy_names)}")

    if len(matched_pairs) < 3:
        print("Error: Need at least 3 matched atoms for alignment")
        sys.exit(1)

    # Extract coordinates for Kabsch
    P = np.array([[atom['x'], atom['y'], atom['z']] for atom, _ in matched_pairs])
    Q = np.array([coords for _, coords in matched_pairs])

    # Compute alignment
    R, t = kabsch_alignment(P, Q)

    # Apply transformation to ALL atoms in mol2 (including H)
    for atom in mol2_atoms:
        coords = np.array([atom['x'], atom['y'], atom['z']])
        transformed = R @ coords + t
        atom['x'] = transformed[0]
        atom['y'] = transformed[1]
        atom['z'] = transformed[2]

    # Compute RMSD to verify
    P_transformed = (R @ P.T).T + t
    rmsd = np.sqrt(np.mean(np.sum((P_transformed - Q)**2, axis=1)))
    print(f"Heavy-atom RMSD after alignment: {rmsd:.4f} Å")

    # Write output
    write_mol2(output_mol2, mol2_atoms, bonds, header)
    print(f"Wrote {output_mol2}")

if __name__ == '__main__':
    main()
