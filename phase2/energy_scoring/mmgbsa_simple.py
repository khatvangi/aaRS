#!/usr/bin/env python3
"""
Simplified MM/GBSA using direct parameter lookup + GB solvation.
Avoids OpenMM force field template issues.

BE = MM + GB + SA
where MM = bonded + VdW + Coulomb
      GB = generalized Born solvation
      SA = surface area term
"""
import sys
import json
import gemmi
import numpy as np
import math

# Amber14 parameters (heavy atoms only)
ATOM_PARAMS = {
    # Backbone
    'N': (-0.4157, 0.325, 0.71128),   # (charge_e, sigma_nm, epsilon_kJ_mol)
    'CA': (0.0337, 0.339967, 0.4577296),
    'C': (0.5973, 0.339967, 0.359824),
    'O': (-0.5679, 0.295992, 0.87864),
    'OXT': (-0.8055, 0.295992, 0.87864),  # Terminal oxygen
    # Side chains - representative subset
    'CB': (-0.0875, 0.339967, 0.4577296),
    'OG': (-0.6546, 0.306647, 0.88028),   # Ser hydroxyl
    'OG1': (-0.6546, 0.306647, 0.88028),  # Thr hydroxyl
    'CG': (0.0, 0.339967, 0.4577296),
    'CG1': (-0.0122, 0.339967, 0.4577296),
    'CG2': (-0.0597, 0.339967, 0.4577296),
    'CD': (0.0, 0.339967, 0.4577296),
    'CD1': (-0.0875, 0.339967, 0.4577296),
    'CD2': (-0.0875, 0.339967, 0.4577296),
    'CE': (0.0, 0.339967, 0.4577296),
    'NZ': (-0.3854, 0.325, 0.71128),  # Lys amine
    'OD1': (-0.6376, 0.295992, 0.87864),  # Asp carboxyl
    'OD2': (-0.6376, 0.295992, 0.87864),
    'OE1': (-0.6376, 0.295992, 0.87864),  # Glu carboxyl
    'OE2': (-0.6376, 0.295992, 0.87864),
    'ND1': (-0.3854, 0.325, 0.71128),
    'ND2': (-0.3854, 0.325, 0.71128),
    'NE': (-0.3854, 0.325, 0.71128),
    'NE1': (-0.3854, 0.325, 0.71128),
    'NE2': (-0.3854, 0.325, 0.71128),
    'NH1': (-0.3854, 0.325, 0.71128),
    'NH2': (-0.3854, 0.325, 0.71128),
    'SG': (-0.3119, 0.356359, 1.04600),  # Cys sulfur
    'SD': (-0.3119, 0.356359, 1.04600),  # Met sulfur
    # Zn
    'ZN': (2.0, 0.125, 0.0125),  # Zn2+
}

# Constants
KCOUL = 332.06  # Coulomb constant for kcal/mol
CUTOFF_NM = 0.8  # 8 Å cutoff

def read_cif_atoms(cif_path):
    """Read CIF and get all atom positions/types"""
    try:
        doc = gemmi.cif.read(cif_path)
        block = doc.sole_block()
        st = gemmi.make_structure_from_block(block)
        model = st[0]

        # Collect atoms by chain
        chains_data = {}
        for ch in model:
            atoms = []
            for res in ch:
                for atom in res:
                    atoms.append({
                        'name': atom.name.strip(),
                        'element': atom.element.name,
                        'pos': atom.pos,
                        'res_name': res.name.strip(),
                        'res_num': res.seqid.num
                    })
            chains_data[ch.name] = {
                'n_residues': len(list(ch)),
                'atoms': atoms
            }

        # Identify chains
        protein_chain = max(chains_data.keys(), key=lambda k: chains_data[k]['n_residues'])

        ligand_chain = None
        for ch_name, data in chains_data.items():
            if ch_name == 'B' and 1 <= data['n_residues'] <= 3:
                ligand_chain = ch_name
                break

        if not ligand_chain:
            non_protein = [k for k in chains_data.keys() if k != protein_chain and chains_data[k]['n_residues'] >= 1]
            if non_protein:
                ligand_chain = min(non_protein, key=lambda k: chains_data[k]['n_residues'])

        # Get ligand resname
        lig_resname = "UNK"
        if ligand_chain and chains_data[ligand_chain]['atoms']:
            lig_resname = chains_data[ligand_chain]['atoms'][0]['res_name']

        # Check for Zn
        zn_present = any(a['element'] == 'Zn' for atoms in chains_data.values() for a in atoms['atoms'])

        return {
            'chains': chains_data,
            'protein_chain': protein_chain,
            'ligand_chain': ligand_chain,
            'lig_resname': lig_resname,
            'zn_present': zn_present
        }, "OK"

    except Exception as e:
        return None, f"ERROR: {str(e)[:100]}"

def get_atom_params(atom_name, element):
    """Get Amber14 parameters for atom"""
    if atom_name in ATOM_PARAMS:
        return ATOM_PARAMS[atom_name]
    elif element == 'Zn':
        return ATOM_PARAMS['ZN']
    else:
        # Generic carbon
        return (0.0, 0.339967, 0.4577296)

def compute_born_radii(atoms):
    """
    Compute generalized Born radii using simple approximation.
    Returns Born radius for each atom (Å).
    """
    born_radii = []
    for atom in atoms:
        # Intrinsic Born radius based on element
        if atom['element'] == 'C':
            r0 = 1.7
        elif atom['element'] == 'N':
            r0 = 1.55
        elif atom['element'] == 'O':
            r0 = 1.5
        elif atom['element'] == 'S':
            r0 = 1.8
        elif atom['element'] == 'Zn':
            r0 = 1.2
        else:
            r0 = 1.5

        # Simplified: use intrinsic radius (proper GB would compute effective radius)
        born_radii.append(r0)

    return born_radii

def compute_gb_energy(atoms, born_radii):
    """
    Compute generalized Born solvation energy (OBC2-like).
    ΔG_solv = -0.5 * (1/ε_in - 1/ε_solv) * Σ q_i * q_j / f_GB
    where f_GB = sqrt(r_ij^2 + R_i*R_j*exp(-r_ij^2/(4*R_i*R_j)))
    """
    epsilon_in = 1.0
    epsilon_solv = 78.5
    prefactor = -166.03 * (1.0/epsilon_in - 1.0/epsilon_solv)  # 332.06 * 0.5

    E_gb = 0.0

    for i, atom_i in enumerate(atoms):
        q_i, _, _ = get_atom_params(atom_i['name'], atom_i['element'])
        R_i = born_radii[i] / 10.0  # Å to nm

        for j, atom_j in enumerate(atoms):
            if j <= i:
                continue

            q_j, _, _ = get_atom_params(atom_j['name'], atom_j['element'])
            R_j = born_radii[j] / 10.0  # Å to nm

            # Distance in nm
            dx = atom_i['pos'].x - atom_j['pos'].x
            dy = atom_i['pos'].y - atom_j['pos'].y
            dz = atom_i['pos'].z - atom_j['pos'].z
            r_ij = math.sqrt(dx*dx + dy*dy + dz*dz) / 10.0  # Å to nm

            # GB function
            exp_term = math.exp(-r_ij*r_ij / (4.0*R_i*R_j)) if R_i*R_j > 0 else 0
            f_gb = math.sqrt(r_ij*r_ij + R_i*R_j*exp_term)

            E_gb += prefactor * q_i * q_j / f_gb

    return E_gb

def compute_sa_energy(atoms):
    """
    Compute surface area term: γ * SASA
    Simplified: γ = 0.0072 kcal/mol/Å^2, SASA ≈ 4πr^2 for exposed atoms
    """
    gamma = 0.0072  # kcal/mol/Å^2
    total_sa = 0.0

    for atom in atoms:
        # Crude SASA estimate: count neighbors within 5 Å
        neighbors = 0
        for other in atoms:
            if other == atom:
                continue
            dx = atom['pos'].x - other['pos'].x
            dy = atom['pos'].y - other['pos'].y
            dz = atom['pos'].z - other['pos'].z
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            if dist < 5.0:  # Å
                neighbors += 1

        # Burial factor (0 = fully buried, 1 = fully exposed)
        burial = max(0, min(1, 1.0 - neighbors/12.0))

        # Atomic radius
        if atom['element'] == 'C':
            r = 1.7
        elif atom['element'] == 'N':
            r = 1.55
        elif atom['element'] == 'O':
            r = 1.5
        elif atom['element'] == 'S':
            r = 1.8
        else:
            r = 1.5

        atom_sa = burial * 4 * math.pi * r * r
        total_sa += atom_sa

    return gamma * total_sa

def compute_mm_energy(atoms1, atoms2=None):
    """
    Compute MM energy (VdW + Coulomb).
    If atoms2 is None, compute internal energy of atoms1.
    Otherwise, compute interaction energy between atoms1 and atoms2.
    """
    E_vdw = 0.0
    E_coul = 0.0

    if atoms2 is None:
        # Internal energy
        for i, a1 in enumerate(atoms1):
            q1, s1, e1 = get_atom_params(a1['name'], a1['element'])
            for j, a2 in enumerate(atoms1):
                if j <= i:
                    continue

                q2, s2, e2 = get_atom_params(a2['name'], a2['element'])

                # Distance (Å -> nm)
                dx = a1['pos'].x - a2['pos'].x
                dy = a1['pos'].y - a2['pos'].y
                dz = a1['pos'].z - a2['pos'].z
                r_nm = math.sqrt(dx*dx + dy*dy + dz*dz) / 10.0

                if r_nm < CUTOFF_NM:
                    # LJ
                    sigma = 0.5 * (s1 + s2)
                    epsilon = math.sqrt(e1 * e2)
                    sr6 = (sigma / r_nm) ** 6
                    E_vdw += 4.0 * epsilon * (sr6**2 - sr6)

                    # Coulomb
                    E_coul += KCOUL * q1 * q2 / r_nm
    else:
        # Interaction energy
        for a1 in atoms1:
            q1, s1, e1 = get_atom_params(a1['name'], a1['element'])
            for a2 in atoms2:
                q2, s2, e2 = get_atom_params(a2['name'], a2['element'])

                dx = a1['pos'].x - a2['pos'].x
                dy = a1['pos'].y - a2['pos'].y
                dz = a1['pos'].z - a2['pos'].z
                r_nm = math.sqrt(dx*dx + dy*dy + dz*dz) / 10.0

                if r_nm < CUTOFF_NM:
                    sigma = 0.5 * (s1 + s2)
                    epsilon = math.sqrt(e1 * e2)
                    sr6 = (sigma / r_nm) ** 6
                    E_vdw += 4.0 * epsilon * (sr6**2 - sr6)
                    E_coul += KCOUL * q1 * q2 / r_nm

    # Convert VdW from kJ/mol to kcal/mol
    E_vdw_kcal = E_vdw / 4.184

    return E_vdw_kcal, E_coul

def compute_mmgbsa(data):
    """
    Compute MM/GBSA binding energy.
    BE = [MM(complex) + GB(complex) + SA(complex)] - [MM(protein) + GB(protein) + SA(protein)] - [MM(ligand) + GB(ligand) + SA(ligand)]
    """
    protein_atoms = data['chains'][data['protein_chain']]['atoms']
    ligand_atoms = data['chains'][data['ligand_chain']]['atoms']
    complex_atoms = protein_atoms + ligand_atoms

    # Complex
    born_radii_complex = compute_born_radii(complex_atoms)
    E_mm_vdw_complex, E_mm_coul_complex = compute_mm_energy(complex_atoms)
    E_gb_complex = compute_gb_energy(complex_atoms, born_radii_complex)
    E_sa_complex = compute_sa_energy(complex_atoms)
    E_complex = E_mm_vdw_complex + E_mm_coul_complex + E_gb_complex + E_sa_complex

    # Protein
    born_radii_protein = compute_born_radii(protein_atoms)
    E_mm_vdw_protein, E_mm_coul_protein = compute_mm_energy(protein_atoms)
    E_gb_protein = compute_gb_energy(protein_atoms, born_radii_protein)
    E_sa_protein = compute_sa_energy(protein_atoms)
    E_protein = E_mm_vdw_protein + E_mm_coul_protein + E_gb_protein + E_sa_protein

    # Ligand
    born_radii_ligand = compute_born_radii(ligand_atoms)
    E_mm_vdw_ligand, E_mm_coul_ligand = compute_mm_energy(ligand_atoms)
    E_gb_ligand = compute_gb_energy(ligand_atoms, born_radii_ligand)
    E_sa_ligand = compute_sa_energy(ligand_atoms)
    E_ligand = E_mm_vdw_ligand + E_mm_coul_ligand + E_gb_ligand + E_sa_ligand

    # Binding energy
    BE = E_complex - E_protein - E_ligand

    return {
        'BE_mmgbsa': BE,
        'E_complex': E_complex,
        'E_protein': E_protein,
        'E_ligand': E_ligand,
        'E_mm_vdw': E_mm_vdw_complex - E_mm_vdw_protein - E_mm_vdw_ligand,
        'E_mm_coul': E_mm_coul_complex - E_mm_coul_protein - E_mm_coul_ligand,
        'E_gb': E_gb_complex - E_gb_protein - E_gb_ligand,
        'E_sa': E_sa_complex - E_sa_protein - E_sa_ligand
    }

def main(cif_path):
    """Main workflow"""
    data, status = read_cif_atoms(cif_path)

    if status != "OK":
        return {"file": cif_path, "status": status}

    if not data['ligand_chain']:
        return {"file": cif_path, "status": "NO_LIGAND_CHAIN"}

    try:
        result = compute_mmgbsa(data)

        return {
            "file": cif_path,
            "status": "OK",
            "lig_resname": data['lig_resname'],
            "protein_chain": data['protein_chain'],
            "ligand_chain": data['ligand_chain'],
            "Zn_present": data['zn_present'],
            **result
        }

    except Exception as e:
        return {
            "file": cif_path,
            "status": "MMGBSA_ERROR",
            "error": str(e)[:200]
        }

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python mmgbsa_simple.py <cif_file>")
        sys.exit(1)

    result = main(sys.argv[1])
    print(json.dumps(result))
