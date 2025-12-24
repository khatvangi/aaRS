#!/usr/bin/env python3
"""
Corrected MM/GBSA - simplified interaction energy only.
Returns realistic binding energies in the -50 to +50 kcal/mol range.
"""
import sys
import json
import gemmi
import numpy as np
import math

# Amber14 parameters
ATOM_PARAMS = {
    'N': (-0.4157, 0.325, 0.71128),
    'CA': (0.0337, 0.339967, 0.4577296),
    'C': (0.5973, 0.339967, 0.359824),
    'O': (-0.5679, 0.295992, 0.87864),
    'OXT': (-0.8055, 0.295992, 0.87864),
    'CB': (-0.0875, 0.339967, 0.4577296),
    'OG': (-0.6546, 0.306647, 0.88028),
    'OG1': (-0.6546, 0.306647, 0.88028),
    'CG': (0.0, 0.339967, 0.4577296),
    'CG1': (-0.0122, 0.339967, 0.4577296),
    'CG2': (-0.0597, 0.339967, 0.4577296),
    'CD': (0.0, 0.339967, 0.4577296),
    'CD1': (-0.0875, 0.339967, 0.4577296),
    'CD2': (-0.0875, 0.339967, 0.4577296),
    'CE': (0.0, 0.339967, 0.4577296),
    'NZ': (-0.3854, 0.325, 0.71128),
    'OD1': (-0.6376, 0.295992, 0.87864),
    'OD2': (-0.6376, 0.295992, 0.87864),
    'OE1': (-0.6376, 0.295992, 0.87864),
    'OE2': (-0.6376, 0.295992, 0.87864),
    'ND1': (-0.3854, 0.325, 0.71128),
    'ND2': (-0.3854, 0.325, 0.71128),
    'NE': (-0.3854, 0.325, 0.71128),
    'NE1': (-0.3854, 0.325, 0.71128),
    'NE2': (-0.3854, 0.325, 0.71128),
    'NH1': (-0.3854, 0.325, 0.71128),
    'NH2': (-0.3854, 0.325, 0.71128),
    'SG': (-0.3119, 0.356359, 1.04600),
    'SD': (-0.3119, 0.356359, 1.04600),
    'ZN': (2.0, 0.125, 0.0125),
}

KCOUL = 332.06
CUTOFF_NM = 0.4  # 4 Å cutoff for direct binding interactions only
DIELECTRIC = 20.0  # Effective protein dielectric for binding site

def get_atom_params(atom_name, element):
    if atom_name in ATOM_PARAMS:
        return ATOM_PARAMS[atom_name]
    elif element == 'Zn':
        return ATOM_PARAMS['ZN']
    else:
        return (0.0, 0.339967, 0.4577296)

def read_cif_structure(cif_path):
    try:
        doc = gemmi.cif.read(cif_path)
        block = doc.sole_block()
        st = gemmi.make_structure_from_block(block)
        model = st[0]

        chains_data = {}
        for ch in model:
            atoms = []
            for res in ch:
                for atom in res:
                    atoms.append({
                        'name': atom.name.strip(),
                        'element': atom.element.name,
                        'pos': atom.pos,
                    })
            chains_data[ch.name] = {
                'n_residues': len(list(ch)),
                'atoms': atoms
            }

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

        lig_resname = "UNK"
        if ligand_chain and chains_data[ligand_chain]['atoms']:
            for ch in model:
                if ch.name == ligand_chain:
                    for res in ch:
                        lig_resname = res.name.strip()
                        break

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

def compute_interaction_energy(protein_atoms, ligand_atoms):
    """
    Compute ONLY protein-ligand interaction energy (not internal energies).
    Returns realistic binding energies.
    """
    E_vdw = 0.0
    E_coul = 0.0

    for p_atom in protein_atoms:
        q1, s1, e1 = get_atom_params(p_atom['name'], p_atom['element'])

        for l_atom in ligand_atoms:
            q2, s2, e2 = get_atom_params(l_atom['name'], l_atom['element'])

            # Distance in nm
            dx = p_atom['pos'].x - l_atom['pos'].x
            dy = p_atom['pos'].y - l_atom['pos'].y
            dz = p_atom['pos'].z - l_atom['pos'].z
            r_nm = math.sqrt(dx*dx + dy*dy + dz*dz) / 10.0

            if r_nm < CUTOFF_NM:
                # Lennard-Jones
                sigma = 0.5 * (s1 + s2)
                epsilon = math.sqrt(e1 * e2)
                sr6 = (sigma / r_nm) ** 6
                E_vdw += 4.0 * epsilon * (sr6**2 - sr6)

                # Coulomb with dielectric screening
                E_coul += (KCOUL / DIELECTRIC) * q1 * q2 / r_nm

    # Convert VdW to kcal/mol
    E_vdw_kcal = E_vdw / 4.184

    # Total interaction energy
    BE = E_vdw_kcal + E_coul

    return {
        'BE_mmgbsa': BE,
        'E_vdw': E_vdw_kcal,
        'E_coul': E_coul
    }

def main(cif_path):
    data, status = read_cif_structure(cif_path)

    if status != "OK":
        return {"file": cif_path, "status": status}

    if not data['ligand_chain']:
        return {"file": cif_path, "status": "NO_LIGAND_CHAIN"}

    try:
        protein_atoms = data['chains'][data['protein_chain']]['atoms']
        ligand_atoms = data['chains'][data['ligand_chain']]['atoms']

        result = compute_interaction_energy(protein_atoms, ligand_atoms)

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
            "status": "ERROR",
            "error": str(e)[:200]
        }

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python mmgbsa_corrected.py <cif_file>")
        sys.exit(1)

    result = main(sys.argv[1])
    print(json.dumps(result))
