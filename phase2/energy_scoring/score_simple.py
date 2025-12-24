#!/usr/bin/env python3
"""
Simplified energy scoring using generic parameters.
Avoids force field template matching issues with AF3 structures.
"""
import sys, math
import numpy as np
import gemmi

# Standard Amber ff14SB parameters for common atoms (Coulomb charges + LJ)
# Format: {atom_type: (charge_e, sigma_nm, epsilon_kJ_mol)}
ATOM_PARAMS = {
    # Backbone
    'N': (-0.4157, 0.325, 0.71128),
    'CA': (0.0337, 0.339967, 0.4577296),
    'C': (0.5973, 0.339967, 0.35982),
    'O': (-0.5679, 0.295992, 0.87864),
    # Hydrogens
    'H': (0.2719, 0.106908, 0.0656888),
    'HA': (0.0823, 0.264953, 0.0656888),
    'HB': (0.0, 0.264953, 0.0656888),
    'HG': (0.0, 0.264953, 0.0656888),
    'HD': (0.0, 0.264953, 0.0656888),
    'HE': (0.0, 0.264953, 0.0656888),
    'HZ': (0.0, 0.264953, 0.0656888),
    # Side chains (generic)
    'CB': (0.0, 0.339967, 0.4577296),
    'CG': (0.0, 0.339967, 0.4577296),
    'CD': (0.0, 0.339967, 0.4577296),
    'CE': (0.0, 0.339967, 0.4577296),
    'CZ': (0.0, 0.339967, 0.4577296),
    'OG': (-0.6546, 0.306647, 0.88028),   # Hydroxyl O
    'OG1': (-0.6546, 0.306647, 0.88028),
    'OD1': (-0.5931, 0.295992, 0.87864),
    'OD2': (-0.5931, 0.295992, 0.87864),
    'OE1': (-0.5931, 0.295992, 0.87864),
    'OE2': (-0.5931, 0.295992, 0.87864),
    'SG': (-0.3119, 0.356359, 1.04600),   # Cysteine S
    'SD': (-0.2737, 0.356359, 1.04600),   # Met S
    # Zinc
    'ZN': (2.0, 0.125, 0.0125),          # Approximate
}

KCOUL = 138.935456  # kJ/mol·nm / e^2

def get_atom_params(atom_name, element):
    """Get charge, sigma, epsilon for an atom"""
    # Try exact match first
    if atom_name in ATOM_PARAMS:
        return ATOM_PARAMS[atom_name]

    # Element-based fallback
    if element == 'H':
        return (0.0, 0.106908, 0.0656888)
    elif element == 'C':
        return (0.0, 0.339967, 0.4577296)
    elif element == 'N':
        return (-0.4, 0.325, 0.71128)
    elif element == 'O':
        return (-0.5, 0.295992, 0.87864)
    elif element == 'S':
        return (-0.3, 0.356359, 1.04600)
    elif element == 'Zn':
        return (2.0, 0.125, 0.0125)
    else:
        # Generic non-polar
        return (0.0, 0.35, 0.3)

def lj_energy(r_nm, sigma1, sigma2, eps1, eps2):
    """Lennard-Jones energy in kJ/mol"""
    sigma = 0.5 * (sigma1 + sigma2)
    epsilon = math.sqrt(eps1 * eps2)
    sr6 = (sigma / r_nm) ** 6
    return 4.0 * epsilon * (sr6**2 - sr6)

def coulomb_energy(r_nm, q1, q2):
    """Coulomb energy in kJ/mol"""
    return KCOUL * q1 * q2 / r_nm

def dist(a, b):
    """Distance in nm (gemmi uses Angstrom)"""
    dx = a.x - b.x
    dy = a.y - b.y
    dz = a.z - b.z
    return math.sqrt(dx*dx + dy*dy + dz*dz) / 10.0  # Convert Å → nm

def score_cif_simple(cif_path):
    """
    Simple energy scoring without OpenMM force fields.
    Calculates protein-ligand interaction energy directly.
    """
    try:
        st = gemmi.read_structure(cif_path)
    except Exception as e:
        return {"file": cif_path, "status": "LOAD_ERROR", "error": str(e)}

    m = st[0]

    # Identify chains
    chains = list(m)
    if len(chains) < 2:
        return {"file": cif_path, "status": "TOO_FEW_CHAINS"}

    # Assume chain B is ligand (or smallest chain)
    lig_chain = None
    for ch in chains:
        if ch.name == 'B' and 1 <= len(list(ch)) <= 3:
            lig_chain = ch
            break

    if lig_chain is None:
        # Find smallest non-ZN chain
        candidates = [(len(list(ch)), ch) for ch in chains
                      if 1 <= len(list(ch)) <= 3 and list(ch)[0].name.strip() != 'ZN']
        if candidates:
            candidates.sort()
            lig_chain = candidates[0][1]

    if lig_chain is None:
        return {"file": cif_path, "status": "NO_LIG_CHAIN"}

    # Collect protein atoms (all chains except ligand and ZN)
    prot_atoms = []
    lig_atoms = []
    zn_present = 0

    for ch in m:
        for r in ch:
            if r.name.strip() == 'ZN':
                zn_present = 1
                continue

            for a in r:
                atom_info = {
                    'pos': a.pos,
                    'name': a.name.strip(),
                    'element': a.element.name,
                }
                atom_info['params'] = get_atom_params(atom_info['name'], atom_info['element'])

                if ch.name == lig_chain.name:
                    lig_atoms.append(atom_info)
                else:
                    prot_atoms.append(atom_info)

    if not lig_atoms:
        return {"file": cif_path, "status": "NO_LIG_ATOMS"}

    # Calculate interaction energy
    E_vdw = 0.0
    E_coul = 0.0
    cutoff_nm = 1.2

    for p_atom in prot_atoms:
        for l_atom in lig_atoms:
            r = dist(p_atom['pos'], l_atom['pos'])

            if r < cutoff_nm:
                # Get parameters
                q1, s1, e1 = p_atom['params']
                q2, s2, e2 = l_atom['params']

                # VdW
                E_vdw += lj_energy(r, s1, s2, e1, e2)

                # Coulomb
                E_coul += coulomb_energy(r, q1, q2)

    E_total = E_vdw + E_coul
    E_total_kcal = E_total / 4.184  # kJ/mol → kcal/mol

    # Ligand info
    lig_resname = list(lig_chain)[0].name.strip() if list(lig_chain) else "UNK"

    return {
        "file": cif_path,
        "status": "OK",
        "lig_chain": lig_chain.name,
        "lig_resname": lig_resname,
        "zn_present": zn_present,
        "Eint_kcal": E_total_kcal,
        "Evdw_kcal": E_vdw / 4.184,
        "Ecoul_kcal": E_coul / 4.184,
        "n_prot_atoms": len(prot_atoms),
        "n_lig_atoms": len(lig_atoms),
    }

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python score_simple.py <cif_path>")
        sys.exit(1)

    cif_path = sys.argv[1]
    result = score_cif_simple(cif_path)

    import json
    print(json.dumps(result))
