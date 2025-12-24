#!/usr/bin/env python3
"""
Extract protein+ligand PDBs from AF3 CIF for AmberTools MMPBSA.
Handles multiple ligands (COMPETITION jobs) and Zn.
"""
import sys
import gemmi

def extract_complex(cif_path, output_prefix):
    """
    Extract from CIF:
    - complex.pdb (protein + first ligand + Zn if present)
    - receptor.pdb (protein + Zn if present, no ligand)
    - ligand.pdb (first ligand only, residue renamed to LIG)

    Returns ligand residue name (e.g., 'THR', 'SER')
    """
    st = gemmi.read_structure(cif_path)
    m = st[0]

    # Identify chains
    chains = []
    for c in m:
        nres = len(list(c))
        chains.append((c.name, nres))

    # Protein = largest chain (typically A)
    prot_chain = max(chains, key=lambda x: x[1])[0]

    # Ligand = chain B (1-3 residues)
    lig_chain = None
    for name, nres in chains:
        if name == 'B' and 1 <= nres <= 3:
            lig_chain = name
            break

    if lig_chain is None:
        raise ValueError(f"No ligand chain B found in {cif_path}")

    # Get ligand residue name
    lig_resname = None
    for c in m:
        if c.name == lig_chain:
            lig_resname = list(c)[0].name.strip()
            break

    # Check for Zn (typically chain D or last chain)
    has_zn = False
    zn_chain = None
    for c in m:
        for r in c:
            if r.name.strip().upper() == 'ZN':
                has_zn = True
                zn_chain = c.name
                break

    # Write complex (protein + ligand + Zn)
    st_complex = st.clone()
    m_complex = st_complex[0]
    for c in list(m_complex):
        if c.name not in [prot_chain, lig_chain] + ([zn_chain] if has_zn else []):
            m_complex.remove_chain(c.name)

    # Rename ligand residue to LIG in complex
    for c in m_complex:
        if c.name == lig_chain:
            for r in c:
                r.name = 'LIG'

    st_complex.write_pdb(f"{output_prefix}_complex.pdb")

    # Write receptor (protein + Zn, no ligand)
    st_receptor = st.clone()
    m_receptor = st_receptor[0]
    for c in list(m_receptor):
        if c.name not in [prot_chain] + ([zn_chain] if has_zn else []):
            m_receptor.remove_chain(c.name)
    st_receptor.write_pdb(f"{output_prefix}_receptor.pdb")

    # Write ligand (chain B only, renamed to LIG)
    st_ligand = st.clone()
    m_ligand = st_ligand[0]
    for c in list(m_ligand):
        if c.name != lig_chain:
            m_ligand.remove_chain(c.name)

    # Rename to LIG
    for c in m_ligand:
        for r in c:
            r.name = 'LIG'

    st_ligand.write_pdb(f"{output_prefix}_ligand.pdb")

    return lig_resname, has_zn

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python prep_from_cif.py input.cif output_prefix")
        sys.exit(1)

    cif_path = sys.argv[1]
    output_prefix = sys.argv[2]

    lig_resname, has_zn = extract_complex(cif_path, output_prefix)
    print(f"Extracted {output_prefix}_*.pdb")
    print(f"  Ligand: {lig_resname}")
    print(f"  Zn: {'Yes' if has_zn else 'No'}")
