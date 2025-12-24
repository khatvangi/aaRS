#!/usr/bin/env python3
import sys, json
import gemmi

AA3 = {
  "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
  "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
}

def pick_protein_and_ligand(model: gemmi.Model):
    chains = []
    for ch in model:
        nres = len(list(ch))
        nat = sum(1 for r in ch for _ in r)
        chains.append((ch.name, nres, nat))
    protein = max(chains, key=lambda x: x[1])[0]

    ligand = None
    # Prefer B if tiny
    for name, nres, _ in chains:
        if name == "B" and 1 <= nres <= 3:
            ligand = name
            break
    # Else smallest non-protein chain that has residues
    if ligand is None:
        other = [x for x in chains if x[0] != protein and x[1] >= 1]
        if other:
            ligand = min(other, key=lambda x: x[1])[0]
    return protein, ligand, chains

def iter_atoms(chain: gemmi.Chain):
    for res in chain:
        for atom in res:
            yield res, atom

def pdb_atom_line(
    rec, serial, atom_name, resname, chain_id, resseq, x, y, z, occ=1.00, b=50.0, element="C"
):
    # PDB fixed columns; atom_name needs right-justification rules
    an = atom_name.strip()
    if len(an) < 4:
        # If atom name starts with a letter, pad left to column 14
        if an[0].isalpha():
            an = f"{an:>4}"
        else:
            an = f"{an:<4}"
    else:
        an = an[:4]
    rn = f"{resname:>3}"[:3]
    el = element.strip().upper()
    return (
        f"{rec:<6}{serial:>5} {an} {rn:>3} {chain_id:1}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{occ:>6.2f}{b:>6.2f}          {el:>2}"
    )

def write_pdb(path, atom_records):
    with open(path, "w") as f:
        serial = 1
        last_chain = None
        for rec in atom_records:
            (rectype, atom_name, resname, chain_id, resseq, pos, bfac, element) = rec
            if last_chain is not None and chain_id != last_chain:
                f.write("TER\n")
            f.write(
                pdb_atom_line(
                    rectype, serial, atom_name, resname, chain_id, resseq,
                    pos.x, pos.y, pos.z, 1.00, bfac, element
                )
                + "\n"
            )
            serial += 1
            last_chain = chain_id
        f.write("END\n")

def main(cif_path, outdir):
    st = gemmi.read_structure(cif_path)
    st.setup_entities()
    model = st[0]

    protein_chain, ligand_chain, chains = pick_protein_and_ligand(model)
    if ligand_chain is None:
        raise SystemExit(f"No ligand chain detected in {cif_path}")

    # collect protein atoms (ATOM records only for standard residues; keep heteros like MSE as HETATM)
    protein_atoms = []
    zn_atoms = []
    ligand_atoms = []

    lig_original_resname = "UNK"

    # Identify Zn: any atom whose element is Zn OR residue name is ZN
    def is_zn(res, atom):
        return (atom.element.name.upper() == "ZN") or (res.name.strip().upper() == "ZN")

    for ch in model:
        if ch.name == protein_chain:
            for res, atom in iter_atoms(ch):
                rn = res.name.strip().upper()
                if is_zn(res, atom):
                    zn_atoms.append((res, atom))
                    continue
                rectype = "ATOM" if rn in AA3 else "HETATM"
                protein_atoms.append((rectype, atom.name.strip(), rn, protein_chain, res.seqid.num,
                                     atom.pos, float(atom.b_iso), atom.element.name))
        elif ch.name == ligand_chain:
            # ligand: force to LIG chain B resseq 1
            first_res = True
            for res, atom in iter_atoms(ch):
                if first_res:
                    lig_original_resname = res.name.strip().upper()
                    first_res = False
                ligand_atoms.append(("HETATM", atom.name.strip(), "LIG", "B", 1,
                                     atom.pos, float(atom.b_iso), atom.element.name))
        else:
            # other chains: if they contain Zn, keep as Zn; otherwise ignore
            for res, atom in iter_atoms(ch):
                if is_zn(res, atom):
                    zn_atoms.append((res, atom))

    # normalize Zn output: single-residue ZN in chain D, resseq 1, atom name ZN
    zn_records = []
    for res, atom in zn_atoms:
        zn_records.append(("HETATM", "ZN", "ZN", "D", 1, atom.pos, float(atom.b_iso), "ZN"))

    # Build outputs
    complex_records = []
    # protein as chain A (keep original chain id)
    complex_records.extend(protein_atoms)
    # Zn as chain D
    complex_records.extend(zn_records)
    # ligand as chain B
    complex_records.extend(ligand_atoms)

    receptor_records = []
    receptor_records.extend(protein_atoms)
    receptor_records.extend(zn_records)

    ligand_records = []
    ligand_records.extend(ligand_atoms)

    # write
    import os
    os.makedirs(outdir, exist_ok=True)
    write_pdb(os.path.join(outdir, "complex.pdb"), complex_records)
    write_pdb(os.path.join(outdir, "receptor.pdb"), receptor_records)
    write_pdb(os.path.join(outdir, "ligand.pdb"), ligand_records)

    meta = {
        "input_cif": cif_path,
        "protein_chain": protein_chain,
        "ligand_chain": ligand_chain,
        "ligand_original_resname": lig_original_resname,
        "zn_count": len(zn_records),
        "chains": chains,
    }
    with open(os.path.join(outdir, "meta.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(json.dumps(meta))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: prep_from_cif.py input.cif OUTDIR")
        sys.exit(2)
    main(sys.argv[1], sys.argv[2])
