#!/usr/bin/env python3
"""
Scan all AF3 CIF files and build manifest with chain info, Zn presence, ligand chain.
"""
import glob, os, math
import pandas as pd
import gemmi

AA = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','SEC','PYL'}
WATER = {'HOH','WAT','DOD'}

def dist(a,b):
    dx=a.x-b.x; dy=a.y-b.y; dz=a.z-b.z
    return math.sqrt(dx*dx+dy*dy+dz*dz)

def detect_ligand_chain(model):
    """Detect ligand chain - prefer chain B if 1-3 residues"""
    chains = list(model)
    # Prefer chain B
    for ch in chains:
        res = list(ch)
        if ch.name == 'B' and 1 <= len(res) <= 3:
            return ch.name
    # Fallback: smallest non-water chain with <=3 residues
    cand=[]
    for ch in chains:
        res = list(ch)
        if 1 <= len(res) <= 3:
            names = {r.name.strip() for r in res}
            if not names.issubset(WATER):
                cand.append((len(res), ch.name, sorted(list(names))))
    cand.sort()
    return cand[0][1] if cand else None

def scan_one(path):
    """Scan a single CIF file"""
    try:
        st = gemmi.read_structure(path)
        m = st[0]
        chain_names = [c.name for c in m]
        chain_nres = {c.name: len(list(c)) for c in m}

        # Collect Zn atoms and ligand atoms
        zn_atoms = []
        lig_atoms = []
        lig_chain = detect_ligand_chain(m)

        for ch in m:
            for r in ch:
                for a in r:
                    el = a.element.name
                    if el == 'Zn' or r.name.strip() == 'ZN':
                        zn_atoms.append(a.pos)
                    if lig_chain and ch.name == lig_chain:
                        lig_atoms.append(a.pos)

        zn_present = int(len(zn_atoms) > 0)

        # Zn-ligand min distance (Angstrom)
        zn_lig_min = None
        if zn_atoms and lig_atoms:
            zn_lig_min = min(dist(z, l) for z in zn_atoms for l in lig_atoms)

        # Ligand residue name
        lig_resname = None
        if lig_chain:
            res = list(next(c for c in m if c.name == lig_chain))
            if res:
                lig_resname = res[0].name.strip()

        return {
            "file": path,
            "chains": ",".join(chain_names),
            "chain_nres": ";".join([f"{k}:{v}" for k,v in chain_nres.items()]),
            "lig_chain": lig_chain,
            "lig_resname": lig_resname,
            "zn_present": zn_present,
            "zn_lig_min_dist_A": zn_lig_min if zn_lig_min is not None else None,
            "zn_atoms": len(zn_atoms),
        }
    except Exception as e:
        return {"file": path, "error": str(e)}

def main():
    # Find all CIF files
    base = "/storage/kiran-stuff/aaRS/phase2"
    os.chdir(base)

    files = sorted(glob.glob("**/*_model.cif", recursive=True))
    print(f"Found {len(files)} CIF files")

    rows=[]
    for i, f in enumerate(files):
        if (i+1) % 100 == 0:
            print(f"  Scanned {i+1}/{len(files)}")
        rows.append(scan_one(f))

    df = pd.DataFrame(rows)
    df.to_csv("energy_scoring/manifest_cifs.csv", index=False)

    print(f"\n✓ Wrote energy_scoring/manifest_cifs.csv ({len(df)} rows)")
    print(f"\nZn present: {df['zn_present'].sum()} / {len(df)}")
    print(f"Ligand detected: {df['lig_chain'].notna().sum()} / {len(df)}")
    print("\nTop ligands:")
    print(df['lig_resname'].value_counts().head(10))

if __name__ == "__main__":
    main()
