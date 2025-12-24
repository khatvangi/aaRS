import glob, pandas as pd, gemmi

def detect(st):
    m=st[0]
    chains=[]
    for c in m:
        nres=len(list(c))
        nat=sum(1 for r in c for _ in r)
        chains.append((c.name, nres, nat))
    # protein = max residues
    prot=max(chains, key=lambda x:x[1])[0]
    # ligand: prefer chain B with 1-3 residues; else smallest non-protein chain
    lig=None
    for name,nres,_ in chains:
        if name=="B" and 1<=nres<=3:
            lig=name; break
    if lig is None:
        other=[x for x in chains if x[0]!=prot and x[1]>=1]
        if other:
            lig=min(other, key=lambda x:x[1])[0]
    # ligand resname
    lig_res="UNK"
    if lig:
        for c in m:
            if c.name==lig:
                lig_res=list(c)[0].name.strip()
                break
    zn=False
    for c in m:
        for r in c:
            for a in r:
                if a.element.name.upper()=="ZN" or r.name.strip().upper()=="ZN":
                    zn=True; break
    return prot, lig, lig_res, zn, chains

def main():
    files=sorted(glob.glob("**/*.cif", recursive=True))
    rows=[]
    for f in files:
        try:
            st=gemmi.read_structure(f)
            prot, lig, lig_res, zn, chains = detect(st)
            rows.append({
                "file": f,
                "protein_chain": prot,
                "ligand_chain": lig,
                "lig_resname": lig_res,
                "zn_present": int(zn),
                "chains": ";".join([f"{n}:{r}:{a}" for n,r,a in chains])
            })
        except Exception as e:
            rows.append({"file": f, "status": "ERROR", "error": str(e)[:200]})
    df=pd.DataFrame(rows)
    df.to_csv("manifest.csv", index=False)
    print("Wrote manifest.csv", len(df))
    print(df["zn_present"].value_counts(dropna=False))
    print(df["ligand_chain"].value_counts(dropna=False).head(10))
if __name__=="__main__":
    main()
