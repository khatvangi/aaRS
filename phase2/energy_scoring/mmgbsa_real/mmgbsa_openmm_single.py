import sys, json, math
import numpy as np
import gemmi
import openmm as mm
import openmm.app as app
from openmm import unit

# --- Settings (keep fixed for comparability)
SOLUTE_DIELECTRIC = 1.0
SOLVENT_DIELECTRIC = 78.5
IONIC_STRENGTH = 0.15  # molar (GBSAOBC in OpenMM doesn't explicitly use salt; kept for record)
DO_MINIMIZE = False     # set True if you want restrained minimization
MAX_MIN_ITERS = 300
RESTRAINT_K = 200.0     # kJ/mol/nm^2 on protein heavy atoms (+ Zn heavy atom)

# cutoff: GBSA is nonbonded with cutoff recommended; use NoCutoff for simplicity in single-frame
NONBONDED_METHOD = app.NoCutoff

# Forcefield: amber14 protein + ions. Try ions.xml if present.
FF_CANDIDATES = [
    ("amber14-all.xml", "amber14/tip3p.xml", "amber14/ions.xml"),
    ("amber14-all.xml",),
]

def read_structure(cif_path):
    doc = gemmi.cif.read(cif_path)
    block = doc.sole_block()
    st = gemmi.make_structure_from_block(block)
    return st

def pick_chains(model):
    chains=[]
    for c in model:
        nres=len(list(c))
        nat=sum(1 for r in c for _ in r)
        chains.append((c.name, nres, nat))
    prot=max(chains, key=lambda x:x[1])[0]
    lig=None
    for name,nres,_ in chains:
        if name=="B" and 1<=nres<=3:
            lig=name; break
    if lig is None:
        other=[x for x in chains if x[0]!=prot and x[1]>=1]
        if other:
            lig=min(other, key=lambda x:x[1])[0]
    return prot, lig

def write_temp_pdb(st, prot_chain, lig_chain, out_pdb):
    # Keep only prot+lig chains; preserve ZN if it is on protein chain (common)
    st2=st.clone()
    m=st2[0]
    keep=set([prot_chain])
    if lig_chain: keep.add(lig_chain)
    for c in list(m):
        if c.name not in keep:
            m.remove_chain(c.name)
    st2.setup_entities()
    st2.write_pdb(out_pdb)

def build_ff():
    last_err=None
    for tup in FF_CANDIDATES:
        try:
            ff=app.ForceField(*tup)
            return ff
        except Exception as e:
            last_err=e
    raise RuntimeError(f"ForceField load failed: {last_err}")

def add_gbsa(system, topology):
    # Add OpenMM GBSA-OBC2. This is the "real GBSA" part.
    gb = mm.GBSAOBCForce()
    gb.setSoluteDielectric(SOLUTE_DIELECTRIC)
    gb.setSolventDielectric(SOLVENT_DIELECTRIC)
    # Populate GB parameters from the NonbondedForce charges + default radii
    # We'll use OpenMM's helper by cloning charges into GB with standard radii.
    # Radii/scale factors: use common OBC2 defaults based on element.
    nb=None
    for f in system.getForces():
        if isinstance(f, mm.NonbondedForce):
            nb=f; break
    if nb is None:
        raise RuntimeError("No NonbondedForce found")

    def elem_from_atom(atom):
        # element symbol
        return atom.element.symbol if atom.element is not None else "C"

    # crude but stable OBC2 radii map (Å)
    RAD = {"H":1.2, "C":1.7, "N":1.55, "O":1.52, "S":1.8, "P":1.8, "ZN":1.39}
    SCALE = {"H":0.85, "C":0.72, "N":0.79, "O":0.85, "S":0.96, "P":0.86, "ZN":0.8}

    atoms=list(topology.atoms())
    for i,atom in enumerate(atoms):
        q, sig, eps = nb.getParticleParameters(i)
        sym=elem_from_atom(atom).upper()
        r=RAD.get(sym, 1.7) * unit.angstrom
        s=SCALE.get(sym, 0.8)
        gb.addParticle(q, r, s)

    system.addForce(gb)

def add_restraints(system, topology, positions, restrain_atom_indices):
    # harmonic position restraints
    force = mm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addPerParticleParameter("k")
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    for i in restrain_atom_indices:
        p = positions[i]
        force.addParticle(i, [RESTRAINT_K, p.x, p.y, p.z])
    system.addForce(force)

def compute_G(topology, positions, ff, keep_ligand_chain=None, drop_ligand=False):
    # optionally drop ligand (receptor energy) or drop protein (ligand energy)
    modeller = app.Modeller(topology, positions)
    if drop_ligand and keep_ligand_chain:
        # delete ligand chain(s)
        dels=[]
        for ch in modeller.topology.chains():
            if ch.id in keep_ligand_chain:
                for r in ch.residues():
                    dels.append(r)
        if dels: modeller.delete(dels)

    if (not drop_ligand) and keep_ligand_chain is None:
        pass

    # add hydrogens based on ff
    modeller.addHydrogens(ff)

    system = ff.createSystem(modeller.topology, nonbondedMethod=NONBONDED_METHOD, constraints=app.HBonds)
    add_gbsa(system, modeller.topology)

    # Restrain protein heavy atoms (and Zn) if minimizing
    if DO_MINIMIZE:
        lig_set=set()
        if keep_ligand_chain:
            for a in modeller.topology.atoms():
                if a.residue.chain.id in keep_ligand_chain:
                    lig_set.add(a.index)
        restr=[]
        for a in modeller.topology.atoms():
            sym = a.element.symbol.upper() if a.element is not None else ""
            if a.index not in lig_set and sym != "H":
                restr.append(a.index)
        add_restraints(system, modeller.topology, modeller.positions, restr)

    integrator = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
    platform = mm.Platform.getPlatformByName("CPU")  # stable + scalable on 64 cores
    ctx = mm.Context(system, integrator, platform)
    ctx.setPositions(modeller.positions)

    if DO_MINIMIZE:
        mm.LocalEnergyMinimizer.minimize(ctx, maxIterations=MAX_MIN_ITERS)

    state = ctx.getState(getEnergy=True)
    G = state.getPotentialEnergy()  # includes MM + GBSA (since GBSA force was added)
    return G.value_in_unit(unit.kilocalories_per_mole)

def main(cif_path):
    st = read_structure(cif_path)
    m = st[0]
    prot, lig = pick_chains(m)
    if lig is None:
        return {"file": cif_path, "status": "NO_LIGAND_CHAIN"}

    tmp="__tmp_complex.pdb"
    write_temp_pdb(st, prot, lig, tmp)

    pdb = app.PDBFile(tmp)
    ff = build_ff()

    lig_ids=[lig]  # ligand chain id (likely 'B')
    try:
        # complex
        Gcx = compute_G(pdb.topology, pdb.positions, ff, keep_ligand_chain=lig_ids, drop_ligand=False)
        # receptor (protein only)
        Gr  = compute_G(pdb.topology, pdb.positions, ff, keep_ligand_chain=lig_ids, drop_ligand=True)

        # ligand only: easiest is to delete all non-ligand residues via Modeller trick
        modeller = app.Modeller(pdb.topology, pdb.positions)
        dels=[]
        for ch in modeller.topology.chains():
            if ch.id not in lig_ids:
                for r in ch.residues():
                    dels.append(r)
        if dels: modeller.delete(dels)
        Gl = compute_G(modeller.topology, modeller.positions, ff, keep_ligand_chain=None, drop_ligand=False)

        dG = Gcx - Gr - Gl

        # ligand residue name
        lig_res="UNK"
        for ch in pdb.topology.chains():
            if ch.id==lig:
                lig_res=list(ch.residues())[0].name
                break

        # Zn present?
        zn_present=0
        for a in pdb.topology.atoms():
            if a.element is not None and a.element.symbol.upper()=="ZN":
                zn_present=1; break

        return {
            "file": cif_path,
            "status": "OK",
            "protein_chain": prot,
            "ligand_chain": lig,
            "lig_resname": lig_res,
            "zn_present": zn_present,
            "G_complex": Gcx,
            "G_receptor": Gr,
            "G_ligand": Gl,
            "BE_dG_bind": dG,
            "minimized": int(DO_MINIMIZE),
        }

    except Exception as e:
        return {"file": cif_path, "status": "ERROR", "error": str(e)[:200]}

if __name__=="__main__":
    if len(sys.argv)<2:
        print("Usage: python mmgbsa_openmm_single.py model.cif")
        sys.exit(1)
    print(json.dumps(main(sys.argv[1])))
