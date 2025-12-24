#!/usr/bin/env python3
"""
Score a single CIF file - protein-ligand interaction energy.
Can run NO-MIN (single point) or MIN (restrained minimization).
"""
import sys, math
import numpy as np
import openmm as mm
import openmm.app as app
from openmm import unit

KCOUL = 138.935456  # kJ/mol·nm / e^2

def kcal(x): return x.value_in_unit(unit.kilocalories_per_mole)
def ang(x): return x.value_in_unit(unit.angstrom)

def pick_lig_chain(top):
    """Prefer chain B if 1-3 residues"""
    for ch in top.chains():
        if ch.id == "B" and 1 <= len(list(ch.residues())) <= 3:
            return ["B"]
    # Fallback: smallest non-water chain
    small=[]
    for ch in top.chains():
        res=list(ch.residues())
        if 1 <= len(res) <= 3:
            names={r.name.strip() for r in res}
            if not names.issubset({"HOH","WAT","DOD"}):
                small.append((len(res), ch.id))
    small.sort()
    return [small[0][1]] if small else []

def atom_set_for_chain(top, chain_ids):
    """Get atom indices for given chain IDs"""
    s=set()
    for a in top.atoms():
        if a.residue.chain.id in chain_ids:
            s.add(a.index)
    return s

def delete_resname(modeller, resname="ZN"):
    """Delete all residues with given name"""
    dels=[]
    for ch in modeller.topology.chains():
        for r in ch.residues():
            if r.name.strip() == resname:
                dels.append(r)
    if dels:
        modeller.delete(dels)

def interaction_force_from_system(system, prot_set, lig_set, cutoff_nm=1.2):
    """Create custom force for protein-ligand interactions only"""
    nb=None
    for f in system.getForces():
        if isinstance(f, mm.NonbondedForce):
            nb=f; break
    if nb is None:
        raise RuntimeError("No NonbondedForce found")

    expr = (
      "4*sqrt(epsilon1*epsilon2)*((0.5*(sigma1+sigma2)/r)^12 - (0.5*(sigma1+sigma2)/r)^6)"
      f" + {KCOUL}*charge1*charge2/r"
    )
    cn=mm.CustomNonbondedForce(expr)
    cn.addPerParticleParameter("charge")
    cn.addPerParticleParameter("sigma")
    cn.addPerParticleParameter("epsilon")
    cn.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
    cn.setCutoffDistance(cutoff_nm*unit.nanometer)

    for i in range(system.getNumParticles()):
        q,sig,eps = nb.getParticleParameters(i)
        cn.addParticle([q.value_in_unit(unit.elementary_charge),
                        sig.value_in_unit(unit.nanometer),
                        eps.value_in_unit(unit.kilojoule_per_mole)])

    cn.addInteractionGroup(prot_set, lig_set)

    # Add exclusions from exceptions
    for i in range(nb.getNumExceptions()):
        a,b,*_ = nb.getExceptionParameters(i)
        cn.addExclusion(int(a), int(b))
    return cn

def add_heavy_atom_restraints(system, topology, positions, prot_set, k_prot=200.0):
    """Add harmonic restraints on protein heavy atoms"""
    force = mm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addPerParticleParameter("k")
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    for atom in topology.atoms():
        i = atom.index
        if i in prot_set and atom.element.symbol != 'H':
            p = positions[i]
            force.addParticle(i, [k_prot, p.x, p.y, p.z])
    system.addForce(force)

def ligand_rmsd(pos0, pos1, lig_idx):
    """Calculate ligand heavy-atom RMSD in Angstrom"""
    pts0=[]; pts1=[]
    for i in lig_idx:
        pts0.append(pos0[i].value_in_unit(unit.angstrom))
        pts1.append(pos1[i].value_in_unit(unit.angstrom))
    if not pts0: return np.nan
    pts0=np.array(pts0); pts1=np.array(pts1)
    d=pts0-pts1
    return float(np.sqrt((d*d).sum(axis=1).mean()))

def score_cif(cif_path, mode="nomin", platform_name="CPU", remove_zn_if_fail=True, maxiter=300):
    """
    Score a CIF file.
    mode: "nomin" (single point) or "min" (restrained minimization)
    """
    try:
        cif = app.PDBxFile(cif_path)
    except Exception as e:
        return {"file":cif_path, "status":"LOAD_ERROR", "error":str(e)}

    top = cif.topology
    pos = cif.positions
    lig_ids = pick_lig_chain(top)
    if not lig_ids:
        return {"file":cif_path, "status":"NO_LIG_CHAIN"}

    ff = app.ForceField("amber14-all.xml")
    mod = app.Modeller(top, pos)

    def run_scoring(modeller, zn_removed_flag):
        """Core scoring function"""
        modeller.addHydrogens(ff)

        if mode == "min":
            system = ff.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
        else:
            system = ff.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff, constraints=None)

        lig_set = atom_set_for_chain(modeller.topology, lig_ids)
        prot_set = set(range(modeller.topology.getNumAtoms())) - lig_set

        # Add protein restraints if minimizing
        if mode == "min":
            add_heavy_atom_restraints(system, modeller.topology, modeller.positions, prot_set, k_prot=200.0)

        # Add interaction-only force
        cn = interaction_force_from_system(system, prot_set, lig_set)
        cn.setForceGroup(1)
        system.addForce(cn)

        # Create context
        integ = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
        platform = mm.Platform.getPlatformByName(platform_name)
        ctx = mm.Context(system, integ, platform)
        ctx.setPositions(modeller.positions)

        # Pre-minimization energy
        Epre = ctx.getState(getEnergy=True, groups=1<<1).getPotentialEnergy()

        result = {
            "file": cif_path,
            "status": "OK",
            "lig_chain": ",".join(lig_ids),
            "zn_removed": zn_removed_flag,
            "mode": mode,
        }

        if mode == "min":
            pos_pre = ctx.getState(getPositions=True).getPositions()

            # Minimize
            mm.LocalEnergyMinimizer.minimize(ctx, maxIterations=maxiter)

            # Post-minimization
            Epost = ctx.getState(getEnergy=True, groups=1<<1).getPotentialEnergy()
            pos_post = ctx.getState(getPositions=True).getPositions()

            # Ligand RMSD
            lig_idx=[]
            for atom in modeller.topology.atoms():
                if atom.residue.chain.id in lig_ids and atom.element.symbol != 'H':
                    lig_idx.append(atom.index)
            rmsd = ligand_rmsd(pos_pre, pos_post, lig_idx)

            result.update({
                "Eint_kcal_premin": kcal(Epre),
                "Eint_kcal_postmin": kcal(Epost),
                "ligand_RMSD_A": rmsd,
            })
        else:
            result["Eint_kcal"] = kcal(Epre)

        del ctx, integ, system
        return result

    # Try with Zn first
    try:
        return run_scoring(mod, 0)
    except Exception as e:
        if not remove_zn_if_fail:
            return {"file":cif_path, "status":"ERROR", "error":str(e)}

        # Retry after removing Zn
        try:
            mod2 = app.Modeller(top, pos)
            delete_resname(mod2, "ZN")
            result = run_scoring(mod2, 1)
            result["status"] = "OK_ZN_REMOVED"
            result["orig_error"] = str(e)
            return result
        except Exception as e2:
            return {"file":cif_path, "status":"ERROR_RETRY", "error":str(e2), "orig_error":str(e)}

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python score_single.py <cif_path> <mode>")
        print("  mode: nomin or min")
        sys.exit(1)

    cif_path = sys.argv[1]
    mode = sys.argv[2]

    result = score_cif(cif_path, mode=mode, platform_name="CPU")

    # Print result as JSON for parsing
    import json
    print(json.dumps(result))
