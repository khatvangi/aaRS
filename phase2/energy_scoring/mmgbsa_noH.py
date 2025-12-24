#!/usr/bin/env python3
"""
Single-structure MM/GBSA with restrained minimization (heavy atoms only).
Uses OpenMM with OBC2 implicit solvent (GBSA).
Works without adding hydrogens - consistent heavy-atom treatment.
"""
import sys
import json
import gemmi
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np
import tempfile
import os

def read_cif_structure(cif_path):
    """Read CIF and identify chains"""
    try:
        doc = gemmi.cif.read(cif_path)
        block = doc.sole_block()
        st = gemmi.make_structure_from_block(block)

        if len(st) == 0:
            return None, "NO_MODELS"

        model = st[0]

        # Identify chains
        chains = []
        for ch in model:
            residues = list(ch)
            chains.append({
                'name': ch.name,
                'n_residues': len(residues),
                'residues': residues
            })

        # Protein chain: largest chain
        protein_chain = max(chains, key=lambda x: x['n_residues'])

        # Ligand chain: prefer chain B if small (1-3 residues), else smallest
        ligand_chain = None
        for ch in chains:
            if ch['name'] == 'B' and 1 <= ch['n_residues'] <= 3:
                ligand_chain = ch
                break

        if not ligand_chain:
            # Fallback: smallest non-water chain != protein
            non_protein = [c for c in chains if c['name'] != protein_chain['name'] and c['n_residues'] >= 1]
            if non_protein:
                ligand_chain = min(non_protein, key=lambda x: x['n_residues'])

        # Find Zn atoms
        zn_atoms = []
        for ch in model:
            for res in ch:
                for atom in res:
                    if atom.element.name == 'Zn':
                        zn_atoms.append({'chain': ch.name, 'res': res.name, 'atom': atom.name})

        return {
            'structure': st,
            'model': model,
            'protein_chain': protein_chain,
            'ligand_chain': ligand_chain,
            'zn_atoms': zn_atoms
        }, "OK"
    except Exception as e:
        return None, f"CIF_READ_ERROR: {str(e)[:100]}"

def minimize_and_score(pdb_file, protein_chain_id, ligand_chain_id, zn_present):
    """
    Restrained minimization + MM/GBSA scoring (heavy atoms only).
    Returns: BE, ligand_RMSD, component energies
    """
    # Load structure
    pdb = app.PDBFile(pdb_file)

    # Force field with implicit solvent (OBC2 GBSA)
    forcefield = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')

    # Create system WITHOUT hydrogens (implicit/obc2 will work on heavy atoms)
    # Use HCT or OBC2 GB model
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None,
        rigidWater=False,
        implicitSolvent=app.OBC2,  # Explicitly set GB model
        soluteDielectric=1.0,
        solventDielectric=78.5
    )

    # Identify atom indices
    protein_heavy = []
    zn_indices = []
    ligand_indices = []
    ligand_positions_orig = []

    for atom_idx, atom in enumerate(pdb.topology.atoms()):
        if atom.residue.chain.id == protein_chain_id and atom.element.symbol != 'H':
            protein_heavy.append(atom_idx)
        elif atom.element.symbol == 'Zn':
            zn_indices.append(atom_idx)
        elif atom.residue.chain.id == ligand_chain_id:
            ligand_indices.append(atom_idx)
            ligand_positions_orig.append(pdb.positions[atom_idx])

    # Add restraints
    restraint_force = mm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint_force.addPerParticleParameter("k")
    restraint_force.addPerParticleParameter("x0")
    restraint_force.addPerParticleParameter("y0")
    restraint_force.addPerParticleParameter("z0")

    # Strong restraints on protein heavy atoms (1000 kcal/mol/Å^2)
    k_strong = 1000.0 * unit.kilocalories_per_mole / unit.angstrom**2
    k_strong_kj = k_strong.value_in_unit(unit.kilojoules_per_mole / unit.nanometer**2)

    for idx in protein_heavy:
        pos = pdb.positions[idx]
        restraint_force.addParticle(idx, [
            k_strong_kj,
            pos.value_in_unit(unit.nanometer)[0],
            pos.value_in_unit(unit.nanometer)[1],
            pos.value_in_unit(unit.nanometer)[2]
        ])

    # Strong restraints on Zn if present
    if zn_present and len(zn_indices) > 0:
        for idx in zn_indices:
            pos = pdb.positions[idx]
            restraint_force.addParticle(idx, [
                k_strong_kj,
                pos.value_in_unit(unit.nanometer)[0],
                pos.value_in_unit(unit.nanometer)[1],
                pos.value_in_unit(unit.nanometer)[2]
            ])

    # Light restraints on ligand (10 kcal/mol/Å^2)
    k_light = 10.0 * unit.kilocalories_per_mole / unit.angstrom**2
    k_light_kj = k_light.value_in_unit(unit.kilojoules_per_mole / unit.nanometer**2)

    for idx in ligand_indices:
        pos = pdb.positions[idx]
        restraint_force.addParticle(idx, [
            k_light_kj,
            pos.value_in_unit(unit.nanometer)[0],
            pos.value_in_unit(unit.nanometer)[1],
            pos.value_in_unit(unit.nanometer)[2]
        ])

    system.addForce(restraint_force)

    # Minimize
    integrator = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
    platform = mm.Platform.getPlatformByName('CPU')
    context = mm.Context(system, integrator, platform)
    context.setPositions(pdb.positions)

    mm.LocalEnergyMinimizer.minimize(context, tolerance=10.0, maxIterations=500)

    state = context.getState(getEnergy=True, getPositions=True)
    positions_min = state.getPositions()
    E_minimized = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)

    del context, integrator

    # Calculate ligand RMSD after minimization
    ligand_rmsd = 0.0
    if len(ligand_indices) > 0:
        rmsd_sum = 0.0
        for i, idx in enumerate(ligand_indices):
            pos_old = ligand_positions_orig[i].value_in_unit(unit.nanometer)
            pos_new = positions_min[idx].value_in_unit(unit.nanometer)
            dx = pos_new[0] - pos_old[0]
            dy = pos_new[1] - pos_old[1]
            dz = pos_new[2] - pos_old[2]
            rmsd_sum += dx**2 + dy**2 + dz**2
        ligand_rmsd = np.sqrt(rmsd_sum / len(ligand_indices)) * 10.0  # nm to Å

    # Write minimized structure
    with tempfile.NamedTemporaryFile(mode='w', suffix='_min.pdb', delete=False) as f:
        temp_pdb_min = f.name
        app.PDBFile.writeFile(pdb.topology, positions_min, open(temp_pdb_min, 'w'))

    # Compute MM/GBSA on minimized structure
    mmgbsa_result = compute_mmgbsa_from_pdb(temp_pdb_min, protein_chain_id, ligand_chain_id)

    os.unlink(temp_pdb_min)

    return {
        'BE_mmgbsa': mmgbsa_result['BE'],
        'E_complex': mmgbsa_result['E_complex'],
        'E_protein': mmgbsa_result['E_protein'],
        'E_ligand': mmgbsa_result['E_ligand'],
        'ligand_RMSD_after_min': ligand_rmsd
    }

def compute_mmgbsa_from_pdb(pdb_file, protein_chain_id, ligand_chain_id):
    """
    Compute MM/GBSA: BE = E(complex) - E(protein) - E(ligand)
    """
    forcefield = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')

    # Complex energy
    pdb = app.PDBFile(pdb_file)
    system_complex = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.NoCutoff,
        implicitSolvent=app.OBC2,
        soluteDielectric=1.0,
        solventDielectric=78.5
    )
    integrator = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
    context = mm.Context(system_complex, integrator, mm.Platform.getPlatformByName('CPU'))
    context.setPositions(pdb.positions)
    E_complex = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    del context, integrator

    # Protein-only energy
    modeller_protein = app.Modeller(pdb.topology, pdb.positions)
    to_delete = [atom for atom in modeller_protein.topology.atoms() if atom.residue.chain.id == ligand_chain_id]
    modeller_protein.delete(to_delete)

    system_protein = forcefield.createSystem(
        modeller_protein.topology,
        nonbondedMethod=app.NoCutoff,
        implicitSolvent=app.OBC2,
        soluteDielectric=1.0,
        solventDielectric=78.5
    )
    integrator = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
    context = mm.Context(system_protein, integrator, mm.Platform.getPlatformByName('CPU'))
    context.setPositions(modeller_protein.positions)
    E_protein = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    del context, integrator

    # Ligand-only energy
    modeller_ligand = app.Modeller(pdb.topology, pdb.positions)
    to_delete = [atom for atom in modeller_ligand.topology.atoms() if atom.residue.chain.id != ligand_chain_id]
    modeller_ligand.delete(to_delete)

    system_ligand = forcefield.createSystem(
        modeller_ligand.topology,
        nonbondedMethod=app.NoCutoff,
        implicitSolvent=app.OBC2,
        soluteDielectric=1.0,
        solventDielectric=78.5
    )
    integrator = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
    context = mm.Context(system_ligand, integrator, mm.Platform.getPlatformByName('CPU'))
    context.setPositions(modeller_ligand.positions)
    E_ligand = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    del context, integrator

    BE = E_complex - E_protein - E_ligand

    return {
        'BE': BE,
        'E_complex': E_complex,
        'E_protein': E_protein,
        'E_ligand': E_ligand
    }

def main(cif_path):
    """Main MM/GBSA workflow"""
    try:
        # Read structure
        info, status = read_cif_structure(cif_path)
        if status != "OK":
            return {"file": cif_path, "status": status}

        if not info['ligand_chain']:
            return {"file": cif_path, "status": "NO_LIGAND_CHAIN"}

        # Get ligand resname
        lig_resname = "UNK"
        if len(info['ligand_chain']['residues']) > 0:
            lig_resname = info['ligand_chain']['residues'][0].name.strip()

        zn_present = len(info['zn_atoms']) > 0
        protein_chain_id = info['protein_chain']['name']
        ligand_chain_id = info['ligand_chain']['name']

        # Convert to PDB
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            temp_pdb = f.name
            info['structure'].write_pdb(temp_pdb)

        try:
            result = minimize_and_score(temp_pdb, protein_chain_id, ligand_chain_id, zn_present)
            os.unlink(temp_pdb)

            return {
                "file": cif_path,
                "status": "OK",
                "lig_resname": lig_resname,
                "protein_chain": protein_chain_id,
                "ligand_chain": ligand_chain_id,
                "Zn_present": zn_present,
                **result
            }

        except Exception as e:
            if os.path.exists(temp_pdb):
                os.unlink(temp_pdb)
            return {
                "file": cif_path,
                "status": "OPENMM_ERROR",
                "error": str(e)[:200]
            }

    except Exception as e:
        return {
            "file": cif_path,
            "status": "ERROR",
            "error": str(e)[:200]
        }

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python mmgbsa_noH.py <cif_file>")
        sys.exit(1)

    result = main(sys.argv[1])
    print(json.dumps(result))
