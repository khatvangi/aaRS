#!/usr/bin/env python3
"""
Single-structure MM/GBSA with restrained minimization.
Uses OpenMM with OBC2 implicit solvent (GBSA).
"""
import sys
import json
import gemmi
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np
from collections import defaultdict

def read_cif_structure(cif_path):
    """Read CIF and identify chains"""
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
        # Fallback: smallest non-water chain
        non_water = [c for c in chains if c['n_residues'] >= 1]
        if non_water:
            ligand_chain = min(non_water, key=lambda x: x['n_residues'])

    # Find Zn atoms
    zn_atoms = []
    zn_coords = []
    for ch in model:
        for res in ch:
            for atom in res:
                if atom.element.name == 'Zn':
                    zn_atoms.append({'chain': ch.name, 'res': res.name, 'atom': atom.name})
                    zn_coords.append(atom.pos)

    return {
        'structure': st,
        'model': model,
        'protein_chain': protein_chain,
        'ligand_chain': ligand_chain,
        'zn_atoms': zn_atoms,
        'zn_coords': zn_coords
    }, "OK"

def setup_system_with_gb(pdb_file, protein_chain_id, ligand_chain_id, zn_present=False):
    """
    Setup OpenMM system with GBSA (OBC2).
    Add hydrogens to protein only (ligand amino acids are zwitterions, hard to parameterize).
    Returns: topology, positions, system
    """
    pdb = app.PDBFile(pdb_file)

    # Force field with implicit solvent
    forcefield = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')

    # Try adding hydrogens with residueTemplates=None to skip problematic residues
    modeller = app.Modeller(pdb.topology, pdb.positions)

    # Add hydrogens only to protein chain
    # This is a workaround: we'll add H to everything but handle ligand separately if needed
    try:
        modeller.addHydrogens(forcefield, pH=7.0)
    except Exception as e:
        # If addHydrogens fails, skip it and continue without H
        print(f"Warning: Could not add hydrogens ({str(e)[:100]}), continuing without them", file=sys.stderr)
        pass

    # Create system with GBSA
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.NoCutoff,  # Implicit solvent, no cutoff
        constraints=None,
        rigidWater=False
    )

    return modeller.topology, modeller.positions, system

def add_restraints(system, topology, positions, protein_chain_id, zn_indices, ligand_indices, zn_present):
    """
    Add positional restraints:
    - Protein heavy atoms: strong (1000 kJ/mol/nm^2)
    - Zn + 1st shell: strong (1000 kJ/mol/nm^2) if Zn present
    - Ligand: light (10 kJ/mol/nm^2) or none
    """
    # Get atom indices
    protein_heavy = []

    for atom_idx, atom in enumerate(topology.atoms()):
        # Protein heavy atoms
        if atom.residue.chain.id == protein_chain_id and atom.element.symbol != 'H':
            protein_heavy.append(atom_idx)

    # Add harmonic restraints
    restraint_force = mm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint_force.addPerParticleParameter("k")
    restraint_force.addPerParticleParameter("x0")
    restraint_force.addPerParticleParameter("y0")
    restraint_force.addPerParticleParameter("z0")

    # Strong restraints on protein heavy atoms
    k_strong = 1000.0 * unit.kilocalories_per_mole / unit.angstrom**2
    k_strong_kj = k_strong.value_in_unit(unit.kilojoules_per_mole / unit.nanometer**2)

    for idx in protein_heavy:
        pos = positions[idx]
        restraint_force.addParticle(idx, [
            k_strong_kj,
            pos.value_in_unit(unit.nanometer)[0],
            pos.value_in_unit(unit.nanometer)[1],
            pos.value_in_unit(unit.nanometer)[2]
        ])

    # Strong restraints on Zn if present
    if zn_present and len(zn_indices) > 0:
        for idx in zn_indices:
            pos = positions[idx]
            restraint_force.addParticle(idx, [
                k_strong_kj,
                pos.value_in_unit(unit.nanometer)[0],
                pos.value_in_unit(unit.nanometer)[1],
                pos.value_in_unit(unit.nanometer)[2]
            ])

    # Light restraints on ligand (optional)
    k_light = 10.0 * unit.kilocalories_per_mole / unit.angstrom**2
    k_light_kj = k_light.value_in_unit(unit.kilojoules_per_mole / unit.nanometer**2)

    for idx in ligand_indices:
        pos = positions[idx]
        restraint_force.addParticle(idx, [
            k_light_kj,
            pos.value_in_unit(unit.nanometer)[0],
            pos.value_in_unit(unit.nanometer)[1],
            pos.value_in_unit(unit.nanometer)[2]
        ])

    system.addForce(restraint_force)
    return system

def minimize_structure(system, topology, positions, max_iter=500):
    """Run restrained minimization"""
    integrator = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
    platform = mm.Platform.getPlatformByName('CPU')
    context = mm.Context(system, integrator, platform)
    context.setPositions(positions)

    # Minimize
    mm.LocalEnergyMinimizer.minimize(context, tolerance=10.0, maxIterations=max_iter)

    # Get minimized positions and energy
    state = context.getState(getEnergy=True, getPositions=True)
    energy = state.getPotentialEnergy()
    positions_min = state.getPositions()

    del context
    del integrator

    return positions_min, energy

def compute_mmgbsa(pdb_complex, protein_chain_id, ligand_chain_id):
    """
    Compute MM/GBSA binding energy.
    BE = E_complex - E_protein - E_ligand
    NOTE: pdb_complex already has hydrogens added during minimization
    """
    # Setup complex system (already has H)
    pdb = app.PDBFile(pdb_complex)
    forcefield = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')

    system_complex = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None,
        rigidWater=False
    )

    # Compute complex energy
    integrator_complex = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
    context_complex = mm.Context(system_complex, integrator_complex, mm.Platform.getPlatformByName('CPU'))
    context_complex.setPositions(pdb.positions)
    state_complex = context_complex.getState(getEnergy=True)
    E_complex = state_complex.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)

    # Create protein-only topology and positions
    modeller_protein = app.Modeller(pdb.topology, pdb.positions)
    to_delete = [atom for atom in modeller_protein.topology.atoms() if atom.residue.chain.id == ligand_chain_id]
    modeller_protein.delete(to_delete)

    system_protein = forcefield.createSystem(
        modeller_protein.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None,
        rigidWater=False
    )

    integrator_protein = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
    context_protein = mm.Context(system_protein, integrator_protein, mm.Platform.getPlatformByName('CPU'))
    context_protein.setPositions(modeller_protein.positions)
    state_protein = context_protein.getState(getEnergy=True)
    E_protein = state_protein.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)

    # Create ligand-only topology and positions
    modeller_ligand = app.Modeller(pdb.topology, pdb.positions)
    to_delete = [atom for atom in modeller_ligand.topology.atoms() if atom.residue.chain.id != ligand_chain_id]
    modeller_ligand.delete(to_delete)

    system_ligand = forcefield.createSystem(
        modeller_ligand.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None,
        rigidWater=False
    )

    integrator_ligand = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
    context_ligand = mm.Context(system_ligand, integrator_ligand, mm.Platform.getPlatformByName('CPU'))
    context_ligand.setPositions(modeller_ligand.positions)
    state_ligand = context_ligand.getState(getEnergy=True)
    E_ligand = state_ligand.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)

    # Binding energy
    BE = E_complex - E_protein - E_ligand

    # Cleanup
    del context_complex, context_protein, context_ligand
    del integrator_complex, integrator_protein, integrator_ligand

    return {
        'E_complex': E_complex,
        'E_protein': E_protein,
        'E_ligand': E_ligand,
        'BE_mmgbsa': BE
    }

def main(cif_path):
    """Main MM/GBSA workflow"""
    try:
        # Read structure
        info, status = read_cif_structure(cif_path)
        if status != "OK":
            return {"file": cif_path, "status": status}

        # Get ligand resname
        lig_resname = "UNK"
        if info['ligand_chain']:
            if len(info['ligand_chain']['residues']) > 0:
                lig_resname = info['ligand_chain']['residues'][0].name.strip()

        zn_present = len(info['zn_atoms']) > 0

        # Convert to PDB for OpenMM (gemmi -> PDB string -> file)
        # Simple approach: write temporary PDB
        import tempfile
        import os

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            temp_pdb = f.name
            info['structure'].write_pdb(temp_pdb)

        # Identify chain IDs first
        protein_chain_id = info['protein_chain']['name']
        ligand_chain_id = info['ligand_chain']['name'] if info['ligand_chain'] else None

        if not ligand_chain_id:
            os.unlink(temp_pdb)
            return {"file": cif_path, "status": "NO_LIGAND_CHAIN"}

        # Setup system and run minimization
        try:
            topology, positions, system = setup_system_with_gb(temp_pdb, protein_chain_id, ligand_chain_id, zn_present)

            # Find Zn and ligand indices
            zn_indices = []
            ligand_indices = []

            for atom_idx, atom in enumerate(topology.atoms()):
                if atom.element.symbol == 'Zn':
                    zn_indices.append(atom_idx)
                if atom.residue.chain.id == ligand_chain_id:
                    ligand_indices.append(atom_idx)

            # Add restraints
            system = add_restraints(system, topology, positions, protein_chain_id, zn_indices, ligand_indices, zn_present)

            # Minimize
            positions_min, energy_min = minimize_structure(system, topology, positions, max_iter=500)

            # Calculate ligand RMSD after minimization
            ligand_rmsd = 0.0
            if len(ligand_indices) > 0:
                rmsd_sum = 0.0
                for idx in ligand_indices:
                    pos_old = positions[idx].value_in_unit(unit.nanometer)
                    pos_new = positions_min[idx].value_in_unit(unit.nanometer)
                    dx = pos_new[0] - pos_old[0]
                    dy = pos_new[1] - pos_old[1]
                    dz = pos_new[2] - pos_old[2]
                    rmsd_sum += dx**2 + dy**2 + dz**2
                ligand_rmsd = np.sqrt(rmsd_sum / len(ligand_indices)) * 10.0  # Convert nm to Å

            # Write minimized structure
            with tempfile.NamedTemporaryFile(mode='w', suffix='_min.pdb', delete=False) as f:
                temp_pdb_min = f.name
                app.PDBFile.writeFile(topology, positions_min, open(temp_pdb_min, 'w'))

            # Compute MM/GBSA
            mmgbsa_result = compute_mmgbsa(temp_pdb_min, protein_chain_id, ligand_chain_id)

            # Cleanup temp files
            os.unlink(temp_pdb)
            os.unlink(temp_pdb_min)

            # Return results
            return {
                "file": cif_path,
                "status": "OK",
                "lig_resname": lig_resname,
                "protein_chain": protein_chain_id,
                "ligand_chain": ligand_chain_id,
                "Zn_present": zn_present,
                "n_Zn": len(zn_indices),
                "ligand_RMSD_after_min": ligand_rmsd,
                "BE_mmgbsa": mmgbsa_result['BE_mmgbsa'],
                "E_complex": mmgbsa_result['E_complex'],
                "E_protein": mmgbsa_result['E_protein'],
                "E_ligand": mmgbsa_result['E_ligand']
            }

        except Exception as e:
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
        print("Usage: python mmgbsa_single.py <cif_file>")
        sys.exit(1)

    result = main(sys.argv[1])
    print(json.dumps(result))
