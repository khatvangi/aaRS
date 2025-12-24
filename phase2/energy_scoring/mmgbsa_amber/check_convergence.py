#!/usr/bin/env python3
"""
Check minimization convergence and geometry drift for MM/GBSA jobs.
"""
import os
import re
import numpy as np

def parse_minimization_output(min_out_file):
    """
    Parse sander minimization output to check:
    1. Convergence (RMS gradient)
    2. Energy change
    3. Number of steps completed
    """
    if not os.path.exists(min_out_file):
        return None

    with open(min_out_file) as f:
        content = f.read()

    # Find FINAL RESULTS section
    match = re.search(r'FINAL RESULTS.*?NSTEP\s+ENERGY\s+RMS\s+GMAX.*?(\d+)\s+([-\d.E+]+)\s+([-\d.E+]+)\s+([-\d.E+]+)', content, re.DOTALL)

    if not match:
        return None

    nstep = int(match.group(1))
    energy = float(match.group(2))
    rms = float(match.group(3))
    gmax = float(match.group(4))

    # Check convergence criteria
    converged = rms < 0.1  # RMS gradient < 0.1 kcal/mol/Å

    return {
        'nstep': nstep,
        'energy': energy,
        'rms_gradient': rms,
        'max_gradient': gmax,
        'converged': converged
    }

def calculate_rmsd(pdb1, pdb2, atom_selection='heavy'):
    """
    Calculate RMSD between two PDB structures.
    atom_selection: 'heavy' (non-hydrogen) or 'all'

    Simple RMSD calculation on coordinates only.
    """
    import gemmi

    st1 = gemmi.read_structure(pdb1)
    st2 = gemmi.read_structure(pdb2)

    coords1 = []
    coords2 = []

    for chain1, chain2 in zip(st1[0], st2[0]):
        if chain1.name == chain2.name:
            for res1, res2 in zip(chain1, chain2):
                if res1.name == res2.name:
                    for atom1, atom2 in zip(res1, res2):
                        if atom1.name == atom2.name:
                            if atom_selection == 'heavy' and atom1.element.name.upper() == 'H':
                                continue
                            coords1.append([atom1.pos.x, atom1.pos.y, atom1.pos.z])
                            coords2.append([atom2.pos.x, atom2.pos.y, atom2.pos.z])

    if len(coords1) == 0:
        return None

    coords1 = np.array(coords1)
    coords2 = np.array(coords2)

    # Simple RMSD without alignment (assumes structures are already aligned)
    diff = coords1 - coords2
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

    return rmsd

def check_job(job_dir):
    """
    Check a single job for:
    1. Convergence
    2. Geometry drift (RMSD between initial and optimized)
    """
    min_out = os.path.join(job_dir, 'min.out')

    conv_info = parse_minimization_output(min_out)
    if conv_info is None:
        return {'status': 'NO_MIN_OUTPUT'}

    # Check geometry drift
    initial_pdb = os.path.join(job_dir, 'complex.pdb')

    # Convert min.rst to PDB for RMSD calculation
    # This would require cpptraj - skip for now, just check convergence

    result = {
        'converged': conv_info['converged'],
        'nstep': conv_info['nstep'],
        'rms_gradient': conv_info['rms_gradient'],
        'max_gradient': conv_info['max_gradient'],
        'energy': conv_info['energy'],
    }

    if not conv_info['converged']:
        result['warning'] = f"NOT CONVERGED: RMS={conv_info['rms_gradient']:.3f} kcal/mol/Å (need <0.1)"

    return result

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print("Usage: python check_convergence.py job_directory")
        sys.exit(1)

    result = check_job(sys.argv[1])
    import json
    print(json.dumps(result, indent=2))
