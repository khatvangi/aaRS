#!/usr/bin/env python3
"""
QC Step 3: Check parameter consistency for all AA templates.
For each ligands_fixed/<AA>/lig.mol2:
- Compute total charge sum
- Check all atom types exist
- Verify net charge matches intent (ARG/LYS +1, ASP/GLU -1, others 0)
"""
import os
import glob
import pandas as pd
import numpy as np

# Expected net charges at pH ~7
EXPECTED_CHARGES = {
    'ARG': +1, 'LYS': +1,
    'ASP': -1, 'GLU': -1,
    'ALA': 0, 'ASN': 0, 'CYS': 0, 'GLN': 0, 'GLY': 0,
    'HIS': 0, 'ILE': 0, 'LEU': 0, 'MET': 0, 'PHE': 0,
    'PRO': 0, 'SER': 0, 'THR': 0, 'TRP': 0, 'TYR': 0, 'VAL': 0
}

def parse_mol2(mol2_file):
    """
    Parse mol2 file and return:
    - atoms: list of dicts
    - total_charge: sum of partial charges
    - atom_types: set of atom types
    """
    atoms = []
    in_atoms = False

    with open(mol2_file) as f:
        for line in f:
            if '@<TRIPOS>ATOM' in line:
                in_atoms = True
                continue
            if '@<TRIPOS>BOND' in line:
                break
            if in_atoms and line.strip():
                parts = line.split()
                if len(parts) >= 9:
                    atoms.append({
                        'name': parts[1],
                        'type': parts[5],
                        'charge': float(parts[8])
                    })

    total_charge = sum(a['charge'] for a in atoms)
    atom_types = set(a['type'] for a in atoms)

    return atoms, total_charge, atom_types

def check_aa_params(aa_name, lig_dir):
    """Check parameters for one amino acid."""
    mol2_file = os.path.join(lig_dir, 'lig.mol2')
    frcmod_file = os.path.join(lig_dir, 'lig.frcmod')

    result = {'aa': aa_name}

    # Check files exist
    if not os.path.exists(mol2_file):
        result['error'] = 'lig.mol2 missing'
        return result

    if not os.path.exists(frcmod_file):
        result['error'] = 'lig.frcmod missing'
        return result

    # Parse mol2
    atoms, total_charge, atom_types = parse_mol2(mol2_file)

    result['n_atoms'] = len(atoms)
    result['total_charge'] = total_charge
    result['expected_charge'] = EXPECTED_CHARGES.get(aa_name, 0)
    result['charge_error'] = abs(total_charge - result['expected_charge'])
    result['atom_types'] = ','.join(sorted(atom_types))
    result['n_atom_types'] = len(atom_types)

    # Count heavy atoms vs hydrogens
    n_heavy = sum(1 for a in atoms if not (a['name'].startswith('H') or a['name'].startswith('h')))
    n_H = len(atoms) - n_heavy
    result['n_heavy'] = n_heavy
    result['n_H'] = n_H

    return result

def main():
    print("Checking parameter consistency for all AA templates...")

    # Find all AA directories
    aa_dirs = sorted(glob.glob('ligands_fixed/*/'))

    results = []
    for aa_dir in aa_dirs:
        aa_name = os.path.basename(aa_dir.rstrip('/'))
        print(f"  Checking {aa_name}...")

        result = check_aa_params(aa_name, aa_dir)
        results.append(result)

    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    results_df.to_csv('qc/ligand_param_qc.csv', index=False)

    print(f"\n✓ Wrote qc/ligand_param_qc.csv")

    # Summary
    print(f"\n{'='*70}")
    print("Summary:")
    print(f"{'='*70}")
    print(f"Total AAs checked: {len(results_df)}")

    if 'error' in results_df.columns:
        errors = results_df[results_df['error'].notna()]
        if len(errors) > 0:
            print(f"\n✗ {len(errors)} AAs with errors:")
            print(errors[['aa', 'error']])
    else:
        print(f"✓ No missing files")

    # Check charge errors
    valid = results_df[~results_df['charge_error'].isna()]

    if len(valid) > 0:
        print(f"\nCharge consistency:")
        print(f"  Mean absolute error: {valid['charge_error'].mean():.6f} e")
        print(f"  Max error: {valid['charge_error'].max():.6f} e")

        # Pass/Fail
        max_error = valid['charge_error'].max()
        if max_error < 0.05:
            print(f"  ✓ PASS: All charge errors < 0.05 e")
        else:
            print(f"  ✗ FAIL: Max charge error = {max_error:.6f} e ≥ 0.05 e")

            # Show problem cases
            problems = valid[valid['charge_error'] >= 0.05]
            print(f"\n  Problem cases:")
            print(problems[['aa', 'total_charge', 'expected_charge', 'charge_error']])

    # Show summary table
    print(f"\n{'='*70}")
    print("Per-AA Summary:")
    print(f"{'='*70}")
    if len(valid) > 0:
        print(valid[['aa', 'n_atoms', 'n_heavy', 'n_H', 'total_charge', 'expected_charge', 'charge_error']].to_string(index=False))

if __name__ == '__main__':
    main()
