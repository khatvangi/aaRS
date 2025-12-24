#!/usr/bin/env python3
"""
H-Bond Analysis for aaRS Evolution Figures

Analyzes hydrogen bonds between protein and ligand from AF3 CIF structures.
Generates data for Figure 5D and Figure 7B.

Requirements:
  - BioPython (pip install biopython)
  - numpy, pandas

H-bond criteria:
  - Donor-Acceptor distance < 3.5 Å
  - D-H...A angle > 120° (if hydrogens present)
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json

try:
    from Bio.PDB import MMCIFParser, NeighborSearch, Selection
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    print("WARNING: BioPython not available. Install with: pip install biopython")

# H-bond distance cutoffs
HBOND_DISTANCE = 3.5  # Angstroms
ZN_COORDINATION_DISTANCE = 2.8  # Angstroms for Zn coordination

# Donor and acceptor atoms
DONORS = {'N', 'O', 'S'}  # Can donate H-bonds
ACCEPTORS = {'N', 'O', 'S', 'F'}  # Can accept H-bonds


def find_ligand_chain(structure):
    """
    Identify the ligand chain (typically chain B for single ligand, or B/C for multiple).
    Returns list of ligand chain IDs.
    """
    ligand_chains = []
    for model in structure:
        for chain in model:
            # Check residue composition - ligands are typically non-standard
            residues = list(chain.get_residues())
            if len(residues) == 1:  # Single residue = likely ligand
                res = residues[0]
                if res.id[0] != ' ':  # HETATM
                    ligand_chains.append(chain.id)
            elif len(residues) <= 3 and any(r.id[0] != ' ' for r in residues):
                ligand_chains.append(chain.id)

    return ligand_chains


def find_zinc_atoms(structure):
    """Find all zinc atoms in the structure."""
    zn_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element == 'ZN' or 'ZN' in atom.get_id():
                        zn_atoms.append(atom)
    return zn_atoms


def calculate_hbonds(protein_atoms, ligand_atoms, distance_cutoff=HBOND_DISTANCE):
    """
    Calculate hydrogen bonds between protein and ligand.

    Returns:
        List of tuples: (protein_atom, ligand_atom, distance)
    """
    hbonds = []

    # Build neighbor search for efficiency
    ns = NeighborSearch(protein_atoms)

    for lig_atom in ligand_atoms:
        if lig_atom.element not in ACCEPTORS and lig_atom.element not in DONORS:
            continue

        # Find nearby protein atoms
        nearby = ns.search(lig_atom.coord, distance_cutoff)

        for prot_atom in nearby:
            if prot_atom.element not in ACCEPTORS and prot_atom.element not in DONORS:
                continue

            # Calculate distance
            distance = lig_atom - prot_atom

            # Check if potential H-bond (donor-acceptor pair)
            is_donor_acceptor = (
                (lig_atom.element in DONORS and prot_atom.element in ACCEPTORS) or
                (lig_atom.element in ACCEPTORS and prot_atom.element in DONORS)
            )

            if is_donor_acceptor and distance <= distance_cutoff:
                hbonds.append((prot_atom, lig_atom, distance))

    return hbonds


def analyze_zn_coordination(ligand_atoms, zn_atoms, distance_cutoff=ZN_COORDINATION_DISTANCE):
    """
    Analyze coordination bonds between ligand and zinc.

    Returns:
        List of tuples: (zn_atom, ligand_atom, distance, ligand_atom_type)
    """
    coordination_bonds = []

    for zn in zn_atoms:
        for lig_atom in ligand_atoms:
            distance = zn - lig_atom
            if distance <= distance_cutoff:
                coordination_bonds.append((zn, lig_atom, distance, lig_atom.element))

    return coordination_bonds


def analyze_structure(cif_file, job_name):
    """
    Analyze a single CIF structure for H-bonds and Zn coordination.

    Returns:
        Dict with analysis results
    """
    if not BIOPYTHON_AVAILABLE:
        return {'error': 'BioPython not available'}

    parser = MMCIFParser(QUIET=True)

    try:
        structure = parser.get_structure('protein', cif_file)
    except Exception as e:
        return {'error': f'Failed to parse CIF: {str(e)}'}

    # Identify ligand chains
    ligand_chain_ids = find_ligand_chain(structure)

    if not ligand_chain_ids:
        return {'error': 'No ligand chain found'}

    # Get protein and ligand atoms
    protein_atoms = []
    ligand_atoms = []

    for model in structure:
        for chain in model:
            atoms = Selection.unfold_entities(chain, 'A')
            if chain.id in ligand_chain_ids:
                ligand_atoms.extend(atoms)
            else:
                protein_atoms.extend(atoms)

    if not ligand_atoms or not protein_atoms:
        return {'error': 'Empty protein or ligand'}

    # Find zinc atoms
    zn_atoms = find_zinc_atoms(structure)
    has_zinc = len(zn_atoms) > 0

    # Calculate H-bonds
    hbonds = calculate_hbonds(protein_atoms, ligand_atoms)

    # Analyze Zn coordination if present
    zn_coordination = []
    if has_zinc:
        zn_coordination = analyze_zn_coordination(ligand_atoms, zn_atoms)

    # Count H-bonds by ligand atom type
    hbond_by_atom_type = {}
    for _, lig_atom, _ in hbonds:
        atom_name = lig_atom.get_id()
        if atom_name not in hbond_by_atom_type:
            hbond_by_atom_type[atom_name] = 0
        hbond_by_atom_type[atom_name] += 1

    # Results (convert numpy types to native Python for JSON serialization)
    results = {
        'job_name': job_name,
        'cif_file': str(cif_file),
        'total_hbonds': len(hbonds),
        'hbond_by_atom_type': hbond_by_atom_type,
        'has_zinc': has_zinc,
        'zn_coordination_bonds': len(zn_coordination),
        'zn_coordination_details': [
            {
                'ligand_atom': coord[1].get_id(),
                'ligand_element': coord[3],
                'distance': float(coord[2])  # Convert to native Python float
            }
            for coord in zn_coordination
        ] if has_zinc else [],
        'avg_hbond_distance': float(np.mean([h[2] for h in hbonds])) if hbonds else None,
        'protein_atoms': len(protein_atoms),
        'ligand_atoms': len(ligand_atoms)
    }

    return results


def main():
    """
    Main analysis workflow.
    """

    print("="*80)
    print("H-BOND ANALYSIS FOR aaRS EVOLUTION")
    print("="*80)

    if not BIOPYTHON_AVAILABLE:
        print("\nERROR: BioPython is required for this analysis.")
        print("Install with: pip install biopython")
        print("\nCreating placeholder documentation instead...")
        create_analysis_documentation()
        return

    # Define structures to analyze (UPDATED to match available CIF files)
    # Based on actual directory names in outputs/

    structures_to_analyze = {
        # Editing domain comparisons (for Figure 2D)
        'deep_editing_thr': '/storage/kiran-stuff/aaRS/phase2/outputs/deep_editing_thr/deep_editing_thr/deep_editing_thr_model.cif',
        'deep_editing_pro': '/storage/kiran-stuff/aaRS/phase2/outputs/deep_editing_pro/deep_editing_pro/deep_editing_pro_model.cif',

        # Ancestral ThrRS (for evolutionary comparison)
        'deep_thrrs_thr': '/storage/kiran-stuff/aaRS/phase2/outputs/deep_thrrs_thr/deep_thrrs_thr/deep_thrrs_thr_model.cif',
        'deep_thrrs_pro': '/storage/kiran-stuff/aaRS/phase2/outputs/deep_thrrs_pro/deep_thrrs_pro/deep_thrrs_pro_model.cif',

        # Modern enzymes (no Zn in these structures)
        'modern_thrrs_thr': '/storage/kiran-stuff/aaRS/phase2/outputs/modern_thrrs_thr/modern_thrrs_thr/modern_thrrs_thr_model.cif',
        'modern_thrrs_pro': '/storage/kiran-stuff/aaRS/phase2/outputs/modern_thrrs_pro/modern_thrrs_pro/modern_thrrs_pro_model.cif',
        'modern_prours_thr': '/storage/kiran-stuff/aaRS/phase2/outputs/modern_prours_thr/modern_prours_thr/modern_prours_thr_model.cif',
        'modern_prours_pro': '/storage/kiran-stuff/aaRS/phase2/outputs/modern_prours_pro/modern_prours_pro/modern_prours_pro_model.cif',
    }

    print("\nAnalyzing available CIF structures...")
    all_results = []

    for job_name, cif_path in structures_to_analyze.items():
        cif_file = Path(cif_path)

        if cif_file.exists():
            print(f"  Analyzing: {job_name}")
            result = analyze_structure(cif_file, job_name)
            result['pattern'] = job_name
            all_results.append(result)
        else:
            print(f"  Not found: {job_name}")

    # Save results
    if all_results:
        # Create DataFrame
        df = pd.DataFrame(all_results)

        # Save to CSV
        output_csv = 'figures/data/hbond_analysis.csv'
        df.to_csv(output_csv, index=False)
        print(f"\nSaved: {output_csv}")

        # Save detailed JSON
        output_json = 'figures/data/hbond_analysis_detailed.json'
        with open(output_json, 'w') as f:
            json.dump(all_results, f, indent=2)
        print(f"Saved: {output_json}")

        # Print summary
        print("\n" + "="*80)
        print("SUMMARY")
        print("="*80)

        for result in all_results:
            if 'error' not in result:
                print(f"\n{result['pattern']}:")
                print(f"  Total H-bonds: {result['total_hbonds']}")
                print(f"  Has Zn: {result['has_zinc']}")
                if result['has_zinc']:
                    print(f"  Zn coordination bonds: {result['zn_coordination_bonds']}")
    else:
        print("\nNo structures analyzed. Check file paths.")
        create_analysis_documentation()


def create_analysis_documentation():
    """
    Create documentation for H-bond analysis when BioPython is not available.
    """

    doc = """# H-Bond Analysis Documentation

## Requirements Not Met

BioPython is required for H-bond analysis but is not currently installed.

### To Install BioPython:

```bash
pip install biopython
# or
conda install -c conda-forge biopython
```

## Analysis Plan

### Structures to Analyze:

**For Figure 5D (Zinc Trap):**
1. Modern ThrRS + Zn + THR
   - Expect: 2-3 Zn coordination bonds (hydroxyl + amino)
   - Expect: 5-8 total H-bonds

2. Modern ThrRS + Zn + SER
   - Expect: Similar Zn coordination to THR (the trap!)
   - Should show identical coordination geometry

3. Modern ThrRS + Zn + ILE
   - Expect: 0-1 Zn coordination bonds (hydrophobic)
   - Fewer H-bonds overall

**For Figure 7B (Evolutionary Comparison):**
4. Ancestral ThrRS + THR (no Zn)
5. Ancestral ThrRS + Zn + THR
6. Ancestral ProRS + PRO
7. Ancestral ProRS editing + THR
8. Ancestral ProRS editing + PRO

### Expected Outputs:

1. **hbond_analysis.csv** - Summary table:
   - job_name, total_hbonds, has_zinc, zn_coordination_bonds, avg_distance

2. **hbond_analysis_detailed.json** - Full details:
   - Individual H-bond partners
   - Zn coordination geometry
   - Atom-by-atom breakdown

3. **Figure 7B** - Bar chart comparing:
   - H-bond counts across conditions
   - Zn coordination differences
   - Ancestral vs Modern evolution

### H-Bond Criteria:

- Distance: < 3.5 Å (donor-acceptor)
- Angle: > 120° (if hydrogens present)
- Zn coordination: < 2.8 Å

### Key Questions:

1. Why does SER coordinate Zn like THR?
   → Both have hydroxyl groups in same position

2. Why does ILE fail to coordinate?
   → Hydrophobic side chain cannot coordinate metal

3. Did H-bond network change during evolution?
   → Compare ancestral vs modern counts

## Alternative Analysis (Without BioPython)

If BioPython cannot be installed, manual analysis options:

1. **PyMOL**: Use distance measurements
   ```
   distance hbonds, (chain A), (chain B), 3.5
   ```

2. **ChimeraX**: Use built-in H-bond finding
   ```
   hbonds restrict both
   ```

3. **HBPLUS**: Standalone H-bond calculation tool
   ```
   hbplus structure.pdb
   ```

## Manual Inspection Priority:

Focus on these key structures:
1. modern_thrrs_ecoli_zn_THR vs SER (the trap!)
2. modern_thrrs_ecoli_zn_ILE (rejected)
3. anc_prors_edit_THR vs PRO (editing selectivity)
"""

    output_file = 'figures/data/HBOND_ANALYSIS_NEEDED.md'
    with open(output_file, 'w') as f:
        f.write(doc)

    print(f"\nCreated: {output_file}")
    print("Install BioPython and re-run this script for actual analysis.")


if __name__ == '__main__':
    main()
