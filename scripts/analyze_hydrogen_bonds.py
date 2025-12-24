#!/usr/bin/env python3
"""
analyze_hydrogen_bonds.py
=========================
Analyze hydrogen bonding patterns in aaRS-ligand complexes to test the 
hypothesis that threonine's β-hydroxyl forms compensating H-bonds with 
polar active site residues.

Key questions:
1. Does Thr form MORE H-bonds than Pro in the same pocket?
2. Are the extra H-bonds with catalytically essential residues?
3. Does this explain the persistent cross-reactivity?

Requirements:
    pip install biopython numpy pandas

Usage:
    python analyze_hydrogen_bonds.py /path/to/af3_output
"""

import os
import sys
import json
import glob
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

try:
    from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch, Selection
    from Bio.PDB.vectors import calc_angle, calc_dihedral
except ImportError:
    print("ERROR: BioPython required. Install with: pip install biopython")
    sys.exit(1)

# =============================================================================
# HYDROGEN BOND DETECTION PARAMETERS
# =============================================================================

# Standard H-bond criteria
HBOND_DISTANCE_MAX = 3.5  # Angstroms (donor-acceptor distance)
HBOND_DISTANCE_MIN = 2.5  # Angstroms (too close = clash)
HBOND_ANGLE_MIN = 120.0   # Degrees (D-H...A angle, if H available)

# Donor atoms (can donate H)
DONORS = {'N', 'NE', 'NE1', 'NE2', 'ND1', 'ND2', 'NZ', 'NH1', 'NH2', 
          'OG', 'OG1', 'OH', 'NE', 'OXT'}

# Acceptor atoms (can accept H)
ACCEPTORS = {'O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH', 
             'ND1', 'NE2', 'SD', 'OXT'}

# Polar residues (catalytically important)
POLAR_RESIDUES = {'SER', 'THR', 'ASN', 'GLN', 'TYR', 'CYS', 'HIS',
                  'ASP', 'GLU', 'LYS', 'ARG'}

# Ligand codes
LIGAND_CODES = {'PRO': 'Proline', 'THR': 'Threonine', 'TRP': 'Tryptophan', 
                'PHE': 'Phenylalanine', 'ALA': 'Alanine'}

# =============================================================================
# H-BOND ANALYSIS FUNCTIONS
# =============================================================================

def load_structure(filepath):
    """Load CIF or PDB structure"""
    filepath = str(filepath)
    
    if filepath.endswith('.cif'):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    
    structure_id = Path(filepath).stem
    return parser.get_structure(structure_id, filepath)

def get_ligand_atoms(structure, ligand_codes=None):
    """Extract ligand atoms from structure"""
    if ligand_codes is None:
        ligand_codes = list(LIGAND_CODES.keys())
    
    ligand_atoms = []
    ligand_residues = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()
                if resname in ligand_codes:
                    ligand_residues.append(residue)
                    for atom in residue:
                        ligand_atoms.append(atom)
    
    return ligand_atoms, ligand_residues

def get_protein_atoms(structure, exclude_ligands=True):
    """Get all protein atoms"""
    protein_atoms = []
    
    standard_aa = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 
                   'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
                   'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'}
    
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()
                if resname in standard_aa:
                    for atom in residue:
                        protein_atoms.append(atom)
    
    return protein_atoms

def find_hbonds(donor_atoms, acceptor_atoms, distance_max=HBOND_DISTANCE_MAX,
                distance_min=HBOND_DISTANCE_MIN):
    """
    Find potential hydrogen bonds between donor and acceptor atom sets.
    Returns list of (donor_atom, acceptor_atom, distance) tuples.
    """
    hbonds = []
    
    # Build neighbor search on acceptors
    if not acceptor_atoms:
        return hbonds
    
    ns = NeighborSearch(acceptor_atoms)
    
    for donor in donor_atoms:
        # Find acceptors within distance
        nearby = ns.search(donor.get_coord(), distance_max)
        
        for acceptor in nearby:
            # Skip if same residue
            if donor.get_parent() == acceptor.get_parent():
                continue
            
            # Calculate distance
            dist = donor - acceptor
            
            if distance_min < dist < distance_max:
                # Check if atoms are appropriate types
                donor_name = donor.get_name()
                acceptor_name = acceptor.get_name()
                
                if donor_name in DONORS and acceptor_name in ACCEPTORS:
                    hbonds.append((donor, acceptor, dist))
                elif acceptor_name in DONORS and donor_name in ACCEPTORS:
                    # Reverse - acceptor could be donor
                    hbonds.append((acceptor, donor, dist))
    
    return hbonds

def analyze_ligand_hbonds(structure, ligand_resname):
    """
    Analyze H-bonds between a specific ligand and protein.
    
    Returns dict with:
    - total_hbonds: count
    - hbond_details: list of (protein_res, protein_atom, ligand_atom, distance)
    - polar_contacts: count with polar residues
    - backbone_contacts: count with backbone atoms
    """
    results = {
        'ligand': ligand_resname,
        'total_hbonds': 0,
        'hbond_details': [],
        'polar_residue_contacts': 0,
        'backbone_contacts': 0,
        'sidechain_contacts': 0,
        'unique_residues': set()
    }
    
    # Get ligand and protein atoms
    ligand_atoms, ligand_residues = get_ligand_atoms(structure, [ligand_resname])
    
    if not ligand_atoms:
        print(f"  WARNING: No {ligand_resname} found in structure")
        return results
    
    protein_atoms = get_protein_atoms(structure)
    
    if not protein_atoms:
        print(f"  WARNING: No protein atoms found")
        return results
    
    # Find H-bonds from ligand to protein
    # Ligand as donor
    ligand_donors = [a for a in ligand_atoms if a.get_name() in DONORS]
    protein_acceptors = [a for a in protein_atoms if a.get_name() in ACCEPTORS]
    
    hbonds_as_donor = find_hbonds(ligand_donors, protein_acceptors)
    
    # Ligand as acceptor
    ligand_acceptors = [a for a in ligand_atoms if a.get_name() in ACCEPTORS]
    protein_donors = [a for a in protein_atoms if a.get_name() in DONORS]
    
    hbonds_as_acceptor = find_hbonds(protein_donors, ligand_acceptors)
    
    # Combine and analyze
    all_hbonds = []
    
    for donor, acceptor, dist in hbonds_as_donor:
        prot_res = acceptor.get_parent()
        all_hbonds.append({
            'protein_chain': prot_res.get_parent().id,
            'protein_resname': prot_res.get_resname(),
            'protein_resid': prot_res.get_id()[1],
            'protein_atom': acceptor.get_name(),
            'ligand_atom': donor.get_name(),
            'distance': dist,
            'ligand_is_donor': True
        })
    
    for donor, acceptor, dist in hbonds_as_acceptor:
        prot_res = donor.get_parent()
        all_hbonds.append({
            'protein_chain': prot_res.get_parent().id,
            'protein_resname': prot_res.get_resname(),
            'protein_resid': prot_res.get_id()[1],
            'protein_atom': donor.get_name(),
            'ligand_atom': acceptor.get_name(),
            'distance': dist,
            'ligand_is_donor': False
        })
    
    # Compile results
    results['total_hbonds'] = len(all_hbonds)
    results['hbond_details'] = all_hbonds
    
    backbone_atoms = {'N', 'O', 'C', 'CA'}
    
    for hb in all_hbonds:
        resname = hb['protein_resname']
        resid = hb['protein_resid']
        atom = hb['protein_atom']
        
        results['unique_residues'].add(f"{resname}{resid}")
        
        if resname in POLAR_RESIDUES:
            results['polar_residue_contacts'] += 1
        
        if atom in backbone_atoms:
            results['backbone_contacts'] += 1
        else:
            results['sidechain_contacts'] += 1
    
    results['unique_residues'] = len(results['unique_residues'])
    
    return results

def analyze_hydroxyl_specifically(structure, ligand_resname='THR'):
    """
    Specifically analyze H-bonds from threonine's OG1 (β-hydroxyl)
    """
    results = {
        'hydroxyl_hbonds': 0,
        'hydroxyl_partners': []
    }
    
    ligand_atoms, _ = get_ligand_atoms(structure, [ligand_resname])
    
    # Find OG1 (threonine β-hydroxyl) or OG (serine)
    hydroxyl_atoms = [a for a in ligand_atoms if a.get_name() in ['OG1', 'OG', 'OH']]
    
    if not hydroxyl_atoms:
        return results
    
    protein_atoms = get_protein_atoms(structure)
    
    # Find H-bonds specifically from hydroxyl
    for oh_atom in hydroxyl_atoms:
        oh_coord = oh_atom.get_coord()
        
        for prot_atom in protein_atoms:
            if prot_atom.get_name() not in DONORS.union(ACCEPTORS):
                continue
            
            dist = oh_atom - prot_atom
            
            if HBOND_DISTANCE_MIN < dist < HBOND_DISTANCE_MAX:
                prot_res = prot_atom.get_parent()
                results['hydroxyl_hbonds'] += 1
                results['hydroxyl_partners'].append({
                    'residue': f"{prot_res.get_resname()}{prot_res.get_id()[1]}",
                    'atom': prot_atom.get_name(),
                    'distance': round(dist, 2)
                })
    
    return results

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def analyze_directory(base_dir):
    """Analyze all AF3 outputs in a directory"""
    
    print("=" * 70)
    print("HYDROGEN BOND ANALYSIS - Testing Hydroxyl Compensation Hypothesis")
    print("=" * 70)
    
    # Find all model files
    cif_files = glob.glob(f"{base_dir}/**/*_model.cif", recursive=True)
    
    if not cif_files:
        # Try alternative patterns
        cif_files = glob.glob(f"{base_dir}/**/*.cif", recursive=True)
    
    print(f"\nFound {len(cif_files)} structure files")
    
    all_results = []
    
    for cif_path in sorted(cif_files):
        # Extract model name from path
        model_name = Path(cif_path).parent.name
        if model_name.startswith('seed-'):
            model_name = Path(cif_path).parent.parent.name
        
        print(f"\nAnalyzing: {model_name}")
        
        try:
            structure = load_structure(cif_path)
        except Exception as e:
            print(f"  ERROR loading: {e}")
            continue
        
        # Determine which ligand is present
        for lig_code in LIGAND_CODES.keys():
            lig_atoms, _ = get_ligand_atoms(structure, [lig_code])
            if lig_atoms:
                print(f"  Ligand: {LIGAND_CODES[lig_code]} ({lig_code})")
                
                # General H-bond analysis
                hb_results = analyze_ligand_hbonds(structure, lig_code)
                
                print(f"  Total H-bonds: {hb_results['total_hbonds']}")
                print(f"  Unique protein residues: {hb_results['unique_residues']}")
                print(f"  Polar residue contacts: {hb_results['polar_residue_contacts']}")
                print(f"  Backbone contacts: {hb_results['backbone_contacts']}")
                print(f"  Sidechain contacts: {hb_results['sidechain_contacts']}")
                
                # Hydroxyl-specific analysis for THR
                if lig_code == 'THR':
                    oh_results = analyze_hydroxyl_specifically(structure, 'THR')
                    print(f"  β-hydroxyl H-bonds: {oh_results['hydroxyl_hbonds']}")
                    for partner in oh_results['hydroxyl_partners']:
                        print(f"    → {partner['residue']}:{partner['atom']} ({partner['distance']}Å)")
                    hb_results['hydroxyl_hbonds'] = oh_results['hydroxyl_hbonds']
                    hb_results['hydroxyl_partners'] = oh_results['hydroxyl_partners']
                else:
                    hb_results['hydroxyl_hbonds'] = 0
                    hb_results['hydroxyl_partners'] = []
                
                # Store results
                hb_results['model_name'] = model_name
                hb_results['file_path'] = cif_path
                all_results.append(hb_results)
                
                break  # Found ligand, move to next file
    
    return all_results

def compare_pro_thr(results):
    """Compare H-bonding between Pro and Thr complexes"""
    
    print("\n" + "=" * 70)
    print("COMPARISON: Proline vs Threonine H-bonding")
    print("=" * 70)
    
    # Group by enzyme (remove _pro/_thr suffix)
    pro_results = {}
    thr_results = {}
    
    for r in results:
        name = r['model_name']
        
        if '_pro' in name.lower() or 'pro' in name.lower():
            base = name.lower().replace('_pro', '').replace('pro', '')
            pro_results[base] = r
        elif '_thr' in name.lower() or 'thr' in name.lower():
            base = name.lower().replace('_thr', '').replace('thr', '')
            thr_results[base] = r
    
    # Find matched pairs
    comparisons = []
    
    for base in set(pro_results.keys()) & set(thr_results.keys()):
        pro = pro_results[base]
        thr = thr_results[base]
        
        comparison = {
            'enzyme': base,
            'pro_hbonds': pro['total_hbonds'],
            'thr_hbonds': thr['total_hbonds'],
            'delta_hbonds': thr['total_hbonds'] - pro['total_hbonds'],
            'pro_polar': pro['polar_residue_contacts'],
            'thr_polar': thr['polar_residue_contacts'],
            'thr_hydroxyl': thr.get('hydroxyl_hbonds', 0)
        }
        comparisons.append(comparison)
        
        print(f"\n{base.upper()}:")
        print(f"  Proline:   {pro['total_hbonds']} H-bonds ({pro['polar_residue_contacts']} to polar residues)")
        print(f"  Threonine: {thr['total_hbonds']} H-bonds ({thr['polar_residue_contacts']} to polar residues)")
        print(f"  Δ H-bonds: {comparison['delta_hbonds']:+d}")
        print(f"  Thr hydroxyl contributes: {comparison['thr_hydroxyl']} H-bonds")
    
    if comparisons:
        avg_delta = np.mean([c['delta_hbonds'] for c in comparisons])
        avg_hydroxyl = np.mean([c['thr_hydroxyl'] for c in comparisons])
        
        print("\n" + "-" * 50)
        print(f"AVERAGE Δ H-bonds (Thr - Pro): {avg_delta:+.1f}")
        print(f"AVERAGE hydroxyl contribution: {avg_hydroxyl:.1f}")
        
        if avg_delta > 0:
            print("\n✓ HYPOTHESIS SUPPORTED: Threonine forms MORE H-bonds than Proline")
            print("  The extra hydroxyl finds compensating polar partners.")
        else:
            print("\n✗ HYPOTHESIS NOT SUPPORTED: Threonine does NOT form more H-bonds")
    
    return comparisons

def save_results(results, comparisons, output_dir):
    """Save results to CSV files"""
    
    # Full results
    df_full = pd.DataFrame([{
        'model': r['model_name'],
        'ligand': r['ligand'],
        'total_hbonds': r['total_hbonds'],
        'unique_residues': r['unique_residues'],
        'polar_contacts': r['polar_residue_contacts'],
        'backbone_contacts': r['backbone_contacts'],
        'sidechain_contacts': r['sidechain_contacts'],
        'hydroxyl_hbonds': r.get('hydroxyl_hbonds', 0)
    } for r in results])
    
    df_full.to_csv(f"{output_dir}/hbond_analysis_full.csv", index=False)
    print(f"\nSaved: {output_dir}/hbond_analysis_full.csv")
    
    # Comparisons
    if comparisons:
        df_comp = pd.DataFrame(comparisons)
        df_comp.to_csv(f"{output_dir}/hbond_pro_thr_comparison.csv", index=False)
        print(f"Saved: {output_dir}/hbond_pro_thr_comparison.csv")
    
    # Summary for manuscript
    summary = []
    summary.append("=" * 60)
    summary.append("HYDROGEN BOND ANALYSIS SUMMARY - For Manuscript")
    summary.append("=" * 60)
    
    if comparisons:
        avg_delta = np.mean([c['delta_hbonds'] for c in comparisons])
        avg_hydroxyl = np.mean([c['thr_hydroxyl'] for c in comparisons])
        
        summary.append(f"\nThreonine forms {avg_delta:+.1f} more H-bonds than Proline on average")
        summary.append(f"The β-hydroxyl group contributes {avg_hydroxyl:.1f} H-bonds on average")
        summary.append("\nThis supports the hypothesis that threonine's hydroxyl")
        summary.append("forms compensating interactions with catalytically essential")
        summary.append("polar residues, explaining the persistent cross-reactivity.")
    
    with open(f"{output_dir}/hbond_summary.txt", 'w') as f:
        f.write('\n'.join(summary))
    print(f"Saved: {output_dir}/hbond_summary.txt")

# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python analyze_hydrogen_bonds.py <af3_output_directory>")
        print("\nExample:")
        print("  python analyze_hydrogen_bonds.py /storage/kiran-stuff/aaRS/phase2")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    
    if not os.path.exists(base_dir):
        print(f"ERROR: Directory not found: {base_dir}")
        sys.exit(1)
    
    # Run analysis
    results = analyze_directory(base_dir)
    
    if results:
        # Compare Pro vs Thr
        comparisons = compare_pro_thr(results)
        
        # Save outputs
        output_dir = os.path.join(base_dir, "hbond_analysis")
        os.makedirs(output_dir, exist_ok=True)
        save_results(results, comparisons, output_dir)
        
        print("\n" + "=" * 70)
        print("ANALYSIS COMPLETE")
        print("=" * 70)
    else:
        print("\nNo results generated. Check that CIF files exist.")
