#!/usr/bin/env python3
"""
Publication-Quality PyMOL Figure Generation Script
For ancestral aaRS promiscuity paper

This script generates high-quality molecular figures showing:
- Figure A: Binding pocket comparison (LUCA ProRS + PRO vs THR)
- Figure B: Cognate vs non-cognate overlay (zoomed binding pocket)
- Figure C: Editing domain inverted specificity
- Figure D: Evolutionary comparison (LUCA vs Modern)
"""

import pymol
from pymol import cmd
import os

# Structure file paths
BASE_DIR = "/storage/kiran-stuff/aaRS/phase2/outputs"
OUTPUT_DIR = "/storage/kiran-stuff/structural_figures"

STRUCTURES = {
    'luca_cat_pro': f"{BASE_DIR}/deep_domain_pro/deep_domain_pro_model.cif",
    'luca_cat_thr': f"{BASE_DIR}/deep_domain_thr/deep_domain_thr_model.cif",
    'luca_edit_pro': f"{BASE_DIR}/deep_editing_pro/deep_editing_pro/deep_editing_pro_model.cif",
    'luca_edit_thr': f"{BASE_DIR}/deep_editing_thr/deep_editing_thr/deep_editing_thr_model.cif",
    'modern_cat_pro': f"{BASE_DIR}/shallow_domain_pro/shallow_domain_pro_model.cif",
    'modern_cat_thr': f"{BASE_DIR}/shallow_domain_thr/shallow_domain_thr_model.cif",
}

def setup_pymol_publication_settings():
    """Configure PyMOL for publication-quality rendering"""
    cmd.set('ray_trace_mode', 1)
    cmd.set('antialias', 2)
    cmd.bg_color('white')
    cmd.set('ray_shadows', 0)
    cmd.set('depth_cue', 0)
    cmd.set('specular', 0.2)
    cmd.set('cartoon_fancy_helices', 1)
    cmd.set('sphere_scale', 0.25)
    cmd.set('stick_radius', 0.15)
    cmd.set('cartoon_transparency', 0.5)

def find_ligand_residues(selection):
    """Identify ligand residues (PRO or THR) in a structure"""
    # Look for common ligand residue names
    ligands = []
    cmd.iterate(f"{selection} and resn PRO+THR+PRO+THR and not polymer.protein",
                "ligands.append((resi, resn))", space={'ligands': ligands})

    if not ligands:
        # Try finding by chain or hetatm
        cmd.iterate(f"{selection} and hetatm",
                    "ligands.append((resi, resn))", space={'ligands': ligands})

    return ligands

def figure_a_binding_pocket_comparison():
    """
    FIGURE A: Binding Pocket Comparison
    LUCA ProRS + PRO (green) vs LUCA ProRS + THR (red)
    Shows active site as transparent surface with ligands as sticks
    """
    print("\n" + "="*60)
    print("GENERATING FIGURE A: Binding Pocket Comparison")
    print("="*60)

    cmd.reinitialize()
    setup_pymol_publication_settings()

    # Load structures
    cmd.load(STRUCTURES['luca_cat_pro'], 'luca_pro')
    cmd.load(STRUCTURES['luca_cat_thr'], 'luca_thr')

    # Align structures
    cmd.align('luca_thr', 'luca_pro')

    # Show protein as cartoon
    cmd.show('cartoon', 'all')
    cmd.color('palegreen', 'luca_pro')
    cmd.color('lightblue', 'luca_thr')

    # Find and display ligands
    print("Finding ligands...")
    ligands_pro = find_ligand_residues('luca_pro')
    ligands_thr = find_ligand_residues('luca_thr')

    print(f"PRO ligands found: {ligands_pro}")
    print(f"THR ligands found: {ligands_thr}")

    # If ligands found, show them
    if ligands_pro:
        cmd.select('lig_pro', f'luca_pro and resn PRO and not polymer.protein')
        cmd.show('sticks', 'lig_pro')
        cmd.color('green', 'lig_pro')
        cmd.set('stick_radius', 0.2, 'lig_pro')

    if ligands_thr:
        cmd.select('lig_thr', f'luca_thr and resn THR and not polymer.protein')
        cmd.show('sticks', 'lig_thr')
        cmd.color('red', 'lig_thr')
        cmd.set('stick_radius', 0.2, 'lig_thr')

    # Show active site residues (common binding site residues)
    # These are typically within 5Å of the ligand
    if ligands_pro:
        cmd.select('pocket', 'luca_pro within 5 of lig_pro')
        cmd.show('surface', 'pocket')
        cmd.set('transparency', 0.5, 'pocket')
        cmd.color('palecyan', 'pocket')

    # Orient view
    cmd.orient('lig_pro or lig_thr')
    cmd.zoom('lig_pro or lig_thr', 10)

    # Ray trace and save
    output_file = f"{OUTPUT_DIR}/figure_a_binding_pocket_comparison.png"
    cmd.png(output_file, width=2000, height=2000, dpi=300, ray=1)
    print(f"✓ Figure A saved: {output_file}")

def figure_b_cognate_noncognate_overlay():
    """
    FIGURE B: Cognate vs Non-cognate Overlay
    Zoomed view of binding pocket showing PRO and THR occupying same space
    """
    print("\n" + "="*60)
    print("GENERATING FIGURE B: Cognate vs Non-cognate Overlay")
    print("="*60)

    cmd.reinitialize()
    setup_pymol_publication_settings()

    # Load structures
    cmd.load(STRUCTURES['luca_cat_pro'], 'luca_pro')
    cmd.load(STRUCTURES['luca_cat_thr'], 'luca_thr')

    # Align structures
    cmd.align('luca_thr', 'luca_pro')

    # Hide protein, show only binding pocket residues
    cmd.hide('everything')

    # Find ligands
    ligands_pro = find_ligand_residues('luca_pro')
    ligands_thr = find_ligand_residues('luca_thr')

    if ligands_pro:
        cmd.select('lig_pro', f'luca_pro and resn PRO and not polymer.protein')
        cmd.show('sticks', 'lig_pro')
        cmd.color('green', 'lig_pro')
        cmd.set('stick_radius', 0.25, 'lig_pro')

    if ligands_thr:
        cmd.select('lig_thr', f'luca_thr and resn THR and not polymer.protein')
        cmd.show('sticks', 'lig_thr')
        cmd.color('red', 'lig_thr')
        cmd.set('stick_radius', 0.25, 'lig_thr')

    # Show key binding residues as lines
    if ligands_pro:
        cmd.select('binding_res', 'luca_pro within 4 of lig_pro')
        cmd.show('lines', 'binding_res')
        cmd.color('gray70', 'binding_res')

    # Add distance measurements if both ligands present
    if ligands_pro and ligands_thr:
        # Measure distance between ligand centers
        cmd.distance('dist', 'lig_pro', 'lig_thr')
        cmd.hide('labels', 'dist')

    # Orient and zoom
    cmd.orient('lig_pro or lig_thr')
    cmd.zoom('lig_pro or lig_thr', 3)

    # Ray trace and save
    output_file = f"{OUTPUT_DIR}/figure_b_cognate_noncognate_overlay.png"
    cmd.png(output_file, width=2000, height=2000, dpi=300, ray=1)
    print(f"✓ Figure B saved: {output_file}")

def figure_c_editing_domain_specificity():
    """
    FIGURE C: Editing Domain Inverted Specificity
    Shows THR fitting better than PRO in editing domain
    """
    print("\n" + "="*60)
    print("GENERATING FIGURE C: Editing Domain Inverted Specificity")
    print("="*60)

    cmd.reinitialize()
    setup_pymol_publication_settings()

    # Load editing domain structures
    cmd.load(STRUCTURES['luca_edit_pro'], 'edit_pro')
    cmd.load(STRUCTURES['luca_edit_thr'], 'edit_thr')

    # Align structures
    cmd.align('edit_thr', 'edit_pro')

    # Show protein as cartoon
    cmd.show('cartoon', 'all')
    cmd.color('wheat', 'edit_pro')
    cmd.color('lightblue', 'edit_thr')

    # Find and display ligands
    ligands_pro = find_ligand_residues('edit_pro')
    ligands_thr = find_ligand_residues('edit_thr')

    print(f"Editing domain - PRO ligands: {ligands_pro}")
    print(f"Editing domain - THR ligands: {ligands_thr}")

    if ligands_pro:
        cmd.select('lig_pro', f'edit_pro and resn PRO and not polymer.protein')
        cmd.show('sticks', 'lig_pro')
        cmd.color('orange', 'lig_pro')
        cmd.set('stick_radius', 0.2, 'lig_pro')

    if ligands_thr:
        cmd.select('lig_thr', f'edit_thr and resn THR and not polymer.protein')
        cmd.show('sticks', 'lig_thr')
        cmd.color('marine', 'lig_thr')
        cmd.set('stick_radius', 0.2, 'lig_thr')

    # Show binding pocket
    if ligands_thr:
        cmd.select('edit_pocket', 'edit_thr within 5 of lig_thr')
        cmd.show('surface', 'edit_pocket')
        cmd.set('transparency', 0.6, 'edit_pocket')
        cmd.color('lightpink', 'edit_pocket')

    # Orient view
    cmd.orient('lig_pro or lig_thr')
    cmd.zoom('lig_pro or lig_thr', 10)

    # Ray trace and save
    output_file = f"{OUTPUT_DIR}/figure_c_editing_domain_specificity.png"
    cmd.png(output_file, width=2000, height=2000, dpi=300, ray=1)
    print(f"✓ Figure C saved: {output_file}")

def figure_d_evolutionary_comparison():
    """
    FIGURE D: Evolutionary Comparison
    Side-by-side: LUCA ProRS pocket vs Modern ProRS pocket
    """
    print("\n" + "="*60)
    print("GENERATING FIGURE D: Evolutionary Comparison")
    print("="*60)

    cmd.reinitialize()
    setup_pymol_publication_settings()

    # Load LUCA and Modern structures
    cmd.load(STRUCTURES['luca_cat_pro'], 'luca')
    cmd.load(STRUCTURES['modern_cat_pro'], 'modern')

    # Align structures
    cmd.align('modern', 'luca')

    # Show protein as cartoon
    cmd.show('cartoon', 'all')
    cmd.color('forest', 'luca')
    cmd.color('slate', 'modern')

    # Find ligands
    ligands_luca = find_ligand_residues('luca')
    ligands_modern = find_ligand_residues('modern')

    print(f"LUCA ligands: {ligands_luca}")
    print(f"Modern ligands: {ligands_modern}")

    if ligands_luca:
        cmd.select('lig_luca', f'luca and resn PRO and not polymer.protein')
        cmd.show('sticks', 'lig_luca')
        cmd.color('green', 'lig_luca')
        cmd.set('stick_radius', 0.2, 'lig_luca')

    if ligands_modern:
        cmd.select('lig_modern', f'modern and resn PRO and not polymer.protein')
        cmd.show('sticks', 'lig_modern')
        cmd.color('cyan', 'lig_modern')
        cmd.set('stick_radius', 0.2, 'lig_modern')

    # Show binding pockets with different colors
    if ligands_luca:
        cmd.select('pocket_luca', 'luca within 5 of lig_luca')
        cmd.show('surface', 'pocket_luca')
        cmd.set('transparency', 0.5, 'pocket_luca')
        cmd.color('palegreen', 'pocket_luca')

    if ligands_modern:
        cmd.select('pocket_modern', 'modern within 5 of lig_modern')
        cmd.show('surface', 'pocket_modern')
        cmd.set('transparency', 0.5, 'pocket_modern')
        cmd.color('lightblue', 'pocket_modern')

    # Orient view
    cmd.orient('lig_luca or lig_modern')
    cmd.zoom('lig_luca or lig_modern', 10)

    # Ray trace and save
    output_file = f"{OUTPUT_DIR}/figure_d_evolutionary_comparison.png"
    cmd.png(output_file, width=2000, height=2000, dpi=300, ray=1)
    print(f"✓ Figure D saved: {output_file}")

def main():
    """Generate all publication figures"""
    print("\n" + "#"*60)
    print("# PUBLICATION FIGURE GENERATION FOR ANCESTRAL aaRS PAPER")
    print("#"*60)

    # Verify structure files exist
    print("\nVerifying structure files...")
    missing = []
    for name, path in STRUCTURES.items():
        if os.path.exists(path):
            print(f"✓ {name}: {path}")
        else:
            print(f"✗ {name}: NOT FOUND - {path}")
            missing.append(name)

    if missing:
        print(f"\n⚠ WARNING: {len(missing)} structure file(s) not found!")
        print("Please check the file paths.")
        return

    # Generate all figures
    try:
        figure_a_binding_pocket_comparison()
        figure_b_cognate_noncognate_overlay()
        figure_c_editing_domain_specificity()
        figure_d_evolutionary_comparison()

        print("\n" + "#"*60)
        print("# ALL FIGURES GENERATED SUCCESSFULLY!")
        print("#"*60)
        print(f"\nOutput directory: {OUTPUT_DIR}")
        print("\nGenerated files:")
        for f in os.listdir(OUTPUT_DIR):
            if f.endswith('.png'):
                print(f"  - {f}")

    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    main()
