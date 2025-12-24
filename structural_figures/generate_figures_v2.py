#!/usr/bin/env python3
"""
Publication-Quality PyMOL Figure Generation Script - Version 2
For ancestral aaRS promiscuity paper

Updated to correctly identify ligands in Chain C
"""

import pymol
from pymol import cmd
import os
import sys

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

def figure_a_binding_pocket_comparison():
    """
    FIGURE A: Binding Pocket Comparison
    LUCA ProRS + PRO (green) vs LUCA ProRS + THR (red)
    """
    print("\n" + "="*60)
    print("GENERATING FIGURE A: Binding Pocket Comparison")
    print("="*60)

    cmd.reinitialize()
    setup_pymol_publication_settings()

    # Load structures
    print("Loading structures...")
    cmd.load(STRUCTURES['luca_cat_pro'], 'luca_pro')
    cmd.load(STRUCTURES['luca_cat_thr'], 'luca_thr')

    # Align protein chains
    print("Aligning structures...")
    cmd.align('luca_thr and chain A', 'luca_pro and chain A')

    # Show protein as cartoon (Chain A only)
    cmd.hide('everything')
    cmd.show('cartoon', 'chain A')
    cmd.color('palegreen', 'luca_pro and chain A')
    cmd.color('lightblue', 'luca_thr and chain A')

    # Show ligands (Chain C) as sticks
    print("Displaying ligands...")
    cmd.select('lig_pro', 'luca_pro and chain C')
    cmd.show('sticks', 'lig_pro')
    cmd.color('green', 'lig_pro')
    cmd.set('stick_radius', 0.25, 'lig_pro')

    cmd.select('lig_thr', 'luca_thr and chain C')
    cmd.show('sticks', 'lig_thr')
    cmd.color('red', 'lig_thr')
    cmd.set('stick_radius', 0.25, 'lig_thr')

    # Show binding pocket residues
    print("Showing binding pocket...")
    cmd.select('pocket', 'luca_pro and chain A within 5 of lig_pro')
    cmd.show('surface', 'pocket')
    cmd.set('transparency', 0.5, 'pocket')
    cmd.color('palecyan', 'pocket')

    # Orient view
    print("Orienting view...")
    cmd.orient('lig_pro or lig_thr')
    cmd.zoom('lig_pro or lig_thr', 8)

    # Ray trace and save
    print("Rendering figure (this may take a minute)...")
    output_file = f"{OUTPUT_DIR}/figure_a_binding_pocket.png"
    cmd.png(output_file, width=2000, height=2000, dpi=300, ray=1)
    print(f"✓ Figure A saved: {output_file}")

def figure_b_cognate_noncognate_overlay():
    """
    FIGURE B: Cognate vs Non-cognate Overlay
    Zoomed view of binding pocket
    """
    print("\n" + "="*60)
    print("GENERATING FIGURE B: Cognate vs Non-cognate Overlay")
    print("="*60)

    cmd.reinitialize()
    setup_pymol_publication_settings()

    # Load structures
    print("Loading structures...")
    cmd.load(STRUCTURES['luca_cat_pro'], 'luca_pro')
    cmd.load(STRUCTURES['luca_cat_thr'], 'luca_thr')

    # Align
    cmd.align('luca_thr and chain A', 'luca_pro and chain A')

    # Hide everything, then show only binding region
    cmd.hide('everything')

    # Show ligands as sticks
    cmd.select('lig_pro', 'luca_pro and chain C')
    cmd.show('sticks', 'lig_pro')
    cmd.color('green', 'lig_pro')
    cmd.set('stick_radius', 0.3, 'lig_pro')

    cmd.select('lig_thr', 'luca_thr and chain C')
    cmd.show('sticks', 'lig_thr')
    cmd.color('red', 'lig_thr')
    cmd.set('stick_radius', 0.3, 'lig_thr')

    # Show key binding residues as lines
    cmd.select('binding_res', 'luca_pro and chain A within 4 of lig_pro')
    cmd.show('lines', 'binding_res')
    cmd.color('gray70', 'binding_res')

    # Labels for ligands
    cmd.label('lig_pro and name CA', '"PRO"')
    cmd.label('lig_thr and name CA', '"THR"')

    # Orient and zoom close
    cmd.orient('lig_pro or lig_thr')
    cmd.zoom('lig_pro or lig_thr', 2)

    # Ray trace and save
    print("Rendering figure...")
    output_file = f"{OUTPUT_DIR}/figure_b_overlay.png"
    cmd.png(output_file, width=2000, height=2000, dpi=300, ray=1)
    print(f"✓ Figure B saved: {output_file}")

def figure_c_editing_domain_specificity():
    """
    FIGURE C: Editing Domain Inverted Specificity
    """
    print("\n" + "="*60)
    print("GENERATING FIGURE C: Editing Domain Inverted Specificity")
    print("="*60)

    cmd.reinitialize()
    setup_pymol_publication_settings()

    # Load editing domain structures
    print("Loading structures...")
    cmd.load(STRUCTURES['luca_edit_pro'], 'edit_pro')
    cmd.load(STRUCTURES['luca_edit_thr'], 'edit_thr')

    # Align
    cmd.align('edit_thr and chain A', 'edit_pro and chain A')

    # Show protein as cartoon
    cmd.hide('everything')
    cmd.show('cartoon', 'chain A')
    cmd.color('wheat', 'edit_pro and chain A')
    cmd.color('lightblue', 'edit_thr and chain A')

    # Show ligands
    cmd.select('lig_pro', 'edit_pro and chain C')
    cmd.show('sticks', 'lig_pro')
    cmd.color('orange', 'lig_pro')
    cmd.set('stick_radius', 0.25, 'lig_pro')

    cmd.select('lig_thr', 'edit_thr and chain C')
    cmd.show('sticks', 'lig_thr')
    cmd.color('marine', 'lig_thr')
    cmd.set('stick_radius', 0.25, 'lig_thr')

    # Show binding pocket
    cmd.select('edit_pocket', 'edit_thr and chain A within 5 of lig_thr')
    cmd.show('surface', 'edit_pocket')
    cmd.set('transparency', 0.6, 'edit_pocket')
    cmd.color('lightpink', 'edit_pocket')

    # Orient view
    cmd.orient('lig_pro or lig_thr')
    cmd.zoom('lig_pro or lig_thr', 8)

    # Ray trace and save
    print("Rendering figure...")
    output_file = f"{OUTPUT_DIR}/figure_c_editing_domain.png"
    cmd.png(output_file, width=2000, height=2000, dpi=300, ray=1)
    print(f"✓ Figure C saved: {output_file}")

def figure_d_evolutionary_comparison():
    """
    FIGURE D: Evolutionary Comparison
    LUCA vs Modern ProRS binding pocket
    """
    print("\n" + "="*60)
    print("GENERATING FIGURE D: Evolutionary Comparison")
    print("="*60)

    cmd.reinitialize()
    setup_pymol_publication_settings()

    # Load LUCA and Modern structures
    print("Loading structures...")
    cmd.load(STRUCTURES['luca_cat_pro'], 'luca')
    cmd.load(STRUCTURES['modern_cat_pro'], 'modern')

    # Align
    cmd.align('modern and chain A', 'luca and chain A')

    # Show protein as cartoon
    cmd.hide('everything')
    cmd.show('cartoon', 'chain A')
    cmd.color('forest', 'luca and chain A')
    cmd.color('slate', 'modern and chain A')

    # Show ligands
    cmd.select('lig_luca', 'luca and chain C')
    cmd.show('sticks', 'lig_luca')
    cmd.color('green', 'lig_luca')
    cmd.set('stick_radius', 0.25, 'lig_luca')

    cmd.select('lig_modern', 'modern and chain C')
    cmd.show('sticks', 'lig_modern')
    cmd.color('cyan', 'lig_modern')
    cmd.set('stick_radius', 0.25, 'lig_modern')

    # Show binding pockets with different colors
    cmd.select('pocket_luca', 'luca and chain A within 5 of lig_luca')
    cmd.show('surface', 'pocket_luca')
    cmd.set('transparency', 0.5, 'pocket_luca')
    cmd.color('palegreen', 'pocket_luca')

    cmd.select('pocket_modern', 'modern and chain A within 5 of lig_modern')
    cmd.show('surface', 'pocket_modern')
    cmd.set('transparency', 0.5, 'pocket_modern')
    cmd.color('lightblue', 'pocket_modern')

    # Orient view
    cmd.orient('lig_luca or lig_modern')
    cmd.zoom('lig_luca or lig_modern', 8)

    # Ray trace and save
    print("Rendering figure...")
    output_file = f"{OUTPUT_DIR}/figure_d_evolutionary.png"
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
            print(f"✓ {name}: Found")
        else:
            print(f"✗ {name}: NOT FOUND - {path}")
            missing.append(name)

    if missing:
        print(f"\n⚠ WARNING: {len(missing)} structure file(s) not found!")
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
        for f in sorted(os.listdir(OUTPUT_DIR)):
            if f.endswith('.png'):
                fsize = os.path.getsize(os.path.join(OUTPUT_DIR, f)) / (1024*1024)
                print(f"  - {f} ({fsize:.1f} MB)")

    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
