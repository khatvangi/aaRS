#!/usr/bin/env python3
"""
Fixed Figure Generation - Correct Protein Alignment
Fixes Figures B and C to properly superimpose ligands in the binding pocket
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

def figure_b_cognate_noncognate_overlay_FIXED():
    """
    FIGURE B (FIXED): Cognate vs Non-cognate Overlay
    Properly aligned proteins so both ligands appear in the SAME pocket
    """
    print("\n" + "="*60)
    print("GENERATING FIGURE B (FIXED): Cognate vs Non-cognate Overlay")
    print("="*60)

    cmd.reinitialize()
    setup_pymol_publication_settings()

    # Load structures
    print("Loading structures...")
    cmd.load(STRUCTURES['luca_cat_pro'], 'luca_pro')
    cmd.load(STRUCTURES['luca_cat_thr'], 'luca_thr')

    # CRITICAL: Align the PROTEINS first (chain A only)
    # This ensures both ligands will be in the SAME binding pocket
    print("Aligning proteins (chain A)...")
    cmd.align('luca_thr and chain A', 'luca_pro and chain A')
    print("Alignment complete - ligands should now be in same pocket")

    # Hide everything first
    cmd.hide('everything')

    # Show ligands as sticks (Chain C contains the ligands)
    print("Displaying ligands...")
    cmd.select('lig_pro', 'luca_pro and chain C')
    cmd.show('sticks', 'lig_pro')
    cmd.color('green', 'lig_pro')
    cmd.set('stick_radius', 0.3, 'lig_pro')

    cmd.select('lig_thr', 'luca_thr and chain C')
    cmd.show('sticks', 'lig_thr')
    cmd.color('red', 'lig_thr')
    cmd.set('stick_radius', 0.3, 'lig_thr')

    # Show binding pocket residues within 5Å of EITHER ligand
    print("Showing binding pocket residues...")
    cmd.select('pocket_res', '(luca_pro and chain A within 5 of (lig_pro or lig_thr))')
    cmd.show('sticks', 'pocket_res')
    cmd.color('wheat', 'pocket_res')
    cmd.set('stick_radius', 0.15, 'pocket_res')

    # Add labels to ligands
    cmd.label('lig_pro and name CA', '"PRO"')
    cmd.label('lig_thr and name CA', '"THR"')
    cmd.set('label_size', 20)
    cmd.set('label_color', 'black')

    # Orient and zoom to the binding pocket
    print("Orienting view...")
    cmd.orient('lig_pro or lig_thr')
    cmd.zoom('lig_pro or lig_thr', 3)

    # Ray trace and save
    print("Rendering figure (this may take a minute)...")
    output_file = f"{OUTPUT_DIR}/figure_b_overlay_FIXED.png"
    cmd.png(output_file, width=2000, height=2000, dpi=300, ray=1)
    print(f"✓ Figure B (FIXED) saved: {output_file}")

def figure_c_editing_domain_specificity_FIXED():
    """
    FIGURE C (FIXED): Editing Domain Inverted Specificity
    Properly aligned proteins so both ligands appear in the SAME pocket
    """
    print("\n" + "="*60)
    print("GENERATING FIGURE C (FIXED): Editing Domain Inverted Specificity")
    print("="*60)

    cmd.reinitialize()
    setup_pymol_publication_settings()

    # Load editing domain structures
    print("Loading structures...")
    cmd.load(STRUCTURES['luca_edit_pro'], 'edit_pro')
    cmd.load(STRUCTURES['luca_edit_thr'], 'edit_thr')

    # CRITICAL: Align the PROTEINS first (chain A only)
    print("Aligning proteins (chain A)...")
    cmd.align('edit_thr and chain A', 'edit_pro and chain A')
    print("Alignment complete - ligands should now be in same pocket")

    # Hide everything first
    cmd.hide('everything')

    # Show protein as cartoon
    cmd.show('cartoon', 'chain A')
    cmd.color('wheat', 'edit_pro and chain A')
    cmd.color('lightblue', 'edit_thr and chain A')
    cmd.set('cartoon_transparency', 0.3)

    # Show ligands
    print("Displaying ligands...")
    cmd.select('lig_pro', 'edit_pro and chain C')
    cmd.show('sticks', 'lig_pro')
    cmd.color('orange', 'lig_pro')
    cmd.set('stick_radius', 0.3, 'lig_pro')

    cmd.select('lig_thr', 'edit_thr and chain C')
    cmd.show('sticks', 'lig_thr')
    cmd.color('marine', 'lig_thr')
    cmd.set('stick_radius', 0.3, 'lig_thr')

    # Show binding pocket surface (around BOTH ligands)
    print("Showing binding pocket surface...")
    cmd.select('edit_pocket', 'edit_thr and chain A within 5 of (lig_pro or lig_thr)')
    cmd.show('surface', 'edit_pocket')
    cmd.set('transparency', 0.5, 'edit_pocket')
    cmd.color('lightpink', 'edit_pocket')

    # Show key binding residues as sticks
    cmd.select('pocket_sticks', 'edit_thr and chain A within 4 of (lig_pro or lig_thr)')
    cmd.show('sticks', 'pocket_sticks')
    cmd.color('gray70', 'pocket_sticks')
    cmd.set('stick_radius', 0.15, 'pocket_sticks')

    # Add labels
    cmd.label('lig_pro and name CA', '"PRO"')
    cmd.label('lig_thr and name CA', '"THR"')
    cmd.set('label_size', 20)
    cmd.set('label_color', 'black')

    # Orient view
    print("Orienting view...")
    cmd.orient('lig_pro or lig_thr')
    cmd.zoom('lig_pro or lig_thr', 5)

    # Ray trace and save
    print("Rendering figure (this may take a minute)...")
    output_file = f"{OUTPUT_DIR}/figure_c_editing_domain_FIXED.png"
    cmd.png(output_file, width=2000, height=2000, dpi=300, ray=1)
    print(f"✓ Figure C (FIXED) saved: {output_file}")

def main():
    """Generate fixed figures B and C"""
    print("\n" + "#"*60)
    print("# FIXING FIGURES B AND C - CORRECT ALIGNMENT")
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

    # Generate fixed figures
    try:
        figure_b_cognate_noncognate_overlay_FIXED()
        figure_c_editing_domain_specificity_FIXED()

        print("\n" + "#"*60)
        print("# FIXED FIGURES GENERATED SUCCESSFULLY!")
        print("#"*60)
        print(f"\nOutput directory: {OUTPUT_DIR}")
        print("\nFixed files:")
        for f in sorted(os.listdir(OUTPUT_DIR)):
            if 'FIXED' in f and f.endswith('.png'):
                fsize = os.path.getsize(os.path.join(OUTPUT_DIR, f)) / (1024*1024)
                print(f"  - {f} ({fsize:.1f} MB)")

        print("\n" + "="*60)
        print("KEY FIX APPLIED:")
        print("- Aligned PROTEINS first (chain A to chain A)")
        print("- Ligands now appear in the SAME binding pocket")
        print("- Both PRO and THR are spatially overlapping as expected")
        print("="*60)

    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    main()
