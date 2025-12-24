#!/usr/bin/env python3
"""
generate_figure5_receiver_first.py
==================================
PyMOL script to generate Figure 5: The "Receiver-First" pattern

Creates side-by-side comparison of:
- Panel A: Modern E. coli ProRS (pocket ipTM 0.95, global ipTM 0.95) - rigid/blue
- Panel B: LUCA ProRS (pocket ipTM 0.78, global ipTM 0.28) - pocket blue, tRNA orange

Run with: pymol -cq generate_figure5_receiver_first.py
Or interactively: run generate_figure5_receiver_first.py
"""

from pymol import cmd, util
import os

# =============================================================================
# CONFIGURATION - UPDATE THESE PATHS
# =============================================================================

# Modern E. coli ProRS + Pro (the "gold standard")
MODERN_CIF = "/storage/kiran-stuff/aaRS/phase2/af3_gaps/af3_output/modern_ecoli_full_pro/modern_ecoli_full_pro_model.cif"

# LUCA ProRS + Pro (ancestral with uncoupled pocket/global)
LUCA_CIF = "/storage/kiran-stuff/aaRS/phase2/af3_output_full/fulllength_deep_pro/fulllength_deep_pro_model.cif"

# Output directory
OUTPUT_DIR = "/storage/kiran-stuff/aaRS/figures"

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def setup_visualization():
    """Set up PyMOL visualization parameters"""
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)
    cmd.set("ray_shadows", 0)

def color_by_plddt(obj_name):
    """Color structure by pLDDT (b-factor) using AF3 color scheme"""
    # AF3 pLDDT color scheme:
    # >90: blue (high confidence)
    # 70-90: cyan
    # 50-70: yellow (low confidence)
    # <50: orange (very low)
    cmd.spectrum("b", "orange_yellow_cyan_blue", obj_name, minimum=50, maximum=90)

def highlight_ligand(obj_name, ligand_resn="PRO"):
    """Show ligand as green spheres"""
    lig_sel = f"{obj_name} and resn {ligand_resn}"
    cmd.show("spheres", lig_sel)
    cmd.color("green", lig_sel)
    cmd.set("sphere_scale", 0.8, lig_sel)

def hide_tRNA(obj_name):
    """Hide tRNA chain (usually chain B) to focus on protein-ligand"""
    # Try common tRNA chain IDs
    for chain in ["B", "C"]:
        sel = f"{obj_name} and chain {chain} and (resn A or resn U or resn G or resn C)"
        if cmd.count_atoms(sel) > 50:  # Likely tRNA
            cmd.hide("everything", f"{obj_name} and chain {chain}")
            print(f"  Hidden tRNA on chain {chain}")
            break

# =============================================================================
# MAIN FIGURE GENERATION
# =============================================================================

def generate_figure5():
    """Generate the Receiver-First comparison figure"""

    print("=" * 60)
    print("Generating Figure 5: Receiver-First Pattern")
    print("=" * 60)

    # Clear everything
    cmd.reinitialize()
    setup_visualization()

    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # =========================================================================
    # Panel A: Modern E. coli ProRS (Rigid/Coupled)
    # =========================================================================
    print("\n[Panel A] Loading Modern E. coli ProRS...")

    if os.path.exists(MODERN_CIF):
        cmd.load(MODERN_CIF, "modern")
        cmd.show("cartoon", "modern")
        color_by_plddt("modern")
        highlight_ligand("modern", "PRO")

        # Optional: hide tRNA to focus on protein
        # hide_tRNA("modern")

        # Orient and save
        cmd.orient("modern")
        cmd.zoom("modern", buffer=5)

        # Save Panel A
        cmd.ray(2400, 2400)
        cmd.png(f"{OUTPUT_DIR}/figure5a_modern_ecoli.png", dpi=300)
        print(f"  Saved: figure5a_modern_ecoli.png")

        # Also save session
        cmd.save(f"{OUTPUT_DIR}/figure5a_modern.pse")
    else:
        print(f"  WARNING: Modern CIF not found at {MODERN_CIF}")
        print("  Please update MODERN_CIF path in script")

    # =========================================================================
    # Panel B: LUCA ProRS (Pocket rigid, Global flexible)
    # =========================================================================
    print("\n[Panel B] Loading LUCA ProRS...")

    cmd.reinitialize()
    setup_visualization()

    if os.path.exists(LUCA_CIF):
        cmd.load(LUCA_CIF, "luca")
        cmd.show("cartoon", "luca")
        color_by_plddt("luca")
        highlight_ligand("luca", "PRO")

        # Orient and save
        cmd.orient("luca")
        cmd.zoom("luca", buffer=5)

        # Save Panel B
        cmd.ray(2400, 2400)
        cmd.png(f"{OUTPUT_DIR}/figure5b_luca.png", dpi=300)
        print(f"  Saved: figure5b_luca.png")

        cmd.save(f"{OUTPUT_DIR}/figure5b_luca.pse")
    else:
        print(f"  WARNING: LUCA CIF not found at {LUCA_CIF}")
        print("  Please update LUCA_CIF path in script")

    # =========================================================================
    # Combined view (both structures)
    # =========================================================================
    print("\n[Combined] Creating side-by-side view...")

    cmd.reinitialize()
    setup_visualization()

    if os.path.exists(MODERN_CIF) and os.path.exists(LUCA_CIF):
        cmd.load(MODERN_CIF, "modern")
        cmd.load(LUCA_CIF, "luca")

        # Show as cartoon
        cmd.show("cartoon", "all")

        # Color by pLDDT
        color_by_plddt("modern")
        color_by_plddt("luca")

        # Highlight ligands
        highlight_ligand("modern", "PRO")
        highlight_ligand("luca", "PRO")

        # Translate LUCA to the right for side-by-side
        cmd.translate([150, 0, 0], "luca")

        # Orient and save
        cmd.zoom("all", buffer=10)
        cmd.ray(3600, 1800)
        cmd.png(f"{OUTPUT_DIR}/figure5_combined.png", dpi=300)
        print(f"  Saved: figure5_combined.png")

        cmd.save(f"{OUTPUT_DIR}/figure5_combined.pse")

    print("\n" + "=" * 60)
    print("Figure 5 generation complete!")
    print("=" * 60)
    print(f"\nOutput files in: {OUTPUT_DIR}")
    print("  - figure5a_modern_ecoli.png (Panel A)")
    print("  - figure5b_luca.png (Panel B)")
    print("  - figure5_combined.png (Side-by-side)")
    print("  - *.pse files (PyMOL sessions for editing)")

# =============================================================================
# ACTIVE SITE ZOOM VIEW
# =============================================================================

def generate_active_site_zoom():
    """Generate zoomed views of the active sites"""

    print("\n" + "=" * 60)
    print("Generating Active Site Zoom Views")
    print("=" * 60)

    cmd.reinitialize()
    setup_visualization()

    if os.path.exists(MODERN_CIF):
        cmd.load(MODERN_CIF, "modern")
        cmd.show("cartoon", "modern")
        color_by_plddt("modern")

        # Find and highlight ligand
        cmd.show("sticks", "modern and resn PRO")
        cmd.color("green", "modern and resn PRO")

        # Show residues within 5A of ligand
        cmd.select("pocket", "modern and byres (resn PRO around 5)")
        cmd.show("sticks", "pocket")

        # Zoom to pocket
        cmd.zoom("pocket", buffer=3)
        cmd.orient("pocket")

        # Show surface of pocket
        cmd.show("surface", "pocket")
        cmd.set("transparency", 0.5, "pocket")

        cmd.ray(2400, 2400)
        cmd.png(f"{OUTPUT_DIR}/figure5a_pocket_zoom.png", dpi=300)
        print(f"  Saved: figure5a_pocket_zoom.png")

    # Same for LUCA
    cmd.reinitialize()
    setup_visualization()

    if os.path.exists(LUCA_CIF):
        cmd.load(LUCA_CIF, "luca")
        cmd.show("cartoon", "luca")
        color_by_plddt("luca")

        cmd.show("sticks", "luca and resn PRO")
        cmd.color("green", "luca and resn PRO")

        cmd.select("pocket", "luca and byres (resn PRO around 5)")
        cmd.show("sticks", "pocket")
        cmd.zoom("pocket", buffer=3)

        cmd.show("surface", "pocket")
        cmd.set("transparency", 0.5, "pocket")

        cmd.ray(2400, 2400)
        cmd.png(f"{OUTPUT_DIR}/figure5b_pocket_zoom.png", dpi=300)
        print(f"  Saved: figure5b_pocket_zoom.png")

# =============================================================================
# RUN
# =============================================================================

if __name__ == "__main__":
    generate_figure5()
    generate_active_site_zoom()
    print("\nDone! Review the PNG files and PSE sessions.")
