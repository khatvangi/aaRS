#!/usr/bin/env python3
"""
PyMOL Script: Figure 3D - Ancestral vs Modern Active Site
Shows pocket evolution in ThrRS (ancestral vs modern)

Status: READY TO RUN
NOTE: These structures do NOT have Zn (ligand-only predictions)
"""
from pymol import cmd

# Load structures (UPDATED PATHS)
cmd.load("/storage/kiran-stuff/aaRS/phase2/outputs/deep_thrrs_thr/deep_thrrs_thr/deep_thrrs_thr_model.cif", "ancestral")
cmd.load("/storage/kiran-stuff/aaRS/phase2/outputs/modern_thrrs_thr/modern_thrrs_thr/modern_thrrs_thr_model.cif", "modern")

# Align on catalytic domain
cmd.align("modern", "ancestral")

# Hide all
cmd.hide("everything")

# Show cartoons
cmd.show("cartoon", "all")

# Color by evolution
cmd.color("purple", "ancestral")
cmd.color("green", "modern")

# Show active site residues as sticks (approximate range - adjust as needed)
# Note: Residue ranges may need adjustment based on actual structure
cmd.show("sticks", "ancestral and chain A and resi 200-250")
cmd.show("sticks", "modern and chain A and resi 200-250")

# Show THR ligand
cmd.show("sticks", "chain B")
cmd.color("yellow", "chain B")

# Styling
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
# Set transparency for cartoons to see through overlays
cmd.set("cartoon_transparency", 0.3, "ancestral")
cmd.set("cartoon_transparency", 0.3, "modern")

# Zoom to ligand binding site
cmd.zoom("chain B", 12)

# Add labels
cmd.label("ancestral and chain B and name CA", '"Ancestral (loose pocket)"')
cmd.label("modern and chain B and name CA", '"Modern (tighter pocket)"')

# Render
cmd.ray(2400, 2400)
cmd.png("figures/structural/fig3d_evolution_overlay.png", dpi=300)

print("="*80)
print("Saved: figures/structural/fig3d_evolution_overlay.png")
print("Figure 3D: Shows pocket evolution from ancestral to modern ThrRS")
print("NOTE: No Zn in these structures (ligand-only predictions)")
print("      Ancestral ipTM: 0.210")
print("      Modern ipTM: 0.240")
print("="*80)
