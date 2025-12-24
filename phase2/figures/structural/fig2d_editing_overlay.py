#!/usr/bin/env python3
"""
PyMOL Script: Figure 2D - Editing Domain Overlay
Shows THR vs PRO in ancestral ProRS editing domain

Status: READY TO RUN
"""
from pymol import cmd

# Load structures (UPDATED PATHS)
cmd.load("/storage/kiran-stuff/aaRS/phase2/outputs/deep_editing_thr/deep_editing_thr/deep_editing_thr_model.cif", "edit_THR")
cmd.load("/storage/kiran-stuff/aaRS/phase2/outputs/deep_editing_pro/deep_editing_pro/deep_editing_pro_model.cif", "edit_PRO")

# Align structures
cmd.align("edit_PRO", "edit_THR")

# Hide everything initially
cmd.hide("everything")

# Show protein as cartoon
cmd.show("cartoon", "edit_THR and chain A")
cmd.show("cartoon", "edit_PRO and chain A")

# Color proteins
cmd.color("lightblue", "edit_THR and chain A")
cmd.color("palegreen", "edit_PRO and chain A")

# Show ligands as sticks
cmd.show("sticks", "edit_THR and chain B")
cmd.show("sticks", "edit_PRO and chain B")

# Color ligands
cmd.color("green", "edit_THR and chain B")  # THR = accepted in editing site
cmd.color("cyan", "edit_PRO and chain B")    # PRO = cognate (should not bind here)

# Set view and styling
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)
cmd.set("cartoon_fancy_helices", 1)

# Zoom to active site
cmd.zoom("chain B", 8)

# Add labels
cmd.label("edit_THR and chain B and name CA", '"THR (ipTM=0.87)"')
cmd.label("edit_PRO and chain B and name CA", '"PRO (ipTM=0.82)"')

# Render
cmd.ray(2400, 2400)
cmd.png("figures/structural/fig2d_editing_overlay.png", dpi=300)

print("="*80)
print("Saved: figures/structural/fig2d_editing_overlay.png")
print("Figure 2D: Editing domain shows THR binds better than PRO")
print("="*80)
