#!/usr/bin/env python3
"""
PyMOL Script: Figure 5C - SER in Zn Site (The Trap!)
Shows SER coordinating Zn just like THR - the zinc trap problem

Status: READY TO RUN (CIF file exists!)
"""
from pymol import cmd

# Load structure
cif_path = "/storage/kiran-stuff/aaRS/phase2/af3_modern_matrix/thrrs_ecoli_zn_jobs/af3_output/modern_thrrs_ecoli_zn_SER/modern_thrrs_ecoli_zn_SER_model.cif"
cmd.load(cif_path, "modern_thrrs_zn_ser")

# Hide everything initially
cmd.hide("everything")

# Show protein as cartoon
cmd.show("cartoon", "modern_thrrs_zn_ser and chain A")
cmd.color("green", "modern_thrrs_zn_ser and chain A")

# Show active site residues
cmd.select("active_site", "modern_thrrs_zn_ser and chain A and resi 500-550")
cmd.show("sticks", "active_site")
cmd.color("wheat", "active_site")

# Show SER ligand
cmd.show("sticks", "modern_thrrs_zn_ser and chain B")
cmd.color("orange", "modern_thrrs_zn_ser and chain B")  # Orange = trap!

# Show Zn ion
cmd.show("spheres", "modern_thrrs_zn_ser and chain C")
cmd.color("gray", "modern_thrrs_zn_ser and chain C")
cmd.set("sphere_scale", 0.5, "modern_thrrs_zn_ser and chain C")

# Measure coordination bonds (should be similar to THR!)
cmd.distance("zn_coord", "modern_thrrs_zn_ser and chain C",
             "modern_thrrs_zn_ser and chain B", 2.8)
cmd.hide("labels", "zn_coord")
cmd.color("red", "zn_coord")
cmd.set("dash_width", 4, "zn_coord")

# Styling
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)

# Zoom to binding site
cmd.zoom("modern_thrrs_zn_ser and chain B or chain C", 8)

# Add labels
cmd.label("modern_thrrs_zn_ser and chain C", '"Zn"')
cmd.label("modern_thrrs_zn_ser and chain B and name CA", '"SER (TRAP!)"')
cmd.set("label_size", 20)
cmd.set("label_color", "black")

# Add annotation box
cmd.select("trap_warning", "modern_thrrs_zn_ser and chain B")

# Render
cmd.ray(2400, 2400)
cmd.png("figures/structural/fig5c_ser_zinc_trap.png", dpi=300)

print("="*80)
print("Saved: figures/structural/fig5c_ser_zinc_trap.png")
print("Figure 5C: SER escapes Zn filter (coordinates just like THR!)")
print("ipTM: 0.950 (97.9% of THR!), Zn_iptm: 0.980")
print("Only 1.02× discrimination - THE ZINC TRAP!")
print("="*80)
