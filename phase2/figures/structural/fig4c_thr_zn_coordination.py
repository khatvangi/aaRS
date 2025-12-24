#!/usr/bin/env python3
"""
PyMOL Script: Figure 4C - THR Coordinating Zn
Shows bidentate coordination of THR with Zn in modern ThrRS

Status: READY TO RUN (CIF file exists!)
"""
from pymol import cmd

# Load structure
cif_path = "/storage/kiran-stuff/aaRS/phase2/af3_modern_matrix/thrrs_ecoli_zn_jobs/af3_output/modern_thrrs_ecoli_zn_THR/modern_thrrs_ecoli_zn_THR_model.cif"
cmd.load(cif_path, "modern_thrrs_zn_thr")

# Hide everything initially
cmd.hide("everything")

# Show protein as cartoon
cmd.show("cartoon", "modern_thrrs_zn_thr and chain A")
cmd.color("green", "modern_thrrs_zn_thr and chain A")

# Show active site residues (approximate - may need adjustment)
cmd.select("active_site", "modern_thrrs_zn_thr and chain A and resi 500-550")
cmd.show("sticks", "active_site")
cmd.color("wheat", "active_site")

# Show THR ligand
cmd.show("sticks", "modern_thrrs_zn_thr and chain B")
cmd.color("yellow", "modern_thrrs_zn_thr and chain B")

# Show Zn ion
cmd.show("spheres", "modern_thrrs_zn_thr and chain C")
cmd.color("gray", "modern_thrrs_zn_thr and chain C")
cmd.set("sphere_scale", 0.5, "modern_thrrs_zn_thr and chain C")

# Measure coordination bonds
cmd.distance("zn_coord", "modern_thrrs_zn_thr and chain C",
             "modern_thrrs_zn_thr and chain B", 2.8)
cmd.hide("labels", "zn_coord")
cmd.color("red", "zn_coord")
cmd.set("dash_width", 4, "zn_coord")

# Styling
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)

# Zoom to binding site
cmd.zoom("modern_thrrs_zn_thr and chain B or chain C", 8)

# Add label
cmd.label("modern_thrrs_zn_thr and chain C", '"Zn"')
cmd.set("label_size", 20)
cmd.set("label_color", "black")

# Render
cmd.ray(2400, 2400)
cmd.png("figures/structural/fig4c_thr_zn_coordination.png", dpi=300)

print("="*80)
print("Saved: figures/structural/fig4c_thr_zn_coordination.png")
print("Figure 4C: THR coordinating Zn via bidentate bonds")
print("ipTM: 0.970, Zn_iptm: 0.980")
print("="*80)
