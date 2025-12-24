#!/usr/bin/env python3
"""
PyMOL Visualization: THR Bidentate Coordination
Shows N-CA-OG1 coordination with Zn (3 atoms)
"""

import pymol
from pymol import cmd
import sys

# Initialize PyMOL
pymol.finish_launching(['pymol', '-c'])

# Load structure
cif_path = "/storage/kiran-stuff/aaRS/phase2/af3_modern_matrix/thrrs_ecoli_zn_jobs/af3_output/modern_thrrs_ecoli_zn_THR/modern_thrrs_ecoli_zn_THR_model.cif"

cmd.load(cif_path, "thrrs_thr")

# Styling
cmd.bg_color("white")
cmd.set("ray_opaque_background", 0)
cmd.set("antialias", 2)
cmd.set("ray_trace_mode", 1)

# Hide everything initially
cmd.hide("everything")

# Identify chains
cmd.select("protein", "chain A")
cmd.select("ligand", "chain B")
cmd.select("zn", "chain C")

# Show protein as cartoon (semi-transparent)
cmd.show("cartoon", "protein")
cmd.color("marine", "protein")
cmd.set("cartoon_transparency", 0.7, "protein")

# Show ligand as sticks
cmd.show("sticks", "ligand")
cmd.color("green", "ligand and elem C")
cmd.color("red", "ligand and elem O")
cmd.color("blue", "ligand and elem N")
cmd.set("stick_radius", 0.3, "ligand")

# Show Zn as sphere
cmd.show("spheres", "zn")
cmd.color("gray50", "zn")
cmd.set("sphere_scale", 0.5, "zn")

# Highlight coordinating atoms
cmd.select("coord_atoms", "ligand and (name N or name CA or name OG1)")
cmd.show("spheres", "coord_atoms")
cmd.set("sphere_scale", 0.4, "coord_atoms")
cmd.set("sphere_transparency", 0.3, "coord_atoms")

# Measure coordination distances
cmd.distance("coord_N", "ligand and name N", "zn", cutoff=3.0)
cmd.distance("coord_CA", "ligand and name CA", "zn", cutoff=3.0)
cmd.distance("coord_OG1", "ligand and name OG1", "zn", cutoff=3.0)

# Style distance labels
cmd.set("dash_color", "yellow", "coord_*")
cmd.set("dash_width", 4, "coord_*")
cmd.set("dash_gap", 0.3, "coord_*")
cmd.set("label_size", 30)
cmd.set("label_color", "black")

# Add labels to coordinating atoms
cmd.label("ligand and name N", '"N"')
cmd.label("ligand and name CA", '"CA"')
cmd.label("ligand and name OG1", '"OG1 (-OH)"')
cmd.label("zn", '"Zn²⁺"')

# Zoom to active site
cmd.select("active_site", "ligand or zn or (protein within 8 of ligand)")
cmd.zoom("active_site", buffer=3)

# Set view angle
cmd.orient("ligand")
cmd.rotate("x", -20)
cmd.rotate("y", 30)

# Add title
cmd.set("label_position", (0, 0, 10))
cmd.pseudoatom("title", pos=[0, 0, 0])
cmd.hide("everything", "title")

# Render
output_png = "figures/structural/coordination_thr_bidentate.png"
cmd.png(output_png, width=1200, height=1200, dpi=300, ray=1)

print(f"✓ Generated: {output_png}")
print("THR: Bidentate coordination via N-CA-OG1 (3 atoms)")

cmd.quit()
