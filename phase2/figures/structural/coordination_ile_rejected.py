#!/usr/bin/env python3
"""
PyMOL Visualization: ILE Monodentate Coordination (REJECTED)
Shows N-CA coordination with Zn (2 atoms, NO -OH)
"""

import pymol
from pymol import cmd
import sys

# Initialize PyMOL
pymol.finish_launching(['pymol', '-c'])

# Load structure
cif_path = "/storage/kiran-stuff/aaRS/phase2/af3_modern_matrix/thrrs_ecoli_zn_jobs/af3_output/modern_thrrs_ecoli_zn_ILE/modern_thrrs_ecoli_zn_ILE_model.cif"

cmd.load(cif_path, "thrrs_ile")

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
cmd.color("red", "ligand and elem C")  # Red for rejected
cmd.color("red", "ligand and elem O")
cmd.color("blue", "ligand and elem N")
cmd.set("stick_radius", 0.3, "ligand")

# Show Zn as sphere
cmd.show("spheres", "zn")
cmd.color("gray50", "zn")
cmd.set("sphere_scale", 0.5, "zn")

# Highlight coordinating atoms (only N and CA for ILE - no -OH)
cmd.select("coord_atoms", "ligand and (name N or name CA)")
cmd.show("spheres", "coord_atoms")
cmd.set("sphere_scale", 0.4, "coord_atoms")
cmd.set("sphere_transparency", 0.3, "coord_atoms")

# Measure coordination distances
cmd.distance("coord_N", "ligand and name N", "zn", cutoff=3.0)
cmd.distance("coord_CA", "ligand and name CA", "zn", cutoff=3.0)

# Style distance labels
cmd.set("dash_color", "red", "coord_*")
cmd.set("dash_width", 4, "coord_*")
cmd.set("dash_gap", 0.3, "coord_*")
cmd.set("label_size", 30)
cmd.set("label_color", "black")

# Add labels to coordinating atoms
cmd.label("ligand and name N", '"N"')
cmd.label("ligand and name CA", '"CA"')
cmd.label("zn", '"Zn²⁺"')

# Show where -OH would be (highlight the hydrophobic side chain)
cmd.select("sidechain", "ligand and (name CB or name CG1 or name CG2 or name CD1)")
cmd.show("sticks", "sidechain")
cmd.color("gray70", "sidechain")
cmd.label("ligand and name CG1", '"Hydrophobic\\n(NO -OH)"')

# Zoom to active site
cmd.select("active_site", "ligand or zn or (protein within 8 of ligand)")
cmd.zoom("active_site", buffer=3)

# Set view angle
cmd.orient("ligand")
cmd.rotate("x", -20)
cmd.rotate("y", 30)

# Render
output_png = "figures/structural/coordination_ile_rejected.png"
cmd.png(output_png, width=1200, height=1200, dpi=300, ray=1)

print(f"✓ Generated: {output_png}")
print("ILE: Monodentate coordination via N-CA (2 atoms, NO -OH)")
print("REJECTED: Lacks third coordination point → ipTM 0.83")

cmd.quit()
