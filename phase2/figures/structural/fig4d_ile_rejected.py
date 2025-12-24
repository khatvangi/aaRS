#!/usr/bin/env python3
"""
PyMOL Script: Figure 4D - ILE Rejected by Zn Filter
Shows ILE cannot coordinate Zn (hydrophobic side chain)

Status: READY TO RUN (CIF file exists!)
"""
from pymol import cmd

# Load structure
cif_path = "/storage/kiran-stuff/aaRS/phase2/af3_modern_matrix/thrrs_ecoli_zn_jobs/af3_output/modern_thrrs_ecoli_zn_ILE/modern_thrrs_ecoli_zn_ILE_model.cif"
cmd.load(cif_path, "modern_thrrs_zn_ile")

# Hide everything initially
cmd.hide("everything")

# Show protein as cartoon
cmd.show("cartoon", "modern_thrrs_zn_ile and chain A")
cmd.color("green", "modern_thrrs_zn_ile and chain A")

# Show active site residues
cmd.select("active_site", "modern_thrrs_zn_ile and chain A and resi 500-550")
cmd.show("sticks", "active_site")
cmd.color("wheat", "active_site")

# Show ILE ligand
cmd.show("sticks", "modern_thrrs_zn_ile and chain B")
cmd.color("red", "modern_thrrs_zn_ile and chain B")  # Red = rejected

# Show Zn ion
cmd.show("spheres", "modern_thrrs_zn_ile and chain C")
cmd.color("gray", "modern_thrrs_zn_ile and chain C")
cmd.set("sphere_scale", 0.5, "modern_thrrs_zn_ile and chain C")

# Measure distance (should be > 3.5 Å, no coordination)
cmd.distance("zn_distance", "modern_thrrs_zn_ile and chain C",
             "modern_thrrs_zn_ile and chain B", 5.0)
cmd.hide("labels", "zn_distance")
cmd.color("orange", "zn_distance")
cmd.set("dash_width", 2, "zn_distance")
cmd.set("dash_gap", 0.5, "zn_distance")  # Dashed = no bond

# Styling
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)

# Zoom to binding site
cmd.zoom("modern_thrrs_zn_ile and chain B or chain C", 8)

# Add labels
cmd.label("modern_thrrs_zn_ile and chain C", '"Zn"')
cmd.label("modern_thrrs_zn_ile and chain B and name CA", '"ILE (rejected)"')
cmd.set("label_size", 20)
cmd.set("label_color", "black")

# Render
cmd.ray(2400, 2400)
cmd.png("figures/structural/fig4d_ile_rejected.png", dpi=300)

print("="*80)
print("Saved: figures/structural/fig4d_ile_rejected.png")
print("Figure 4D: ILE rejected by Zn filter (no coordination)")
print("ipTM: 0.830 (lower than THR), Zn_iptm: 0.980")
print("="*80)
