#!/usr/bin/env pymol
# PyMOL script to visualize modern E. coli vs LUCA structures
# Shows the evolution from promiscuous to specific binding

# Clear everything
reinitialize

# ===================================================================
# LOAD STRUCTURES
# ===================================================================

# LUCA structures
load /storage/kiran-stuff/aaRS/phase2/af3_output_full/fulllength_deep_pro/fulllength_deep_pro_model.cif, luca_pro
load /storage/kiran-stuff/aaRS/phase2/af3_output_full/fulllength_deep_thr/fulllength_deep_thr_model.cif, luca_thr

# Modern E. coli structures
load /storage/kiran-stuff/aaRS/phase2/af3_gaps/af3_output/modern_ecoli_full_pro/modern_ecoli_full_pro_model.cif, modern_pro
load /storage/kiran-stuff/aaRS/phase2/af3_gaps/af3_output/modern_ecoli_full_thr/modern_ecoli_full_thr_model.cif, modern_thr

# ===================================================================
# ALIGN STRUCTURES
# ===================================================================

# Align all to LUCA ProRS protein chain (chain A)
super luca_thr and chain A, luca_pro and chain A
super modern_pro and chain A, luca_pro and chain A
super modern_thr and chain A, luca_pro and chain A

# ===================================================================
# STYLING - Protein Chains
# ===================================================================

# Hide everything initially
hide everything

# Show cartoon for all proteins
show cartoon, chain A
cartoon automatic

# Color by structure and confidence (b-factor = pLDDT)
# LUCA - Pink tones
spectrum b, red_white_pink, luca_pro and chain A, minimum=50, maximum=90
spectrum b, red_white_pink, luca_thr and chain A, minimum=50, maximum=90

# Modern - Blue tones
spectrum b, blue_white_cyan, modern_pro and chain A, minimum=70, maximum=95
spectrum b, blue_white_cyan, modern_thr and chain A, minimum=70, maximum=95

# ===================================================================
# STYLING - Ligands (Proline and Threonine)
# ===================================================================

# Show ligands as spheres
show spheres, (resn PRO or resn THR) and not chain A
set sphere_scale, 0.4

# Color ligands
color green, resn PRO and not chain A  # Proline = green
color orange, resn THR and not chain A  # Threonine = orange

# ===================================================================
# IDENTIFY AND SHOW BINDING POCKET RESIDUES
# ===================================================================

# Select residues within 5Å of ligands
select pocket_luca_pro, (luca_pro and chain A) within 5 of (luca_pro and resn PRO)
select pocket_luca_thr, (luca_thr and chain A) within 5 of (luca_thr and resn THR)
select pocket_modern_pro, (modern_pro and chain A) within 5 of (modern_pro and resn PRO)
select pocket_modern_thr, (modern_thr and chain A) within 5 of (modern_thr and resn THR)

# Show pocket residues as sticks
show sticks, pocket_luca_pro
show sticks, pocket_luca_thr
show sticks, pocket_modern_pro
show sticks, pocket_modern_thr

# ===================================================================
# CREATE SCENES FOR DIFFERENT VIEWS
# ===================================================================

# Scene 1: LUCA ProRS with both Pro and Thr overlaid
disable all
enable luca_pro
zoom luca_pro and chain A
orient (luca_pro and resn PRO)
scene scene_1_luca_promiscuous, store

# Scene 2: Modern E. coli ProRS (high specificity)
disable all
enable modern_pro
zoom modern_pro and chain A
orient (modern_pro and resn PRO)
scene scene_2_modern_specific, store

# Scene 3: LUCA vs Modern side-by-side overlay
disable all
enable luca_pro
enable modern_pro
zoom chain A
scene scene_3_overlay, store

# Scene 4: Focus on binding pocket comparison
disable all
enable luca_pro
enable modern_pro
zoom (pocket_luca_pro or pocket_modern_pro)
scene scene_4_pocket_zoom, store

# ===================================================================
# RENDER HIGH-QUALITY IMAGES
# ===================================================================

# Set rendering parameters
bg_color white
set ray_trace_mode, 1
set ray_shadows, 1
set ambient, 0.4
set specular, 0.5
set depth_cue, 0

# Output directory
cd /storage/kiran-stuff/aaRS/phase2/pymol_renders

# Render Scene 1: LUCA promiscuous binding
scene scene_1_luca_promiscuous
ray 2400, 2400
png figure4a_luca_promiscuous_pocket.png, dpi=300

# Render Scene 2: Modern specific binding
scene scene_2_modern_specific
ray 2400, 2400
png figure4b_modern_specific_pocket.png, dpi=300

# Render Scene 3: Overlay comparison
scene scene_3_overlay
ray 2400, 2400
png figure4c_luca_vs_modern_overlay.png, dpi=300

# Render Scene 4: Pocket zoom
scene scene_4_pocket_zoom
ray 2400, 2400
png figure4d_pocket_comparison_zoom.png, dpi=300

# ===================================================================
# SAVE SESSION
# ===================================================================

save /storage/kiran-stuff/aaRS/phase2/modern_vs_luca_analysis.pse

print "=============================================="
print "PyMOL rendering complete!"
print "Scenes created:"
print "  1. LUCA promiscuous (Pro + Thr binding)"
print "  2. Modern specific (optimized Pro binding)"
print "  3. LUCA vs Modern overlay"
print "  4. Pocket comparison zoom"
print ""
print "Files saved to: phase2/pymol_renders/"
print "Session saved: modern_vs_luca_analysis.pse"
print "=============================================="

# Exit (comment out if running interactively)
# quit
