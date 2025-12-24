# PyMOL script for aaRS structural analysis
# Generated for publication-quality figures

# Initialize
reinitialize
bg_color white
set ray_opaque_background, on
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, 0
set specular, 0.2
set orthoscopic, on

# Color scheme
set_color luca_green, [0.2, 0.7, 0.3]
set_color luca_red, [0.9, 0.2, 0.2]
set_color modern_blue, [0.3, 0.5, 0.8]
set_color orange_edit, [1.0, 0.6, 0.0]
set_color marine_edit, [0.0, 0.5, 0.7]

print("=" * 60)
print("TASK 1: LUCA CATALYTIC DOMAIN - PRO vs THR SUPERIMPOSITION")
print("=" * 60)

# Load LUCA structures
load /storage/kiran-stuff/aaRS/phase2/outputs/deep_domain_pro/deep_domain_pro_model.cif, luca_pro
load /storage/kiran-stuff/aaRS/phase2/outputs/deep_domain_thr/deep_domain_thr_model.cif, luca_thr

# Align protein chains (chain A is protein, chain B+ are ligands)
super luca_thr and chain A, luca_pro and chain A
print("Alignment RMSD:")
rms_ca luca_thr and chain A and name CA, luca_pro and chain A and name CA

# Create combined object for visualization
create luca_combined, luca_pro or luca_thr

# Identify ligands (non-protein chains)
# PRO ligand (chain B in luca_pro)
select pro_ligand, luca_pro and not chain A
# THR ligand (chain B in luca_thr)
select thr_ligand, luca_thr and not chain A

# Calculate ligand RMSD (critical measurement!)
print("\nLigand RMSD (PRO vs THR position):")
rms pro_ligand, thr_ligand

# Measure distance between ligand centers
print("\nDistance between ligand centers:")
distance lig_center_dist, pro_ligand and name CA, thr_ligand and name CA

# FIGURE 1: LUCA Promiscuity - Both ligands in pocket
print("\nGenerating FIGURE 1: LUCA ProRS Promiscuity...")
hide everything
show surface, luca_pro and chain A
set transparency, 0.3
color gray80, luca_pro and chain A

# Show ligands as spheres
show spheres, pro_ligand
show spheres, thr_ligand
color luca_green, pro_ligand
color luca_red, thr_ligand
set sphere_scale, 0.5

# Find binding pocket residues (within 5A of ligands)
select pocket_residues, (luca_pro and chain A) within 5 of (pro_ligand or thr_ligand)
show sticks, pocket_residues
color gray50, pocket_residues
set stick_radius, 0.15

# Set view
orient pocket_residues
zoom pocket_residues, 8

# Ray trace and save
ray 2000, 2000
png /storage/kiran-stuff/aaRS/structural_figures/v2/figure1_luca_promiscuity_overlay.png

print("  Saved: figure1_luca_promiscuity_overlay.png")

# FIGURE 2: Zoomed Ligand Overlay
print("\nGenerating FIGURE 2: Zoomed Ligand Overlay...")
hide everything
show sticks, pro_ligand or thr_ligand
set stick_radius, 0.3
color luca_green, pro_ligand
color luca_red, thr_ligand

# Show nearby residues
select nearby_residues, (luca_pro and chain A) within 4 of (pro_ligand or thr_ligand)
show lines, nearby_residues
color gray70, nearby_residues

# Transparent surface
show surface, nearby_residues
set transparency, 0.7
color gray90, nearby_residues

# Add distance label
distance ligand_rmsd, pro_ligand, thr_ligand
hide labels, ligand_rmsd

orient pro_ligand or thr_ligand
zoom pro_ligand or thr_ligand, 5

ray 2000, 2000
png /storage/kiran-stuff/aaRS/structural_figures/v2/figure2_ligand_overlay_zoom.png
print("  Saved: figure2_ligand_overlay_zoom.png")

# Calculate pocket volume for LUCA
print("\n" + "=" * 60)
print("POCKET VOLUME CALCULATION - LUCA ProRS")
print("=" * 60)

# Method: Calculate volume of pocket residues
select luca_pocket, (luca_pro and chain A) within 5 of pro_ligand
print("Pocket residues (5Å from PRO ligand):")
iterate luca_pocket, print(f"{resi} {resn}")

# Export pocket for external volume calculation
save /storage/kiran-stuff/aaRS/structural_figures/v2/luca_pocket.pdb, luca_pocket or pro_ligand

print("=" * 60)
print("TASK 2: MODERN ProRS STRUCTURES")
print("=" * 60)

# Load modern structures
delete all
load /storage/kiran-stuff/aaRS/phase2/outputs/shallow_domain_pro/shallow_domain_pro_model.cif, modern_pro
load /storage/kiran-stuff/aaRS/phase2/outputs/shallow_domain_thr/shallow_domain_thr_model.cif, modern_thr

# Align
super modern_thr and chain A, modern_pro and chain A
print("Modern alignment RMSD:")
rms_ca modern_thr and chain A and name CA, modern_pro and chain A and name CA

# Identify ligands
select modern_pro_lig, modern_pro and not chain A
select modern_thr_lig, modern_thr and not chain A

# Calculate modern pocket
select modern_pocket, (modern_pro and chain A) within 5 of modern_pro_lig
save /storage/kiran-stuff/aaRS/structural_figures/v2/modern_pocket.pdb, modern_pocket or modern_pro_lig

print("=" * 60)
print("TASK 3: LUCA vs MODERN COMPARISON")
print("=" * 60)

# Load both for comparison
delete all
load /storage/kiran-stuff/aaRS/phase2/outputs/deep_domain_pro/deep_domain_pro_model.cif, luca_pro
load /storage/kiran-stuff/aaRS/phase2/outputs/shallow_domain_pro/shallow_domain_pro_model.cif, modern_pro

# Align the two structures
super modern_pro and chain A, luca_pro and chain A
print("LUCA vs Modern alignment RMSD:")
rms_ca modern_pro and chain A and name CA, luca_pro and chain A and name CA

# FIGURE 4: LUCA vs Modern pocket comparison
print("\nGenerating FIGURE 4: LUCA vs Modern Pocket Comparison...")

# Select pockets
select luca_pkt, (luca_pro and chain A) within 5 of (luca_pro and not chain A)
select modern_pkt, (modern_pro and chain A) within 5 of (modern_pro and not chain A)

# Show surfaces
hide everything
show surface, luca_pkt
show surface, modern_pkt
set transparency, 0.4
color luca_green, luca_pkt
color modern_blue, modern_pkt

# Show ligands
select luca_lig, luca_pro and not chain A
select modern_lig, modern_pro and not chain A
show spheres, luca_lig or modern_lig
set sphere_scale, 0.4
color tv_green, luca_lig
color tv_blue, modern_lig

orient luca_pkt or modern_pkt
zoom luca_pkt or modern_pkt, 8

ray 2000, 2000
png /storage/kiran-stuff/aaRS/structural_figures/v2/figure4_luca_vs_modern_pocket.png
print("  Saved: figure4_luca_vs_modern_pocket.png")

# Export for volume calculation
save /storage/kiran-stuff/aaRS/structural_figures/v2/luca_modern_comparison.pdb, luca_pkt or modern_pkt or luca_lig or modern_lig

print("=" * 60)
print("TASK 4: EDITING DOMAIN ANALYSIS")
print("=" * 60)

delete all
load /storage/kiran-stuff/aaRS/phase2/outputs/deep_editing_pro/deep_editing_pro/seed-1_sample-0/deep_editing_pro_seed-1_sample-0_model.cif, editing_pro
load /storage/kiran-stuff/aaRS/phase2/outputs/deep_editing_thr/deep_editing_thr/seed-1_sample-0/deep_editing_thr_seed-1_sample-0_model.cif, editing_thr

# Align editing domains
super editing_thr and chain A, editing_pro and chain A
print("Editing domain alignment RMSD:")
rms_ca editing_thr and chain A and name CA, editing_pro and chain A and name CA

# Identify ligands
select edit_pro_lig, editing_pro and not chain A
select edit_thr_lig, editing_thr and not chain A

print("\nEditing domain ligand RMSD:")
rms edit_pro_lig, edit_thr_lig

# FIGURE 3: Editing domain inverted specificity
print("\nGenerating FIGURE 3: Editing Domain Inverted Specificity...")
hide everything
show surface, editing_pro and chain A
set transparency, 0.3
color wheat, editing_pro and chain A

# Show ligands
show spheres, edit_pro_lig or edit_thr_lig
set sphere_scale, 0.5
color orange_edit, edit_pro_lig
color marine_edit, edit_thr_lig

# Show pocket residues
select edit_pocket, (editing_pro and chain A) within 5 of (edit_pro_lig or edit_thr_lig)
show sticks, edit_pocket
color gray60, edit_pocket
set stick_radius, 0.15

orient edit_pocket
zoom edit_pocket, 8

ray 2000, 2000
png /storage/kiran-stuff/aaRS/structural_figures/v2/figure3_editing_domain_inverted.png
print("  Saved: figure3_editing_domain_inverted.png")

print("=" * 60)
print("PyMOL ANALYSIS COMPLETE")
print("=" * 60)
print(f"\nStructural figures saved to: /storage/kiran-stuff/aaRS/structural_figures/v2")
print("\nNext: Run pocket volume calculations with fpocket or manual analysis")

# Save session for manual inspection
save /storage/kiran-stuff/aaRS/structural_figures/v2/structural_analysis.pse

quit
