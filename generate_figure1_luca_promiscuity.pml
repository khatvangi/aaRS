# Figure 1: LUCA ProRS Promiscuity - Protein surface + PRO/THR overlay
# Shows LUCA ProRS catalytic domain can bind both proline and threonine
# Run: pymol -c generate_figure1_luca_promiscuity.pml

# Load structures
load /storage/kiran-stuff/aaRS/phase2/outputs/deep_domain_pro/deep_domain_pro_model.cif, luca_pro
load /storage/kiran-stuff/aaRS/phase2/outputs/deep_domain_thr/deep_domain_thr_model.cif, luca_thr

# Align structures on protein chain
align luca_thr and chain A, luca_pro and chain A

# Hide everything first
hide everything

# Select protein (chain A)
select protein_pro, luca_pro and chain A
select protein_thr, luca_thr and chain A

# Select amino acid ligands
# In AF3 outputs, ligands are typically in separate chains or as HETATM records
# Try multiple selection strategies
select proline_lig, luca_pro and (resn PRO and not polymer)
select threonine_lig, luca_thr and (resn THR and not polymer)

# If above doesn't work, try selecting by chain (typically chain C for ligands)
if proline_lig == 0:
    select proline_lig, luca_pro and chain C
    select threonine_lig, luca_thr and chain C

# Show protein as surface (from PRO structure - they're aligned)
show surface, protein_pro
color gray80, protein_pro
set transparency, 0.3, protein_pro

# Show amino acids as spheres
show spheres, proline_lig
show spheres, threonine_lig
color green, proline_lig
color firebrick, threonine_lig

# Also show as sticks for better visibility
show sticks, proline_lig
show sticks, threonine_lig
set stick_radius, 0.3

# Zoom to binding pocket region (15 Angstrom around ligands)
zoom proline_lig, 15

# Publication-quality rendering settings
bg_color white
set ray_shadows, 0
set antialias, 2
set ray_trace_mode, 1
set surface_quality, 2
set sphere_scale, 0.5
set depth_cue, 0
set specular, 0.2

# Set view angle for best pocket visualization
set_view (\
     0.891877115,    0.365063906,   -0.266199440,\
     0.139847383,   -0.778408408,   -0.611991167,\
    -0.430310041,    0.509316742,   -0.744924605,\
     0.000000000,    0.000000000, -150.000000000,\
    25.123456001,   30.234567001,   15.345678001,\
   100.000000000,  200.000000000,  -20.000000000 )

# Add labels
pseudoatom pro_label, pos=[proline_lig and elem C and name CA]
pseudoatom thr_label, pos=[threonine_lig and elem C and name CA]
label pro_label, "PRO (cognate)"
label thr_label, "THR (non-cognate)"
set label_size, 20
set label_color, black

# Render at high resolution
viewport 2400, 2400
ray 2400, 2400
png /storage/kiran-stuff/aaRS/figure1_luca_promiscuity_clean.png, dpi=300

# Save session for later editing
save /storage/kiran-stuff/aaRS/figure1_session.pse

print "Figure 1 completed: figure1_luca_promiscuity_clean.png"

quit
