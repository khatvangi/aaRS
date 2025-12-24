# Figure 2: LUCA Editing Domain - Inverted Specificity
# Shows editing domain preferentially binds THR over PRO (opposite of catalytic)
# Run: pymol -c generate_figure2_editing_inverted.pml

# Load editing domain structures
load /storage/kiran-stuff/aaRS/phase2/outputs/deep_editing_pro/deep_editing_pro/deep_editing_pro_model.cif, edit_pro
load /storage/kiran-stuff/aaRS/phase2/outputs/deep_editing_thr/deep_editing_thr/deep_editing_thr_model.cif, edit_thr

# Align structures
align edit_thr and chain A, edit_pro and chain A

# Hide everything
hide everything

# Select protein and ligands
select protein_edit, edit_pro and chain A
select pro_edit_lig, edit_pro and (resn PRO and not polymer)
select thr_edit_lig, edit_thr and (resn THR and not polymer)

# Fallback to chain selection if needed
if pro_edit_lig == 0:
    select pro_edit_lig, edit_pro and chain C
    select thr_edit_lig, edit_thr and chain C

# Show protein surface
show surface, protein_edit
color wheat, protein_edit
set transparency, 0.3, protein_edit

# Show ligands as spheres and sticks
show spheres, pro_edit_lig
show spheres, thr_edit_lig
show sticks, pro_edit_lig
show sticks, thr_edit_lig
color orange, pro_edit_lig
color marine, thr_edit_lig
set stick_radius, 0.3
set sphere_scale, 0.5

# Zoom to editing site
zoom pro_edit_lig, 15

# Publication rendering settings
bg_color white
set ray_shadows, 0
set antialias, 2
set ray_trace_mode, 1
set surface_quality, 2
set depth_cue, 0
set specular, 0.2

# Add labels with ipTM scores
pseudoatom pro_label, pos=[pro_edit_lig and elem C and name CA]
pseudoatom thr_label, pos=[thr_edit_lig and elem C and name CA]
label pro_label, "PRO (ipTM=0.14)"
label thr_label, "THR (ipTM=0.45)"
set label_size, 20
set label_color, black

# Render at high resolution
viewport 2400, 2400
ray 2400, 2400
png /storage/kiran-stuff/aaRS/figure2_editing_inverted_clean.png, dpi=300

# Save session
save /storage/kiran-stuff/aaRS/figure2_session.pse

print "Figure 2 completed: figure2_editing_inverted_clean.png"

quit
