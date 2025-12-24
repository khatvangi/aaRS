# Figure 3: LUCA vs Modern ProRS - Contact Paradox
# Side-by-side comparison showing more contacts but worse specificity
# Run: pymol -c generate_figure3_luca_vs_modern.pml

# Load LUCA and Modern ProRS structures (with proline)
load /storage/kiran-stuff/aaRS/phase2/outputs/deep_domain_pro/deep_domain_pro_model.cif, luca
load /storage/kiran-stuff/aaRS/phase2/outputs/shallow_domain_pro/shallow_domain_pro_model.cif, modern

# Align modern to LUCA (on protein chain A)
align modern and chain A, luca and chain A

# Hide everything
hide everything

# Select proteins
select prot_luca, luca and chain A
select prot_modern, modern and chain A

# Select ligands (proline in both)
select lig_luca, luca and (resn PRO and not polymer)
select lig_modern, modern and (resn PRO and not polymer)

# Fallback to chain C if needed
if lig_luca == 0:
    select lig_luca, luca and chain C
    select lig_modern, modern and chain C

# Color scheme: LUCA = light gray, Modern = light blue
show surface, prot_luca
show surface, prot_modern
color gray80, prot_luca
color lightblue, prot_modern
set transparency, 0.3

# Show ligands as sticks and spheres
show spheres, lig_luca
show spheres, lig_modern
show sticks, lig_luca
show sticks, lig_modern
color green, lig_luca
color green, lig_modern
set stick_radius, 0.3
set sphere_scale, 0.5

# For side-by-side: translate modern structure 50 Angstroms
translate [50, 0, 0], modern

# Zoom to show both structures
zoom all
center

# Publication rendering
bg_color white
set ray_shadows, 0
set antialias, 2
set ray_trace_mode, 1
set surface_quality, 2
set depth_cue, 0
set specular, 0.2

# Add labels
pseudoatom luca_label, pos=[lig_luca and elem C and name CA]
pseudoatom modern_label, pos=[lig_modern and elem C and name CA]
label luca_label, "LUCA ProRS (3.5 Gya)"
label modern_label, "Modern ProRS"
set label_size, 24
set label_color, black
set label_bg_color, white
set label_bg_transparency, 0.3

# Adjust view to show both side by side
turn y, 20

# Render at high resolution (wider for side-by-side)
viewport 3600, 1800
ray 3600, 1800
png /storage/kiran-stuff/aaRS/figure3_luca_vs_modern_clean.png, dpi=300

# Save session
save /storage/kiran-stuff/aaRS/figure3_session.pse

print "Figure 3 completed: figure3_luca_vs_modern_clean.png"

quit
