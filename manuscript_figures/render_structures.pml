# PyMOL Script for Manuscript Structure Rendering
# Generated automatically

# Setup
bg_color white
set ray_shadows, 0
set antialias, 2
set hash_max, 300

# LUCA ProRS + PRO (cognate)
load /storage/kiran-stuff/aaRS/phase2/outputs/deep_domain_pro/deep_domain_pro_model.cif, deep_domain_pro

# Style for deep_domain_pro
hide everything, deep_domain_pro
show cartoon, deep_domain_pro
color marine, deep_domain_pro and chain A

# Show ligand as spheres
select ligand_deep_domain_pro, deep_domain_pro and organic
show spheres, ligand_deep_domain_pro
color tv_red, ligand_deep_domain_pro and elem C
color tv_blue, ligand_deep_domain_pro and elem N
color red, ligand_deep_domain_pro and elem O

# Show binding site residues
select binding_site_deep_domain_pro, deep_domain_pro and (byres (ligand_deep_domain_pro around 5))
show sticks, binding_site_deep_domain_pro
color yellow, binding_site_deep_domain_pro and elem C

# Center and orient
zoom ligand_deep_domain_pro, 8
orient ligand_deep_domain_pro

# Add label
set label_color, black
set label_size, 20

# Render
ray 1200, 1200
png /storage/kiran-stuff/aaRS/manuscript_figures/structures/deep_domain_pro_render.png, dpi=300

# Clean up for next structure
delete deep_domain_pro
delete ligand_deep_domain_pro
delete binding_site_deep_domain_pro

# LUCA ProRS + THR (non-cognate)
load /storage/kiran-stuff/aaRS/phase2/outputs/deep_domain_thr/deep_domain_thr_model.cif, deep_domain_thr

# Style for deep_domain_thr
hide everything, deep_domain_thr
show cartoon, deep_domain_thr
color marine, deep_domain_thr and chain A

# Show ligand as spheres
select ligand_deep_domain_thr, deep_domain_thr and organic
show spheres, ligand_deep_domain_thr
color tv_red, ligand_deep_domain_thr and elem C
color tv_blue, ligand_deep_domain_thr and elem N
color red, ligand_deep_domain_thr and elem O

# Show binding site residues
select binding_site_deep_domain_thr, deep_domain_thr and (byres (ligand_deep_domain_thr around 5))
show sticks, binding_site_deep_domain_thr
color yellow, binding_site_deep_domain_thr and elem C

# Center and orient
zoom ligand_deep_domain_thr, 8
orient ligand_deep_domain_thr

# Add label
set label_color, black
set label_size, 20

# Render
ray 1200, 1200
png /storage/kiran-stuff/aaRS/manuscript_figures/structures/deep_domain_thr_render.png, dpi=300

# Clean up for next structure
delete deep_domain_thr
delete ligand_deep_domain_thr
delete binding_site_deep_domain_thr

# Shallow ProThrRS + PRO
load /storage/kiran-stuff/aaRS/phase2/outputs/shallow_domain_pro/shallow_domain_pro_model.cif, shallow_domain_pro

# Style for shallow_domain_pro
hide everything, shallow_domain_pro
show cartoon, shallow_domain_pro
color marine, shallow_domain_pro and chain A

# Show ligand as spheres
select ligand_shallow_domain_pro, shallow_domain_pro and organic
show spheres, ligand_shallow_domain_pro
color tv_red, ligand_shallow_domain_pro and elem C
color tv_blue, ligand_shallow_domain_pro and elem N
color red, ligand_shallow_domain_pro and elem O

# Show binding site residues
select binding_site_shallow_domain_pro, shallow_domain_pro and (byres (ligand_shallow_domain_pro around 5))
show sticks, binding_site_shallow_domain_pro
color yellow, binding_site_shallow_domain_pro and elem C

# Center and orient
zoom ligand_shallow_domain_pro, 8
orient ligand_shallow_domain_pro

# Add label
set label_color, black
set label_size, 20

# Render
ray 1200, 1200
png /storage/kiran-stuff/aaRS/manuscript_figures/structures/shallow_domain_pro_render.png, dpi=300

# Clean up for next structure
delete shallow_domain_pro
delete ligand_shallow_domain_pro
delete binding_site_shallow_domain_pro

# Shallow ProThrRS + THR
load /storage/kiran-stuff/aaRS/phase2/outputs/shallow_domain_thr/shallow_domain_thr_model.cif, shallow_domain_thr

# Style for shallow_domain_thr
hide everything, shallow_domain_thr
show cartoon, shallow_domain_thr
color marine, shallow_domain_thr and chain A

# Show ligand as spheres
select ligand_shallow_domain_thr, shallow_domain_thr and organic
show spheres, ligand_shallow_domain_thr
color tv_red, ligand_shallow_domain_thr and elem C
color tv_blue, ligand_shallow_domain_thr and elem N
color red, ligand_shallow_domain_thr and elem O

# Show binding site residues
select binding_site_shallow_domain_thr, shallow_domain_thr and (byres (ligand_shallow_domain_thr around 5))
show sticks, binding_site_shallow_domain_thr
color yellow, binding_site_shallow_domain_thr and elem C

# Center and orient
zoom ligand_shallow_domain_thr, 8
orient ligand_shallow_domain_thr

# Add label
set label_color, black
set label_size, 20

# Render
ray 1200, 1200
png /storage/kiran-stuff/aaRS/manuscript_figures/structures/shallow_domain_thr_render.png, dpi=300

# Clean up for next structure
delete shallow_domain_thr
delete ligand_shallow_domain_thr
delete binding_site_shallow_domain_thr

# Modern ProRS + PRO
load /storage/kiran-stuff/aaRS/phase2/outputs/modern_prours_pro/modern_prours_pro/modern_prours_pro_model.cif, modern_prours_pro

# Style for modern_prours_pro
hide everything, modern_prours_pro
show cartoon, modern_prours_pro
color marine, modern_prours_pro and chain A

# Show ligand as spheres
select ligand_modern_prours_pro, modern_prours_pro and organic
show spheres, ligand_modern_prours_pro
color tv_red, ligand_modern_prours_pro and elem C
color tv_blue, ligand_modern_prours_pro and elem N
color red, ligand_modern_prours_pro and elem O

# Show binding site residues
select binding_site_modern_prours_pro, modern_prours_pro and (byres (ligand_modern_prours_pro around 5))
show sticks, binding_site_modern_prours_pro
color yellow, binding_site_modern_prours_pro and elem C

# Center and orient
zoom ligand_modern_prours_pro, 8
orient ligand_modern_prours_pro

# Add label
set label_color, black
set label_size, 20

# Render
ray 1200, 1200
png /storage/kiran-stuff/aaRS/manuscript_figures/structures/modern_prours_pro_render.png, dpi=300

# Clean up for next structure
delete modern_prours_pro
delete ligand_modern_prours_pro
delete binding_site_modern_prours_pro

# Modern ProRS + THR
load /storage/kiran-stuff/aaRS/phase2/outputs/modern_prours_thr/modern_prours_thr/modern_prours_thr_model.cif, modern_prours_thr

# Style for modern_prours_thr
hide everything, modern_prours_thr
show cartoon, modern_prours_thr
color marine, modern_prours_thr and chain A

# Show ligand as spheres
select ligand_modern_prours_thr, modern_prours_thr and organic
show spheres, ligand_modern_prours_thr
color tv_red, ligand_modern_prours_thr and elem C
color tv_blue, ligand_modern_prours_thr and elem N
color red, ligand_modern_prours_thr and elem O

# Show binding site residues
select binding_site_modern_prours_thr, modern_prours_thr and (byres (ligand_modern_prours_thr around 5))
show sticks, binding_site_modern_prours_thr
color yellow, binding_site_modern_prours_thr and elem C

# Center and orient
zoom ligand_modern_prours_thr, 8
orient ligand_modern_prours_thr

# Add label
set label_color, black
set label_size, 20

# Render
ray 1200, 1200
png /storage/kiran-stuff/aaRS/manuscript_figures/structures/modern_prours_thr_render.png, dpi=300

# Clean up for next structure
delete modern_prours_thr
delete ligand_modern_prours_thr
delete binding_site_modern_prours_thr

# Done!
print 'Structure rendering complete'
quit