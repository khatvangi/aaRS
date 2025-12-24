
# Manual PyMOL Rendering Guide

If the automated script doesn't work, here's how to render manually:

## 1. Load Structure
```
load /path/to/deep_domain_pro_model.cif, pro_cognate
```

## 2. Style the Protein
```
hide everything
show cartoon, pro_cognate
color marine, pro_cognate and chain A
```

## 3. Show the Ligand
```
select ligand, pro_cognate and organic
show spheres, ligand
color tv_red, ligand and elem C
color tv_blue, ligand and elem N
color red, ligand and elem O
```

## 4. Show Binding Site
```
select binding_site, pro_cognate and (byres (ligand around 5))
show sticks, binding_site
color yellow, binding_site and elem C
```

## 5. Orient and Center
```
zoom ligand, 8
orient ligand
```

## 6. Adjust View
Use mouse to rotate to best angle showing:
- Ligand in binding pocket
- Key interacting residues
- Clear view of catalytic site

## 7. Render High Quality
```
bg_color white
set ray_shadows, 0
set antialias, 2
ray 1200, 1200
png /path/to/output/pro_cognate_render.png, dpi=300
```

## 8. Repeat for Each Structure

Do this for all 6 key structures:
1. deep_domain_pro (LUCA ProRS + PRO)
2. deep_domain_thr (LUCA ProRS + THR)
3. shallow_domain_pro (Shallow + PRO)
4. shallow_domain_thr (Shallow + THR)
5. modern_prours_pro (Modern + PRO)
6. modern_prours_thr (Modern + THR)

## Tips for Publication Quality

1. **Consistent orientation**: Use same view angle for comparable structures
2. **Clear labeling**: Make sure ligand and key residues are visible
3. **Resolution**: Always use 1200x1200 at 300 DPI minimum
4. **Background**: White background for publications
5. **Colors**: Keep consistent across all panels

## Alternative: Use Sessions

Save your styled view as a session:
```
save my_view.pse
```

Then load others and apply same settings:
```
load new_structure.cif
@my_view.pse  # Apply saved view
ray 1200, 1200
png output.png
```
