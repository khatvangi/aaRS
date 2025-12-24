#!/usr/bin/env python3
"""
PyMOL Script Generator for Structure Rendering
Creates PyMOL commands to render your AF3 structures

Usage:
    python generate_pymol_script.py
    
This creates render_structures.pml that you can run in PyMOL:
    pymol render_structures.pml
"""

import os

BASE_DIR = '/storage/kiran-stuff/aaRS'
PHASE2_OUTPUTS = os.path.join(BASE_DIR, 'phase2/outputs')
OUTPUT_DIR = os.path.join(BASE_DIR, 'manuscript_figures/structures')
PYMOL_SCRIPT = os.path.join(BASE_DIR, 'manuscript_figures/render_structures.pml')

def find_cif_files():
    """Find all CIF structure files from AF3 outputs"""
    structures = {}
    
    for model_dir in os.listdir(PHASE2_OUTPUTS):
        model_path = os.path.join(PHASE2_OUTPUTS, model_dir)
        if not os.path.isdir(model_path):
            continue
        
        # Look for model CIF files
        possible_cifs = [
            os.path.join(model_path, f"{model_dir}_model.cif"),
            os.path.join(model_path, model_dir, f"{model_dir}_model.cif"),
        ]
        
        for cif in possible_cifs:
            if os.path.exists(cif):
                structures[model_dir] = cif
                break
    
    return structures

def generate_pymol_script(structures):
    """Generate PyMOL script for rendering"""
    
    script_lines = [
        "# PyMOL Script for Manuscript Structure Rendering",
        "# Generated automatically",
        "",
        "# Setup",
        "bg_color white",
        "set ray_shadows, 0",
        "set antialias, 2",
        "set hash_max, 300",
        "",
    ]
    
    # Define structures we want for Figure 1 Panel D
    target_structures = {
        'deep_domain_pro': 'LUCA ProRS + PRO (cognate)',
        'deep_domain_thr': 'LUCA ProRS + THR (non-cognate)',
        'shallow_domain_pro': 'Shallow ProThrRS + PRO',
        'shallow_domain_thr': 'Shallow ProThrRS + THR',
        'modern_prours_pro': 'Modern ProRS + PRO',
        'modern_prours_thr': 'Modern ProRS + THR',
    }
    
    for i, (model_name, description) in enumerate(target_structures.items()):
        if model_name not in structures:
            print(f"Warning: {model_name} not found in outputs")
            continue
        
        cif_file = structures[model_name]
        obj_name = model_name
        output_file = os.path.join(OUTPUT_DIR, f"{model_name}_render.png")
        
        script_lines.extend([
            f"# {description}",
            f"load {cif_file}, {obj_name}",
            "",
            f"# Style for {obj_name}",
            f"hide everything, {obj_name}",
            f"show cartoon, {obj_name}",
            f"color marine, {obj_name} and chain A",  # Protein
            "",
            f"# Show ligand as spheres",
            f"select ligand_{obj_name}, {obj_name} and organic",
            f"show spheres, ligand_{obj_name}",
            f"color tv_red, ligand_{obj_name} and elem C",  # Carbon atoms
            f"color tv_blue, ligand_{obj_name} and elem N",  # Nitrogen
            f"color red, ligand_{obj_name} and elem O",      # Oxygen
            "",
            f"# Show binding site residues",
            f"select binding_site_{obj_name}, {obj_name} and (byres (ligand_{obj_name} around 5))",
            f"show sticks, binding_site_{obj_name}",
            f"color yellow, binding_site_{obj_name} and elem C",
            "",
            f"# Center and orient",
            f"zoom ligand_{obj_name}, 8",
            f"orient ligand_{obj_name}",
            "",
            f"# Add label",
            f"set label_color, black",
            f"set label_size, 20",
            "",
            f"# Render",
            f"ray 1200, 1200",
            f"png {output_file}, dpi=300",
            "",
            f"# Clean up for next structure",
            f"delete {obj_name}",
            f"delete ligand_{obj_name}",
            f"delete binding_site_{obj_name}",
            "",
        ])
    
    script_lines.extend([
        "# Done!",
        "print 'Structure rendering complete'",
        "quit",
    ])
    
    return '\n'.join(script_lines)

def generate_manual_pymol_guide():
    """Generate manual PyMOL rendering guide"""
    
    guide = """
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
"""
    
    return guide

def main():
    """Generate PyMOL rendering scripts and guides"""
    print("=" * 60)
    print("PYMOL STRUCTURE RENDERING SCRIPT GENERATOR")
    print("=" * 60)
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Find structures
    print("\nSearching for AF3 structure files...")
    structures = find_cif_files()
    
    if not structures:
        print("\n❌ No structure files found!")
        print(f"Looked in: {PHASE2_OUTPUTS}")
        print("\nMake sure your AF3 outputs are in phase2/outputs/")
        return
    
    print(f"\n✓ Found {len(structures)} structures:")
    for name in sorted(structures.keys())[:10]:
        print(f"  - {name}")
    if len(structures) > 10:
        print(f"  ... and {len(structures) - 10} more")
    
    # Generate PyMOL script
    print("\nGenerating PyMOL script...")
    script = generate_pymol_script(structures)
    
    with open(PYMOL_SCRIPT, 'w') as f:
        f.write(script)
    
    print(f"✓ PyMOL script saved: {PYMOL_SCRIPT}")
    
    # Generate manual guide
    guide_file = os.path.join(os.path.dirname(PYMOL_SCRIPT), 'PyMOL_manual_guide.md')
    guide = generate_manual_pymol_guide()
    
    with open(guide_file, 'w') as f:
        f.write(guide)
    
    print(f"✓ Manual guide saved: {guide_file}")
    
    # Instructions
    print("\n" + "=" * 60)
    print("HOW TO USE")
    print("=" * 60)
    print("\nOption 1 - Automated (if PyMOL command line works):")
    print(f"  pymol -c {PYMOL_SCRIPT}")
    print(f"  # Renders will be saved to: {OUTPUT_DIR}")
    
    print("\nOption 2 - Manual (recommended for publication quality):")
    print(f"  1. Open PyMOL GUI")
    print(f"  2. Follow instructions in: {guide_file}")
    print(f"  3. Manually adjust each view for best presentation")
    print(f"  4. Save PNGs to: {OUTPUT_DIR}")
    
    print("\n" + "=" * 60)
    print("NEXT STEPS")
    print("=" * 60)
    print("1. Render all 6 key structures")
    print("2. Place PNGs in manuscript_figures/structures/")
    print("3. Re-run generate_all_figures.py to insert them into Figure 1")
    
    print("\n✓ Done!")

if __name__ == '__main__':
    main()
