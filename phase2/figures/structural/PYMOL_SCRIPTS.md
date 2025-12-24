# PyMOL Scripts for Structural Renders
## aaRS Evolution - Publication Quality Structures

**Generated:** 2025-12-18
**Target:** 5 structural panels for Figures 2-5

---

## 🔬 Prerequisites

### Software
- **PyMOL** (version 2.0+)
  - Open source: https://pymol.org
  - Commercial: https://pymol.org/2/
  - Alternative: **ChimeraX** (https://www.cgl.ucsf.edu/chimerax/)

### CIF Files Needed

Based on available outputs, we have these structures:
- `deep_domain_pro` - Ancestral ProRS catalytic
- `deep_domain_thr` - Ancestral ThrRS catalytic
- `deep_editing_pro` - Ancestral ProRS editing
- `modern_thrrs_pro` - Modern ThrRS with PRO
- `modern_thrrs_thr` - Modern ThrRS with THR
- `modern_prours_pro` - Modern ProRS with PRO
- `modern_prours_thr` - Modern ProRS with THR

**Note:** Full matrix with all 20 AAs + Zn may require additional AF3 runs.

---

## 📜 PyMOL Script 1: Figure 2D - Editing Domain Overlay

**Purpose:** Show THR vs PRO in ProRS editing domain

### Files Needed
- `deep_editing_pro/deep_editing_pro_model.cif` (or equivalent with THR)
- Need separate runs with THR and PRO ligands

### PyMOL Script

```python
# fig2d_editing_overlay.py
from pymol import cmd

# Load structures
cmd.load("path/to/editing_THR.cif", "edit_THR")
cmd.load("path/to/editing_PRO.cif", "edit_PRO")

# Align structures
cmd.align("edit_PRO", "edit_THR")

# Hide everything initially
cmd.hide("everything")

# Show protein as cartoon
cmd.show("cartoon", "edit_THR and chain A")
cmd.show("cartoon", "edit_PRO and chain A")

# Color proteins
cmd.color("lightblue", "edit_THR and chain A")
cmd.color("palegreen", "edit_PRO and chain A")

# Show ligands as sticks
cmd.show("sticks", "edit_THR and chain B")
cmd.show("sticks", "edit_PRO and chain B")

# Color ligands
cmd.color("green", "edit_THR and chain B")  # THR = target
cmd.color("cyan", "edit_PRO and chain B")    # PRO = excluded

# Set view and styling
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)
cmd.set("cartoon_fancy_helices", 1)

# Zoom to active site
cmd.zoom("chain B", 8)

# Add labels
cmd.label("edit_THR and chain B and name CA", '"THR (binds better)"')
cmd.label("edit_PRO and chain B and name CA", '"PRO (excluded)"')

# Render
cmd.ray(2400, 2400)
cmd.png("figures/structural/fig2d_editing_overlay.png", dpi=300)

print("Saved: fig2d_editing_overlay.png")
```

### Expected Output
- Overlay showing both ligands
- THR in green (binds better, 0.87)
- PRO in cyan (excluded, 0.82)
- Protein surface/cartoon showing binding pocket

---

## 📜 PyMOL Script 2: Figure 3D - Ancestral vs Modern Active Site

**Purpose:** Show pocket tightening around Zn during evolution

### Files Needed
- `deep_domain_thr/deep_domain_thr_model.cif` (ancestral)
- `modern_thrrs_thr/modern_thrrs_thr_model.cif` (modern)

### PyMOL Script

```python
# fig3d_evolution_overlay.py
from pymol import cmd

# Load structures
cmd.load("path/to/anc_thrrs_THR.cif", "ancestral")
cmd.load("path/to/mod_thrrs_THR.cif", "modern")

# Align on catalytic domain
cmd.align("modern", "ancestral")

# Hide all
cmd.hide("everything")

# Show cartoons
cmd.show("cartoon", "all")

# Color by evolution
cmd.color("purple", "ancestral")
cmd.color("green", "modern")

# Show active site residues as sticks
cmd.show("sticks", "ancestral and (resi 100-150)")  # Adjust residue range
cmd.show("sticks", "modern and (resi 100-150)")

# Show Zn if present
cmd.show("sphere", "name ZN")
cmd.set("sphere_scale", 0.5, "name ZN")
cmd.color("gray", "name ZN")

# Show THR ligand
cmd.show("sticks", "chain B")
cmd.color("yellow", "chain B")

# Styling
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("transparency", 0.3, "cartoon")

# Zoom to active site
cmd.zoom("name ZN", 12)

# Render
cmd.ray(2400, 2400)
cmd.png("figures/structural/fig3d_evolution_overlay.png", dpi=300)

print("Saved: fig3d_evolution_overlay.png")
```

### Expected Output
- Ancestral (purple) with looser pocket
- Modern (green) with tighter Zn coordination
- THR ligand shown
- Zn sphere at center

---

## 📜 PyMOL Script 3: Figure 4C - THR Coordinating Zn

**Purpose:** Close-up of THR bidentate coordination to Zn

### Files Needed
- `modern_thrrs_thr/modern_thrrs_thr_model.cif` (with Zn)

### PyMOL Script

```python
# fig4c_thr_zn_coordination.py
from pymol import cmd

# Load structure
cmd.load("path/to/modern_thrrs_zn_THR.cif", "structure")

# Hide all
cmd.hide("everything")

# Show active site residues as sticks
cmd.show("sticks", "(resi 100-150) or (chain B)")  # Adjust range

# Show Zn
cmd.show("sphere", "name ZN")
cmd.set("sphere_scale", 0.6, "name ZN")
cmd.color("gray", "name ZN")

# Color THR ligand
cmd.color("green", "chain B")
cmd.set("stick_radius", 0.15, "chain B")

# Color protein residues
cmd.color("wheat", "protein")

# Highlight Zn-coordinating atoms
cmd.select("zn_coord", "(chain B) within 3 of (name ZN)")
cmd.show("spheres", "zn_coord and (name O or name N)")
cmd.set("sphere_scale", 0.3, "zn_coord")
cmd.color("red", "zn_coord and name O")
cmd.color("blue", "zn_coord and name N")

# Distance measurements for Zn coordination
cmd.distance("zn_O", "name ZN", "chain B and name O", cutoff=3.0)
cmd.distance("zn_N", "name ZN", "chain B and name N", cutoff=3.0)

# Hide distance labels, just show dashes
cmd.hide("labels", "zn_*")
cmd.set("dash_color", "yellow")
cmd.set("dash_width", 3)

# Styling
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)

# Zoom to Zn
cmd.zoom("name ZN", 5)

# Add label
cmd.label("chain B and name CA", '"THR (bidentate)"')

# Render
cmd.ray(2400, 2400)
cmd.png("figures/structural/fig4c_thr_zn_coordination.png", dpi=300)

print("Saved: fig4c_thr_zn_coordination.png")
```

### Expected Output
- Close-up of Zn (gray sphere)
- THR (green sticks) with hydroxyl and amino groups coordinating
- Yellow dashes showing coordination bonds (~2.1-2.5 Å)
- Protein residues in background

---

## 📜 PyMOL Script 4: Figure 4D - ILE Rejected

**Purpose:** Show ILE cannot coordinate Zn (hydrophobic)

### Files Needed
- `modern_thrrs_ile/modern_thrrs_ile_model.cif` (with Zn, if available)

### PyMOL Script

```python
# fig4d_ile_rejected.py
from pymol import cmd

# Load structure
cmd.load("path/to/modern_thrrs_zn_ILE.cif", "structure")

# Hide all
cmd.hide("everything")

# Show active site
cmd.show("sticks", "(resi 100-150) or (chain B)")

# Show Zn
cmd.show("sphere", "name ZN")
cmd.set("sphere_scale", 0.6, "name ZN")
cmd.color("gray", "name ZN")

# Color ILE ligand (red = rejected)
cmd.color("red", "chain B")
cmd.set("stick_radius", 0.15, "chain B")

# Color protein
cmd.color("wheat", "protein")

# Show distance to Zn (should be longer, no coordination)
cmd.distance("gap", "name ZN", "chain B and name C*", cutoff=4.0)
cmd.set("dash_color", "red")
cmd.set("dash_width", 2)
cmd.hide("labels", "gap")

# Styling
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)

# Same view as 4C for comparison
cmd.zoom("name ZN", 5)

# Add label
cmd.label("chain B and name CA", '"ILE (no coordination)"')

# Render
cmd.ray(2400, 2400)
cmd.png("figures/structural/fig4d_ile_rejected.png", dpi=300)

print("Saved: fig4d_ile_rejected.png")
```

### Expected Output
- Same view as Figure 4C
- ILE (red sticks) with hydrophobic side chain
- Gap between ILE and Zn (no coordination)
- Distance > 3.5 Å (compared to < 2.5 Å for THR)

---

## 📜 PyMOL Script 5: Figure 5C - SER in Zn Site (The Trap)

**Purpose:** Show SER coordinates Zn identically to THR

### Files Needed
- `modern_thrrs_ser/modern_thrrs_ser_model.cif` (with Zn, if available)

### PyMOL Script

```python
# fig5c_ser_zinc_trap.py
from pymol import cmd

# Load structure
cmd.load("path/to/modern_thrrs_zn_SER.cif", "structure")

# Hide all
cmd.hide("everything")

# Show active site
cmd.show("sticks", "(resi 100-150) or (chain B)")

# Show Zn
cmd.show("sphere", "name ZN")
cmd.set("sphere_scale", 0.6, "name ZN")
cmd.color("gray", "name ZN")

# Color SER ligand (orange = trapped!)
cmd.color("orange", "chain B")
cmd.set("stick_radius", 0.15, "chain B")

# Color protein
cmd.color("wheat", "protein")

# Highlight Zn-coordinating atoms
cmd.select("zn_coord", "(chain B) within 3 of (name ZN)")
cmd.show("spheres", "zn_coord and (name O or name N)")
cmd.set("sphere_scale", 0.3, "zn_coord")
cmd.color("red", "zn_coord and name O")
cmd.color("blue", "zn_coord and name N")

# Distance measurements
cmd.distance("zn_O", "name ZN", "chain B and name O", cutoff=3.0)
cmd.distance("zn_N", "name ZN", "chain B and name N", cutoff=3.0)
cmd.hide("labels", "zn_*")
cmd.set("dash_color", "orange")
cmd.set("dash_width", 3)

# Styling
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)

# Same view as 4C
cmd.zoom("name ZN", 5)

# Add label with warning
cmd.label("chain B and name CA", '"SER (TRAPPED - looks like THR!)"')

# Add annotation
cmd.pseudoatom("warning", pos=[10, 10, 10])
cmd.label("warning", '"Editing domain required!"')

# Render
cmd.ray(2400, 2400)
cmd.png("figures/structural/fig5c_ser_zinc_trap.png", dpi=300)

print("Saved: fig5c_ser_zinc_trap.png")
```

### Expected Output
- Same view and coordination as Figure 4C (THR)
- SER (orange sticks) with identical coordination pattern
- Orange dashes showing coordination (~2.1-2.5 Å like THR)
- Warning annotation about editing requirement

---

## 🚀 Running the Scripts

### Option 1: Command Line

```bash
# Navigate to phase2 directory
cd /storage/kiran-stuff/aaRS/phase2/

# Run each script
pymol -c -r figures/structural/fig2d_editing_overlay.py
pymol -c -r figures/structural/fig3d_evolution_overlay.py
pymol -c -r figures/structural/fig4c_thr_zn_coordination.py
pymol -c -r figures/structural/fig4d_ile_rejected.py
pymol -c -r figures/structural/fig5c_ser_zinc_trap.py
```

### Option 2: Interactive PyMOL

```bash
# Open PyMOL
pymol

# Then in PyMOL command line:
run figures/structural/fig2d_editing_overlay.py
```

### Option 3: ChimeraX Alternative

If PyMOL is not available, ChimeraX can be used:

```bash
# Example ChimeraX command
chimerax --nogui structure.cif --script render_script.cxc
```

---

## 📋 Checklist Before Rendering

- [ ] All required CIF files are available
- [ ] File paths updated in scripts
- [ ] PyMOL/ChimeraX installed
- [ ] Output directory exists: `figures/structural/`
- [ ] Resolution set to 300 DPI
- [ ] Background is white
- [ ] Ray tracing enabled for quality

---

## 🎨 Styling Consistency

All renders should use:
- **Background**: White
- **Resolution**: 2400x2400 px (300 DPI at 8 inches)
- **Protein cartoon**: Smooth, alpha helix emphasis
- **Ligand sticks**: Radius 0.15
- **Zn sphere**: Scale 0.5-0.6, gray color
- **Ray tracing**: On, no shadows
- **Antialiasing**: Level 2

---

## 🔧 Troubleshooting

### Problem: CIF files not found
**Solution:** Check actual directory names in `outputs/` and update paths

### Problem: PyMOL not installed
**Solution:** Install PyMOL or use ChimeraX as alternative

### Problem: Zn not visible
**Solution:** Check if structure actually contains Zn (some may not)

### Problem: Alignment fails
**Solution:** Use `super` instead of `align` for more flexible alignment

---

## 📊 Expected File Sizes

- PNG files: ~1-2 MB each (2400x2400 px, 300 DPI)
- 5 renders total: ~5-10 MB

---

## ⏱️ Estimated Time

- Script setup: 10-15 minutes per render
- Rendering time: 1-2 minutes per image
- Total: ~1 hour for all 5 renders

---

**Note:** These scripts provide templates. Actual residue numbers and file paths need to be adjusted based on your specific CIF files.

---

**Generated:** 2025-12-18
**For:** aaRS Evolution manuscript structural figures
