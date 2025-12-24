# Alignment Fix Applied to Figures B and C

## Problem Identified
The initial figures B and C had alignment issues where ligands appeared to be floating outside the protein instead of being positioned in the binding pocket.

## Root Cause
When superimposing two structures with different ligands:
- Initial approach: Used `cmd.align('structure1', 'structure2')` which aligned everything
- Problem: This moves the ligands independently of their proteins
- Result: Ligands appeared disconnected from the binding pocket

## Solution Applied

### Correct Alignment Strategy:
1. **Align only the PROTEINS** (chain A to chain A)
2. **Keep ligands attached** to their respective proteins
3. **Result:** Both ligands now appear in the SAME binding pocket space

### Code Fix:
```python
# WRONG (old approach):
cmd.align('luca_thr', 'luca_pro')  # Aligns everything, breaks protein-ligand relationship

# CORRECT (new approach):
cmd.align('luca_thr and chain A', 'luca_pro and chain A')  # Aligns only proteins
```

## Fixed Figures

### Figure B: Cognate vs Non-cognate Overlay (FIXED)
- File: `figure_b_overlay.png` (206 KB)
- PRO (green sticks) and THR (red sticks) now **overlap in the same pocket**
- Pocket residues (wheat) shown within 5Å of both ligands
- Demonstrates spatial overlap of cognate and non-cognate amino acids

**Key visualization improvements:**
- Both ligands clearly visible in binding pocket
- Pocket residues shown as context
- Labels added (PRO, THR)
- Zoomed to 3Å around ligands for detail

### Figure C: Editing Domain Inverted Specificity (FIXED)
- File: `figure_c_editing_domain.png` (732 KB)
- PRO (orange) and THR (marine blue) now **both inside editing pocket**
- Pink transparent surface shows the editing domain pocket
- Gray sticks show key binding residues

**Key visualization improvements:**
- Surface encompasses both ligands
- Cartoon shows protein context (30% transparent)
- Pocket residues visible
- Labels added

## Verification

Both figures now correctly show:
1. ✅ Ligands positioned inside the binding pocket
2. ✅ PRO and THR spatially overlapping
3. ✅ Pocket residues/surface surrounding both ligands
4. ✅ Clear spatial relationship demonstrating promiscuity

## Technical Details

### Structure Architecture:
- Chain A: Protein (ProRS catalytic or editing domain)
- Chain B: tRNA
- Chain C: Ligand (PRO or THR as HETATM)

### Alignment Process:
```python
# Load both structures
cmd.load(STRUCTURES['luca_cat_pro'], 'luca_pro')
cmd.load(STRUCTURES['luca_cat_thr'], 'luca_thr')

# Align ONLY the protein chains
cmd.align('luca_thr and chain A', 'luca_pro and chain A')

# Now ligands are automatically in the same pocket space
# because they maintain their position relative to their protein
```

### Pocket Display Strategy:
```python
# Select pocket residues within 5Å of EITHER ligand
cmd.select('pocket_res',
           '(luca_pro and chain A within 5 of (lig_pro or lig_thr))')

# This ensures the pocket encompasses both ligands
```

## Files Generated

**Fixed scripts:**
- `fix_alignment_figures.py` - Script that generates correctly aligned figures

**Fixed figures:**
- `figure_b_overlay_FIXED.png` - Original fixed version
- `figure_b_overlay.png` - Replaced with fixed version
- `figure_c_editing_domain_FIXED.png` - Original fixed version
- `figure_c_editing_domain.png` - Replaced with fixed version

**Unchanged figures (already correct):**
- `figure_a_binding_pocket.png` - Didn't need fixing
- `figure_d_evolutionary.png` - Didn't need fixing

## Date Fixed
November 25, 2025

## Summary
The alignment fix ensures that the figures accurately represent the biological reality: both PRO and THR bind in the same pocket, demonstrating the ancestral promiscuity of the enzyme. This is critical for the paper's main thesis.
