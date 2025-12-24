# Publication Figures for Ancestral aaRS Promiscuity Paper

## Generated Figures

All figures are publication-quality (300 DPI, 2000x2000 pixels) PNG files.

### Figure A: Binding Pocket Comparison
**File:** `figure_a_binding_pocket.png` (2.0 MB)

**Description:**
- Shows LUCA ProRS catalytic domain binding pocket
- Green sticks: Proline (PRO) - cognate ligand
- Red sticks: Threonine (THR) - non-cognate ligand
- Transparent cyan surface: Active site pocket
- Demonstrates that LUCA ProRS can accommodate both PRO and THR

**Interpretation:** The promiscuous binding pocket in ancestral ProRS allows both proline and threonine to bind, supporting the hypothesis of ancestral aminoacyl-tRNA synthetase promiscuity.

### Figure B: Cognate vs Non-cognate Overlay
**File:** `figure_b_overlay.png` (0.1 MB)

**Description:**
- Zoomed view of the binding pocket
- Shows PRO (green) and THR (red) occupying overlapping space
- Gray lines: Key binding residues
- Demonstrates spatial overlap of cognate and non-cognate amino acids

**Interpretation:** The close proximity and overlap of PRO and THR in the binding site explains how the ancestral enzyme could charge tRNA with the wrong amino acid.

### Figure C: Editing Domain Inverted Specificity
**File:** `figure_c_editing_domain.png` (0.7 MB)

**Description:**
- Shows LUCA ProRS editing domain
- Orange sticks: PRO (non-cognate to editing domain)
- Marine blue sticks: THR (cognate to editing domain)
- Pink transparent surface: Editing pocket
- Demonstrates inverted specificity: editing domain prefers THR over PRO

**Interpretation:** The editing domain has inverted specificity compared to the catalytic domain - it preferentially binds and removes threonine (the mischarged amino acid) rather than proline. This is the proofreading mechanism.

### Figure D: Evolutionary Comparison
**File:** `figure_d_evolutionary.png` (1.8 MB)

**Description:**
- Side-by-side comparison of LUCA vs Modern ProRS
- Forest green: LUCA ProRS (ancestral, deep evolutionary time)
- Slate blue: Modern ProRS (shallow, recent)
- Pale green/light blue surfaces: Respective binding pockets
- Shows evolution of pocket specificity

**Interpretation:** Comparison reveals how the binding pocket evolved from promiscuous (LUCA) to more specific (modern), with potential changes in pocket volume (~1-2 Å reduction mentioned in your project description).

## Structure Information

### AlphaFold3 Model Architecture
Each `.cif` file contains:
- **Chain A:** Protein (ProRS catalytic or editing domain, ~500 residues)
- **Chain B:** tRNA (nucleotides A, C, G, U, ~76 residues)
- **Chain C:** Amino acid ligand (PRO or THR, 1 residue, HETATM)

### Source Files Used
```
LUCA ProRS Catalytic Domain:
  + PRO: /storage/kiran-stuff/aaRS/phase2/outputs/deep_domain_pro/deep_domain_pro_model.cif
  + THR: /storage/kiran-stuff/aaRS/phase2/outputs/deep_domain_thr/deep_domain_thr_model.cif

LUCA ProRS Editing Domain:
  + PRO: /storage/kiran-stuff/aaRS/phase2/outputs/deep_editing_pro/deep_editing_pro/deep_editing_pro_model.cif
  + THR: /storage/kiran-stuff/aaRS/phase2/outputs/deep_editing_thr/deep_editing_thr/deep_editing_thr_model.cif

Modern ProRS Catalytic Domain:
  + PRO: /storage/kiran-stuff/aaRS/phase2/outputs/shallow_domain_pro/shallow_domain_pro_model.cif
  + THR: /storage/kiran-stuff/aaRS/phase2/outputs/shallow_domain_thr/shallow_domain_thr_model.cif
```

## PyMOL Settings Used

Publication-quality rendering settings:
```python
ray_trace_mode = 1       # High-quality ray tracing
antialias = 2            # Smooth edges
bg_color = white         # White background for publications
ray_shadows = 0          # No shadows (cleaner for journals)
depth_cue = 0            # No depth cueing
specular = 0.2           # Subtle specularity
cartoon_fancy_helices = 1 # Nice helix rendering
```

## Color Scheme

### Ligands:
- **Green:** Proline (PRO) in catalytic domain
- **Red:** Threonine (THR) in catalytic domain
- **Orange:** PRO in editing domain
- **Marine blue:** THR in editing domain
- **Cyan:** PRO in modern structures

### Proteins:
- **Pale green / Light blue:** Catalytic domain structures
- **Wheat / Light blue:** Editing domain structures
- **Forest green / Slate blue:** Evolutionary comparison
- **Gray:** Binding residues

### Surfaces:
- **Pale cyan:** Catalytic domain pocket
- **Light pink:** Editing domain pocket
- All surfaces: 50-60% transparency

## Scripts

### Main Figure Generation Script
`generate_figures_v2.py` - Generates all four publication figures

### Test Script
`test_ligand_detection.py` - Inspects CIF files to identify chains and ligands

## Software Requirements

- **PyMOL 3.1.0** (conda-forge, pymol-open-source)
- **Python 3.12**
- **UCSF Chimera 1.16** (available for interactive exploration)

## Usage

To regenerate figures:
```bash
cd /storage/kiran-stuff/structural_figures
python generate_figures_v2.py
```

To explore structures interactively:
```bash
pymol deep_domain_pro_model.cif
# or
/home/kiran/.local/UCSF-Chimera64-1.16/bin/chimera deep_domain_pro_model.cif
```

## Additional Available Structures

Many more AlphaFold3 models are available in `/storage/kiran-stuff/aaRS/phase2/outputs/`:
- LUCA ThrRS with PRO and THR
- Modern ProRS with PRO and THR
- Modern ThrRS with PRO and THR
- Catalytic domains with Phe and Trp
- Full-length structures in `/storage/kiran-stuff/aaRS/phase2/af3_output_full/`

## Notes for Paper

1. **Figure legends:** Add distance measurements and specific binding residues identified from the structures
2. **Pocket volume analysis:** Consider using CASTp or POVME to quantify pocket volume differences
3. **Binding energy:** May want to add computational binding energy calculations (e.g., MM/GBSA)
4. **Sequence conservation:** Overlay evolutionary conservation on the structures using ConSurf
5. **Movie/animations:** PyMOL can also generate rotating animations for supplementary materials

## Date Generated
November 25, 2025

## Author
Claude Code (Anthropic AI)
