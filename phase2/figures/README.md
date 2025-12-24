# aaRS Evolution Figures - Publication Quality

**Generated:** 2025-12-18
**Target:** Cell/NSMB
**Resolution:** 300 DPI PNG + vector PDF

---

## 📊 Figure Inventory

### Figure 1: The Ancestral "Bucket" Problem
**Message:** Ancestral Class IIa aaRS was promiscuous and could not discriminate cognate from non-cognate.

| Panel | Description | File | Status |
|-------|-------------|------|--------|
| **A** | Phylogenetic tree | TBD | ⏳ Need IQ-TREE output |
| **B** | Domain architecture | TBD | ⏳ Need BioRender |
| **C** | Ancestral ThrRS landscape | `figure1/panel_c_anc_thrrs.png` | ✅ DONE |
| **D** | Ancestral ProRS landscape | `figure1/panel_d_anc_prors.png` | ✅ DONE |

**Key Finding:**
- Ancestral ThrRS: THR ranks **#8 out of 20** (ARG and ILE are better at 0.87 vs 0.85)
- Ancestral ProRS: PRO ranks **#3 out of 19** (GLU is best at 0.89 vs 0.85)

---

### Figure 2: ProRS - The Double Sieve Solution
**Message:** ProRS retained catalytic promiscuity and evolved editing domain for post-transfer correction.

| Panel | Description | File | Status |
|-------|-------------|------|--------|
| **A** | Double sieve mechanism | TBD | ⏳ Need BioRender |
| **B** | Modern ProRS catalytic | `figure2/panel_b_mod_prors_catalytic.png` | ✅ DONE |
| **C** | ProRS editing domain | `figure2/panel_c_prors_editing.png` | ✅ DONE |
| **D** | Structure: THR vs PRO in editing | TBD | ⏳ Need PyMOL/ChimeraX |

**Key Finding:**
- Modern ProRS catalytic: ALA = 98%, VAL = 97%, LEU = 95% of PRO
- ProRS editing: THR ranks #1 (0.87), PRO ranks #2 (0.82)
- Editing domain preferentially binds misacylation products

---

### Figure 3: The Zinc Disconnect
**Message:** Zn-binding evolved BEFORE Zn-mediated discrimination. Tool before function.

| Panel | Description | File | Status |
|-------|-------------|------|--------|
| **A** | Zn coordination chemistry | TBD | ⏳ Need BioRender |
| **B** | Competition experiments | `figure3/panel_b_competitions.png` | ✅ DONE |
| **C** | Structural overlay | TBD | ⏳ Need PyMOL/ChimeraX |

**Key Finding:**
- Ancestral ThrRS + Zn: THR 0.89 vs ILE 0.88 = **1.03x** (dead heat)
- Modern ThrRS + Zn: THR 0.96 vs ILE 0.79 = **1.22x** (discriminates!)
- Zn ipTM is high (0.92-0.98) in both, but only modern couples Zn to substrate selection

---

### Figure 4: The Zinc Filter
**Message:** Modern ThrRS uses Zn to reject hydrophobic AAs via steric exclusion.

| Panel | Description | File | Status |
|-------|-------------|------|--------|
| **A** | Zn filter mechanism | TBD | ⏳ Need BioRender |
| **B** | Heatmap: All 20 AAs | `figure4/panel_b_zinc_filter_heatmap.png` | ✅ DONE |
| **C** | Structure: THR coordinating Zn | TBD | ⏳ Need PyMOL/ChimeraX |
| **D** | Structure: ILE rejected | TBD | ⏳ Need PyMOL/ChimeraX |

**Key Finding:**
- THR: 0.97 (100%)
- Hydrophobics rejected: ILE 85.6%, LEU 86.6%, PRO 84.5%
- But SER escapes: 97.9% (the trap!)

---

### Figure 5: The Zinc Trap
**Message:** Zn filter fails against SER because both THR and SER coordinate Zn via hydroxyl.

| Panel | Description | File | Status |
|-------|-------------|------|--------|
| **A** | Schematic: SER coordinates Zn | TBD | ⏳ Need BioRender |
| **B** | Bar chart: THR vs SER vs ILE | `figure5/panel_ab_zinc_trap.png` | ✅ DONE |
| **C** | Structure: SER in Zn site | TBD | ⏳ Need PyMOL/ChimeraX |
| **D** | H-bond network comparison | TBD | ⏳ Need CIF analysis |

**Key Finding:**
- THR/SER discrimination: **1.02x** (too close!)
- THR/ILE discrimination: **1.17x** (good)
- Both THR and SER show Zn ipTM ~0.98
- Editing domain is REQUIRED, not optional

---

### Figure 6: Evolutionary Synthesis
**Message:** Two divergent solutions to ancestral promiscuity, constrained by chemistry.

| Panel | Description | File | Status |
|-------|-------------|------|--------|
| **ALL** | Comprehensive comparison | `figure6/comprehensive_synthesis.png` | ✅ DONE |

**Includes:**
- Row 1: Ancestral state (both promiscuous)
- Row 2: Modern state (ThrRS fixed, ProRS unfixed)
- Row 3: Comparison table

**Key Finding:**
- ThrRS: Structural solution (Zn filter) + secondary editing
- ProRS: Kinetic solution (editing domain) + retained promiscuity
- Chemical constraints drove divergent evolution

---

### Figure 7: Validation
**Message:** AF3 predictions validated by experimental data and structural analysis.

| Panel | Description | File | Status |
|-------|-------------|------|--------|
| **A** | AF3 vs experimental | TBD | ⏳ Need Tawfik data |
| **B** | H-bond analysis | TBD | ⏳ Need CIF analysis |
| **C** | Contact count | TBD | ⏳ Need CIF analysis |
| **D** | Validation summary | TBD | ⏳ TBD |

---

## 📁 Directory Structure

```
figures/
├── figure1/          # Ancestral bucket
│   ├── panel_c_anc_thrrs.png/pdf ✅
│   └── panel_d_anc_prors.png/pdf ✅
├── figure2/          # ProRS double sieve
│   ├── panel_b_mod_prors_catalytic.png/pdf ✅
│   └── panel_c_prors_editing.png/pdf ✅
├── figure3/          # Zinc disconnect
│   └── panel_b_competitions.png/pdf ✅
├── figure4/          # Zinc filter
│   └── panel_b_zinc_filter_heatmap.png/pdf ✅
├── figure5/          # Zinc trap
│   └── panel_ab_zinc_trap.png/pdf ✅
├── figure6/          # Evolutionary synthesis
│   └── comprehensive_synthesis.png/pdf ✅
├── figure7/          # Validation (TBD)
├── data/             # Processed data files
│   ├── categorized_predictions.csv
│   ├── fig1c_anc_thrrs_no_zn.csv
│   ├── fig1d_anc_prors.csv
│   ├── fig2b_mod_prors_catalytic.csv
│   ├── fig2c_prors_editing.csv
│   ├── fig3_competitions.csv
│   ├── fig4b_mod_thrrs_zn_all.csv
│   ├── fig5b_zinc_trap.csv
│   └── catalog_summary.json
└── scripts/          # Generation scripts
    ├── 01_catalog_structures.py ✅
    ├── 02_generate_figure1.py ✅
    ├── 03_generate_figure2.py ✅
    ├── 04_generate_zinc_figures.py ✅
    └── 05_generate_figure6_synthesis.py ✅
```

---

## 🎨 Color Scheme (Standard)

### By Role
- **Cognate**: `#2ecc71` (Green)
- **Near-cognate/Trapped**: `#f39c12` (Orange)
- **Rejected**: `#3498db` (Blue)
- **Better than cognate**: `#e74c3c` (Red)
- **Neutral/Low**: `#95a5a6` (Gray)

### By Enzyme
- **ThrRS**: Blues (#3498db → #2c3e50)
- **ProRS**: Oranges (#f39c12 → #e67e22)

### By Era
- **Ancestral**: Light colors
- **Modern**: Dark colors

---

## 📈 Data Summary

### Total Predictions Analyzed
- **133 high-quality predictions** (pTM ≥ 0.40, no tRNA runs)
- Ancestral ThrRS: 43 predictions
- Modern ThrRS: 27 predictions
- Ancestral ProRS: 38 predictions
- Modern ProRS: 21 predictions
- With Zn: 44 predictions
- Competition experiments: 3

### Key Metrics

| Metric | Ancestral ThrRS | Modern ThrRS | Ancestral ProRS | Modern ProRS |
|--------|----------------|--------------|-----------------|--------------|
| Cognate rank | #8/20 | #1/21 | #3/19 | #1/20 |
| Cognate score | 0.850 | 0.970 | 0.850 | 0.950 |
| Top non-cognate | ARG (0.870) | SER (0.950) | GLU (0.890) | ALA (0.930) |
| Discrimination | None | THR/ILE: 1.17x | None | Retained promiscuity |

---

## 🔧 Still Needed

### Structural Renders (PyMOL/ChimeraX)
1. **Figure 1B**: Domain architecture schematic
2. **Figure 2D**: THR vs PRO in editing domain (overlay)
3. **Figure 3C**: Ancestral vs Modern ThrRS active site (overlay)
4. **Figure 4C**: THR coordinating Zn (bidentate)
5. **Figure 4D**: ILE rejected by Zn pocket
6. **Figure 5C**: SER in Zn site (looks like THR!)

### Schematics (BioRender)
1. **Figure 1B**: Domain architecture (Cat + Edit for ProRS, Cat for ThrRS)
2. **Figure 2A**: Double sieve mechanism
3. **Figure 3A**: Zn coordination chemistry
4. **Figure 4A**: Zn filter mechanism
5. **Figure 5A**: Zinc trap concept
6. **Figure 6A**: Evolutionary tree/divergence

### Analysis (From CIF Files)
1. **Figure 5D**: H-bond network analysis
2. **Figure 7B**: Complete H-bond comparison
3. **Figure 7C**: Contact count analysis
4. **Figure 7A**: Validation against experimental data

---

## 🚀 How to Regenerate Figures

```bash
# Navigate to phase2 directory
cd /storage/kiran-stuff/aaRS/phase2/

# Run all scripts in order
python3 figures/scripts/01_catalog_structures.py    # Organize data
python3 figures/scripts/02_generate_figure1.py       # Figure 1 (C, D)
python3 figures/scripts/03_generate_figure2.py       # Figure 2 (B, C)
python3 figures/scripts/04_generate_zinc_figures.py  # Figures 3, 4, 5
python3 figures/scripts/05_generate_figure6_synthesis.py  # Figure 6

# All figures saved as PNG (300 dpi) + PDF (vector)
```

---

## 📝 Notes for Publication

1. **Font fallback**: Scripts use Arial, but fall back to DejaVu Sans (publication acceptable)
2. **Resolution**: All PNG outputs are 300 DPI (print quality)
3. **Vector formats**: PDF files provided for all figures (scalable)
4. **Color blindness**: Consider running through Coblis for accessibility
5. **Figure legends**: Will need to be written separately for manuscript

---

## 🎯 Figure Status Summary

| Figure | Panels Complete | Panels Needed | Progress |
|--------|----------------|---------------|----------|
| Fig 1 | 2/4 (C, D) | A, B | 50% |
| Fig 2 | 2/4 (B, C) | A, D | 50% |
| Fig 3 | 1/3 (B) | A, C | 33% |
| Fig 4 | 1/4 (B) | A, C, D | 25% |
| Fig 5 | 1/4 (AB) | A, C, D | 25% |
| Fig 6 | 1/1 (ALL) | None | 100% |
| Fig 7 | 0/4 | A, B, C, D | 0% |

**Overall:** 8/24 panels complete (33%)

**Data visualization panels:** 8/8 (100%) ✅
**Structural renders:** 0/9 (0%) ⏳
**BioRender schematics:** 0/6 (0%) ⏳
**Validation analysis:** 0/1 (0%) ⏳

---

## 📧 Contact

For questions about figure generation or data:
- See `AF3_RESULTS_CORRECTED.csv` for master data
- See `AF3_EVOLUTIONARY_NARRATIVE_FULL.txt` for complete analysis
- See `AF3_KEY_FINDINGS.md` for summary

**Generated by Claude Code on 2025-12-18**
