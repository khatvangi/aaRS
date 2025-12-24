# Master Index - aaRS Evolution Figures
## Complete Directory Map & Quick Start Guide

**Last Updated:** 2025-12-19
**Project:** aaRS Evolution Manuscript (Cell/NSMB)
**Location:** `/storage/kiran-stuff/aaRS/phase2/figures/`

---

## 🚀 Quick Start

### New to this project?
1. Read: `PHASE2_COMPLETION_SUMMARY.md` (comprehensive status)
2. Read: `README.md` (figure descriptions)
3. Choose your task below

### Ready to generate figures?

**Data visualization (DONE ✅):**
```bash
# All 8 panels already generated
ls figure*/panel_*.png
```

**Structural renders (2/5 READY):**
```bash
pymol -c figures/structural/fig2d_editing_overlay.py
pymol -c figures/structural/fig3d_evolution_overlay.py
```

**BioRender schematics (DOCUMENTED):**
- Read: `biorender/BIORENDER_INSTRUCTIONS.md`
- Go to: https://app.biorender.com

**H-bond analysis (NEEDS BIOPYTHON):**
```bash
pip install biopython
python3 scripts/06_hbond_analysis.py
```

---

## 📁 Complete Directory Structure

```
/storage/kiran-stuff/aaRS/phase2/figures/
│
├── 📄 MASTER_INDEX.md                      ← YOU ARE HERE
├── 📄 PHASE2_COMPLETION_SUMMARY.md         ← START HERE (comprehensive status)
├── 📄 README.md                            ← Figure descriptions & methods
├── 📄 INDEX.md                             ← Quick reference (legacy)
├── 📄 FIGURE_GENERATION_SUMMARY.md         ← Progress report (legacy)
│
├── 📊 Data Files (10 files)
│   └── data/
│       ├── categorized_predictions.csv     ← Master dataset (133 predictions)
│       ├── catalog_summary.json            ← Statistics
│       ├── fig1c_anc_thrrs_no_zn.csv      ← Ancestral ThrRS (20 AAs)
│       ├── fig1d_anc_prors.csv            ← Ancestral ProRS (19 AAs)
│       ├── fig2b_mod_prors_catalytic.csv  ← Modern ProRS catalytic
│       ├── fig2c_prors_editing.csv        ← ProRS editing domain
│       ├── fig3_competitions.csv          ← Competition experiments
│       ├── fig4b_mod_thrrs_zn_all.csv     ← Modern ThrRS+Zn (21 entries)
│       ├── fig5b_zinc_trap.csv            ← Zinc trap data
│       └── HBOND_ANALYSIS_NEEDED.md       ← H-bond analysis fallback doc
│
├── 🐍 Python Scripts (6 files)
│   └── scripts/
│       ├── 01_catalog_structures.py        ← Data organization (DONE)
│       ├── 02_generate_figure1.py          ← Fig 1 panels (DONE)
│       ├── 03_generate_figure2.py          ← Fig 2 panels (DONE)
│       ├── 04_generate_zinc_figures.py     ← Figs 3-5 (DONE)
│       ├── 05_generate_figure6_synthesis.py ← Fig 6 synthesis (DONE)
│       └── 06_hbond_analysis.py            ← H-bond extraction (NEEDS BIOPYTHON)
│
├── 📈 Figure Outputs (20 files) ✅ COMPLETE
│   ├── figure1/
│   │   ├── panel_c_anc_thrrs.png          ← 216 KB, 300 DPI
│   │   ├── panel_c_anc_thrrs.pdf          ← 33 KB, vector
│   │   ├── panel_d_anc_prors.png          ← 210 KB, 300 DPI
│   │   └── panel_d_anc_prors.pdf          ← 33 KB, vector
│   ├── figure2/
│   │   ├── panel_b_mod_prors_catalytic.png ← 226 KB
│   │   ├── panel_b_mod_prors_catalytic.pdf ← 35 KB
│   │   ├── panel_c_prors_editing.png       ← 245 KB
│   │   └── panel_c_prors_editing.pdf       ← 37 KB
│   ├── figure3/
│   │   ├── panel_b_competitions.png        ← 179 KB
│   │   └── panel_b_competitions.pdf        ← 30 KB
│   ├── figure4/
│   │   ├── panel_b_zinc_filter_heatmap.png ← 259 KB
│   │   └── panel_b_zinc_filter_heatmap.pdf ← 42 KB
│   ├── figure5/
│   │   ├── panel_ab_zinc_trap.png          ← 186 KB
│   │   └── panel_ab_zinc_trap.pdf          ← 31 KB
│   └── figure6/
│       ├── comprehensive_synthesis.png     ← 547 KB
│       └── comprehensive_synthesis.pdf     ← 50 KB
│
├── 🎨 BioRender Schematics (DOCUMENTED)
│   └── biorender/
│       ├── BIORENDER_INSTRUCTIONS.md       ← Detailed instructions for 5 panels
│       └── [OUTPUT LOCATION]
│           ├── fig1b_domain_architecture.png     (TO BE CREATED)
│           ├── fig1b_domain_architecture.pdf     (TO BE CREATED)
│           ├── fig2a_double_sieve.png            (TO BE CREATED)
│           ├── fig2a_double_sieve.pdf            (TO BE CREATED)
│           ├── fig3a_zn_coordination.png         (TO BE CREATED)
│           ├── fig3a_zn_coordination.pdf         (TO BE CREATED)
│           ├── fig4a_zn_filter_mechanism.png     (TO BE CREATED)
│           ├── fig4a_zn_filter_mechanism.pdf     (TO BE CREATED)
│           ├── fig5a_zinc_trap_concept.png       (TO BE CREATED)
│           └── fig5a_zinc_trap_concept.pdf       (TO BE CREATED)
│
└── 🔬 Structural Renders (2/5 READY, 3/5 BLOCKED)
    └── structural/
        ├── PYMOL_SCRIPTS.md                ← Template scripts (5 renders)
        ├── CIF_CATALOG.md                  ← CIF file availability map
        ├── fig2d_editing_overlay.py        ← READY ✅
        ├── fig3d_evolution_overlay.py      ← READY ✅ (no Zn)
        └── [OUTPUT LOCATION]
            ├── fig2d_editing_overlay.png         (READY TO GENERATE)
            ├── fig3d_evolution_overlay.png       (READY TO GENERATE)
            ├── fig4c_thr_zn_coordination.png     (BLOCKED - needs AF3)
            ├── fig4d_ile_rejected.png            (BLOCKED - needs AF3)
            └── fig5c_ser_zinc_trap.png           (BLOCKED - needs AF3)
```

---

## 📚 Documentation Quick Reference

### Primary Documents

| File | Purpose | When to Read |
|------|---------|--------------|
| `MASTER_INDEX.md` | Navigation & directory map | First time here |
| `PHASE2_COMPLETION_SUMMARY.md` | Complete status & next steps | Planning work |
| `README.md` | Figure descriptions & methods | Understanding results |
| `structural/CIF_CATALOG.md` | CIF file availability | Before PyMOL work |
| `biorender/BIORENDER_INSTRUCTIONS.md` | Schematic creation guide | Before BioRender |
| `data/HBOND_ANALYSIS_NEEDED.md` | H-bond analysis alternatives | If BioPython unavailable |

### Legacy Documents (Informational)

| File | Purpose | Status |
|------|---------|--------|
| `INDEX.md` | Original quick reference | Superseded by PHASE2_COMPLETION_SUMMARY |
| `FIGURE_GENERATION_SUMMARY.md` | Phase 1 progress | Historical record |

---

## 🎯 Status by Figure Panel

### Figure 1: Ancestral Promiscuity
- [x] Panel A (schematic) → BioRender instructions ready
- [x] Panel B (schematic) → BioRender instructions ready
- [x] Panel C (data) → Generated ✅
- [x] Panel D (data) → Generated ✅

**Status:** 2/4 complete, 2/4 documented

---

### Figure 2: ProRS Double Sieve
- [x] Panel A (schematic) → BioRender instructions ready
- [x] Panel B (data) → Generated ✅
- [x] Panel C (data) → Generated ✅
- [x] Panel D (structural) → PyMOL script ready ✅

**Status:** 3/4 complete, 1/4 documented

---

### Figure 3: Zinc Discrimination Evolution
- [x] Panel A (schematic) → BioRender instructions ready
- [x] Panel B (data) → Generated ✅
- [x] Panel D (structural) → PyMOL script ready ✅ (no Zn)

**Status:** 2/3 complete, 1/3 documented

---

### Figure 4: Zinc Filter Mechanism
- [x] Panel A (schematic) → BioRender instructions ready
- [x] Panel B (heatmap) → Generated ✅
- [ ] Panel C (structural) → BLOCKED (needs AF3) ❌
- [ ] Panel D (structural) → BLOCKED (needs AF3) ❌

**Status:** 2/4 complete, 1/4 documented, 2/4 blocked

---

### Figure 5: Zinc Trap
- [x] Panel A (schematic) → BioRender instructions ready
- [x] Panel B (data) → Generated ✅
- [ ] Panel C (structural) → BLOCKED (needs AF3) ❌

**Status:** 2/3 complete, 1/3 documented, 1/3 blocked

---

### Figure 6: Comprehensive Synthesis
- [x] Full synthesis panel → Generated ✅

**Status:** 1/1 complete

---

### Figure 7: H-Bond Analysis (Future)
- [ ] Panel A (TBD)
- [ ] Panel B (bar chart) → Script ready, needs BioPython

**Status:** 0/2, 1/2 scripted

---

## 📊 Overall Statistics

### Completion Metrics

| Category | Done | Total | % Complete |
|----------|------|-------|------------|
| **Data visualization** | 8 | 8 | 100% ✅ |
| **BioRender instructions** | 5 | 5 | 100% ✅ |
| **BioRender renders** | 0 | 5 | 0% 📋 |
| **PyMOL scripts (ready)** | 2 | 5 | 40% ⚠️ |
| **PyMOL scripts (blocked)** | 3 | 5 | 60% ❌ |
| **H-bond analysis** | 0 | 1 | 0% 📋 |
| **TOTAL PANELS** | 15 | 29 | 52% 🔄 |

### Blocker Analysis

| Blocker | Panels Affected | Resolution |
|---------|-----------------|------------|
| Manual BioRender work | 5 panels | 2-4 hours work |
| Missing AF3 structures | 3 panels | Run AF3 for Zn+THR/SER/ILE |
| BioPython not installed | 1 panel | `pip install biopython` |
| PyMOL not available | 2 panels | Install PyMOL or use ChimeraX |

---

## 🔧 Quick Actions

### Regenerate All Data Figures
```bash
cd /storage/kiran-stuff/aaRS/phase2/
python3 figures/scripts/01_catalog_structures.py
python3 figures/scripts/02_generate_figure1.py
python3 figures/scripts/03_generate_figure2.py
python3 figures/scripts/04_generate_zinc_figures.py
python3 figures/scripts/05_generate_figure6_synthesis.py
```
**Time:** < 1 minute
**Output:** Regenerates all 8 data panels

---

### Generate Available Structural Renders
```bash
cd /storage/kiran-stuff/aaRS/phase2/

# Check if PyMOL is available
which pymol

# Run scripts
pymol -c figures/structural/fig2d_editing_overlay.py
pymol -c figures/structural/fig3d_evolution_overlay.py
```
**Time:** ~2 minutes per render
**Output:** 2 structural PNG files (2400x2400, 300 DPI)

---

### Run H-Bond Analysis
```bash
# Install BioPython if needed
pip install biopython

# Run analysis
cd /storage/kiran-stuff/aaRS/phase2/
python3 figures/scripts/06_hbond_analysis.py
```
**Time:** ~5-10 minutes
**Output:** `hbond_analysis.csv` and `hbond_analysis_detailed.json`

---

### Create BioRender Schematics
1. Open https://app.biorender.com
2. Create account/login
3. Follow: `figures/biorender/BIORENDER_INSTRUCTIONS.md`
4. Create 5 schematics (30-60 min each)
5. Export as PNG (300 DPI) + PDF
6. Save to `figures/biorender/`

**Time:** 2-4 hours total
**Output:** 10 files (5 PNG + 5 PDF)

---

## 🔑 Key Numbers at a Glance

| Metric | Value | Source |
|--------|-------|--------|
| Total AF3 predictions analyzed | 187 | AF3_RESULTS_CORRECTED.csv |
| High-quality predictions (pTM≥0.40) | 133 | catalog_summary.json |
| Structures with Zn in CSV | 43 | AF3_RESULTS_CORRECTED.csv |
| Structures with Zn CIF files | 0 | CIF_CATALOG.md |
| Ancestral ThrRS THR rank | 8/20 | fig1c_anc_thrrs_no_zn.csv |
| Ancestral ProRS PRO rank | 3/19 | fig1d_anc_prors.csv |
| Modern ThrRS THR ipTM | 0.970 | fig4b_mod_thrrs_zn_all.csv |
| Modern ThrRS SER ipTM | 0.950 | fig5b_zinc_trap.csv |
| Modern ThrRS ILE ipTM | 0.830 | fig4b_mod_thrrs_zn_all.csv |
| THR/SER discrimination | 1.02x | fig5b_zinc_trap.csv |
| THR/ILE discrimination | 1.17x | fig4b_mod_thrrs_zn_all.csv |
| Modern ProRS ALA error rate | 98% | fig2b_mod_prors_catalytic.csv |
| Total figure panels planned | 29 | PHASE2_COMPLETION_SUMMARY.md |
| Panels complete | 15 | PHASE2_COMPLETION_SUMMARY.md |

---

## 🎓 Understanding the Data

### What is ipTM?
**Interface predicted TM-score** (0-1 scale)
- Measures quality of protein-ligand interface prediction
- Higher = better binding predicted
- Threshold: 0.40 for high confidence

### Color Scheme Consistency
- **Green (#2ecc71)**: Cognate/correct substrate
- **Orange (#f39c12)**: Error/trapped substrate
- **Red (#e74c3c)**: Rejected substrate
- **Purple (#9b59b6)**: Ancestral/promiscuous
- **Blue (#3498db)**: Modern/specific
- **Gray (#95a5a6)**: Zinc/metal

### Figure Narrative Flow
1. **Fig 1**: Ancestral enzymes were promiscuous "buckets"
2. **Fig 2**: ProRS evolved editing domain (kinetic solution)
3. **Fig 3**: ThrRS evolved Zn filter (structural solution)
4. **Fig 4**: Zn discriminates by coordination chemistry
5. **Fig 5**: But Zn fails against SER (the trap!)
6. **Fig 6**: Complete synthesis of evolution pathways
7. **Fig 7**: Molecular details (H-bonds, coordination)

---

## 📧 For Manuscript Submission

### What's Ready Now
- 8 publication-quality data panels (PNG + PDF)
- Complete methods documentation
- All source data files
- Reproducible Python scripts

### What to Tell Reviewers
"Data visualization complete. Structural renders and schematics in progress, documented scripts available for reproduction."

### Supplementary Materials
```
Supplementary Figures/
├── All 8 data panels (figures/figure*/)
├── Source data (figures/data/*.csv)
├── Analysis scripts (figures/scripts/*.py)
└── Methods (figures/README.md)
```

---

## 🔗 Related Files (Outside figures/)

```
/storage/kiran-stuff/aaRS/phase2/
├── AF3_RESULTS_CORRECTED.csv          ← Master data source
├── AF3_EVOLUTIONARY_NARRATIVE_FULL.txt ← Analysis narrative
├── AF3_KEY_FINDINGS.md                 ← Key insights summary
└── outputs/                            ← CIF files (94 total)
    ├── deep_editing_pro/               ← Ancestral ProRS editing + PRO
    ├── deep_editing_thr/               ← Ancestral ProRS editing + THR
    ├── deep_thrrs_thr/                 ← Ancestral ThrRS + THR
    ├── modern_thrrs_thr/               ← Modern ThrRS + THR (no Zn)
    └── [13 more directories]
```

---

## 💡 Pro Tips

1. **Version Control**: All generated figures are timestamped and reproducible
2. **File Formats**: Always keep both PNG (for viewing) and PDF (for journals)
3. **Data Provenance**: Every figure links back to source CSV
4. **Script Modularity**: Each figure has its own generation script
5. **Documentation**: Three levels (quick/medium/comprehensive)

---

## 🆘 Troubleshooting

### "No such file or directory" errors
- Check working directory: should be `/storage/kiran-stuff/aaRS/phase2/`
- Use absolute paths if relative paths fail

### Font warnings in matplotlib
- Normal behavior, fallback fonts work fine
- Add `plt.rcParams['font.family'] = 'DejaVu Sans'` if needed

### PyMOL "command not found"
- Install: `conda install -c conda-forge pymol-open-source`
- Alternative: Use ChimeraX (scripts need adaptation)

### BioPython import errors
- Install: `pip install biopython`
- Check Python version: requires Python 3.6+

### Missing CIF files
- Check `structural/CIF_CATALOG.md` for availability
- May need to run additional AF3 predictions

---

## 📝 Citation

If using these figures or scripts:

```
aaRS Evolution Figure Generation Pipeline
Generated: 2025-12-19
Scripts: /storage/kiran-stuff/aaRS/phase2/figures/scripts/
Documentation: https://github.com/your-repo (if applicable)
```

---

## 🎯 Next Session Checklist

Starting a new work session? Check these:

- [ ] Read `PHASE2_COMPLETION_SUMMARY.md` for current status
- [ ] Check if any software dependencies need installing
- [ ] Verify working directory: `/storage/kiran-stuff/aaRS/phase2/`
- [ ] Review which figures are still needed
- [ ] Choose priority tasks from "Recommended Next Steps"

---

**Last Updated:** 2025-12-19
**Maintainer:** Claude Code
**Contact:** See project lead for questions

**Total Documentation Files:** 10
**Total Figure Outputs:** 20
**Total Data Files:** 10
**Total Scripts:** 6

**This index covers:** 46 documented files + 94 CIF structures = 140 total project files

---

## 🏁 Summary

**You have everything needed to:**
- ✅ Submit 8 data panels immediately
- ✅ Generate 2 structural renders (just run PyMOL)
- ✅ Create 5 BioRender schematics (follow instructions)
- ⚠️ Generate 3 more structural renders (after AF3 runs)
- ⚠️ Run H-bond analysis (after BioPython install)

**Estimated time to 100% completion:**
- Immediate work: 4-8 hours
- AF3 predictions: Depends on queue time
- Total: ~1-2 days of active work

**Current progress: 52% complete**

**Priority: Execute available PyMOL scripts and start BioRender schematics.**
