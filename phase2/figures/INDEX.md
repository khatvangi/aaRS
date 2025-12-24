# Figure Index - Quick Reference

## 📁 Directory: `/storage/kiran-stuff/aaRS/phase2/figures/`

**Total Size:** 2.5 MB
**Generated:** 2025-12-18

---

## 🎯 Quick Access

### Complete Figures Ready for Manuscript

| Figure | File | Size | Description |
|--------|------|------|-------------|
| **Fig 1C** | `figure1/panel_c_anc_thrrs.png` | 216 KB | Ancestral ThrRS - THR ranks #8 |
| **Fig 1D** | `figure1/panel_d_anc_prors.png` | 210 KB | Ancestral ProRS - PRO ranks #3 |
| **Fig 2B** | `figure2/panel_b_mod_prors_catalytic.png` | 226 KB | Modern ProRS - ALA 98% |
| **Fig 2C** | `figure2/panel_c_prors_editing.png` | 245 KB | Editing domain - THR #1 |
| **Fig 3B** | `figure3/panel_b_competitions.png` | 179 KB | Competitions - 1.03x → 1.22x |
| **Fig 4B** | `figure4/panel_b_zinc_filter_heatmap.png` | 259 KB | Heatmap - all 20 AAs |
| **Fig 5AB** | `figure5/panel_ab_zinc_trap.png` | 186 KB | Zinc trap - SER escapes |
| **Fig 6** | `figure6/comprehensive_synthesis.png` | 547 KB | Complete synthesis |

### Vector Formats (for journals)

All figures also available as PDF:
- `figure*/panel_*.pdf` (33-50 KB each, scalable)

---

## 📊 Data Files

| File | Purpose | Rows |
|------|---------|------|
| `data/categorized_predictions.csv` | Master dataset | 133 |
| `data/fig1c_anc_thrrs_no_zn.csv` | Anc ThrRS | 20 |
| `data/fig1d_anc_prors.csv` | Anc ProRS | 19 |
| `data/fig2b_mod_prors_catalytic.csv` | Mod ProRS cat | 20 |
| `data/fig2c_prors_editing.csv` | ProRS edit | 19 |
| `data/fig3_competitions.csv` | Competitions | 3 |
| `data/fig4b_mod_thrrs_zn_all.csv` | Mod ThrRS+Zn | 21 |
| `data/fig5b_zinc_trap.csv` | Zinc trap | 3 |
| `data/catalog_summary.json` | Stats | - |

---

## 🔧 Scripts (Reusable)

| Script | Function | Runtime |
|--------|----------|---------|
| `scripts/01_catalog_structures.py` | Organize data | <1 min |
| `scripts/02_generate_figure1.py` | Fig 1 panels | <10 sec |
| `scripts/03_generate_figure2.py` | Fig 2 panels | <10 sec |
| `scripts/04_generate_zinc_figures.py` | Figs 3,4,5 | <15 sec |
| `scripts/05_generate_figure6_synthesis.py` | Fig 6 | <10 sec |

**To regenerate all figures:**
```bash
cd /storage/kiran-stuff/aaRS/phase2/
python3 figures/scripts/01_catalog_structures.py
python3 figures/scripts/02_generate_figure1.py
python3 figures/scripts/03_generate_figure2.py
python3 figures/scripts/04_generate_zinc_figures.py
python3 figures/scripts/05_generate_figure6_synthesis.py
```

---

## 📖 Documentation

- **README.md** - Comprehensive guide with panel descriptions
- **FIGURE_GENERATION_SUMMARY.md** - Detailed progress report
- **INDEX.md** - This file (quick reference)

---

## ✅ Completion Status

**Data Visualization:** 8/8 (100%) ✅
**Total Panels:** 8/24 (33%)

**Still Needed:**
- 5 structural renders (PyMOL/ChimeraX)
- 5 BioRender schematics
- 3 validation analyses

---

## 🔑 Key Numbers at a Glance

| Metric | Value |
|--------|-------|
| Total AF3 predictions analyzed | 187 |
| High-quality predictions (pTM≥0.40) | 133 |
| Ancestral ThrRS THR rank | 8/20 |
| Ancestral ProRS PRO rank | 3/19 |
| Modern ThrRS THR score | 0.970 |
| Modern ThrRS SER score | 0.950 (97.9%) |
| Modern ProRS ALA score | 0.930 (97.9%) |
| Anc ThrRS+Zn discrimination | 1.03x (none) |
| Mod ThrRS+Zn discrimination | 1.22x (works!) |

---

## 📧 Quick Links

- Parent directory: `/storage/kiran-stuff/aaRS/phase2/`
- Master data: `../AF3_RESULTS_CORRECTED.csv`
- Analysis output: `../AF3_EVOLUTIONARY_NARRATIVE_FULL.txt`
- Key findings: `../AF3_KEY_FINDINGS.md`

---

**Last updated:** 2025-12-18
