# ProRS Manuscript Figures - Quick Reference

## Generated Files

```
/storage/kiran-stuff/aaRS/
├── figure1_luca_promiscuity_clean.png       617 KB  [LUCA + PRO/THR overlay]
├── figure2_editing_inverted_clean.png        85 KB  [Editing domain]
├── figure3_luca_vs_modern_clean.png         613 KB  [Evolution comparison]
├── figure5_iptm_bars.png                    205 KB  [Bar charts 3-panel]
├── figure5_iptm_comprehensive.png           235 KB  [Bar charts 6-panel]
├── figure5_iptm_bars.pdf                     23 KB  [Vector]
└── figure5_iptm_bars.svg                     74 KB  [Vector]
```

## Key Results at a Glance

### Cross-Reactivity (THR/PRO ratio)
- **LUCA ProRS:** 82.7% (PRO: 0.75, THR: 0.62)
- **Modern ProRS:** 89.2% (PRO: 0.83, THR: 0.74) ← WORSE!
- **LUCA Editing:** 321% (PRO: 0.14, THR: 0.45) ← INVERTED!

### Structural Conservation
- **LUCA vs Modern:** RMSD = 0.63 Å (highly conserved)
- **PRO vs THR binding:** RMSD = 0.94 Å (nearly identical)

## Quick Commands

```bash
# View all figures
bash view_all_figures.sh

# Regenerate all figures
bash regenerate_all_figures.sh

# View in PyMOL (interactive)
/home/kiran/miniforge3/bin/pymol figure1_session.pse

# View documentation
less ADDITIONAL_FIGURES_README.md

# View summary
cat FIGURE_GENERATION_SUMMARY.txt
```

## Color Scheme

**Structural Figures:**
- PRO (catalytic): Green
- THR (catalytic): Red/Firebrick
- PRO (editing): Orange  
- THR (editing): Marine Blue
- Protein surface: Gray80 or Light Blue

**Bar Charts:**
- PRO bars: Blue (#2E86AB)
- THR bars: Purple/Magenta (#A23B72)

## Figure Captions (Draft)

**Figure 1:** LUCA ProRS catalytic domain bound to cognate (PRO, green) and non-cognate (THR, red) substrates, demonstrating ancestral promiscuity.

**Figure 2:** LUCA ProRS editing domain exhibits inverted specificity, binding THR (marine, ipTM=0.45) 3.2× stronger than PRO (orange, ipTM=0.14).

**Figure 3:** Side-by-side comparison of LUCA (3.5 Gya, left) and modern (right) ProRS. Despite 3.5 billion years of evolution, binding pocket architecture remains conserved (RMSD=0.63 Å).

**Figure 5:** Quantification of substrate binding specificity showing persistent cross-reactivity across evolutionary time.

## Contact
Generated: November 26, 2025
Documentation: ADDITIONAL_FIGURES_README.md
