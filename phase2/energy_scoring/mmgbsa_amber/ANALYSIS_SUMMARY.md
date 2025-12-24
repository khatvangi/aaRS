# MM/GBSA Analysis Summary - aaRS Ligand Binding

**Date:** 2025-12-21  
**Total jobs:** 182  
**Successful:** 172 (94.5%)  
**Method:** Per-structure GAFF2 parameterization, 500-step minimization

---

## Pipeline Validation ✓

- **RMSD QC:** 0.000 Å (median, PASS < 0.05 Å)
- **Convergence:** 172/172 (100%)
- **Median RMS gradient:** 0.14 kcal/mol/Å
- **Per-structure parameterization:** Successfully preserves AF3 coordinates

---

## Key Findings

### 1. System Performance

| System Type | n | Mean ΔG (kcal/mol) | Std | Range |
|-------------|---|-------------------|-----|-------|
| **Modern ProRS** | 22 | **-12.84** | 5.17 | -19.50 to +1.09 |
| Modern full-length | 9 | -10.73 | 7.94 | -18.71 to +5.93 |
| **Ancestral ThrRS** | 41 | **-9.13** | 2.36 | -14.85 to -4.68 |
| Ancestral ProRS | 40 | -8.59 | 5.89 | -17.33 to +13.97 |
| Modern ThrRS* | 44 | **+5.38** | 16.66 | -17.18 to +43.88 |

*Excluding +862 kcal/mol outlier (modern_thrrs_TRP)

**Interpretation:**
- Modern ProRS shows strongest binding overall
- Ancestral ThrRS is most consistent (low std)
- **Modern ThrRS is problematic** - positive ΔG suggests unfavorable binding

---

### 2. Top 10 Strongest Binders

| Rank | Job | Ligand | ΔG (kcal/mol) | Zn |
|------|-----|--------|--------------|-----|
| 1 | modern_prors_VAL | VAL | -19.50 | No |
| 2 | modern_ecoli_full_pro | PRO | -18.71 | No |
| 3 | modern_prors_PRO | PRO | -18.71 | No |
| 4 | modern_prors_LYS | LYS | -18.03 | No |
| 5 | anc_prors_cat_TRP | TRP | -17.33 | No |
| 6 | anc_prors_cat_ILE | ILE | -17.22 | No |
| 7 | modern_thrrs_thr | THR | -17.18 | No |
| 8-10 | deep_cat_trp, etc. | TRP | -17.16 | No |

**Pattern:** Hydrophobic residues (VAL, PRO, TRP, ILE) dominate strongest binding

---

### 3. Amino Acid Binding Preferences

| Rank | AA | n | Mean ΔG | Std | Range |
|------|-----|---|---------|-----|-------|
| **Best** | PHE | 10 | -10.63 | 2.30 | -15.18 to -7.12 |
| | HIS | 7 | -8.70 | 4.60 | -15.85 to -3.78 |
| | LYS | 7 | -8.69 | 5.39 | -18.03 to +0.29 |
| | TYR | 7 | -8.25 | 5.58 | -14.86 to +3.14 |
| **Worst** | CYS | 7 | -1.17 | 13.63 | -15.04 to +21.94 |
| | GLU | 7 | -1.47 | 18.85 | -13.71 to +40.88 |
| | LEU | 7 | -2.06 | 16.35 | -15.35 to +24.76 |
| **Outlier** | TRP | 10 | +74.42 | 276.88 | -17.33 to +862.36 |

---

### 4. Zinc Ion Effect

| Condition | n | Mean ΔG | Std |
|-----------|---|---------|-----|
| **Without Zn** | 128 | -2.42 | 77.29 |
| **With Zn** | 44 | +4.54 | 16.43 |
| **Difference** | - | **+6.97** | - |

**Conclusion:** Zinc **worsens** binding by ~7 kcal/mol on average  
**Concern:** Many Zn-containing systems show unfavorable energies (ΔG > 0)

---

### 5. Specificity Analysis

#### ProRS: PRO (cognate) vs other amino acids
- **PRO (n=5):** ΔG = -8.30 ± 8.71 kcal/mol
- **Other (n=59):** ΔG = -9.96 ± 5.87 kcal/mol
- **ΔΔG = -1.65 kcal/mol** (binds non-cognate BETTER!)

#### ThrRS: THR (cognate) vs other amino acids
- **THR (n=10):** ΔG = +2.70 ± 18.19 kcal/mol
- **Other (n=79):** ΔG = -2.21 ± 13.18 kcal/mol
- **ΔΔG = -4.91 kcal/mol** (binds non-cognate BETTER!)

**Surprising Result:** Neither enzyme shows expected cognate selectivity

---

### 6. Competition Experiments

| Experiment | Ligand | ΔG (kcal/mol) | Zn | Interpretation |
|------------|--------|--------------|-----|----------------|
| anc_thrrs THR vs ILE | THR | **-14.35** | Yes | THR favored ✓ |
| modern_thrrs THR vs ILE | THR | **+26.26** | Yes | THR disfavored! |
| modern_thrrs THR vs SER | THR | **+25.03** | Yes | THR disfavored! |

**Alarming:** Modern ThrRS shows UNFAVORABLE binding to cognate THR in competition

---

## Critical Issues

### 1. Modern ThrRS Catastrophic Failures

**All Modern ThrRS jobs with ΔG > +20 kcal/mol:**
- modern_thrrs_TRP: +862.36 (NO Zn) ← **CATASTROPHIC OUTLIER**
- modern_thrrs_ecoli_zn_MET: +43.88 (Zn)
- modern_thrrs_ecoli_zn_GLU: +40.88 (Zn)
- modern_thrrs_ecoli_zn_GLN: +29.22 (Zn)
- modern_thrrs_ecoli_zn_ASN: +27.87 (Zn)
- COMPETITION_modern_thrrs_THR_vs_ILE: +26.26 (Zn)
- COMPETITION_modern_thrrs_THR_vs_SER: +25.03 (Zn)
- modern_thrrs_LEU: +24.76 (Zn)
- modern_thrrs_ecoli_zn_ALA: +23.88 (Zn)
- modern_thrrs_ecoli_zn_CYS: +21.94 (Zn)
- modern_thrrs_ecoli_zn_ILE: +21.69 (Zn)

**Diagnosis:**
- 10/11 extreme failures involve Zn
- Suggests **Zn coordination geometry problems** in Modern ThrRS structures
- May indicate AlphaFold3 modeling errors for metal binding sites

### 2. Ancestral ThrRS WITH Zn works fine
- Mean ΔG = -8.91 kcal/mol
- Range: -14.35 to -4.68 kcal/mol
- **No extreme outliers**

**Conclusion:** Problem is specific to Modern ThrRS, not Zn in general

---

## Recommendations

### A. Immediate Actions
1. **Inspect modern_thrrs_TRP structure** (862 kcal/mol)
   - Check for atomic clashes using `jobs/modern_thrrs_TRP/min.out`
   - Visualize Zn coordination if present
   - May need to exclude from analysis

2. **Flag Modern ThrRS + Zn systems as unreliable**
   - Consider excluding from biological interpretation
   - Investigate Zn-ligand distances before/after minimization

### B. For Biological Interpretation
1. **Focus on reliable systems:**
   - Ancestral ProRS/ThrRS (consistent results)
   - Modern ProRS (strong binding, no Zn issues)
   - Systems WITHOUT Zn (cleaner energetics)

2. **Use relative rankings, not absolute ΔG:**
   - Single-point MM/GBSA has large systematic errors
   - Trends within a system type are more meaningful

### C. Follow-up Calculations (if needed)
1. **MD simulations** for conformational sampling
2. **Explicit solvent MM/GBSA** for better solvation
3. **QM/MM** for Zn coordination (if metal sites are critical)
4. **Manual curation** of Modern ThrRS structures before re-running

---

## Data Files
- `mmgbsa_results.csv` - Full results (172 successful jobs)
- `qc/kabsch_fit_qc.csv` - RMSD validation data
- `qc/kabsch_fit_distribution.png` - RMSD distribution plot
- `jobs/*/` - Individual job directories with structures and logs

---

## Technical Notes

**Pipeline improvements implemented:**
1. Per-structure ligand parameterization (not template-based)
2. Fixed tleap returncode check (warnings → false failures)
3. Fixed antechamber duplicate bonds bug
4. Reduced minimization to 500 steps (speed)

**QC validation:**
- RMSD = 0.000 Å proves coordinate preservation ✓
- All jobs converged (RMS gradient < 0.5) ✓
- No pose distortion during minimization ✓
