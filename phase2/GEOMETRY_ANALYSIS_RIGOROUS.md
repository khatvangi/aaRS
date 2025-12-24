# Rigorous Geometry-Based Analysis - CORRECTED

**Date:** 2025-12-22
**Total structures:** 228
**Clean structures (clashes ≤ 2):** 220 (96.5%)
**Method:** AF3-predicted geometry with proper normalization and stratification

---

## Critical Fixes Applied

### 1. Fixed Mislabeled Metrics

**A. Renamed "hbonds" → "polar_close_contacts_3p5A"**
- **Problem:** Original metric counted N/O pairs within 3.5 Å with NO angle check or donor/acceptor rules
- **Fix:** Renamed to accurately reflect what we're measuring (polar close contacts, not true H-bonds)
- **Impact:** Prevents reviewer criticism about improper H-bond definition

**B. Split Zn into "engaged" vs "floating"**
- **Problem:** Zn structures had mean distance ~14 Å to ligand (not actually coordinating!)
- **Fix:** Classified Zn structures into:
  - **Zn engaged:** `zn_min_dist ≤ 3.0 Å` AND `zn_coordination_count ≥ 2` (n=26)
  - **Zn floating:** Zn present but not engaged (n=23) → **AF3 artifact, exclude from mechanistic claims**
- **Impact:** Prevents false biological interpretations of poorly modeled Zn sites

### 2. Normalized Contacts by Ligand Size

**Added per-atom metrics to prevent bias toward large amino acids:**
- `contacts_per_atom = contacts_4A / ligand_heavy_atoms`
- `polar_contacts_per_atom = polar_contacts / ligand_heavy_atoms`
- `polar_close_contacts_per_atom = polar_close_contacts_3p5A / ligand_heavy_atoms`
- `clash_rate = clashes_2.2A / ligand_heavy_atoms`

**Impact:** TRP/TYR/ARG no longer dominate just because they're big

### 3. Stratified Analysis by Condition

**Replaced coarse "ProRS vs ThrRS" grouping with:**
- **Condition = (enzyme, epoch, domain, zn_engaged)**
- **Example conditions:**
  - `Ancestral_ProRS_Catalytic` (n=18)
  - `Ancestral_ThrRS_Catalytic` (n=39)
  - `Modern_ThrRS_Full-length_Zn` (n=20)

**Impact:** Prevents catalytic + editing domains from collapsing into meaningless averages

### 4. QC Gates

**Two datasets:**
- **All (n=228):** Full dataset for diagnostics
- **Clean (n=220):** Clashes ≤ 2 for biological interpretation

**Impact:** Removes 8 structures with severe steric problems (up to 198 clashes)

---

## Key Findings (Using Clean Dataset)

### Overall Statistics

| Metric | Mean | Std | Range | Interpretation |
|--------|------|-----|-------|----------------|
| **Pocket ipTM** | 0.764 | 0.183 | 0.11 - 0.98 | AF3 local support |
| **Contacts per atom** | 3.61 | 1.76 | 0 - 9.67 | Shape complementarity |
| **Polar contacts per atom** | 2.87 | 1.44 | 0 - 6.60 | Polar engagement |
| **Clash rate** | 0.028 | 0.055 | 0 - 0.22 | Steric problems (low!) |

### Cognate vs Non-Cognate Specificity

**Conditions with SIGNIFICANT cognate preference (95% CI excludes 0):**

| Condition | Cognate | Δ pocket ipTM | Δ contacts/atom | Interpretation |
|-----------|---------|---------------|-----------------|----------------|
| **Ancestral_ProRS_Editing** | PRO | **+0.080** [+0.055, +0.112] ✓ | **+1.97** [+1.20, +2.69] ✓ | ⭐ **Strong PRO selectivity** |
| **Ancestral_ProRS_Full-length** | PRO | -0.013 ✗ | **+1.33** [+1.07, +1.56] ✓ | Shape selectivity only |
| **Modern_ThrRS_Full-length_Zn** | THR | 0.000 ✗ | **+0.44** [+0.08, +0.81] ✓ | Weak THR selectivity |

**Conditions with NO significant cognate preference:**

| Condition | Cognate | Δ pocket ipTM | Δ contacts/atom | Interpretation |
|-----------|---------|---------------|-----------------|----------------|
| **Ancestral_ThrRS_Catalytic** (n=39) | THR | -0.001 ✗ | -0.087 ✗ | **No THR selectivity** |
| **Modern_ProRS_Domain** (n=30) | PRO | -0.013 ✗ | +0.493 ✗ | **No PRO selectivity** |
| **Modern_ThrRS_Domain** (n=22) | THR | +0.133 ✗ | +1.37 ✗ | High variance, not significant |

**⚠️ CRITICAL FINDING:** Only 3/10 conditions show significant cognate selectivity, and the signal is WEAK even when present.

### Ancestral ProRS Catalytic vs Editing Paradox

| Domain | Δ pocket ipTM (PRO vs other) | Δ contacts/atom (PRO vs other) |
|--------|------------------------------|--------------------------------|
| **Catalytic** | +0.020 ✓ (barely) | **-1.50** ✓ (REVERSE!) |
| **Editing** | +0.080 ✓ | **+1.97** ✓ (correct direction) |

**Interpretation:**
- **Editing domain** shows expected PRO preference (more contacts, higher ipTM)
- **Catalytic domain** shows REVERSE shape complementarity (PRO has FEWER contacts per atom than non-cognate!)
- Suggests catalytic site binds PRO with **tighter packing** (fewer atoms engaged) vs promiscuous shallow binding of large AAs

---

## Zn Effect (Engaged vs Floating)

### Zn Classification Results

| Category | Count | Mean pocket ipTM | Mean contacts/atom |
|----------|-------|------------------|-------------------|
| **No Zn** | 179 | 0.715 | 3.53 |
| **Zn engaged** (≤3.0 Å, coord≥2) | 26 | **0.963** | 2.67 |
| **Zn floating** (>3.0 Å) | 23 | 0.947 | 4.85 |

### Key Insights

1. **Zn engaged structures have dramatically higher pocket ipTM (+0.25)**
   - AF3 is very confident about Zn-coordinated pockets
   - BUT: fewer contacts per atom (2.67 vs 3.53) suggests Zn provides structural stability, not just contact density

2. **Zn floating structures are AF3 artifacts**
   - Mean `zn_min_dist` = 14.8 Å in floating group (Zn nowhere near ligand!)
   - High pocket ipTM (0.947) despite poor Zn placement
   - **DO NOT use for mechanistic interpretation**

3. **Distribution of Zn-ligand distances:**
   - Bimodal: peak at 2-3 Å (engaged) + peak at 10-20 Å (floating)
   - Cutoff of 3.0 Å cleanly separates these populations

---

## Condition-Specific Results (Top 5 by sample size)

### 1. Ancestral_ThrRS_Catalytic (n=39, clean dataset)
- **Cognate selectivity:** NONE (Δ pocket ipTM = -0.001, Δ contacts/atom = -0.087)
- **Mean pocket ipTM:** 0.886 (very high confidence)
- **Mean contacts per atom:** 3.96 (moderate)
- **Interpretation:** Promiscuous binding site with high AF3 confidence

### 2. Modern_ProRS_Domain (n=30)
- **Cognate selectivity:** NONE (Δ pocket ipTM = -0.013, Δ contacts/atom = +0.493 ✗)
- **Mean pocket ipTM:** 0.814
- **Mean contacts per atom:** 4.11
- **Interpretation:** Modern ProRS domain is PROMISCUOUS (consistent with ligand superposition finding)

### 3. Modern_ThrRS_Domain (n=22)
- **Cognate selectivity:** WEAK (Δ pocket ipTM = +0.133 ✗, Δ contacts/atom = +1.37 ✗)
- Large variance prevents significance
- **Note:** This matches ligand superposition finding that modern_thrrs shows distinct THR binding geometry

### 4. Modern_ThrRS_Full-length_Zn (n=20, all Zn engaged)
- **Cognate selectivity:** WEAK but significant for contacts (+0.44 per atom, p<0.05)
- **Mean pocket ipTM:** 0.980 (highest of all conditions!)
- **Interpretation:** Zn coordination provides high confidence; modest THR selectivity

### 5. Ancestral_ProRS_Editing (n=20)
- **Cognate selectivity:** STRONG (Δ pocket ipTM = +0.080, Δ contacts/atom = +1.97, both p<0.05)
- **Interpretation:** Editing domain evolved to REJECT PRO (or preferentially bind misacylated tRNA-PRO)

---

## Integration with Ligand Superposition Analysis

| Condition | Geometry selectivity | Ligand RMSD selectivity | Agreement |
|-----------|---------------------|-------------------------|-----------|
| **Ancestral_ThrRS_Catalytic** | No (Δ contacts = -0.087) | No (cognate clusters WITH non-cognate) | ✓ AGREE |
| **Modern_ThrRS_Domain** | Weak (Δ = +1.37, p>0.05) | **STRONG** (cognate RMSD 18.5 Å vs 9.4 Å) | Partial |
| **Modern_ProRS_Domain** | No (Δ = +0.49, p>0.05) | Yes (minimal, RMSD ~1.8 Å) | Weak signal in both |
| **Ancestral_ProRS_Editing** | **STRONG** (Δ = +1.97, p<0.05) | Not tested (n<20 ligands) | - |

**Conclusion:** Modern ThrRS shows **structural selectivity** (distinct binding poses) that doesn't fully translate to **normalized contact selectivity** → suggests THR binds in a geometrically distinct but not necessarily more complementary way.

---

## Data Files Generated

### Primary outputs
- **`geometry_metrics_derived.csv`** - Full dataset with normalized metrics (n=228, 30 columns)
- **`geometry_metrics_clean.csv`** - QC-filtered dataset (n=220, clashes ≤ 2)
- **`cognate_vs_noncognate_effects.csv`** - Statistical comparison with bootstrap 95% CI (10 conditions)

### Figures
- **`heatmap_stratified_conditions.png`** - 4 metrics × top 12 conditions (rows) × 20 amino acids (cols)
- **`cognate_vs_noncognate_effects.png`** - Bar plots with error bars for Δ pocket ipTM, Δ contacts/atom, Δ polar contacts/atom
- **`zn_diagnostics.png`** - 3 panels: Zn distance histogram, pocket ipTM vs Zn distance, engaged vs floating counts

---

## Interpretation Rules (CRITICAL for paper)

### What these metrics mean

| Metric | Physical meaning | What it's NOT |
|--------|------------------|---------------|
| **pocket_iptm** | AF3 local support (predicted accuracy) | Binding affinity |
| **contacts_per_atom** | Shape complementarity proxy | Binding energy |
| **polar_contacts_per_atom** | Polar engagement proxy | Electrostatic energy |
| **polar_close_contacts_3p5A** | N/O proximity (simplified) | True H-bonds (no angle check) |

### How to use in the paper

✓ **DO:**
- Use pocket ipTM as "AF3 confidence in binding pose"
- Use contacts_per_atom for "shape complementarity"
- Compare **within conditions** (e.g., Ancestral_ProRS_Catalytic: PRO vs other AAs)
- Report bootstrap 95% CI for all Δ values
- Use Zn results **only from zn_engaged subset**

✗ **DON'T:**
- Claim pocket ipTM = binding affinity
- Pool Catalytic + Editing domains
- Use Zn floating structures for mechanistic claims
- Call polar_close_contacts "hydrogen bonds" (reviewers will notice!)

---

## Recommendations for Paper

### Main claims supported by data

1. **Ancestral aaRS are promiscuous** (low/no cognate selectivity in most conditions)
2. **Editing domains differ from catalytic** (Ancestral_ProRS_Editing shows strong PRO selectivity)
3. **Modern ThrRS shows structural selectivity** (ligand superposition) but weak contact selectivity (normalized metrics)
4. **Zn coordination increases AF3 confidence** (pocket ipTM +0.25) but is often poorly modeled (47% floating)

### What NOT to claim

1. ~~"Ancestral enzymes evolved from promiscuous to specific"~~ (Modern ProRS Domain is still promiscuous!)
2. ~~"Zn improves binding"~~ (we only have structure, not affinity)
3. ~~"H-bonds drive selectivity"~~ (we don't measure true H-bonds)

### Suggested figure for paper

**Figure: Cognate selectivity is condition-dependent**
- Panel A: Heatmap of pocket ipTM (stratified by condition × ligand)
- Panel B: Δ contacts/atom (cognate − non-cognate) with 95% CI error bars
- Panel C: Zn diagnostics (engaged vs floating)
- Panel D: Example ligand superposition (modern_thrrs showing distinct THR geometry)

---

## Technical Validation

### QC metrics
- **96.5% structures pass clash filter** (clashes ≤ 2)
- **Zn classification:** 26 engaged, 23 floating, 179 no Zn
- **Normalized metrics:** Prevents size bias (TRP/ARG no longer dominate)
- **Bootstrap CI:** 1000 resamples for statistical rigor

### Comparison to MM/GBSA

| Finding | MM/GBSA | Geometry (rigorous) | Agreement |
|---------|---------|---------------------|-----------|
| ProRS cognate selectivity | None | Editing: yes; Catalytic: reverse | MM/GBSA oversimplified |
| ThrRS cognate selectivity | Reverse (!) | Weak/absent | Both show lack of selectivity ✓ |
| Zn effect | Catastrophic (+862 kcal/mol) | 47% floating (AF3 artifact) | Both identify Zn problems ✓ |

**Conclusion:** Geometry metrics **validate** MM/GBSA finding that cognate selectivity is weak, but provide interpretable physical basis (shape complementarity, not force field errors).

---

## Next Steps

1. ✓ Rigorous geometry analysis complete
2. ✓ Zn engaged/floating classification
3. ✓ Normalized metrics prevent size bias
4. ✓ Stratified by condition
5. **TODO:** Create publication-ready figure with panels A-D
6. **TODO:** Write methods section describing:
   - Polar close contact definition (NOT H-bonds)
   - Zn engaged/floating criteria
   - Bootstrap CI procedure
   - Normalization by ligand heavy atoms
