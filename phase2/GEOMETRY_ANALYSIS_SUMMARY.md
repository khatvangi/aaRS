# Geometry-Based Analysis Summary

**Replaced MM/GBSA with pure geometry metrics** - no force fields, just physical features AF3 actually predicts.

---

## Metrics Computed (228 structures)

### Primary Metrics
1. **Pocket ipTM** (from AF3) - local confidence in binding site
2. **Contacts (4Å)** - heavy atom protein-ligand contacts
3. **H-bonds** - hydrogen bonds (N/O distance < 3.5 Å)
4. **Polar fraction** - proportion of polar vs nonpolar contacts
5. **Clashes (<2.2Å)** - steric overlap indicating poor fit

### Zn Metrics (49 structures with Zn)
6. **Zn min distance** - closest Zn-ligand distance
7. **Zn coordination count** - atoms within 2.4 Å of Zn
8. **Zn ligand contacts** - ligand atoms near Zn

---

## Key Findings

### 1. Specificity Analysis (Geometry-Based)

#### ProRS: PRO (cognate) vs other amino acids
- **PRO:** pocket ipTM = 0.702, contacts = 40.1, H-bonds = 9.7
- **Other:** pocket ipTM = 0.783, contacts = 29.5, H-bonds = 8.8
- **Result:** PRO binds with LOWER pocket ipTM but MORE contacts
  - Suggests tighter packing but lower local confidence
  - Δ pocket ipTM = -0.081 (PRO worse than non-cognate!)

#### ThrRS: THR (cognate) vs other amino acids
- **THR:** pocket ipTM = 0.848, contacts = 66.7, H-bonds = 27.1
- **Other:** pocket ipTM = 0.817, contacts = 35.1, H-bonds = 9.5
- **Result:** THR binds with HIGHER pocket ipTM and MORE contacts
  - Δ pocket ipTM = +0.031 (THR slightly better)
  - Δ contacts = +31.6 (THR much more engaged)

**Conclusion:** ThrRS shows weak geometry-based selectivity for THR; ProRS does not.

---

### 2. Zinc Effect

| Metric | Without Zn (n=179) | With Zn (n=49) | Difference |
|--------|-------------------|----------------|------------|
| Pocket ipTM | 0.715 ± 0.171 | **0.955 ± 0.041** | **+0.240** ✓ |
| Contacts | 32.4 ± 18.3 | 42.1 ± 44.4 | +9.7 |
| H-bonds | 8.9 ± 5.4 | 14.1 ± 19.2 | +5.2 |
| Clashes | 0.5 ± 1.8 | 8.2 ± 39.6 | +7.7 ✗ |

**Interpretation:**
- Zn presence **increases pocket ipTM** dramatically (+0.24)
- But also **increases clashes** (~16x higher)
- Suggests AF3 models Zn structures with high confidence but poor geometry
- **Consistent with MM/GBSA finding** that Zn systems were problematic

### Zn Coordination Geometry
- Min Zn-ligand distance: **13.85 ± 13.66 Å** (very large!)
- Coordination count: 4.5 ± 3.3 (reasonable tetrahedral)
- Ligand contacts: 1.7 ± 1.8 (low engagement)

**Problem:** Zn is far from ligand (~14 Å on average) suggesting AF3 places Zn correctly for protein but not engaging ligand properly.

---

### 3. Promiscuity Index

**Top 3 Most Promiscuous** (high pocket ipTM across ligands):
1. modern_prors_ALA (ipTM=0.930, 3 ligands)
2. modern_ecoli_full_thr (ipTM=0.875, 4 ligands)
3. shallow_domain_thr (ipTM=0.740, 3 ligands)

**Least Promiscuous:**
- deep_editing_pro (ipTM=0.125, 6 ligands) - very selective/poor binding

---

### 4. Overall Statistics

| Metric | Mean | Std | Range |
|--------|------|-----|-------|
| Pocket ipTM | 0.763 | 0.181 | 0.11 - 0.98 |
| Contacts (4Å) | 34.5 | 26.4 | 0 - 248 |
| H-bonds | 10.0 | 10.3 | 0 - 104 |
| Polar fraction | 0.77 | 0.19 | 0 - 1.0 |
| Clashes | 2.2 | 18.5 | 0 - 198 |

**Note:** High std and max values for clashes (198!) indicate some structures have severe steric problems.

---

## Comparison to MM/GBSA

| Finding | MM/GBSA | Geometry Metrics |
|---------|---------|------------------|
| **ProRS specificity** | No selectivity for PRO | No selectivity for PRO ✓ |
| **ThrRS specificity** | Reverse selectivity | Weak positive selectivity |
| **Zn effect** | Catastrophic failures (+862 kcal/mol) | High confidence but clashes |
| **Method reliability** | Force field errors dominate | Geometry errors visible but interpretable |

**Agreement:** Both methods show lack of strong cognate selectivity.

**Difference:** Geometry metrics identify the problem as **AF3 structural issues** (clashes, poor Zn placement) rather than force field artifacts.

---

## Recommendations

### For the paper (IMMEDIATE)

1. **Use geometry metrics as primary signal:**
   - Pocket ipTM (AF3 native metric)
   - Contacts and H-bonds (shape complementarity)
   - Polar fraction (interaction character)

2. **Abandon MM/GBSA entirely**
   - Keep as "tested and rejected" in methods
   - One sentence: "MM/GBSA was unstable due to Zn parameterization issues"

3. **Focus interpretation on:**
   - **Ancestral systems** (most structures, cleaner)
   - **Non-Zn systems** (fewer clashes, more reliable)
   - **Relative trends** (pocket ipTM rankings) not absolute values

4. **Flag Zn structures as AF3 modeling issue:**
   - High pocket ipTM despite poor geometry suggests AF3 overconfident
   - Zn ~14 Å from ligand = not actually coordinating
   - Recommend: "Zn coordination requires QM/MM or manual curation"

### For immediate analysis

1. **Filter out high-clash structures** (clashes > 10) for clean heatmaps
2. **Separate Zn vs no-Zn analyses** completely
3. **Create ligand preference rankings** using pocket ipTM
4. **Compute promiscuity scores** for each enzyme variant

---

## Output Files

- `geometry_metrics.csv` - **Main data file** (228 structures, 23 columns)
- `heatmap_geometry_metrics.png` - Enzyme × ligand heatmaps for 5 metrics
- `comparison_plots_geometry.png` - 4-panel comparison plots

---

## Next Steps

1. ✓ Replace MM/GBSA analysis with geometry metrics
2. ✓ Generate heatmaps and visualizations
3. **TODO:** Filter clashes > 10, re-run heatmaps
4. **TODO:** Create promiscuity rankings per enzyme
5. **TODO:** Compute statistical significance (t-tests for specificity)
6. **TODO:** Generate publication-ready figures

---

## Technical Notes

**Why this is better than MM/GBSA:**
- No force field parameterization errors
- No Zn coordination chemistry failures  
- Metrics AF3 was actually trained to predict
- Directly interpretable physical quantities
- Scales to 1000+ structures trivially

**Limitations:**
- No energetics (can't predict ΔG)
- No dynamics (single snapshot)
- Clash count doesn't distinguish severity
- H-bond definition is simplified (no angles)

**But:** For ranking ligand preferences and identifying promiscuity patterns, geometry metrics are sufficient and defensible.
