# CORRECTED AF3 SUBSTRATE BINDING ANALYSIS

## Data Verification Results

**USER CONCERN**: Identical ipTM values for ProRS and ThrRS looked suspicious
**FINDING**: The identical values were an ERROR in my initial analysis. Here are the ACTUAL values from raw JSON files:

---

## Complete Corrected Dataset

### LUCA ProRS (Ancestral)
| Substrate | Pocket ipTM | Type |
|-----------|-------------|------|
| Proline   | 0.78        | Cognate |
| Threonine | 0.70        | Non-cognate |
| **Difference** | **11%** | **Mild preference for cognate** |

**Source**:
- `/storage/kiran-stuff/aaRS/phase2/af3_output_full/fulllength_deep_pro/` (chain_pair_iptm[0][2] = 0.78)
- `/storage/kiran-stuff/aaRS/phase2/af3_output_full/fulllength_deep_thr/` (chain_pair_iptm[0][2] = 0.70)

---

### LUCA ThrRS (Ancestral)
| Substrate | Pocket ipTM | Type |
|-----------|-------------|------|
| Proline   | 0.88        | Non-cognate |
| Threonine | 0.89        | Cognate |
| **Difference** | **1%** | **TRUE promiscuity** |

**Source**:
- `/storage/kiran-stuff/aaRS/phase2/outputs/deep_thrrs_pro/` (chain_pair_iptm[0][2] = 0.88)
- `/storage/kiran-stuff/aaRS/phase2/outputs/deep_thrrs_thr/` (chain_pair_iptm[0][2] = 0.89)

---

### Modern E. coli ProRS
| Substrate | Pocket ipTM | Type |
|-----------|-------------|------|
| Proline   | 0.95        | Cognate |
| Threonine | 0.87        | Non-cognate |
| **Difference** | **9%** | **Maintained mild preference** |

**Source**:
- `/storage/kiran-stuff/aaRS/phase2/af3_gaps/af3_output/modern_ecoli_full_pro/` (chain_pair_iptm[0][1] = 0.95)
- `/storage/kiran-stuff/aaRS/phase2/af3_gaps/af3_output/modern_ecoli_full_thr/` (chain_pair_iptm[0][1] = 0.87)

---

### Modern ThrRS
| Substrate | Pocket ipTM | Type |
|-----------|-------------|------|
| Proline   | 0.57        | Non-cognate |
| Threonine | 0.84        | Cognate |
| **Difference** | **47%** | **STRONG specificity evolved!** |

**Source**:
- `/storage/kiran-stuff/aaRS/phase2/outputs/modern_thrrs_pro/` (chain_pair_iptm[0][2] = 0.57)
- `/storage/kiran-stuff/aaRS/phase2/outputs/modern_thrrs_thr/` (chain_pair_iptm[0][2] = 0.84)

---

## Key Findings (CORRECTED)

### 1. Differential Evolution of Specificity

**ProRS trajectory**:
- LUCA: 11% preference → Modern: 9% preference
- **MAINTAINED promiscuity** throughout evolution

**ThrRS trajectory**:
- LUCA: 1% preference → Modern: 47% preference
- **EVOLVED strong specificity** from ancestral promiscuity

### 2. Asymmetric Evolutionary Pressure

- **LUCA ThrRS was truly promiscuous** (1% difference) - the most promiscuous enzyme
- **LUCA ProRS had mild discrimination** (11% difference)
- **Modern ThrRS underwent dramatic specificity evolution** (47% preference)
- **Modern ProRS maintained ancestral promiscuity** (9% preference)

### 3. Why ThrRS Evolved Stronger Specificity

Possible explanations:
1. **Directional selection**: ThrRS faced stronger pressure to avoid Pro misincorporation
2. **Asymmetric consequences**: Thr→Pro errors may be more deleterious than Pro→Thr errors
3. **Editing domain compensates for ProRS**: ProRS has editing domain (3.2× preference for Thr-tRNA), ThrRS may lack this
4. **Structural constraints**: ProRS retained promiscuity due to catalytic/editing domain coordination

---

## Error in Previous Analysis

**CLAIMED (INCORRECT)**:
```
All enzymes show 100% promiscuity with identical ipTM values:
- LUCA ProRS: Pro 0.78 / Thr 0.78
- Modern ProRS: Pro 0.95 / Thr 0.95
```

**ACTUAL (CORRECTED)**:
```
Enzymes show differential specificity:
- LUCA ProRS: Pro 0.78 / Thr 0.70 (11% preference)
- LUCA ThrRS: Pro 0.88 / Thr 0.89 (1% - truly promiscuous)
- Modern ProRS: Pro 0.95 / Thr 0.87 (9% preference)
- Modern ThrRS: Pro 0.57 / Thr 0.84 (47% preference - highly specific)
```

**ROOT CAUSE**: I incorrectly extracted or rounded values in the initial analysis instead of reading the raw JSON files directly.

---

## Biological Interpretation

### The "Promiscuity Paradox" Resolved

**Old interpretation (WRONG)**:
"Evolution did NOT eliminate promiscuity - all enzymes bind both substrates equally"

**New interpretation (CORRECT)**:
"Evolution shows ASYMMETRIC refinement - ThrRS evolved strong specificity (47% preference) while ProRS MAINTAINED ancestral promiscuity (9% preference), likely because ProRS editing domain compensates for catalytic promiscuity"

### Why ProRS Can Afford Promiscuity

1. **Editing domain provides quality control**: 3.2× preference for Thr-tRNA^Pro enables error correction
2. **Two-stage fidelity**: Catalytic promiscuity + editing specificity = overall accuracy
3. **Metabolic advantage**: Pro/Thr structural similarity allowed ancestral promiscuity to persist

### Why ThrRS Needed Specificity

1. **No editing domain** (or less effective editing)
2. **Cannot tolerate Pro misincorporation**
3. **Evolved 47% stronger discrimination** to achieve fidelity

---

## Matrix Indexing Notes

**Modern predictions (2x2 matrix)**:
- Chains: [Protein, Ligand]
- Pocket ipTM = chain_pair_iptm[0][1] (protein-ligand interface)

**LUCA predictions (3x3 matrix)**:
- Chains: [Protein, tRNA, Ligand]
- Pocket ipTM = chain_pair_iptm[0][2] (protein-ligand interface)

---

## Files That Need Correction

1. `/storage/kiran-stuff/aaRS/phase2/generate_figure3_modern_comparison.py` - uses incorrect identical values
2. `/storage/kiran-stuff/aaRS/phase2/af3_gaps/af3_output/comprehensive_af3_results.csv` - may need verification
3. Any manuscript text claiming "100% promiscuity" or identical binding

---

**VERIFIED**: 2025-12-10
**Confidence**: HIGH - directly extracted from raw JSON files
