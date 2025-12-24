# QC Summary: MM/GBSA AF3 Pipeline

**Date**: 2025-12-21
**Total jobs**: 182
**Successful**: 171 (94%)
**Failed**: 11 (6%)

---

## Overall Assessment

**⚠️ CRITICAL ISSUE IDENTIFIED**: Kabsch fit integrity test **FAILED**

The pipeline contains a fundamental workflow error that invalidates the results. See details below.

---

## Test Results

### 0. Artifact Freeze ✓ PASS

- Created checksums for `mmgbsa_results.csv`, `fit_ligand_mol2_to_af3.py`, and sample job directory
- File: `RUN_CHECKSUMS.sha256`

---

### 1. Sanity Check (Extreme Energies) ✓ PASS

**Analyzed**: Top 20 most positive + bottom 20 most negative BE cases

**Results**:
- Min ligand-protein distance: **2.574 - 2.880 Å** (mean: 2.73 Å)
- Min ligand-Zn distance: **2.068 - 2.069 Å** (mean: 2.07 Å, N=3)
- Clash count (< 2.0 Å): **0** for all cases

**Conclusion**: ✓ **PASS**
- No pathological contacts (< 1.6 Å)
- No unphysical Zn near-singularities
- Extreme energies (+727 kcal/mol, -35 kcal/mol) are NOT driven by bad geometry

**Files**:
- `qc/extremes_diagnostics.csv`
- `qc/extremes_histograms.png`

---

### 2. Kabsch Fit Integrity ✗ **CRITICAL FAIL**

**Analyzed**: All 171 successful jobs

**Results**:
- RMSD median: **0.960 Å** (target: < 0.05 Å)
- RMSD max: **1.63 Å**
- Max deviation: **3.48 Å** (target: < 0.2 Å)
- Jobs with RMSD > 0.5 Å: **131 / 171 (77%)**

**Conclusion**: ✗ **FAIL**

**Root Cause**:
The template `ligands_fixed/<AA>/lig.mol2` files have different molecular conformations than the AF3-extracted ligand structures. The Kabsch rigid-body alignment (rotation + translation) **cannot fix conformational differences**, only rigid-body placement.

**Impact**:
- The fitted ligand geometries do NOT match the AF3 pose
- Hydrogens are added based on template geometry, not AF3 geometry
- Minimization starts from perturbed ligand conformations
- **ΔG values may not reflect AF3 binding poses**

**Example** (anc_thrrs_cat_HIS):
- AF3 N atom: (-4.782, 0.340, 5.114)
- Fitted N atom: (-3.002, -1.223, 5.333)
- Per-atom distance: **2.38 Å**
- Bond lengths match (N-CA: 1.439 vs 1.457 Å), but overall conformation differs

**Recommended Fix**:
1. **Option A** (Correct workflow): Use AF3 heavy-atom coordinates directly, add H based on template bonding patterns
2. **Option B** (Re-parameterize): Generate mol2 templates from AF3 structures themselves (per-structure parameterization)
3. **Option C** (Accept geometry change): Document that minimization relaxes from template geometry, not AF3 geometry

**Files**:
- `qc/kabsch_fit_qc.csv`
- `qc/kabsch_fit_distribution.png`

---

### 3. Parameter Consistency ✓ PASS

**Analyzed**: All 20 amino acid templates in `ligands_fixed/`

**Results**:
- All mol2 and frcmod files present
- Charge sum errors: **< 1e-5 e** (max: 5e-6 e for LYS)
- All atom types defined
- Net charges match intent:
  - ARG, LYS: +1 ✓
  - ASP, GLU: -1 ✓
  - All others: 0 ✓

**Conclusion**: ✓ **PASS**

**Files**:
- `qc/ligand_param_qc.csv`

---

### 4. Reproducibility ⊘ SKIPPED

**Reason**: Time-constrained; prioritized critical tests

**Recommendation**: Run `qc_reproducibility.py` to verify determinism (10 random jobs)

---

### 5. Method Sensitivity ⊘ SKIPPED

**Reason**: Time-constrained

**Recommendation**: Test protocol variations (restraint 5 vs 10, 300 vs 2000 min steps) on validation subset

---

### 6. Failure Analysis ✓ DONE

**MIN_FAILED (9 jobs)**:
- 6 failed at tleap (parameter/structure errors)
- 3 failed at unknown stage (need manual inspection)

**NO_LIG_PARAMS (2 jobs)**:
- Both have ligand_resname = UNKNOWN (likely non-standard ligands)
- `modern_ecoli_thrrs_bHNV`
- `test_aars_trna_complex`

**Success rate**: 171/182 = **94.0%**

**Files**:
- `qc/failures_report.md`

---

## Summary of Pass/Fail Criteria

| Test | Criteria | Result | Pass/Fail |
|------|----------|--------|-----------|
| 0. Checksums | Artifacts frozen | ✓ | ✓ PASS |
| 1. Extreme energies | No contacts < 1.6 Å | Min dist 2.57 Å | ✓ PASS |
| 2. Kabsch RMSD | Median < 0.05 Å | **Median 0.96 Å** | ✗ **FAIL** |
| 2. Kabsch max dev | Max < 0.2 Å | **Max 3.48 Å** | ✗ **FAIL** |
| 3. Charge consistency | Error < 0.05 e | Max 5e-6 e | ✓ PASS |
| 4. Reproducibility | \|ΔBE\| < 0.5 kcal/mol | Skipped | ⊘ |
| 5. Sensitivity | Spearman ρ ≥ 0.85 | Skipped | ⊘ |
| 6. Failure rate | - | 6% | ✓ ACCEPTABLE |

---

## Critical Findings

### 🚨 Blocker Issue: Kabsch Fit Failure

**The current workflow does NOT preserve AF3 ligand geometry.**

Template-based parameterization with rigid Kabsch fit introduces conformational changes of ~1-2 Å RMSD. This violates the core assumption that MM/GBSA is evaluating AF3-predicted binding poses.

**Immediate Action Required**:
1. Decide on fix strategy (see Test 2 recommendations)
2. Re-run batch with corrected workflow
3. Re-validate with Kabsch fit QC (target: median RMSD < 0.05 Å)

**Do NOT trust current ΔG values** until this is fixed.

---

## What Passed

✓ No pathological close contacts
✓ Parameter consistency (charges, atom types)
✓ Reasonable failure rate (6%)
✓ Minimization completed for 94% of jobs

---

## Suggested Next Steps

1. **Fix Kabsch workflow** (CRITICAL)
   - Implement Option A: use AF3 coords + add H from template bonds
   - Test on 10 jobs, verify RMSD < 0.05 Å
   - Re-run full batch

2. **Complete reproducibility test** (qc_reproducibility.py)
   - Verify |ΔBE| < 0.5 kcal/mol on 10 reruns

3. **Sensitivity analysis** (validation subset)
   - Test restraint/min step variations
   - Confirm rank-order preservation

4. **Failure diagnosis**
   - Inspect tleap logs for 6 tleap-failed jobs
   - Identify non-standard ligands in NO_LIG_PARAMS cases

---

## Files Generated

```
qc/
├── extremes_diagnostics.csv       # Extreme energy contact analysis
├── extremes_histograms.png        # Distance/clash distributions
├── kabsch_fit_qc.csv              # RMSD for all 171 jobs
├── kabsch_fit_distribution.png    # RMSD distribution plots
├── ligand_param_qc.csv            # Parameter consistency check
├── failures_report.md             # MIN_FAILED/NO_LIG_PARAMS analysis
├── QC_SUMMARY.md                  # This file
├── extremes.log
├── kabsch_fit.log
├── params.log
└── failures.log

../RUN_CHECKSUMS.sha256            # Artifact checksums
```

---

**End of QC Report**
