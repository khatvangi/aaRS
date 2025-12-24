# Critical Patches Applied - 2025-12-22

## Summary

Applied two critical fixes identified by user review:

1. **Bootstrap CI bug** (02_sensitivity_sweep.py, 03_pose_gate.py)
2. **Metric mislabeling** ("hbonds" → "polar_close_contacts_3p5A")

---

## Patch 1: Bootstrap CI Bug Fix

### Problem

Both `02_sensitivity_sweep.py` and `03_pose_gate.py` were computing bootstrap CI on **pooled data** instead of on the **difference** (cognate − non-cognate).

**Old code:**
```python
effect = np.mean(cog_vals) - np.mean(noncog_vals)
ci_low, ci_high = bootstrap_ci(
    np.concatenate([cog_vals, noncog_vals]),  # WRONG: pooled data
    n_bootstrap=10000
)
```

This returned CI for the **pooled mean** (~3-5 for contacts/atom), not the **effect** (usually -1 to +2).

### Fix Applied

Added correct `bootstrap_ci_diff` function that bootstraps the difference:

```python
def bootstrap_ci_diff(x, y, n_bootstrap=10000, ci_level=0.95, seed=0):
    """
    Bootstrap confidence interval for difference in means (x - y).

    This is the CORRECT way to compute CI for cognate vs non-cognate effects.
    The pooled bootstrap was computing CI for the pooled mean, not the difference.
    """
    rng = np.random.default_rng(seed)
    x = np.asarray(x)
    y = np.asarray(y)

    diffs = np.empty(n_bootstrap, dtype=float)
    for i in range(n_bootstrap):
        xb = rng.choice(x, size=len(x), replace=True)
        yb = rng.choice(y, size=len(y), replace=True)
        diffs[i] = xb.mean() - yb.mean()

    alpha = 1 - ci_level
    lo = np.percentile(diffs, 100 * alpha / 2)
    hi = np.percentile(diffs, 100 * (1 - alpha / 2))

    return lo, hi
```

**New code:**
```python
effect = np.mean(cog_vals) - np.mean(noncog_vals)
ci_low, ci_high = bootstrap_ci_diff(  # CORRECT: bootstrap the difference
    cog_vals, noncog_vals,
    n_bootstrap=10000,
    seed=0
)
```

### Files Modified

1. **scripts_out/02_sensitivity_sweep.py**
   - Added `bootstrap_ci_diff()` function (lines 40-61)
   - Updated call site (line 160-164)

2. **scripts_out/03_pose_gate.py**
   - Imported `bootstrap_ci_diff` from sensitivity_sweep module (line 20)
   - Updated call site (line 215)
   - **Reran to regenerate `pose_gate/effects_after_gate.csv`**

### Impact

**Before fix:** CI values in `effects_after_gate.csv` were ~3-5 (pooled mean)

**After fix:** CI values now bracket the effect:
- Example: effect = -0.683, CI = [-0.978, -0.389] ✓

**Conclusion:** The permutation p-values were always correct (not affected). Only the CI reporting was wrong. Since no effects were borderline significant, the "0 survives FDR" conclusion remains **valid**.

---

## Patch 2: Metric Renaming (hbonds → polar_close_contacts_3p5A)

### Problem

Metric was labeled "hbonds" but was actually **N/O distance-only pairs** (≤3.5 Å):
- ❌ No angle check
- ❌ No explicit hydrogens
- ❌ No donor/acceptor geometry

**This is NOT true H-bonding** and reviewers would correctly flag this as incorrect.

### Fix Applied

Renamed everywhere:
- Function: `compute_hbonds()` → `compute_polar_close_contacts()`
- Column: `hbonds` → `polar_close_contacts_3p5A`
- Labels: "H-bonds" → "Polar close contacts (≤3.5Å)"

Updated docstring:
```python
def compute_polar_close_contacts(protein_atoms, ligand_atoms):
    """
    Count polar close contacts between protein and ligand.

    NOTE: This is NOT true H-bonding (no angle check, no explicit hydrogens).
    Counts N/O atom pairs within 3.5 Å as a proxy for polar interactions.

    Distance cutoff: 3.5 Å
    """
```

### Files Modified

1. **compute_geometry_metrics.py**
   - Renamed function (line 108)
   - Updated docstring (lines 109-116)
   - Updated call site (line 253)
   - Updated return dict (line 264)
   - Updated column order (line 355)
   - Updated print statement (line 383)

2. **analyze_geometry_metrics.py**
   - Updated metrics dict (line 36)
   - Updated all print statements for ProRS (lines 95, 100, 105)
   - Updated all print statements for ThrRS (lines 116, 121, 126)
   - Updated Zn analysis (lines 141, 147)
   - Updated violin plot (lines 221-224)

3. **geometry_qc_report.py** (already correct)
   - Already renamed in this script ✓

### Files Unmodified (Already Correct)

- `geometry_metrics_derived.csv` - already has `polar_close_contacts_3p5A` ✓
- `geometry_metrics_clean.csv` - already has `polar_close_contacts_3p5A` ✓
- All scripts in `scripts_out/` - use correct naming ✓

### Original File (Still Has Old Name)

- `geometry_metrics.csv` - still has `hbonds` (legacy, not used for publication)

---

## Verification

### Bootstrap CI Fix

**Before:**
```
Row 11: effect=-0.683, ci_low=3.2, ci_high=4.8  # WRONG (pooled mean)
```

**After:**
```
Row 11: effect=-0.683, ci_low=-0.978, ci_high=-0.389  # CORRECT (brackets effect)
```

✅ Verified: CI now correctly brackets the effect value

### Metric Renaming

**Before:**
```
Column: hbonds
Label: "H-bonds"
Function: compute_hbonds()
```

**After:**
```
Column: polar_close_contacts_3p5A
Label: "Polar close contacts (≤3.5Å)"
Function: compute_polar_close_contacts()
Docstring: "NOTE: This is NOT true H-bonding..."
```

✅ Verified: All references updated consistently

---

## Statistical Impact

### Bootstrap CI Bug

**Impact on conclusions:** **NONE**

- Permutation p-values were always correct
- BH-FDR correction was applied to p-values, not CIs
- "0 survives FDR" conclusion remains **valid**
- CIs are for reporting only, not hypothesis testing

**What changed:** Confidence interval reporting is now correct for publication

### Metric Renaming

**Impact on conclusions:** **NONE**

- Metric computation unchanged (same formula)
- Only the **name** and **label** changed
- No recomputation needed for existing data
- **Prevents reviewer criticism** about incorrect H-bond definition

---

## Files Regenerated

✅ **pose_gate/effects_after_gate.csv** - regenerated with correct CIs

**NOT regenerated** (renaming only, no computation change):
- geometry_metrics.csv (legacy, not used)
- All analysis outputs use derived/clean files which already have correct names

---

## Recommendation for Paper

### Methods Section

Add this language:

> **Polar close contacts:** We computed N/O/S atom pairs within 3.5 Å as a proxy for polar interactions. Note that these are **not true hydrogen bonds**, as AF3 structures lack explicit hydrogens and we did not apply angle criteria. We refer to these as "polar close contacts" to distinguish from geometric H-bond definitions.

> **Bootstrap confidence intervals:** Effect sizes (cognate − non-cognate) were computed with 95% bootstrap CIs by resampling cognate and non-cognate groups independently (10,000 iterations), then computing the difference distribution. This correctly accounts for variance in both groups.

---

## Conclusion

Both patches are **editorial/statistical correctness fixes**:

1. ✅ Bootstrap CI now reports correct values (effect bracketed by CI)
2. ✅ Metric renamed to prevent reviewer confusion (polar contacts, not H-bonds)

**Biological conclusions unchanged:**
- Zero cognate effects survive FDR (p-values were always correct)
- Ancestral aaRS show promiscuous binding
- Zn structures split 50/50 engaged/floating

**Publication readiness:** ✓ Ready with these fixes applied
