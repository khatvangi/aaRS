# Sensitivity Sweep Rerun - Status

**Started:** 2025-12-22 20:15 UTC
**Process ID:** 2461830
**Status:** ⏳ RUNNING
**Estimated completion:** ~4-5 hours (based on previous run: 4h 29m)

---

## What's Being Rerun

**Script:** `scripts_out/02_sensitivity_sweep.py`
**Output:** `sweep/effects_stability_table.csv` (will be regenerated)
**Log:** `sweep/sensitivity_sweep_CORRECTED.log`

**Why:** The old `effects_stability_table.csv` contains **incorrect bootstrap CIs** (pooled means instead of difference CIs). The pose gate script was already fixed and rerun successfully, but the main sensitivity sweep needs to be regenerated.

---

## Progress Monitoring

Check if still running:
```bash
ps aux | grep "02_sensitivity_sweep" | grep -v grep
```

Check progress (will show setting X/81):
```bash
tail -50 sweep/sensitivity_sweep_CORRECTED.log
```

Check CPU usage (should be ~100%):
```bash
top -p 2461830
```

---

## When Complete

### 1. Validate Bootstrap CIs

Run the validation script:
```bash
python3 validate_ci_sanity.py sweep/effects_stability_table.csv
```

**Expected output:**
```
Total comparisons: 3240
Bad CI rows (effect not in [ci_low, ci_high]): 0
Percentage bad: 0.00%

✓ ALL CIs CORRECTLY BRACKET EFFECTS
  Analysis is statistically valid for publication

VALIDATION: PASS ✓
```

### 2. Check Key Statistics

**Old (broken) file:**
- Bad CI rows: 2619/3240 (80.8%) ← CIs were pooled means
- Example: effect=-0.125, CI=[2.50, 2.625] (impossible!)

**New (corrected) file should have:**
- Bad CI rows: ~0 (allowing only rare numerical edge cases)
- Example: effect=-0.125, CI=[-0.5, +0.2] (brackets effect ✓)

### 3. Verify FDR Conclusion Unchanged

The "0 effects survive FDR" conclusion should remain **unchanged** because:
- Permutation p-values were always correct
- BH-FDR correction uses p-values, not CIs
- CIs are for reporting only, not hypothesis testing

Quick check:
```bash
python3 -c "import pandas as pd; df=pd.read_csv('sweep/effects_stability_table.csv'); print('Sig at q<0.05:', (df['q'] < 0.05).sum(), '/', len(df))"
```

Expected: `Sig at q<0.05: 0 / 3240`

---

## What Changed

### Before (Broken)
```python
# scripts_out/02_sensitivity_sweep.py (OLD)
effect = np.mean(cog_vals) - np.mean(noncog_vals)
ci_low, ci_high = bootstrap_ci(
    np.concatenate([cog_vals, noncog_vals]),  # ← WRONG: bootstrapping pooled data
    n_bootstrap=10000
)
```

**Result:** CI represents pooled mean (~3-5), not effect difference

### After (Fixed)
```python
# scripts_out/02_sensitivity_sweep.py (CORRECTED)
effect = np.mean(cog_vals) - np.mean(noncog_vals)
ci_low, ci_high = bootstrap_ci_diff(
    cog_vals, noncog_vals,  # ← CORRECT: bootstrapping difference
    n_bootstrap=10000,
    seed=0
)
```

**Result:** CI correctly brackets effect (e.g., effect=-0.5, CI=[-1.2, +0.1])

---

## Files Status

| File | Status | Notes |
|------|--------|-------|
| `scripts_out/02_sensitivity_sweep.py` | ✅ FIXED | Added bootstrap_ci_diff() |
| `scripts_out/03_pose_gate.py` | ✅ FIXED + RERUN | Used bootstrap_ci_diff() |
| `pose_gate/effects_after_gate.csv` | ✅ REGENERATED | CIs correct |
| `sweep/effects_stability_table.csv` | ⏳ REGENERATING | Running now... |
| `validate_ci_sanity.py` | ✅ CREATED | Use to verify output |

---

## Timeline

| Time | Event |
|------|-------|
| 20:15 | Sensitivity sweep started (PID 2461830) |
| ~00:45 | Estimated completion (4-5 hours) |

---

## After Completion - Publication Checklist

Once `effects_stability_table.csv` is regenerated and validated:

✅ Bootstrap CIs are correct (effect bracketed by CI)
✅ Permutation p-values unchanged (were always correct)
✅ BH-FDR q-values unchanged (depend on p, not CI)
✅ "0 effects survive FDR" conclusion unchanged
✅ All 81 parameter settings analyzed
✅ Statistical rigor validated
✅ Ready for publication

---

## Contact Information for Claude Code Session

If the process completes and you want to verify, run:

```bash
# Quick validation
python3 validate_ci_sanity.py sweep/effects_stability_table.csv

# Check FDR result unchanged
python3 -c "
import pandas as pd
df = pd.read_csv('sweep/effects_stability_table.csv')
print(f'Total comparisons: {len(df)}')
print(f'Significant at q<0.05: {(df.q < 0.05).sum()}')
print(f'Min q-value: {df.q.min():.4f}')
"

# Compare file sizes (new should be same size as old)
ls -lh sweep/effects_stability_table.csv
```

Expected output:
- Validation: PASS ✓
- Significant at q<0.05: 0
- Min q-value: ~0.128 (same as before)
- File size: ~402K (same as old version)

---

## Notes

- The sensitivity sweep processes **81 parameter combinations** × ~10 conditions × 4 metrics
- Each comparison requires 10,000 bootstrap iterations + 10,000 permutations
- Total: ~3,240 comparisons requiring ~65 million bootstrap/permutation iterations
- This is computationally intensive but necessary for publication-quality statistics

**Estimated total CPU time:** 4-5 hours on a single core
