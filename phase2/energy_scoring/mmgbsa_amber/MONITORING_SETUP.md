# Monitoring Setup - Active

## Current Status (2025-12-21 11:12)

**Main Batch**: 88/182 complete (48%)
- 148/182 minimizations done
- 60 sander processes running
- Using Gasteiger charges

**Validation Jobs**: 2/5 ready
- ✓ anc_thrrs_cat_zn_THR
- ✓ anc_thrrs_cat_zn_SER
- ⏳ modern_thrrs_ecoli_zn_THR
- ⏳ modern_thrrs_ecoli_zn_SER
- ⏳ modern_thrrs_ecoli_THR_zinc

## Active Monitoring

### 1. Auto-Validation Monitor (PID: 1939698)
**Script**: `auto_validate_when_ready.sh`
**Log**: `validation_monitor.log`
**Frequency**: Every 2 minutes
**Action**: Automatically runs `validate_gasteiger_vs_am1bcc.py` when all 5 validation jobs complete

```bash
# Check status:
tail -f validation_monitor.log
```

### 2. Periodic Status Checks (PID: 1940705)
**Script**: `quick_status.sh` (called from background loop)
**Log**: `periodic_checks.log`
**Frequency**: Every 10 minutes (10 checks total)

```bash
# Check status:
tail -f periodic_checks.log
```

## Manual Checks

### Quick Status
```bash
./quick_status.sh
```

### Full Batch Status
```bash
./monitor_batch.sh
```

### Validation Job Status
```bash
for job in anc_thrrs_cat_zn_THR anc_thrrs_cat_zn_SER modern_thrrs_ecoli_zn_THR modern_thrrs_ecoli_zn_SER modern_thrrs_ecoli_THR_zinc; do
    ls -lh jobs/$job/min.rst 2>/dev/null || echo "⏳ $job not ready"
done
```

### Main Batch Progress
```bash
find jobs -name "FINAL_RESULTS_MMPBSA.dat" | wc -l
```

## What Happens When Validation Jobs Complete

The auto-validator will:
1. Detect all 5 validation jobs have `min.rst` files
2. Automatically run `validate_gasteiger_vs_am1bcc.py`
3. Output results to `validation_monitor.log`
4. Create `validation_gasteiger_vs_am1bcc.csv`

You'll see one of these verdicts:

### ✓ Gasteiger Acceptable
```
✓ SIGN MATCHES - Gasteiger and AM1-BCC agree on ranking
  Both predict SER is worse than THR
  ΔΔG difference: 1.2 kcal/mol (15%)

✓ GASTEIGER IS ACCEPTABLE
  Proceed with current batch using Gasteiger charges.
```
**Action**: Let main batch finish, results are trustworthy

### ⚠ Large Difference
```
✓ SIGN MATCHES - Gasteiger and AM1-BCC agree on ranking
  Both predict SER is worse than THR
  ΔΔG difference: 4.5 kcal/mol (45%)

⚠ WARNING: Large ΔΔG difference
  Consider re-running with AM1-BCC for higher accuracy.
```
**Action**: Decision point - qualitative ranking OK, but quantitative ΔΔG may be off

### ✗ Sign Mismatch (ABORT!)
```
✗ SIGN MISMATCH - Gasteiger and AM1-BCC DISAGREE on ranking!
  Gasteiger predicts SER is worse than THR
  AM1-BCC predicts SER is better than THR

✗ GASTEIGER IS UNRELIABLE
  ABORT current batch and re-run with AM1-BCC charges!
```
**Action**:
1. `pkill -f sander` (stop batch)
2. Re-parameterize all 20 amino acids with AM1-BCC
3. Restart batch from scratch

## Expected Timeline

**Best estimate**: ~1-2 hours until validation ready
- Main batch: ~50% complete
- Large proteins take longer (FULL jobs ~20-30 min/job)
- Smaller domains faster (~5-10 min/job)

**Validation jobs** are in the "domain" category:
- `anc_thrrs_cat_zn_*`: ancestral ThrRS catalytic domain (278 residues)
- `modern_thrrs_ecoli_zn_*`: modern ThrRS (401 residues)

These should complete within the next 30-60 minutes.

## Files Generated

### Logs
- `batch.log` - Main batch progress (from run_mmpbsa_batch.py)
- `validation_monitor.log` - Auto-validator log
- `periodic_checks.log` - Periodic status checks

### Results (when validation completes)
- `validation_gasteiger_vs_am1bcc.csv` - Comparison data
- `jobs/*/validate_gasteiger/` - Re-run with Gasteiger params
- `jobs/*/validate_am1bcc/` - Re-run with AM1-BCC params

### Scripts
- `auto_validate_when_ready.sh` - Auto-validator
- `quick_status.sh` - Quick status check
- `monitor_batch.sh` - Full batch monitor
- `validate_gasteiger_vs_am1bcc.py` - Validation script

## Troubleshooting

### If auto-validator stops
```bash
# Check if it's still running:
ps aux | grep auto_validate_when_ready.sh

# Restart if needed:
nohup ./auto_validate_when_ready.sh > validation_monitor.log 2>&1 &
```

### If you need to stop monitoring
```bash
# Kill auto-validator:
pkill -f auto_validate_when_ready.sh

# Kill periodic checks:
pkill -f quick_status.sh
```

### If validation fails
```bash
# Run manually:
python validate_gasteiger_vs_am1bcc.py

# Check for errors:
tail -100 validation_monitor.log
```

## Next Steps After Validation

1. **If Gasteiger is acceptable:**
   - Wait for main batch to complete (88/182 → 182/182)
   - Analyze full results in `mmgbsa_results.csv`
   - Calculate ΔΔG for all substrate pairs
   - Generate plots and comparison tables

2. **If Gasteiger is unreliable:**
   - Stop batch: `pkill -f sander`
   - Re-parameterize all 20 AAs with AM1-BCC (script ready)
   - Clean jobs: `rm -rf jobs && mkdir jobs`
   - Restart: `python run_mmpbsa_batch.py`
   - ETA: Another ~3-4 hours for full batch

3. **If unclear/borderline:**
   - Consult with user
   - Consider running a subset with AM1-BCC for comparison
   - Document uncertainty in results
