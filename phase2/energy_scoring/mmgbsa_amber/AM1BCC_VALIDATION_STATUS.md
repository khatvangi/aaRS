# AM1-BCC vs Gasteiger Validation Status

## Critical Correction Applied

**Previous diagnosis (INCORRECT):**
> "AM1-BCC fails for zwitterions because PDB lacks charge assignments"

**Correct diagnosis:**
AM1-BCC failed because:
1. **Missing `-nc 0` parameter** - must explicitly specify net charge
2. **No hydrogens in input PDB** - sqm was running on structure with odd number of electrons (55e)

**Solution implemented:**
```bash
obabel ligand.pdb -O ligand_H.pdb -h -p 7.0  # Add hydrogens first (→ 64e, even!)
antechamber -i ligand_H.pdb -c bcc -nc 0     # Then run AM1-BCC with correct net charge
```

## Current Status

### ✓ Completed

1. **AM1-BCC parameterization of THR and SER** ✓
   - Pre-added hydrogens with `obabel`
   - Ran antechamber with `-c bcc -nc 0`
   - Generated:
     - `ligands/THR/lig_am1bcc.mol2` (AM1-BCC charges)
     - `ligands/THR/lig_am1bcc.frcmod`
     - `ligands/SER/lig_am1bcc.mol2`
     - `ligands/SER/lig_am1bcc.frcmod`

2. **Gasteiger parameterization of all 20 amino acids** ✓
   - All with hydrogens added (`-j 5`)
   - Files: `ligands/*/lig.mol2` and `lig.frcmod`

3. **Main MM/GBSA batch running** ✓
   - Progress: **88/182 complete (48%)**
   - Using Gasteiger charges
   - Proper restraints (protein+Zn restrained, ligand free)
   - Geometry drift tracking enabled

### ⏳ Validation Script Ready (Waiting for Data)

**Script:** `validate_gasteiger_vs_am1bcc.py`

**What it does:**
1. Takes 5 ThrRS+Zn validation jobs:
   - `anc_thrrs_cat_zn_THR` (THR, cognate)
   - `anc_thrrs_cat_zn_SER` (SER, similar)
   - `modern_thrrs_ecoli_zn_THR` (THR)
   - `modern_thrrs_ecoli_zn_SER` (SER)
   - `modern_thrrs_ecoli_THR_zinc` (THR)

2. Runs MM/GBSA with **both** charge sets:
   - Gasteiger charges (`lig.mol2`)
   - AM1-BCC charges (`lig_am1bcc.mol2`)

3. Compares **ΔΔG (SER vs THR)**:
   - If **sign matches** → Gasteiger is acceptable
   - If **sign differs** → ABORT batch, redo with AM1-BCC

4. Outputs:
   - Table of ΔG values for each job × charge method
   - ΔΔG comparison
   - **VERDICT**: Accept Gasteiger or Abort

**Status:**
- ⏳ Waiting for validation jobs to complete in main batch
- Check with: `ls jobs/anc_thrrs_cat_zn_THR/min.rst`

## When to Run Validation

**Option 1: Wait for full batch completion** (safest)
- Main batch finishes (~1-2 hours from now)
- All 5 validation jobs will be complete
- Run: `python validate_gasteiger_vs_am1bcc.py`

**Option 2: Run early (if validation jobs ready)**
- Check if these 5 jobs have `min.rst` files:
  ```bash
  ls jobs/anc_thrrs_cat_zn_THR/min.rst
  ls jobs/anc_thrrs_cat_zn_SER/min.rst
  ls jobs/modern_thrrs_ecoli_zn_THR/min.rst
  ls jobs/modern_thrrs_ecoli_zn_SER/min.rst
  ls jobs/modern_thrrs_ecoli_THR_zinc/min.rst
  ```
- If all exist: `python validate_gasteiger_vs_am1bcc.py`
- Validation takes ~15-20 min (re-running MM/GBSA with AM1-BCC params)

## Possible Outcomes

### Outcome 1: Gasteiger is acceptable ✓
```
✓ SIGN MATCHES - Gasteiger and AM1-BCC agree on ranking
  Both predict SER is worse than THR
  ΔΔG difference: 1.2 kcal/mol (15%)

✓ GASTEIGER IS ACCEPTABLE
  Proceed with current batch using Gasteiger charges.
```

**Action:** Let batch finish, use results as-is

---

### Outcome 2: Gasteiger is questionable ⚠
```
✓ SIGN MATCHES - Gasteiger and AM1-BCC agree on ranking
  Both predict SER is worse than THR
  ΔΔG difference: 4.5 kcal/mol (45%)

⚠ WARNING: Large ΔΔG difference
  Consider re-running with AM1-BCC for higher accuracy.
```

**Action:**
- Could use current batch for **qualitative** ranking
- But should re-run with AM1-BCC for **quantitative** ΔΔG values

---

### Outcome 3: Gasteiger is unreliable ✗
```
✗ SIGN MISMATCH - Gasteiger and AM1-BCC DISAGREE on ranking!
  Gasteiger predicts SER is worse than THR
  AM1-BCC predicts SER is better than THR

✗ GASTEIGER IS UNRELIABLE
  ABORT current batch and re-run with AM1-BCC charges!
```

**Action:**
1. Kill main batch: `pkill -f sander`
2. Re-parameterize all 20 amino acids with AM1-BCC:
   - Update `parameterize_ligands.py` to use obabel + AM1-BCC
   - Or run AM1-BCC version for each AA individually
3. Restart batch from scratch

## Files Created

### AM1-BCC Parameters
- `ligands/THR/lig_am1bcc.mol2` (AM1-BCC charges, 17 atoms with H)
- `ligands/THR/lig_am1bcc.frcmod`
- `ligands/SER/lig_am1bcc.mol2` (AM1-BCC charges, 14 atoms with H)
- `ligands/SER/lig_am1bcc.frcmod`

### Gasteiger Parameters (Current Batch)
- `ligands/*/lig.mol2` (Gasteiger charges, all with H)
- `ligands/*/lig.frcmod`

### Scripts
- `parameterize_am1bcc_correct.py` - Re-parameterize THR/SER with AM1-BCC
- `validate_gasteiger_vs_am1bcc.py` - Run validation comparison
- `monitor_batch.sh` - Monitor main batch progress

### Documentation
- `CORRECTIONS_APPLIED.md` - Full pipeline corrections
- `AM1BCC_VALIDATION_STATUS.md` (this file)

## Next Steps

1. **Monitor main batch**: `./monitor_batch.sh`
   - Currently: 88/182 complete (48%)
   - ETA: ~1-2 hours

2. **When validation jobs complete** (check with `ls jobs/anc_thrrs_cat_zn_THR/min.rst`):
   ```bash
   python validate_gasteiger_vs_am1bcc.py
   ```

3. **Review validation results:**
   - Read output and verdict
   - Check `validation_gasteiger_vs_am1bcc.csv`

4. **Decision:**
   - ✓ Gasteiger acceptable → Continue with main batch
   - ⚠ Large difference → Consider re-running with AM1-BCC
   - ✗ Sign mismatch → ABORT and re-run with AM1-BCC

## Technical Details

### Why AM1-BCC Initially Failed

Original error:
```
QMMM: System specified with odd number of electrons (39)
QMMM: but odd spin (1). You most likely have the charge of
QMMM: QM region (qmcharge) set incorrectly.
```

Root cause:
- PDB with 8 heavy atoms (no H) → 55 electrons (odd!)
- sqm defaults to spin=1 (singlet) → requires EVEN electrons
- `-nc 0` was not enough - needed to add H first

Fix:
- `obabel -h` adds hydrogens → 17 atoms → 64 electrons (even!)
- Then `antechamber -c bcc -nc 0` works correctly

### Charge Assignment Methods

**Gasteiger:**
- Fast empirical method based on electronegativity
- Good for relative comparisons
- Can be less accurate for charged/polar groups
- Used in current batch (all 182 jobs)

**AM1-BCC:**
- QM-derived charges (semi-empirical AM1 + bond charge correction)
- More accurate, especially for polar/charged groups
- Slower (~2 min per molecule with sqm)
- Used for validation (THR and SER only)

### Validation Strategy

Compare **ΔΔG sign** rather than absolute values because:
- We care about **ranking** (THR > SER > ILE) not exact kcal/mol
- Single-snapshot MM/GBSA has ~2-3 kcal/mol uncertainty anyway
- Sign mismatch = wrong ranking = fundamentally wrong conclusion
- Sign match + reasonable magnitude = reliable ranking
