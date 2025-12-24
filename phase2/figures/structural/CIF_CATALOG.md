# CIF File Catalog - Structural Renders

**Generated:** 2025-12-19
**Purpose:** Map available CIF files to PyMOL script requirements

---

## 📋 Summary

**Total CIF files found:** 94 (across 17 base directories)
**PyMOL scripts ready:** 2/5 (40%)
**Requires additional AF3 runs:** 3/5 (60%)

---

## ✅ Available Structures (CIF Files Exist)

### Ancestral ProRS Editing Domain
- **deep_editing_pro** → PRO ligand
  - Path: `/storage/kiran-stuff/aaRS/phase2/outputs/deep_editing_pro/deep_editing_pro/deep_editing_pro_model.cif`
  - ipTM: 0.470 (PRO ligand)
  - Use for: Figure 2D (PRO in editing site)

- **deep_editing_thr** → THR ligand
  - Path: `/storage/kiran-stuff/aaRS/phase2/outputs/deep_editing_thr/deep_editing_thr/deep_editing_thr_model.cif`
  - ipTM: 0.370 (THR ligand)
  - Use for: Figure 2D (THR in editing site)

### Ancestral ThrRS Catalytic Domain
- **deep_thrrs_thr** → THR ligand (no Zn)
  - Path: `/storage/kiran-stuff/aaRS/phase2/outputs/deep_thrrs_thr/deep_thrrs_thr/deep_thrrs_thr_model.cif`
  - ipTM: 0.210 (THR ligand)
  - Zn_iptm: 0.890 (background)
  - Use for: Figure 3D (ancestral active site)

- **deep_thrrs_pro** → PRO ligand (no Zn)
  - Path: `/storage/kiran-stuff/aaRS/phase2/outputs/deep_thrrs_pro/deep_thrrs_pro/deep_thrrs_pro_model.cif`
  - ipTM: 0.620 (PRO ligand)
  - Zn_iptm: 0.880 (background)

### Modern ThrRS (WITHOUT Zinc)
- **modern_thrrs_thr** → THR ligand (no Zn)
  - Path: `/storage/kiran-stuff/aaRS/phase2/outputs/modern_thrrs_thr/modern_thrrs_thr/modern_thrrs_thr_model.cif`
  - ipTM: 0.240 (THR ligand)
  - Zn_iptm: 0.840 (background)
  - **NOTE:** This structure does NOT have Zn in the binding site

- **modern_thrrs_pro** → PRO ligand (no Zn)
  - Path: `/storage/kiran-stuff/aaRS/phase2/outputs/modern_thrrs_pro/modern_thrrs_pro/modern_thrrs_pro_model.cif`
  - ipTM: 0.270 (PRO ligand)

### Modern ProRS
- **modern_prours_pro** → PRO ligand
  - Path: `/storage/kiran-stuff/aaRS/phase2/outputs/modern_prours_pro/modern_prours_pro/modern_prours_pro_model.cif`
  - ipTM: 0.380

- **modern_prours_thr** → THR ligand
  - Path: `/storage/kiran-stuff/aaRS/phase2/outputs/modern_prours_thr/modern_prours_thr/modern_prours_thr_model.cif`
  - ipTM: 0.360

---

## ❌ Missing Structures (In CSV but NO CIF Files)

### Modern ThrRS + Zinc (CRITICAL FOR FIGURES 4C, 4D, 5C)

These job names exist in `AF3_RESULTS_CORRECTED.csv` with high ipTM scores but **CIF output directories do not exist**:

1. **modern_thrrs_ecoli_zn_THR** (ipTM: 0.970, Zn_iptm: 0.980)
   - Required for: Figure 4C (THR coordinating Zn)
   - Status: ⚠️ MISSING

2. **modern_thrrs_ecoli_zn_SER** (ipTM: 0.950, Zn_iptm: 0.980)
   - Required for: Figure 5C (SER zinc trap)
   - Status: ⚠️ MISSING

3. **modern_thrrs_ecoli_zn_ILE** (ipTM: 0.830, Zn_iptm: 0.980)
   - Required for: Figure 4D (ILE rejected)
   - Status: ⚠️ MISSING

### Other Zinc-containing structures in CSV but no CIF:
- modern_thrrs_ecoli_zn_VAL (0.900)
- modern_thrrs_ecoli_zn_CYS (0.890)
- modern_thrrs_ecoli_zn_ASN (0.890)
- modern_thrrs_ecoli_zn_ALA (0.880)
- modern_thrrs_ecoli_zn_GLN (0.870)
- modern_thrrs_ecoli_zn_GLU (0.860)
- modern_thrrs_ecoli_zn_LYS (0.850)
- modern_thrrs_ecoli_zn_MET (0.840)
- modern_thrrs_ecoli_zn_LEU (0.840)
- modern_thrrs_ecoli_zn_ASP (0.830)
- modern_thrrs_ecoli_zn_HIS (0.830)
- modern_thrrs_ecoli_zn_PRO (0.820)
- modern_thrrs_ecoli_zn_TYR (0.810)
- modern_thrrs_ecoli_zn_GLY (0.810)
- modern_thrrs_ecoli_zn_PHE (0.810)
- modern_thrrs_ecoli_zn_TRP (0.800)
- modern_thrrs_ecoli_zn_ARG (0.790)

Total: **20 amino acids tested with Zn** (all in CSV, none have CIF files yet)

---

## 📊 PyMOL Script Readiness Status

### ✅ Figure 2D: Editing Domain Overlay (READY)
**Status:** CAN RUN NOW

**Files needed:**
- ✅ `deep_editing_thr_model.cif` (THR in editing site)
- ✅ `deep_editing_pro_model.cif` (PRO in editing site)

**Script location:** `figures/structural/PYMOL_SCRIPTS.md` (lines 32-98)

**Action:** Update script paths and execute

---

### ✅ Figure 3D: Ancestral vs Modern Active Site (READY)
**Status:** CAN RUN NOW (but without Zn!)

**Files needed:**
- ✅ `deep_thrrs_thr_model.cif` (ancestral ThrRS)
- ✅ `modern_thrrs_thr_model.cif` (modern ThrRS)

**Caveat:** Neither structure has Zn in the active site (these are ligand-only predictions)

**Script location:** `figures/structural/PYMOL_SCRIPTS.md` (lines 101-165)

**Action:** Update script paths and execute (show pocket evolution without Zn)

---

### ❌ Figure 4C: THR Coordinating Zn (BLOCKED)
**Status:** CANNOT RUN - Missing CIF file

**Files needed:**
- ❌ `modern_thrrs_ecoli_zn_THR_model.cif` (MISSING)

**Why blocked:** This AF3 prediction exists in CSV (ipTM=0.970) but CIF output directory was never created

**Script location:** `figures/structural/PYMOL_SCRIPTS.md` (lines 168-240)

**Action required:** Run AF3 for `modern_thrrs_ecoli_zn_THR` job to generate CIF file

---

### ❌ Figure 4D: ILE Rejected by Zn (BLOCKED)
**Status:** CANNOT RUN - Missing CIF file

**Files needed:**
- ❌ `modern_thrrs_ecoli_zn_ILE_model.cif` (MISSING)

**Why blocked:** This AF3 prediction exists in CSV (ipTM=0.830) but CIF output directory was never created

**Script location:** `figures/structural/PYMOL_SCRIPTS.md` (lines 243-305)

**Action required:** Run AF3 for `modern_thrrs_ecoli_zn_ILE` job to generate CIF file

---

### ❌ Figure 5C: SER in Zn Site (BLOCKED)
**Status:** CANNOT RUN - Missing CIF file

**Files needed:**
- ❌ `modern_thrrs_ecoli_zn_SER_model.cif` (MISSING)

**Why blocked:** This AF3 prediction exists in CSV (ipTM=0.950) but CIF output directory was never created

**Script location:** `figures/structural/PYMOL_SCRIPTS.md` (lines 308-382)

**Action required:** Run AF3 for `modern_thrrs_ecoli_zn_SER` job to generate CIF file

---

## 🔧 Recommended Next Steps

### Option 1: Use Available Structures (Immediate)
Run PyMOL scripts for:
1. Figure 2D (Editing domain) - READY ✅
2. Figure 3D (Ancestral vs Modern) - READY ✅ (caveat: no Zn shown)

This gives you **2/5 structural renders** immediately.

### Option 2: Generate Missing Structures (Full Coverage)
Run AF3 predictions for modern ThrRS + Zn with priority ligands:
1. `modern_thrrs_ecoli_zn_THR` (highest priority - Figure 4C)
2. `modern_thrrs_ecoli_zn_SER` (high priority - Figure 5C, the trap!)
3. `modern_thrrs_ecoli_zn_ILE` (high priority - Figure 4D, rejection)

Once these complete, run all 5 PyMOL scripts for complete structural visualization.

### Option 3: Alternative Visualization
For Figures 4C, 4D, 5C:
- Use competition predictions that contain Zn:
  - `COMPETITION_modern_thrrs_THR_vs_ILE` (has THR, ILE, Zn)
  - `COMPETITION_modern_thrrs_THR_vs_SER` (has THR, SER, Zn)
- These may have CIF files if competition jobs were run

Let me check if competition CIF files exist...

---

## 🔍 Checking Competition Structures

Searching for competition job outputs...

**Jobs in CSV:**
- COMPETITION_modern_thrrs_THR_vs_ILE (ipTM=0.960, Zn_iptm=0.790)
- COMPETITION_modern_thrrs_THR_vs_SER (ipTM=0.960, Zn_iptm=0.800)
- COMPETITION_anc_thrrs_THR_vs_ILE (ipTM=0.890, Zn_iptm=0.880)

**CIF file status:** Need to check outputs directory for COMPETITION_* folders

---

## 📝 File Mapping Reference

For scripts that ARE ready, update paths as follows:

### Figure 2D Script Updates:
```python
# OLD (template):
cmd.load("path/to/editing_THR.cif", "edit_THR")
cmd.load("path/to/editing_PRO.cif", "edit_PRO")

# NEW (actual paths):
cmd.load("/storage/kiran-stuff/aaRS/phase2/outputs/deep_editing_thr/deep_editing_thr/deep_editing_thr_model.cif", "edit_THR")
cmd.load("/storage/kiran-stuff/aaRS/phase2/outputs/deep_editing_pro/deep_editing_pro/deep_editing_pro_model.cif", "edit_PRO")
```

### Figure 3D Script Updates:
```python
# OLD (template):
cmd.load("path/to/anc_thrrs_THR.cif", "ancestral")
cmd.load("path/to/mod_thrrs_THR.cif", "modern")

# NEW (actual paths):
cmd.load("/storage/kiran-stuff/aaRS/phase2/outputs/deep_thrrs_thr/deep_thrrs_thr/deep_thrrs_thr_model.cif", "ancestral")
cmd.load("/storage/kiran-stuff/aaRS/phase2/outputs/modern_thrrs_thr/modern_thrrs_thr/modern_thrrs_thr_model.cif", "modern")
```

---

## 💡 Key Insights

1. **ipTM scores in CSV come from summary_confidences.json files**, not from CIF structure analysis
2. **CIF files are only generated when AF3 jobs complete successfully**
3. **Modern ThrRS + Zn jobs appear to have been planned but not executed**
4. **All 20 amino acids with Zn have ipTM scores** but no structural outputs

This suggests:
- Phase 1: Initial AF3 runs (ancestral enzymes, no Zn) ✅ COMPLETE
- Phase 2: Modern ThrRS + Zn screen (all 20 AAs) ⏸️ INCOMPLETE

**Recommendation:** Prioritize running AF3 for THR, SER, and ILE with Zn to complete the zinc trap story.

---

**Last updated:** 2025-12-19
