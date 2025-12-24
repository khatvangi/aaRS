# Sequence Verification Report - batch_local_domains

## Summary

**Status**: ⚠️ **CRITICAL ERROR FOUND**

- ✅ ProRS catalytic domain: **CORRECT**
- ✅ ProRS editing domain: **CORRECT**
- ❌ ThrRS catalytic domain: **INCORRECT - PLACEHOLDER SEQUENCE**

---

## Verification Details

### 1. ProRS Catalytic Domain ✅

**File**: `anc_prors_cat_*.json`

**Current sequence** (first 50 chars):
```
QDVGKFVELPGAEMGKVTTRYPPEASGYLHIGHAKAALLNQYYQLNFKGK...
```

**Reference** (from `phase2/inputs/af3_jsons_domain/deep_domain_pro.json`):
```
QDVGKFVELPGAEMGKVTTRYPPEASGYLHIGHAKAALLNQYYQLNFKGK...
```

**Status**: ✅ MATCH (437 amino acids)

**Verification**:
```bash
# Length: 437 aa
# Starts with: QDVGKFVELPGAEMGKVTTRY
# Ends with: GDIIQLQRRGFFIC
```

---

### 2. ProRS Editing Domain ✅

**File**: `anc_prors_edit_*.json`

**Current sequence** (first 50 chars):
```
RTSEFSLWQEGHTAHATREEAEEEVRRMLHDIYARFAREYLAIPVIKGRK...
```

**Reference** (from `phase2/inputs/af3_jsons_validation/deep_editing_pro.json`):
```
RTSEFSLWQEGHTAHATREEAEEEVRRMLHDIYARFAREYLAIPVIKGRK...
```

**Status**: ✅ MATCH (370 amino acids)

**Verification**:
```bash
# Length: 370 aa
# Starts with: RTSEFSLWQEGHTAHATREEAE
# Ends with: QNFARA
```

---

### 3. ThrRS Catalytic Domain ❌

**File**: `anc_thrrs_cat_*.json`

**Current sequence**:
```
MAAAAA
```

**Expected sequence** (from `phase2/inputs/af3_jsons_validation/deep_thrrs_pro.json`):
```
EEEYLLKPMNCPHHIRIYASLRPRSYRDLPLRLAEFGTVSYRYEQSGELHGLTRVRGFTQDDAHIFCRPLDQIKGEFLGVLDLVLLVFKLTFGLNGDDYRARIGVRDKVNDPKSGKYVVSYVCRNCGYRIEGARGTDIECPVCHSHDLEGDEENWALAERQIKEAVEEYKNEVCGLPGNYTIEPGDAAFYGPKLDFVVKDALGREWQLGTIQVDYNLPERFDLTYVGGLNNNNNNNNNNEEINDNNEDGQELEHRPVMIHRAPFGSIERFIGILIEHY
```

**Status**: ❌ **INCORRECT - PLACEHOLDER SEQUENCE "MAAAAA"**

**Expected**:
```bash
# Length: should be 280 amino acids
# Should start with: EEEYLLKPMNCPHHIRIYASL
# Should end with: FGSIERFIGILIEHY
```

---

## Ancient Reconstruction Files (Reference)

Based on the phase2/inputs directories, here are the 5 ancient reconstruction JSON files:

### ProRS (3 files):
1. **Full-length**: `af3_jsons_fulllength/fulllength_deep_pro.json`
   - Full ProRS with tRNA and ligand
   - ~2200+ amino acids

2. **Catalytic domain**: `af3_jsons_domain/deep_domain_pro.json`
   - 437 amino acids
   - Sequence: QDVGKFVELPGAEMGKVTTRY...

3. **Editing domain**: `af3_jsons_validation/deep_editing_pro.json`
   - 370 amino acids
   - Sequence: RTSEFSLWQEGHTAHATREEAE...

### ThrRS (2 files):
1. **Full-length** (with cognate ligand): `af3_jsons_validation/deep_thrrs_thr.json`
   - 280 amino acids
   - Sequence: EEEYLLKPMNCPHHIRIYASL...
   - This appears to be the catalytic domain

2. **Full-length** (with non-cognate ligand): `af3_jsons_validation/deep_thrrs_pro.json`
   - Same 280 amino acid sequence
   - Different ligand (Pro instead of Thr)

**Note**: There's no separate "ThrRS editing domain" in the ancient reconstructions, which matches the biological expectation that ThrRS may lack an editing domain or have a much smaller one.

---

## Recommended Fix

### Fix ThrRS Catalytic Domain Sequences

All 20 `anc_thrrs_cat_*.json` files need to be updated with the correct ThrRS sequence:

**Correct sequence (280 aa)**:
```
EEEYLLKPMNCPHHIRIYASLRPRSYRDLPLRLAEFGTVSYRYEQSGELHGLTRVRGFTQDDAHIFCRPLDQIKGEFLGVLDLVLLVFKLTFGLNGDDYRARIGVRDKVNDPKSGKYVVSYVCRNCGYRIEGARGTDIECPVCHSHDLEGDEENWALAERQIKEAVEEYKNEVCGLPGNYTIEPGDAAFYGPKLDFVVKDALGREWQLGTIQVDYNLPERFDLTYVGGLNNNNNNNNNNEEINDNNEDGQELEHRPVMIHRAPFGSIERFIGILIEHY
```

### Files to Update

```bash
/storage/kiran-stuff/aaRS/phase2/batch_local_domains/anc_thrrs_cat_ALA.json
/storage/kiran-stuff/aaRS/phase2/batch_local_domains/anc_thrrs_cat_ARG.json
/storage/kiran-stuff/aaRS/phase2/batch_local_domains/anc_thrrs_cat_ASN.json
/storage/kiran-stuff/aaRS/phase2/batch_local_domains/anc_thrrs_cat_ASP.json
/storage/kiran-stuff/aaRS/phase2/batch_local_domains/anc_thrrs_cat_CYS.json
/storage/kiran-stuff/aaRS/phase2/batch_local_domains/anc_thrrs_cat_GLN.json
/storage/kiran-stuff/aaRS/phase2/batch_local_domains/anc_thrrs_cat_GLU.json
/storage/kiran-stuff/aaRS/phase2/batch_local_domains/anc_thrrs_cat_GLY.json
/storage/kiran-stuff/aaRS/phase2/batch_local_domains/anc_thrrs_cat_HIS.json
/storage/kiran-stuff/aaRS/phase2/batch_local_domains/anc_thrrs_cat_ILE.json
/storage/kiran-stuff/aaRS/phase2/batch_local_domains/anc_thrrs_cat_LEU.json
/storage/kiran-stuff/aaRS/phase2/batch_local_domains/anc_thrrs_cat_LYS.json
/storage/kiran-stuff/aaRS/phase2/batch_local_domains/anc_thrrs_cat_MET.json
/storage/kiran-stuff/aaRS/phase2/batch_local_domains/anc_thrrs_cat_PHE.json
/storage/kiran-stuff/aaRS/phase2/batch_local_domains/anc_thrrs_cat_PRO.json
/storage/kiran-stuff/aaRS/phase2/batch_local_domains/anc_thrrs_cat_SER.json
/storage/kiran-stuff/aaRS/phase2/batch_local_domains/anc_thrrs_cat_THR.json
/storage/kiran-stuff/aaRS/phase2/batch_local_domains/anc_thrrs_cat_TRP.json
/storage/kiran-stuff/aaRS/phase2/batch_local_domains/anc_thrrs_cat_TYR.json
/storage/kiran-stuff/aaRS/phase2/batch_local_domains/anc_thrrs_cat_VAL.json
```

---

## Summary Table

| Domain | File Pattern | Length (aa) | Status | Notes |
|--------|-------------|-------------|--------|-------|
| ProRS catalytic | `anc_prors_cat_*.json` | 437 | ✅ Correct | Matches deep_domain_pro.json |
| ProRS editing | `anc_prors_edit_*.json` | 370 | ✅ Correct | Matches deep_editing_pro.json |
| ThrRS catalytic | `anc_thrrs_cat_*.json` | 6 | ❌ **WRONG** | Has "MAAAAA" instead of 280 aa sequence |

---

## Impact

**If not fixed**:
- All AF3 predictions for ThrRS catalytic domain will fail or produce meaningless results
- The 20-amino acid matrix for ThrRS will be invalid
- Cannot compare ProRS vs ThrRS substrate promiscuity

**Urgency**: 🔴 **HIGH** - Must fix before running AF3 calculations

---

## Verification Command

```bash
# Check ProRS catalytic length
cat anc_prors_cat_PRO.json | jq -r '.sequences[0].protein.sequence' | wc -c
# Should output: 438 (437 + newline)

# Check ProRS editing length
cat anc_prors_edit_PRO.json | jq -r '.sequences[0].protein.sequence' | wc -c
# Should output: 371 (370 + newline)

# Check ThrRS catalytic length (CURRENT - WRONG)
cat anc_thrrs_cat_PRO.json | jq -r '.sequences[0].protein.sequence' | wc -c
# Currently outputs: 7 (6 + newline) - WRONG!
# Should output: 281 (280 + newline)
```

---

**Generated**: 2025-12-11
**Verified by**: Claude Code sequence verification
