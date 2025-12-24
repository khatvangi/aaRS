# Failure Analysis

## MIN_FAILED (9 jobs)

| Job Name | Failure Point | Suggested Fix |
|----------|---------------|---------------|
| proRS_modern_Pro | tleap | Check tleap.log for parameter/structure errors |
| fulllength_deep_thr | unknown | Check min.out manually |
| proRS_LUCA_Thr | unknown | Check min.out manually |
| COMPETITION_FULL_anc_prors_PRO_vs_THR | tleap | Check tleap.log for parameter/structure errors |
| fulllength_deep_pro | tleap | Check tleap.log for parameter/structure errors |
| fulllength_shallow_pro | tleap | Check tleap.log for parameter/structure errors |
| fulllength_shallow_thr | unknown | Check min.out manually |
| COMPETITION_FULL_anc_prors_PRO_vs_GLU | tleap | Check tleap.log for parameter/structure errors |
| proRS_LUCA_Pro | tleap | Check tleap.log for parameter/structure errors |

## NO_LIG_PARAMS (2 jobs)

| Job Name | Ligand | Issue | Suggested Fix |
|----------|--------|-------|---------------|
| modern_ecoli_thrrs_bHNV | UNKNOWN | Directory ligands_fixed/UNKNOWN missing | Create parameters for UNKNOWN |
| test_aars_trna_complex | UNKNOWN | Directory ligands_fixed/UNKNOWN missing | Create parameters for UNKNOWN |

## Summary

- Total failures: 11
- Success rate: 94.0%
