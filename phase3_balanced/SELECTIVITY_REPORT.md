# Cognate vs Non-Cognate Selectivity Analysis

## Summary

**Metric:** iptm (interface predicted TM-score from AlphaFold3)
**Total AF3 runs:** 1,101
**FDR correction:** Benjamini-Hochberg across 4 tests

## Results

| Enzyme | Construct | n_cog | n_non | Effect | 95% CI | q-value | Significant |
|--------|-----------|-------|-------|--------|--------|---------|-------------|
| **ProRS** | **ancestral** | 108 | 186 | **+0.050** | [+0.022, +0.080] | **0.006** | **YES** |
| ProRS | modern | 62 | 120 | -0.134 | [-0.203, -0.068] | 0.00004 | YES (wrong direction) |
| ThrRS | ancestral | 48 | 189 | +0.014 | [-0.003, +0.029] | 0.139 | NO |
| ThrRS | modern | 66 | 214 | -0.006 | [-0.060, +0.051] | 0.835 | NO |

## Key Finding

**Ancestral ProRS shows significant cognate selectivity:**
- Cognate (PRO) has higher iptm than non-cognates
- Effect: +0.050 (95% CI: +0.022 to +0.080)
- q-value: 0.006 (survives FDR correction)

## Interpretation

1. **Ancestral ProRS**: Clear evidence of cognate selectivity in AF3 predictions. The enzyme preferentially accommodates its cognate substrate (proline) with higher predicted interface quality.

2. **Modern ProRS**: Shows reversed pattern (cognate LOWER). This may reflect:
   - Different binding mode in modern vs ancestral
   - AF3 limitation in modeling modern enzyme specificity
   - Possible misannotation in construct naming

3. **ThrRS (both constructs)**: No significant cognate selectivity detected. Either:
   - ThrRS selectivity operates through mechanisms AF3 doesn't capture
   - Sample sizes still insufficient for ThrRS effect size
   - ThrRS may rely more on editing domain (not tested here)

## Methods

- **Effect:** mean(cognate iptm) - mean(non-cognate iptm)
- **CI:** Bootstrap (10,000 resamples, independent resampling of both groups)
- **p-value:** Permutation test (100,000 permutations)
- **FDR:** Benjamini-Hochberg across all 4 tests

## Data Files

- Raw scores: `all_af3_iptm_scores.csv`
- Results table: `iptm_selectivity_results.csv`
