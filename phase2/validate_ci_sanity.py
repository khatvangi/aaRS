#!/usr/bin/env python3
"""
Validate that bootstrap CIs correctly bracket effect values.

Usage:
  python validate_ci_sanity.py sweep/effects_stability_table.csv

Expected: ~0 bad rows (CI should bracket effect for nearly all rows)
"""
import sys
import pandas as pd

def validate_ci(csv_file):
    print(f"Loading {csv_file}...")
    df = pd.read_csv(csv_file)

    print(f"Total rows: {len(df)}")

    # Check if CI brackets effect
    bad = df[~((df['ci_low'] <= df['effect']) & (df['effect'] <= df['ci_high']))]

    print(f"\n{'='*70}")
    print(f"CI SANITY CHECK")
    print(f"{'='*70}")
    print(f"Total comparisons: {len(df)}")
    print(f"Bad CI rows (effect not in [ci_low, ci_high]): {len(bad)}")
    print(f"Percentage bad: {100*len(bad)/len(df):.2f}%")

    if len(bad) > 0:
        print(f"\n{'='*70}")
        print(f"FIRST 10 BAD ROWS:")
        print(f"{'='*70}")
        print(bad.head(10)[['setting_id', 'condition', 'metric', 'effect', 'ci_low', 'ci_high']].to_string(index=False))

        print(f"\n{'='*70}")
        print(f"DIAGNOSIS")
        print(f"{'='*70}")

        # Check if CIs look like pooled means
        mean_effect = df['effect'].abs().mean()
        mean_ci_low = df['ci_low'].mean()
        mean_ci_high = df['ci_high'].mean()

        print(f"Mean |effect|: {mean_effect:.3f}")
        print(f"Mean ci_low: {mean_ci_low:.3f}")
        print(f"Mean ci_high: {mean_ci_high:.3f}")

        if mean_ci_low > 1 and mean_ci_high > 2:
            print("\n⚠️  CIs appear to be POOLED MEANS, not difference CIs")
            print("   → Need to rerun with bootstrap_ci_diff()")
        else:
            print("\n✓ CIs are in correct range for difference")
            print("  Bad rows may be numerical edge cases")
    else:
        print(f"\n✓ ALL CIs CORRECTLY BRACKET EFFECTS")
        print(f"  Analysis is statistically valid for publication")

    # Sample size summary
    print(f"\n{'='*70}")
    print(f"SAMPLE SIZE SUMMARY")
    print(f"{'='*70}")

    # Group by n_cog, n_non
    small_n = df[(df['n_cog'] <= 2) | (df['n_non'] <= 3)]
    print(f"Comparisons with small n (n_cog≤2 or n_non≤3): {len(small_n)}/{len(df)} ({100*len(small_n)/len(df):.1f}%)")

    if len(small_n) > 100:
        print(f"\n⚠️  WARNING: Many comparisons have small sample sizes")
        print(f"   Consider flagging these in supplementary materials")

    return len(bad) == 0

if __name__ == '__main__':
    if len(sys.argv) < 2:
        csv_file = 'sweep/effects_stability_table.csv'
    else:
        csv_file = sys.argv[1]

    passed = validate_ci(csv_file)

    print(f"\n{'='*70}")
    print(f"VALIDATION: {'PASS ✓' if passed else 'FAIL ✗'}")
    print(f"{'='*70}")

    sys.exit(0 if passed else 1)
