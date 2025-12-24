#!/usr/bin/env python3
"""
Test reproducibility: run recompute_metrics twice, compare results.
"""
import pandas as pd
import numpy as np
import sys
import os

# Import recompute script
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

# Import with name mangling to handle numeric prefix
import importlib.util
spec = importlib.util.spec_from_file_location("recompute_metrics", os.path.join(script_dir, "00_recompute_metrics.py"))
recompute_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(recompute_module)
recompute_main = recompute_module.main

def compare_dataframes(df1, df2, numeric_columns):
    """Compare numeric columns between two dataframes."""

    differences = {}

    for col in numeric_columns:
        if col not in df1.columns or col not in df2.columns:
            continue

        # Handle NaN
        mask = df1[col].notna() & df2[col].notna()

        if mask.sum() == 0:
            continue

        vals1 = df1.loc[mask, col].values
        vals2 = df2.loc[mask, col].values

        abs_diff = np.abs(vals1 - vals2)

        differences[col] = {
            'max_abs_diff': abs_diff.max(),
            'mean_abs_diff': abs_diff.mean(),
            'n_compared': len(abs_diff)
        }

    return differences

def main():
    print("="*70)
    print("REPRODUCIBILITY TEST")
    print("="*70)

    # Run 1
    print("\nRun 1...")
    df1 = recompute_main(
        input_csv='geometry_metrics_clean.csv',
        output_csv='qc/recomputed_run1.csv'
    )

    # Run 2
    print("\nRun 2...")
    df2 = recompute_main(
        input_csv='geometry_metrics_clean.csv',
        output_csv='qc/recomputed_run2.csv'
    )

    # Compare
    print("\n" + "="*70)
    print("COMPARING RUNS")
    print("="*70)

    numeric_cols = [
        'contacts_4.0A', 'clashes_2.2A', 'polar_close_3.5A',
        'contacts_per_atom', 'polar_contacts_per_atom', 'clash_rate',
        'polar_close_contacts_per_atom', 'zn_min_dist_hetero'
    ]

    diffs = compare_dataframes(df1, df2, numeric_cols)

    # Create report
    report_rows = []
    for col, stats in diffs.items():
        report_rows.append({
            'column': col,
            'max_abs_diff': stats['max_abs_diff'],
            'mean_abs_diff': stats['mean_abs_diff'],
            'n_compared': stats['n_compared']
        })

    df_report = pd.DataFrame(report_rows)
    df_report.to_csv('qc/reproducibility_diff_report.csv', index=False)
    print("\nSaved qc/reproducibility_diff_report.csv")

    # Print summary
    print("\nDifference Summary:")
    print(df_report.to_string(index=False))

    # Check pass/fail
    max_diff = df_report['max_abs_diff'].max()
    threshold = 1e-10

    if max_diff < threshold:
        status = "PASS"
        print(f"\n✓ REPRODUCIBILITY TEST PASSED (max diff = {max_diff:.2e} < {threshold})")
    else:
        status = "FAIL"
        print(f"\n✗ REPRODUCIBILITY TEST FAILED (max diff = {max_diff:.2e} >= {threshold})")

    # Write pass file
    with open('qc/reproducibility_pass.txt', 'w') as f:
        f.write(f"Status: {status}\n")
        f.write(f"Max difference: {max_diff:.2e}\n")
        f.write(f"Threshold: {threshold}\n")

    print("\nSaved qc/reproducibility_pass.txt")

    print("\n" + "="*70)
    print("REPRODUCIBILITY TEST COMPLETE")
    print("="*70)

if __name__ == '__main__':
    main()
