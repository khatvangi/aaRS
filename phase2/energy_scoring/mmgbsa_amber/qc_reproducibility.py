#!/usr/bin/env python3
"""
QC Step 4: Test reproducibility by re-running 10 random jobs.
Pass criteria: |ΔBE| < 0.5 kcal/mol for all 10
"""
import pandas as pd
import numpy as np
import os
import sys
import shutil
import random

# Import run_one_job from batch script
sys.path.insert(0, os.path.dirname(__file__))
from run_mmpbsa_batch import run_one_job

def main():
    # Load original results
    df = pd.read_csv('mmgbsa_results.csv')
    ok_df = df[df['status'] == 'OK'].copy()

    print(f"Total OK jobs: {len(ok_df)}")

    # Select 10 random jobs
    np.random.seed(42)
    sample_jobs = ok_df.sample(n=min(10, len(ok_df)), random_state=42)

    print(f"\nSelected {len(sample_jobs)} jobs for reproducibility test:")
    for job_name in sample_jobs['job_name']:
        print(f"  {job_name}")

    # Re-run each job
    results = []
    for idx, row in sample_jobs.iterrows():
        job_name = row['job_name']
        original_BE = row['BE_dG_bind']

        print(f"\nRe-running {job_name}...")

        # Backup and remove old job directory
        job_dir = os.path.join('jobs', job_name)
        backup_dir = os.path.join('jobs_backup_repro', job_name)

        if os.path.exists(backup_dir):
            shutil.rmtree(backup_dir)

        os.makedirs(os.path.dirname(backup_dir), exist_ok=True)
        if os.path.exists(job_dir):
            shutil.move(job_dir, backup_dir)

        # Re-run
        rerun_result = run_one_job(row.to_dict())

        # Extract new BE
        if rerun_result.get('status') == 'OK':
            new_BE = rerun_result['BE_dG_bind']
            delta_BE = abs(new_BE - original_BE)

            print(f"  Original: {original_BE:.4f} kcal/mol")
            print(f"  Rerun:    {new_BE:.4f} kcal/mol")
            print(f"  |ΔBE|:    {delta_BE:.4f} kcal/mol")

            results.append({
                'job_name': job_name,
                'ligand_resname': row['ligand_resname'],
                'original_BE': original_BE,
                'rerun_BE': new_BE,
                'delta_BE': delta_BE,
                'status': 'OK'
            })
        else:
            print(f"  FAILED: {rerun_result.get('error', 'unknown')}")
            results.append({
                'job_name': job_name,
                'ligand_resname': row['ligand_resname'],
                'original_BE': original_BE,
                'rerun_BE': np.nan,
                'delta_BE': np.nan,
                'status': rerun_result.get('status', 'FAILED')
            })

    # Save results
    results_df = pd.DataFrame(results)
    results_df.to_csv('qc/rerun_delta.csv', index=False)

    print(f"\n✓ Wrote qc/rerun_delta.csv")

    # Summary
    valid = results_df[results_df['status'] == 'OK']

    print(f"\n{'='*70}")
    print("Summary:")
    print(f"{'='*70}")
    print(f"Successful reruns: {len(valid)}/{len(results_df)}")

    if len(valid) > 0:
        print(f"\n|ΔBE| statistics (kcal/mol):")
        print(f"  Mean: {valid['delta_BE'].mean():.4f}")
        print(f"  Max:  {valid['delta_BE'].max():.4f}")
        print(f"  Min:  {valid['delta_BE'].min():.4f}")

        # Pass/Fail
        max_delta = valid['delta_BE'].max()
        if max_delta < 0.5:
            print(f"\n✓ PASS: Max |ΔBE| = {max_delta:.4f} < 0.5 kcal/mol")
        else:
            print(f"\n✗ FAIL: Max |ΔBE| = {max_delta:.4f} ≥ 0.5 kcal/mol")

            # Show problem cases
            problems = valid[valid['delta_BE'] >= 0.5]
            print(f"\n  Problem cases:")
            print(problems[['job_name', 'original_BE', 'rerun_BE', 'delta_BE']])

    # Restore backup directories
    print(f"\nRestoring backup directories...")
    for idx, row in sample_jobs.iterrows():
        job_name = row['job_name']
        job_dir = os.path.join('jobs', job_name)
        backup_dir = os.path.join('jobs_backup_repro', job_name)

        if os.path.exists(job_dir):
            shutil.rmtree(job_dir)

        if os.path.exists(backup_dir):
            shutil.move(backup_dir, job_dir)

    print(f"✓ Restored original job directories")

if __name__ == '__main__':
    main()
