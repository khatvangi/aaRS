#!/usr/bin/env python3
"""
Run MM/GBSA in parallel across all AF3-generated CIF files (60 cores).
"""
import os
import sys
import json
import glob
import pandas as pd
from multiprocessing import Pool, cpu_count
import subprocess

def score_one_wrapper(cif_path):
    """Wrapper to call mmgbsa_simple.py and parse JSON output"""
    try:
        python_exe = "/storage/kiran-stuff/blast_env/bin/python"
        cmd = [python_exe, "energy_scoring/mmgbsa_simple.py", cif_path]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

        if result.returncode == 0:
            return json.loads(result.stdout.strip())
        else:
            return {
                "file": cif_path,
                "status": "SUBPROCESS_ERROR",
                "error": result.stderr[:200]
            }
    except subprocess.TimeoutExpired:
        return {
            "file": cif_path,
            "status": "TIMEOUT",
            "error": "Exceeded 120s timeout"
        }
    except Exception as e:
        return {
            "file": cif_path,
            "status": "WRAPPER_ERROR",
            "error": str(e)
        }

def main(ncores=None):
    """
    Run MM/GBSA scoring in parallel.
    ncores: number of cores (default: 60)
    """
    if ncores is None:
        ncores = 60

    print(f"="*80)
    print(f"PARALLEL MM/GBSA SCORING")
    print(f"="*80)
    print(f"Using {ncores} CPU cores")

    # Change to base directory
    base = "/storage/kiran-stuff/aaRS/phase2"
    os.chdir(base)

    # Find all AF3-generated CIF files (only from af3_output/ and /outputs/ directories)
    files = sorted(glob.glob("**/af3_output/**/*_model.cif", recursive=True))
    files += sorted(glob.glob("outputs/**/*_model.cif", recursive=True))
    files = list(set(files))  # Remove duplicates

    print(f"Found {len(files)} AF3-generated CIF files")

    # Run in parallel
    print(f"\nStarting parallel processing...")
    with Pool(processes=ncores) as pool:
        results = []
        for i, result in enumerate(pool.imap_unordered(score_one_wrapper, files), 1):
            results.append(result)
            if i % 50 == 0:
                print(f"  Completed {i}/{len(files)}")

    # Save results
    df = pd.DataFrame(results)
    output_file = "energy_scoring/mmgbsa_scores.csv"
    df.to_csv(output_file, index=False)

    print(f"\n{'='*80}")
    print(f"COMPLETE!")
    print(f"{'='*80}")
    print(f"✓ Saved: {output_file} ({len(df)} rows)")
    print(f"\nStatus summary:")
    print(df['status'].value_counts())

    # Statistics for successful runs
    ok = df[df['status'] == 'OK']
    if len(ok) > 0:
        print(f"\nMM/GBSA statistics:")
        print(f"  Structures scored: {len(ok)}")
        print(f"  BE range: [{ok['BE_mmgbsa'].min():.1f}, {ok['BE_mmgbsa'].max():.1f}] kcal/mol")
        print(f"  BE median: {ok['BE_mmgbsa'].median():.1f} kcal/mol")

        # Group by Zn presence
        with_zn = ok[ok['Zn_present'] == True]
        no_zn = ok[ok['Zn_present'] == False]
        print(f"\n  With Zn: {len(with_zn)} structures, mean BE = {with_zn['BE_mmgbsa'].mean():.1f}")
        print(f"  No Zn:   {len(no_zn)} structures, mean BE = {no_zn['BE_mmgbsa'].mean():.1f}")

        # Top/bottom 5 by BE
        print(f"\nMost favorable binding (lowest BE):")
        print(ok.nsmallest(5, 'BE_mmgbsa')[['file', 'lig_resname', 'BE_mmgbsa', 'Zn_present']].to_string(index=False))

        print(f"\nLeast favorable binding (highest BE):")
        print(ok.nlargest(5, 'BE_mmgbsa')[['file', 'lig_resname', 'BE_mmgbsa', 'Zn_present']].to_string(index=False))

if __name__ == "__main__":
    ncores = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    main(ncores)
