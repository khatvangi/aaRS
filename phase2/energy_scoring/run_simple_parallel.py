#!/usr/bin/env python3
"""
Run simple energy scoring in parallel across all CPU cores.
"""
import os, sys, json, glob
import pandas as pd
from multiprocessing import Pool, cpu_count
import subprocess

def score_one_wrapper(cif_path):
    """Wrapper to call score_simple.py and parse JSON output"""
    try:
        python_exe = "/storage/kiran-stuff/blast_env/bin/python"
        cmd = [python_exe, "energy_scoring/score_simple.py", cif_path]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)

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
            "error": "Exceeded 60s timeout"
        }
    except Exception as e:
        return {
            "file": cif_path,
            "status": "WRAPPER_ERROR",
            "error": str(e)
        }

def main(ncores=None):
    """
    Run scoring in parallel.
    ncores: number of cores (default: all available)
    """
    if ncores is None:
        ncores = cpu_count()

    print(f"="*80)
    print(f"PARALLEL ENERGY SCORING - SIMPLE METHOD")
    print(f"="*80)
    print(f"Using {ncores} CPU cores")

    # Change to base directory
    base = "/storage/kiran-stuff/aaRS/phase2"
    os.chdir(base)

    # Find all CIF files
    files = sorted(glob.glob("**/*_model.cif", recursive=True))
    print(f"Found {len(files)} CIF files")

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
    output_file = "energy_scoring/scores_simple.csv"
    df.to_csv(output_file, index=False)

    print(f"\n{'='*80}")
    print(f"COMPLETE!")
    print(f"{'='*80}")
    print(f"✓ Saved: {output_file} ({len(df)} rows)")
    print(f"\nStatus summary:")
    print(df['status'].value_counts())

    # Energy statistics
    ok = df[df['status'] == 'OK']
    if len(ok) > 0:
        print(f"\nEnergy statistics (kcal/mol):")
        print(f"  Eint range: [{ok['Eint_kcal'].min():.1f}, {ok['Eint_kcal'].max():.1f}]")
        print(f"  Eint median: {ok['Eint_kcal'].median():.1f}")
        print(f"  Evdw median: {ok['Evdw_kcal'].median():.1f}")
        print(f"  Ecoul median: {ok['Ecoul_kcal'].median():.1f}")

        # Top/bottom 5 by energy
        print(f"\nLowest Eint (most favorable):")
        print(ok.nsmallest(5, 'Eint_kcal')[['file', 'lig_resname', 'Eint_kcal']].to_string(index=False))

        print(f"\nHighest Eint (least favorable):")
        print(ok.nlargest(5, 'Eint_kcal')[['file', 'lig_resname', 'Eint_kcal']].to_string(index=False))

if __name__ == "__main__":
    ncores = int(sys.argv[1]) if len(sys.argv) > 1 else 60  # Use 60 cores by default
    main(ncores)
