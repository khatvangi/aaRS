#!/usr/bin/env python3
"""
Run energy scoring in parallel across all CPU cores.
"""
import os, sys, json, glob
import pandas as pd
from multiprocessing import Pool, cpu_count
import subprocess

def score_one_wrapper(args):
    """Wrapper to call score_single.py and parse JSON output"""
    cif_path, mode = args

    try:
        # Call score_single.py
        python_exe = "/storage/kiran-stuff/blast_env/bin/python"
        cmd = [python_exe, "energy_scoring/score_single.py", cif_path, mode]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        if result.returncode == 0:
            # Parse JSON output
            return json.loads(result.stdout.strip())
        else:
            return {
                "file": cif_path,
                "status": "SUBPROCESS_ERROR",
                "error": result.stderr[:500]
            }
    except subprocess.TimeoutExpired:
        return {
            "file": cif_path,
            "status": "TIMEOUT",
            "error": "Exceeded 300s timeout"
        }
    except Exception as e:
        return {
            "file": cif_path,
            "status": "WRAPPER_ERROR",
            "error": str(e)
        }

def main(mode="nomin", ncores=None):
    """
    Run scoring in parallel.
    mode: "nomin" or "min"
    ncores: number of cores (default: all available)
    """
    if ncores is None:
        ncores = cpu_count()

    print(f"="*80)
    print(f"PARALLEL ENERGY SCORING - MODE: {mode.upper()}")
    print(f"="*80)
    print(f"Using {ncores} CPU cores")

    # Find all CIF files
    base = "/storage/kiran-stuff/aaRS/phase2"
    os.chdir(base)

    files = sorted(glob.glob("**/*_model.cif", recursive=True))
    print(f"Found {len(files)} CIF files")

    # Create argument list
    args_list = [(f, mode) for f in files]

    # Run in parallel
    print(f"\nStarting parallel processing...")
    with Pool(processes=ncores) as pool:
        results = []
        for i, result in enumerate(pool.imap_unordered(score_one_wrapper, args_list)):
            results.append(result)
            if (i+1) % 10 == 0:
                print(f"  Completed {i+1}/{len(files)}")

    # Save results
    df = pd.DataFrame(results)
    output_file = f"energy_scoring/scores_{mode}.csv"
    df.to_csv(output_file, index=False)

    print(f"\n{'='*80}")
    print(f"COMPLETE!")
    print(f"{'='*80}")
    print(f"✓ Saved: {output_file} ({len(df)} rows)")
    print(f"\nStatus summary:")
    print(df['status'].value_counts())

    # Mode-specific stats
    if mode == "min" and "ligand_RMSD_A" in df.columns:
        ok = df[df['status'].str.startswith('OK', na=False)]
        if len(ok) > 0:
            print(f"\nLigand RMSD after minimization:")
            print(f"  Median: {ok['ligand_RMSD_A'].median():.3f} Å")
            print(f"  P95: {ok['ligand_RMSD_A'].quantile(0.95):.3f} Å")
            unstable = ok[ok['ligand_RMSD_A'] > 1.0]
            print(f"  Unstable poses (RMSD > 1.0 Å): {len(unstable)}")

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "nomin"
    ncores = int(sys.argv[2]) if len(sys.argv) > 2 else None

    if mode not in ["nomin", "min"]:
        print("Error: mode must be 'nomin' or 'min'")
        sys.exit(1)

    main(mode, ncores)
