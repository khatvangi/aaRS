#!/usr/bin/env python3
"""
Monitor MM/GBSA batch progress and check for convergence issues.
"""
import os
import time
import glob
import pandas as pd
from check_convergence import parse_minimization_output

def monitor_jobs(jobs_dir="jobs", check_interval=30):
    """
    Monitor running jobs and report convergence status.
    """
    print("Monitoring MM/GBSA jobs...")
    print(f"Jobs directory: {jobs_dir}")
    print(f"Check interval: {check_interval} seconds")
    print("=" * 80)

    last_count = 0

    while True:
        # Find all job directories
        job_dirs = sorted(glob.glob(os.path.join(jobs_dir, "*")))

        if len(job_dirs) == 0:
            print("No job directories found yet...")
            time.sleep(check_interval)
            continue

        # Count completed minimizations
        completed = []
        running = []
        not_started = []

        for job_dir in job_dirs:
            min_out = os.path.join(job_dir, "min.out")
            min_rst = os.path.join(job_dir, "min.rst")

            if os.path.exists(min_rst):
                # Minimization completed
                conv_info = parse_minimization_output(min_out)
                if conv_info:
                    completed.append((job_dir, conv_info))
            elif os.path.exists(min_out):
                # Minimization running
                running.append(job_dir)
            else:
                # Not started
                not_started.append(job_dir)

        # Report status
        print(f"\n[{time.strftime('%H:%M:%S')}] Status:")
        print(f"  Completed: {len(completed)}/{len(job_dirs)}")
        print(f"  Running:   {len(running)}")
        print(f"  Pending:   {len(not_started)}")

        # Check convergence of completed jobs
        if len(completed) > last_count:
            print(f"\nConvergence check (last {len(completed) - last_count} completed):")

            converged_count = 0
            not_converged = []

            for job_dir, conv_info in completed[last_count:]:
                job_name = os.path.basename(job_dir)
                if conv_info['converged']:
                    converged_count += 1
                else:
                    not_converged.append((job_name, conv_info))

            print(f"  ✓ Converged: {converged_count}/{len(completed) - last_count}")

            if not_converged:
                print(f"  ✗ NOT converged: {len(not_converged)}")
                for job_name, conv_info in not_converged[:5]:  # Show first 5
                    print(f"    - {job_name}: RMS={conv_info['rms_gradient']:.3f} kcal/mol/Å (need <0.1)")

            last_count = len(completed)

        # Check if all jobs are done
        if len(completed) == len(job_dirs):
            print("\n" + "=" * 80)
            print("All jobs completed!")
            print("=" * 80)

            # Final summary
            total_converged = sum(1 for _, conv_info in completed if conv_info['converged'])
            print(f"\nFinal convergence summary:")
            print(f"  Total jobs: {len(completed)}")
            print(f"  Converged: {total_converged} ({100*total_converged/len(completed):.1f}%)")
            print(f"  NOT converged: {len(completed) - total_converged}")

            # Save convergence report
            convergence_data = []
            for job_dir, conv_info in completed:
                convergence_data.append({
                    'job_name': os.path.basename(job_dir),
                    'converged': conv_info['converged'],
                    'nstep': conv_info['nstep'],
                    'rms_gradient': conv_info['rms_gradient'],
                    'max_gradient': conv_info['max_gradient'],
                    'energy': conv_info['energy']
                })

            df = pd.DataFrame(convergence_data)
            df.to_csv('convergence_report.csv', index=False)
            print(f"\nSaved: convergence_report.csv")

            break

        time.sleep(check_interval)

if __name__ == '__main__':
    monitor_jobs()
