#!/usr/bin/env python3
"""
QC Step 6: Analyze MIN_FAILED and NO_LIG_PARAMS failures.
"""
import pandas as pd
import os

def analyze_failures():
    df = pd.read_csv('mmgbsa_results.csv')

    # MIN_FAILED
    min_failed = df[df['status'] == 'MIN_FAILED']
    print(f"MIN_FAILED: {len(min_failed)} jobs")

    if len(min_failed) > 0:
        print(f"\nJobs that failed at minimization:")
        for idx, row in min_failed.iterrows():
            job_name = row['job_name']
            job_dir = os.path.join('jobs', job_name)

            # Check which step failed
            has_complex_prmtop = os.path.exists(os.path.join(job_dir, 'complex.prmtop'))
            has_min_out = os.path.exists(os.path.join(job_dir, 'min.out'))

            if not has_complex_prmtop:
                failure_point = 'tleap'
            elif not has_min_out:
                failure_point = 'sander'
            else:
                # Check min.out for error
                min_out = os.path.join(job_dir, 'min.out')
                with open(min_out) as f:
                    content = f.read()
                    if 'ERROR' in content or 'FATAL' in content:
                        failure_point = 'sander_error'
                    else:
                        failure_point = 'unknown'

            print(f"  {job_name}: {failure_point}")

    # NO_LIG_PARAMS
    no_params = df[df['status'] == 'NO_LIG_PARAMS']
    print(f"\nNO_LIG_PARAMS: {len(no_params)} jobs")

    if len(no_params) > 0:
        print(f"\nJobs missing ligand parameters:")
        for idx, row in no_params.iterrows():
            job_name = row['job_name']
            lig_resname = row.get('ligand_resname', row.get('ligand', 'UNKNOWN'))

            # Handle NaN
            if pd.isna(lig_resname):
                lig_resname = 'UNKNOWN'

            lig_dir = os.path.join('ligands_fixed', str(lig_resname))
            has_mol2 = os.path.exists(os.path.join(lig_dir, 'lig.mol2'))
            has_frcmod = os.path.exists(os.path.join(lig_dir, 'lig.frcmod'))

            print(f"  {job_name} (ligand={lig_resname}): mol2={has_mol2}, frcmod={has_frcmod}")

    # Create failure report
    with open('qc/failures_report.md', 'w') as f:
        f.write("# Failure Analysis\n\n")

        f.write(f"## MIN_FAILED ({len(min_failed)} jobs)\n\n")
        if len(min_failed) > 0:
            f.write("| Job Name | Failure Point | Suggested Fix |\n")
            f.write("|----------|---------------|---------------|\n")

            for idx, row in min_failed.iterrows():
                job_name = row['job_name']
                job_dir = os.path.join('jobs', job_name)

                has_complex_prmtop = os.path.exists(os.path.join(job_dir, 'complex.prmtop'))
                has_min_out = os.path.exists(os.path.join(job_dir, 'min.out'))

                if not has_complex_prmtop:
                    failure_point = 'tleap'
                    fix = 'Check tleap.log for parameter/structure errors'
                elif not has_min_out:
                    failure_point = 'sander start'
                    fix = 'Check topology/coordinate compatibility'
                else:
                    min_out = os.path.join(job_dir, 'min.out')
                    with open(min_out) as mf:
                        content = mf.read()
                        if 'vlimit exceeded' in content.lower():
                            failure_point = 'vlimit'
                            fix = 'Bad contacts - check initial structure'
                        elif 'NaN' in content or 'Infinity' in content:
                            failure_point = 'numerical instability'
                            fix = 'Reduce minimization step size or add restraints'
                        else:
                            failure_point = 'unknown'
                            fix = 'Check min.out manually'

                f.write(f"| {job_name} | {failure_point} | {fix} |\n")
        else:
            f.write("No MIN_FAILED jobs.\n")

        f.write(f"\n## NO_LIG_PARAMS ({len(no_params)} jobs)\n\n")
        if len(no_params) > 0:
            f.write("| Job Name | Ligand | Issue | Suggested Fix |\n")
            f.write("|----------|--------|-------|---------------|\n")

            for idx, row in no_params.iterrows():
                job_name = row['job_name']
                lig_resname = row.get('ligand_resname', row.get('ligand', 'UNKNOWN'))

                # Handle NaN
                if pd.isna(lig_resname):
                    lig_resname = 'UNKNOWN'

                lig_dir = os.path.join('ligands_fixed', str(lig_resname))
                if not os.path.exists(lig_dir):
                    issue = f'Directory ligands_fixed/{lig_resname} missing'
                    fix = f'Create parameters for {lig_resname}'
                else:
                    has_mol2 = os.path.exists(os.path.join(lig_dir, 'lig.mol2'))
                    has_frcmod = os.path.exists(os.path.join(lig_dir, 'lig.frcmod'))

                    if not has_mol2:
                        issue = 'lig.mol2 missing'
                        fix = 'Run parameterization for this AA'
                    elif not has_frcmod:
                        issue = 'lig.frcmod missing'
                        fix = 'Run parmchk2 for this AA'
                    else:
                        issue = 'unknown'
                        fix = 'Check file permissions'

                f.write(f"| {job_name} | {lig_resname} | {issue} | {fix} |\n")
        else:
            f.write("No NO_LIG_PARAMS jobs.\n")

        f.write("\n## Summary\n\n")
        f.write(f"- Total failures: {len(min_failed) + len(no_params)}\n")
        f.write(f"- Success rate: {(len(df[df['status'] == 'OK']) / len(df) * 100):.1f}%\n")

    print(f"\n✓ Wrote qc/failures_report.md")

if __name__ == '__main__':
    analyze_failures()
