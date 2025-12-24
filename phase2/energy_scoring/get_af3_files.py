#!/usr/bin/env python3
"""
Get list of CIF files matching AF3_RESULTS_CORRECTED.csv job names.
"""
import pandas as pd
import glob
import re

# Load AF3 jobs
df = pd.read_csv('AF3_RESULTS_CORRECTED.csv')
job_names = set(df['job_name'].unique())

print(f"AF3 jobs to find: {len(job_names)}")

# Find all AF3 CIF files
all_cifs = sorted(glob.glob("**/af3_output/**/*_model.cif", recursive=True))
all_cifs += sorted(glob.glob("outputs/**/*_model.cif", recursive=True))
all_cifs = list(set(all_cifs))

print(f"Total CIF files found: {len(all_cifs)}")

# Extract job name from each CIF
def extract_job_name(filepath):
    filename = filepath.split('/')[-1].replace('_model.cif', '')
    # Remove seed-X_sample-Y pattern if present
    job_name = re.sub(r'_seed-\d+_sample-\d+$', '', filename)
    return job_name

# Filter to matching files (take best sample for each job)
matched_files = {}
for cif in all_cifs:
    job = extract_job_name(cif)
    if job in job_names:
        # Keep first match (sorted, so seed-1_sample-0 comes first)
        if job not in matched_files:
            matched_files[job] = cif

print(f"Matched CIF files: {len(matched_files)}")
print(f"Missing jobs: {len(job_names) - len(matched_files)}")

# Save to file list
with open('energy_scoring/af3_files_to_score.txt', 'w') as f:
    for job, cif in sorted(matched_files.items()):
        f.write(f"{cif}\n")

print(f"\n✓ Saved: energy_scoring/af3_files_to_score.txt ({len(matched_files)} files)")

# Show which jobs are missing
missing = job_names - set(matched_files.keys())
if missing:
    print(f"\nMissing jobs ({len(missing)}):")
    for job in sorted(missing)[:10]:
        print(f"  {job}")
