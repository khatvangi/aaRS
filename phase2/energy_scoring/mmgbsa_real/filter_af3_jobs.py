import pandas as pd
import re

# Load AF3 results
af3_df = pd.read_csv('../../AF3_RESULTS_CORRECTED.csv')
af3_jobs = set(af3_df['job_name'].unique())

print(f"AF3 jobs to score: {len(af3_jobs)}")

# Load manifest
manifest_df = pd.read_csv('../../manifest.csv')

def extract_job_name(filepath):
    filename = filepath.split('/')[-1].replace('_model.cif', '')
    job_name = re.sub(r'_seed-\d+_sample-\d+$', '', filename)
    return job_name

# Filter manifest to AF3 jobs only
manifest_df['job_name'] = manifest_df['file'].apply(extract_job_name)
af3_manifest_all = manifest_df[manifest_df['job_name'].isin(af3_jobs)].copy()

print(f"Total matched CIF files (all samples): {len(af3_manifest_all)}")

# Take only ONE sample per job (first one alphabetically)
af3_manifest = af3_manifest_all.sort_values('file').groupby('job_name').first().reset_index()

print(f"One sample per job: {len(af3_manifest)}")
print(f"Zn present: {af3_manifest['zn_present'].sum()}")
print(f"No Zn: {(af3_manifest['zn_present']==0).sum()}")

# Save filtered manifest
af3_manifest.to_csv('manifest_af3only.csv', index=False)
print("\nSaved: manifest_af3only.csv")
