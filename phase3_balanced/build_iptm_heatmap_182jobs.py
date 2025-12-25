#!/usr/bin/env python3
"""
build_iptm_heatmap_182jobs.py

Build heatmap from the 182 unique AF3 jobs (not 1187 CIF files).
Properly deduplicate by job+ligand, removing timestamp duplicates.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Load data
df = pd.read_csv('/storage/kiran-stuff/aaRS/phase2/figures/data/comprehensive_ligand_analysis.csv')

print(f"Raw CIF rows: {len(df)}")

# Filter OUT seed-level files
df = df[~df['job_name'].str.contains('seed-', na=False)]
df = df[df['job_name'].notna()]

# Remove timestamp duplicates
df['job_base'] = df['job_name'].str.replace(r'_\d{8}_\d{6}$', '', regex=True)

# Remove duplicates and entries with nan iptm
df = df[df['AA_iptm'].notna()]
df = df.drop_duplicates(subset=['job_base', 'ligand'], keep='first')

# Remove non-job entries
df = df[~df['job_name'].str.contains('final_results', case=False, na=False)]

print(f"Unique AF3 jobs: {len(df)}")

# Extract construct name (remove ligand suffix)
def get_construct(job):
    parts = job.rsplit('_', 1)
    ligands = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
               'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
    if len(parts) == 2 and parts[1].upper() in ligands:
        return parts[0]
    return job

df['construct'] = df['job_base'].apply(get_construct)

# Keep only constructs with 15+ ligands for the main heatmap
construct_counts = df.groupby('construct')['ligand'].nunique()
full_constructs = construct_counts[construct_counts >= 15].index.tolist()

print(f"\nConstructs with 15+ ligands: {len(full_constructs)}")
for c in sorted(full_constructs):
    print(f"  {c}: {construct_counts[c]} ligands")

df_full = df[df['construct'].isin(full_constructs)]

# Create pivot table
pivot = df_full.pivot_table(
    index='construct',
    columns='ligand',
    values='AA_iptm',
    aggfunc='first'
)

print(f"\nPivot shape: {pivot.shape}")

# Define cognate amino acids
COGNATES = {'ProRS': 'PRO', 'ThrRS': 'THR'}

def get_enzyme(c):
    cl = c.lower()
    if 'prors' in cl or 'prours' in cl:
        return 'ProRS'
    elif 'thrrs' in cl:
        return 'ThrRS'
    return None

# Order ligands
ligand_order = ['PRO', 'THR', 'SER', 'ALA', 'GLY', 'VAL', 'ILE', 'LEU',
                'MET', 'CYS', 'TYR', 'PHE', 'TRP', 'ASN', 'GLN',
                'ASP', 'GLU', 'LYS', 'ARG', 'HIS']
ligand_order = [l for l in ligand_order if l in pivot.columns]
pivot = pivot[ligand_order]

# Order constructs
def sort_key(c):
    enzyme = get_enzyme(c) or 'ZZZ'
    era = 0 if 'anc' in c.lower() else 1
    return (enzyme, era, c)

pivot = pivot.reindex(sorted(pivot.index, key=sort_key))

# Create figure
fig, ax = plt.subplots(figsize=(14, 6))

im = ax.imshow(pivot.values, cmap='RdYlBu_r', aspect='auto', vmin=0.3, vmax=1.0)

cbar = plt.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label('ipTM Score', fontsize=11)

ax.set_xticks(range(len(pivot.columns)))
ax.set_xticklabels(pivot.columns, fontsize=10, rotation=45, ha='right')
ax.set_yticks(range(len(pivot.index)))
ax.set_yticklabels(pivot.index, fontsize=10)

# Highlight cognate cells
for i, construct in enumerate(pivot.index):
    enzyme = get_enzyme(construct)
    if enzyme and enzyme in COGNATES:
        cognate = COGNATES[enzyme]
        if cognate in pivot.columns:
            j = list(pivot.columns).index(cognate)
            rect = plt.Rectangle((j-0.5, i-0.5), 1, 1,
                                fill=False, edgecolor='black', linewidth=3)
            ax.add_patch(rect)

# Add cell values
for i in range(len(pivot.index)):
    for j in range(len(pivot.columns)):
        val = pivot.iloc[i, j]
        if pd.notna(val):
            color = 'white' if val > 0.7 else 'black'
            ax.text(j, i, f'{val:.2f}', ha='center', va='center',
                   fontsize=8, color=color)

ax.set_xlabel('Amino Acid Ligand', fontsize=12)
ax.set_ylabel('Protein Construct', fontsize=12)
ax.set_title(f'AF3 ipTM Selectivity Matrix (n=182 jobs)\n(Black boxes = Cognate)',
             fontsize=13, fontweight='bold')

plt.tight_layout()

outpath = Path('/storage/kiran-stuff/aaRS/phase3_balanced/FIG_iptm_182jobs.png')
plt.savefig(outpath, dpi=150, bbox_inches='tight')
plt.savefig(outpath.with_suffix('.pdf'), bbox_inches='tight')
print(f"\nSaved: {outpath}")

# Selectivity summary
print("\n" + "="*60)
print("SELECTIVITY SUMMARY (182 jobs)")
print("="*60)

for construct in pivot.index:
    enzyme = get_enzyme(construct)
    if not enzyme or enzyme not in COGNATES:
        continue

    cognate = COGNATES[enzyme]
    if cognate not in pivot.columns:
        continue

    cog_iptm = pivot.loc[construct, cognate]
    non_cog = pivot.loc[construct, [c for c in pivot.columns if c != cognate]]
    non_cog_mean = non_cog.mean()
    delta = cog_iptm - non_cog_mean

    row = pivot.loc[construct].dropna()
    rank = (row > cog_iptm).sum() + 1

    print(f"\n{construct}:")
    print(f"  {cognate}: {cog_iptm:.3f} (rank #{rank}/{len(row)})")
    print(f"  Non-cog mean: {non_cog_mean:.3f}, Δ = {delta:+.3f}")
