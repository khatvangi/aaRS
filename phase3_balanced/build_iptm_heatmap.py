#!/usr/bin/env python3
"""
build_iptm_heatmap.py

Build the "mic drop" heatmap from comprehensive_ligand_analysis.csv
Rows = Protein constructs
Columns = Ligands (all 20 amino acids)
Color = ipTM score (higher = better predicted binding)

Highlight cognate amino acids to show selectivity.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Load data
df = pd.read_csv('/storage/kiran-stuff/aaRS/phase2/figures/data/comprehensive_ligand_analysis.csv')

# Filter to job-level (not seed-level)
df = df[~df['job_name'].str.contains('seed-', na=False)]
df = df[df['job_name'].notna()]

# Extract construct name
def extract_construct(job):
    if pd.isna(job):
        return None
    parts = job.rsplit('_', 1)
    ligands = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
               'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
    if len(parts) == 2 and parts[1].upper() in ligands:
        return parts[0]
    return job

df['construct'] = df['job_name'].apply(extract_construct)

# Keep only constructs with 10+ ligands
construct_counts = df.groupby('construct')['ligand'].nunique()
full_constructs = construct_counts[construct_counts >= 10].index.tolist()

print(f"Constructs with 10+ ligands: {len(full_constructs)}")
for c in full_constructs:
    print(f"  {c}: {construct_counts[c]} ligands")

df = df[df['construct'].isin(full_constructs)]

# Create pivot table
pivot = df.pivot_table(
    index='construct',
    columns='ligand',
    values='AA_iptm',
    aggfunc='mean'
)

print(f"\nPivot shape: {pivot.shape}")
print(pivot)

# Define cognate amino acids for each enzyme type
COGNATES = {
    'ProRS': 'PRO',
    'ThrRS': 'THR',
}

# Determine enzyme type from construct name
def get_enzyme(construct):
    c = construct.lower()
    if 'prors' in c or 'prours' in c:
        return 'ProRS'
    elif 'thrrs' in c:
        return 'ThrRS'
    return None

# Order ligands sensibly
ligand_order = ['PRO', 'THR', 'SER', 'ALA', 'GLY', 'VAL', 'ILE', 'LEU',
                'MET', 'CYS', 'TYR', 'PHE', 'TRP', 'ASN', 'GLN',
                'ASP', 'GLU', 'LYS', 'ARG', 'HIS']
ligand_order = [l for l in ligand_order if l in pivot.columns]
pivot = pivot[ligand_order]

# Order constructs by enzyme then era
def construct_sort_key(c):
    enzyme = get_enzyme(c)
    era = 0 if 'anc' in c.lower() or 'deep' in c.lower() else 1
    return (enzyme or 'ZZZ', era, c)

construct_order = sorted(pivot.index, key=construct_sort_key)
pivot = pivot.reindex(construct_order)

# Create figure
fig, ax = plt.subplots(figsize=(14, 8))

# Plot heatmap
im = ax.imshow(pivot.values, cmap='RdYlBu_r', aspect='auto', vmin=0.3, vmax=1.0)

# Add colorbar
cbar = plt.colorbar(im, ax=ax, shrink=0.7)
cbar.set_label('ipTM Score\n(Interface Predicted TM)', fontsize=11)

# Set axis labels
ax.set_xticks(range(len(pivot.columns)))
ax.set_xticklabels(pivot.columns, fontsize=10, rotation=45, ha='right')
ax.set_yticks(range(len(pivot.index)))
ax.set_yticklabels(pivot.index, fontsize=10)

# Highlight cognate cells with border
for i, construct in enumerate(pivot.index):
    enzyme = get_enzyme(construct)
    if enzyme and enzyme in COGNATES:
        cognate = COGNATES[enzyme]
        if cognate in pivot.columns:
            j = list(pivot.columns).index(cognate)
            # Draw rectangle around cognate
            rect = plt.Rectangle((j-0.5, i-0.5), 1, 1,
                                fill=False, edgecolor='black', linewidth=3)
            ax.add_patch(rect)

# Add values to cells
for i in range(len(pivot.index)):
    for j in range(len(pivot.columns)):
        val = pivot.iloc[i, j]
        if pd.notna(val):
            color = 'white' if val > 0.7 else 'black'
            ax.text(j, i, f'{val:.2f}', ha='center', va='center',
                   fontsize=7, color=color)

ax.set_xlabel('Amino Acid Ligand', fontsize=12)
ax.set_ylabel('Protein Construct', fontsize=12)
ax.set_title('AF3 ipTM Selectivity Matrix\n(Black boxes = Cognate amino acid)',
             fontsize=13, fontweight='bold')

plt.tight_layout()

# Save
outpath = Path('/storage/kiran-stuff/aaRS/phase3_balanced/FIG_iptm_selectivity_heatmap.png')
plt.savefig(outpath, dpi=150, bbox_inches='tight')
plt.savefig(outpath.with_suffix('.pdf'), bbox_inches='tight')
print(f"\nSaved: {outpath}")
print(f"Saved: {outpath.with_suffix('.pdf')}")

# Calculate selectivity statistics
print("\n" + "="*60)
print("SELECTIVITY ANALYSIS")
print("="*60)

for construct in pivot.index:
    enzyme = get_enzyme(construct)
    if not enzyme or enzyme not in COGNATES:
        continue

    cognate = COGNATES[enzyme]
    if cognate not in pivot.columns:
        continue

    cog_iptm = pivot.loc[construct, cognate]
    non_cog_iptm = pivot.loc[construct, [c for c in pivot.columns if c != cognate]].mean()
    delta = cog_iptm - non_cog_iptm

    # Rank of cognate
    row = pivot.loc[construct].dropna()
    rank = (row > cog_iptm).sum() + 1

    print(f"\n{construct}:")
    print(f"  Cognate ({cognate}): {cog_iptm:.3f}")
    print(f"  Non-cognate mean: {non_cog_iptm:.3f}")
    print(f"  Δ = {delta:+.3f}")
    print(f"  Rank: {rank}/{len(row)}")

# Save summary CSV
summary_rows = []
for construct in pivot.index:
    enzyme = get_enzyme(construct)
    cognate = COGNATES.get(enzyme)

    for ligand in pivot.columns:
        val = pivot.loc[construct, ligand]
        if pd.notna(val):
            summary_rows.append({
                'construct': construct,
                'enzyme': enzyme,
                'ligand': ligand,
                'is_cognate': ligand == cognate,
                'iptm': val
            })

summary_df = pd.DataFrame(summary_rows)
summary_path = Path('/storage/kiran-stuff/aaRS/phase3_balanced/iptm_matrix_summary.csv')
summary_df.to_csv(summary_path, index=False)
print(f"\nSaved summary: {summary_path}")
