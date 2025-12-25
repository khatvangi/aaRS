#!/usr/bin/env python3
"""
plot_vina_heatmap.py

Create heatmap from Vina scores for 182 AF3 jobs.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Load Vina results
df = pd.read_csv("/storage/kiran-stuff/aaRS/phase2/energy_scoring/vina_182jobs/vina_scores_182jobs.csv")

print(f"Total rows: {len(df)}")
print(f"Successful: {df['vina_score'].notna().sum()}")
print(f"Failed: {df['vina_score'].isna().sum()}")

# Filter successful
df = df[df['vina_score'].notna()]

# Extract construct from job_name
def get_construct(job):
    if pd.isna(job):
        return None
    # Remove ligand suffix
    parts = job.rsplit('_', 1)
    ligands = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
               'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
    if len(parts) == 2 and parts[1].upper() in ligands:
        return parts[0]
    return job

df['construct'] = df['job_name'].apply(get_construct)

# Keep only constructs with 10+ ligands
construct_counts = df.groupby('construct')['ligand'].nunique()
full_constructs = construct_counts[construct_counts >= 10].index.tolist()

print(f"\nConstructs with 10+ ligands: {len(full_constructs)}")
for c in sorted(full_constructs):
    print(f"  {c}: {construct_counts[c]} ligands")

df_full = df[df['construct'].isin(full_constructs)]

# Create pivot table
pivot = df_full.pivot_table(
    index='construct',
    columns='ligand',
    values='vina_score',
    aggfunc='first'
)

print(f"\nPivot shape: {pivot.shape}")

# Cognate amino acids
COGNATES = {'ProRS': 'PRO', 'ThrRS': 'THR'}

def get_enzyme(c):
    cl = c.lower()
    if 'prors' in cl:
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

# Clip extreme values for visualization
pivot_clipped = pivot.clip(lower=-6, upper=6)

# Create figure
fig, ax = plt.subplots(figsize=(14, 6))

# Note: For Vina, MORE NEGATIVE = BETTER BINDING
# So we reverse the colormap
im = ax.imshow(pivot_clipped.values, cmap='RdYlBu', aspect='auto', vmin=-6, vmax=6)

cbar = plt.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label('Vina Score (kcal/mol)\n(More negative = better)', fontsize=10)

ax.set_xticks(range(len(pivot.columns)))
ax.set_xticklabels(pivot.columns, fontsize=10, rotation=45, ha='right')
ax.set_yticks(range(len(pivot.index)))
ax.set_yticklabels(pivot.index, fontsize=9)

# Highlight cognate cells
for i, construct in enumerate(pivot.index):
    enzyme = get_enzyme(construct)
    if enzyme and enzyme in COGNATES:
        cognate = COGNATES[enzyme]
        if cognate in pivot.columns:
            j = list(pivot.columns).index(cognate)
            rect = plt.Rectangle((j-0.5, i-0.5), 1, 1,
                                fill=False, edgecolor='black', linewidth=2.5)
            ax.add_patch(rect)

# Add cell values
for i in range(len(pivot.index)):
    for j in range(len(pivot.columns)):
        val = pivot.iloc[i, j]
        if pd.notna(val):
            # Color text based on value
            color = 'white' if abs(val) > 3 else 'black'
            ax.text(j, i, f'{val:.1f}', ha='center', va='center',
                   fontsize=6, color=color)

ax.set_xlabel('Amino Acid Ligand', fontsize=11)
ax.set_ylabel('Protein Construct', fontsize=11)
ax.set_title('AutoDock Vina Binding Scores (n=155 jobs)\n(Black boxes = Cognate; Blue = stronger binding)',
             fontsize=12, fontweight='bold')

plt.tight_layout()

# Save
outpath = Path('/storage/kiran-stuff/aaRS/phase2/energy_scoring/vina_182jobs/FIG_vina_heatmap.png')
plt.savefig(outpath, dpi=150, bbox_inches='tight')
plt.savefig(outpath.with_suffix('.pdf'), bbox_inches='tight')
print(f"\nSaved: {outpath}")

# Selectivity analysis
print("\n" + "="*60)
print("SELECTIVITY ANALYSIS (Vina)")
print("="*60)

for construct in pivot.index:
    enzyme = get_enzyme(construct)
    if not enzyme or enzyme not in COGNATES:
        continue

    cognate = COGNATES[enzyme]
    if cognate not in pivot.columns:
        continue

    cog_score = pivot.loc[construct, cognate]
    non_cog = pivot.loc[construct, [c for c in pivot.columns if c != cognate]].dropna()
    non_cog_mean = non_cog.mean()

    # For Vina, lower (more negative) is better
    # So cognate should be MORE NEGATIVE than non-cognate for selectivity
    delta = non_cog_mean - cog_score  # Positive delta = cognate binds better

    # Rank (lower score = better, so rank 1 = most negative)
    row = pivot.loc[construct].dropna()
    rank = (row < cog_score).sum() + 1  # Count how many are better (more negative)

    print(f"\n{construct}:")
    print(f"  {cognate}: {cog_score:.2f} kcal/mol (rank #{rank}/{len(row)})")
    print(f"  Non-cog mean: {non_cog_mean:.2f}")
    print(f"  Δ = {delta:+.2f} (positive = cognate binds better)")

# Save summary
summary_path = Path('/storage/kiran-stuff/aaRS/phase2/energy_scoring/vina_182jobs/vina_selectivity_summary.csv')
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
                'vina_score': val
            })

pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
print(f"\nSaved: {summary_path}")
