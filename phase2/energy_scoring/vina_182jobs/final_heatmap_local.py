#!/usr/bin/env python3
"""
Final Vina Heatmap - ALL jobs with local minimization (consistent methodology)
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from matplotlib.colors import LinearSegmentedColormap

OUTPUT_DIR = Path("/storage/kiran-stuff/aaRS/phase2/energy_scoring/vina_182jobs")

# Load LOCAL minimized data
df = pd.read_csv(OUTPUT_DIR / "vina_scores_LOCAL.csv")

print("=" * 60)
print("FINAL HEATMAP - Local Minimization (Consistent)")
print("=" * 60)
print(f"Total jobs: {len(df)}")
print(f"Score range: {df['vina_local'].min():.2f} to {df['vina_local'].max():.2f}")
print(f"Positive scores: {(df['vina_local'] > 0).sum()}")

# Extract construct
def get_construct(job):
    if pd.isna(job):
        return None
    parts = job.rsplit('_', 1)
    ligands = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
               'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
    if len(parts) == 2 and parts[1].upper() in ligands:
        return parts[0]
    return job

df['construct'] = df['job_name'].apply(get_construct)

# Filter to constructs with 10+ ligands
construct_counts = df.groupby('construct')['ligand'].nunique()
full_constructs = construct_counts[construct_counts >= 10].index.tolist()
df_full = df[df['construct'].isin(full_constructs)]

print(f"Constructs with 10+ ligands: {len(full_constructs)}")

# Create pivot using LOCAL scores
pivot = df_full.pivot_table(index='construct', columns='ligand', values='vina_local', aggfunc='first')

# Order ligands
ligand_order = ['PRO', 'THR', 'SER', 'ALA', 'GLY', 'VAL', 'ILE', 'LEU',
                'MET', 'CYS', 'TYR', 'PHE', 'TRP', 'ASN', 'GLN',
                'ASP', 'GLU', 'LYS', 'ARG', 'HIS']
ligand_order = [l for l in ligand_order if l in pivot.columns]
pivot = pivot[ligand_order]

# Order constructs by enzyme type
def get_enzyme(c):
    cl = c.lower()
    if 'prors' in cl: return 'ProRS'
    elif 'thrrs' in cl: return 'ThrRS'
    return 'Other'

def sort_key(c):
    enzyme = get_enzyme(c)
    order = {'ProRS': 0, 'ThrRS': 1, 'Other': 2}
    era = 0 if 'anc' in c.lower() else 1
    return (order.get(enzyme, 2), era, c)

pivot = pivot.reindex(sorted(pivot.index, key=sort_key))

COGNATES = {'ProRS': 'PRO', 'ThrRS': 'THR'}

# === CREATE FIGURE ===
fig, ax = plt.subplots(figsize=(14, 6))

# Cap values: -6 to +2
pivot_capped = pivot.clip(lower=-6, upper=2)

# Custom colormap: blue (binding) -> white -> red (rejection)
colors = ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0', '#f7f7f7',
          '#fddbc7', '#f4a582', '#d6604d', '#b2182b']
cmap = LinearSegmentedColormap.from_list('binding', colors, N=256)

im = ax.imshow(pivot_capped.values, cmap=cmap, aspect='auto', vmin=-6, vmax=2)

cbar = plt.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label('Vina Score (kcal/mol)\nBlue=Binding | Red=Rejection', fontsize=10)
cbar.set_ticks([-6, -4, -2, 0, 2])
cbar.set_ticklabels(['-6', '-4', '-2', '0', '>0'])

# Axis labels
ax.set_xticks(range(len(pivot.columns)))
ax.set_xticklabels(pivot.columns, fontsize=10, rotation=45, ha='right')
ax.set_yticks(range(len(pivot.index)))
ax.set_yticklabels(pivot.index, fontsize=9)

# Highlight cognate cells
for i, construct in enumerate(pivot.index):
    enzyme = get_enzyme(construct)
    if enzyme in COGNATES:
        cognate = COGNATES[enzyme]
        if cognate in pivot.columns:
            j = list(pivot.columns).index(cognate)
            rect = plt.Rectangle((j-0.5, i-0.5), 1, 1,
                                fill=False, edgecolor='gold', linewidth=3)
            ax.add_patch(rect)

# Add cell values
for i in range(len(pivot.index)):
    for j in range(len(pivot.columns)):
        val = pivot.iloc[i, j]
        if pd.notna(val):
            if val < -3:
                color = 'white'
            elif val < 0:
                color = 'black'
            else:
                color = 'white'
            if val > 0:
                txt = f'+{val:.0f}' if val >= 1 else f'+{val:.1f}'
            else:
                txt = f'{val:.1f}'
            ax.text(j, i, txt, ha='center', va='center', fontsize=6, color=color, fontweight='bold')

ax.set_ylabel('Protein Construct', fontsize=11)
ax.set_xlabel('Amino Acid Ligand', fontsize=11)
ax.set_title('AutoDock Vina Binding Energies (Local Minimization)\n(Gold = Cognate substrate)',
             fontsize=13, fontweight='bold')

plt.tight_layout()

outpath = OUTPUT_DIR / "FIG_vina_LOCAL_final.png"
plt.savefig(outpath, dpi=200, bbox_inches='tight', facecolor='white')
plt.savefig(outpath.with_suffix('.pdf'), bbox_inches='tight', facecolor='white')
print(f"\nSaved: {outpath}")
plt.close()

# === SELECTIVITY SUMMARY ===
print("\n" + "=" * 70)
print("SELECTIVITY SUMMARY (Local Minimization)")
print("=" * 70)
print(f"\n{'Construct':<30} {'Cognate':<6} {'Score':>8} {'Non-cog':>10} {'Δ':>8}")
print("-" * 70)

for construct in pivot.index:
    enzyme = get_enzyme(construct)
    if enzyme not in COGNATES:
        continue
    cognate = COGNATES[enzyme]
    if cognate not in pivot.columns:
        continue
    cog_score = pivot.loc[construct, cognate]
    non_cog = pivot.loc[construct, [c for c in pivot.columns if c != cognate]].dropna()
    non_cog_mean = non_cog.mean()
    delta = non_cog_mean - cog_score
    print(f"{construct:<30} {cognate:<6} {cog_score:>8.2f} {non_cog_mean:>10.2f} {delta:>+8.2f}")

print("\n" + "=" * 70)
print("INTERPRETATION")
print("=" * 70)
print("Δ > 0: Cognate binds BETTER than average non-cognate (selective)")
print("Δ < 0: Cognate binds WORSE than average non-cognate (promiscuous)")
