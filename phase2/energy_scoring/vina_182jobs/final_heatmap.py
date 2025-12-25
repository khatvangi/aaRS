#!/usr/bin/env python3
"""
Final Vina Heatmap for Publication

- Cognate clashes fixed with local minimization
- Decoy clashes treated as "Steric Rejection" (capped at +2 for visualization)
- Proper color scale showing binding (blue) vs rejection (red)
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

OUTPUT_DIR = Path("/storage/kiran-stuff/aaRS/phase2/energy_scoring/vina_182jobs")

# Load data
df = pd.read_csv(OUTPUT_DIR / "vina_scores_182jobs_fixed.csv")

# Apply cognate fixes from local minimization
cognate_fixes = {
    ('proRS_modern_Pro', 'PRO'): -2.939,
    ('COMPETITION_FULL_anc_prors_PRO_vs_GLU', 'PRO'): -4.073,
    ('deep_thrrs_thr', 'THR'): -4.001,
}

for (job, lig), new_score in cognate_fixes.items():
    mask = (df['job_name'] == job) & (df['ligand'] == lig)
    old_score = df.loc[mask, 'vina_score'].values[0] if mask.any() else None
    df.loc[mask, 'vina_score'] = new_score
    print(f"Fixed cognate: {job} + {lig}: {old_score:.2f} → {new_score:.2f}")

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

# Create pivot
pivot = df_full.pivot_table(index='construct', columns='ligand', values='vina_score', aggfunc='first')

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

# Cognate definitions
COGNATES = {'ProRS': 'PRO', 'ThrRS': 'THR'}

print(f"\nPivot shape: {pivot.shape}")
print(f"Constructs: {list(pivot.index)}")

# === CREATE FIGURE ===
fig, ax = plt.subplots(figsize=(14, 6))

# Cap values for visualization: binding (-6 to 0 blue) vs rejection (0 to +2 red)
pivot_capped = pivot.clip(lower=-6, upper=2)

# Custom colormap: blue (binding) -> white (neutral) -> red (rejection)
from matplotlib.colors import LinearSegmentedColormap
colors = ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0', '#f7f7f7',
          '#fddbc7', '#f4a582', '#d6604d', '#b2182b']
cmap = LinearSegmentedColormap.from_list('binding', colors, N=256)

im = ax.imshow(pivot_capped.values, cmap=cmap, aspect='auto', vmin=-6, vmax=2)

cbar = plt.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label('Vina Score (kcal/mol)\nBlue=Binding | Red=Rejection', fontsize=10)
cbar.set_ticks([-6, -4, -2, 0, 2])
cbar.set_ticklabels(['-6', '-4', '-2', '0', '>0\n(Steric\nRejection)'])

# Axis labels
ax.set_xticks(range(len(pivot.columns)))
ax.set_xticklabels(pivot.columns, fontsize=10, rotation=45, ha='right')
ax.set_yticks(range(len(pivot.index)))
ax.set_yticklabels(pivot.index, fontsize=9)

# Highlight cognate cells with gold border
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
            # Color based on value
            if val < -3:
                color = 'white'
            elif val < 0:
                color = 'black'
            else:
                color = 'white'

            # Show actual value or "+" for rejections
            if val > 0:
                txt = f'+{val:.0f}' if val >= 1 else f'+{val:.1f}'
            else:
                txt = f'{val:.1f}'

            ax.text(j, i, txt, ha='center', va='center', fontsize=6, color=color, fontweight='bold')

# Add enzyme type labels on right
ax.set_ylabel('Protein Construct', fontsize=11)
ax.set_xlabel('Amino Acid Ligand', fontsize=11)

# Title
ax.set_title('AutoDock Vina Binding Energy Landscape\n(Gold boxes = Cognate substrates)',
             fontsize=13, fontweight='bold')

plt.tight_layout()

# Save
outpath = OUTPUT_DIR / "FIG_vina_final.png"
plt.savefig(outpath, dpi=200, bbox_inches='tight', facecolor='white')
plt.savefig(outpath.with_suffix('.pdf'), bbox_inches='tight', facecolor='white')
print(f"\nSaved: {outpath}")
plt.close()

# === SELECTIVITY SUMMARY TABLE ===
print("\n" + "=" * 70)
print("FINAL SELECTIVITY SUMMARY")
print("=" * 70)

print(f"\n{'Construct':<30} {'Cognate':<6} {'Score':>8} {'Non-cog Mean':>12} {'Δ':>8}")
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

    # For non-cognate mean, use all values (including positive = rejections)
    non_cog_mean = non_cog.mean()
    delta = non_cog_mean - cog_score

    era = 'anc' if 'anc' in construct.lower() else 'mod'
    print(f"{construct:<30} {cognate:<6} {cog_score:>8.2f} {non_cog_mean:>12.2f} {delta:>+8.2f}")

# Save updated CSV
df.to_csv(OUTPUT_DIR / "vina_scores_FINAL.csv", index=False)
print(f"\nSaved: {OUTPUT_DIR / 'vina_scores_FINAL.csv'}")
