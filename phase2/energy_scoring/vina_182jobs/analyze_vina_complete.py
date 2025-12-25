#!/usr/bin/env python3
"""
analyze_vina_complete.py

Complete analysis of Vina scores:
1. Filtered heatmap (exclude positive scores / clashes)
2. Correlation with AF3 ipTM
3. Ancestral vs modern selectivity comparison
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from scipy import stats

OUTPUT_DIR = Path("/storage/kiran-stuff/aaRS/phase2/energy_scoring/vina_182jobs")

# Load fixed Vina results
df = pd.read_csv(OUTPUT_DIR / "vina_scores_182jobs_fixed.csv")

print("=" * 60)
print("VINA ANALYSIS - COMPLETE")
print("=" * 60)
print(f"Total jobs: {len(df)}")
print(f"Successful: {df['vina_score'].notna().sum()}")

# ============================================================
# 1. FILTERED HEATMAP
# ============================================================
print("\n" + "=" * 60)
print("1. GENERATING FILTERED HEATMAP")
print("=" * 60)

# Extract construct from job_name
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

# Keep only constructs with 10+ ligands
construct_counts = df.groupby('construct')['ligand'].nunique()
full_constructs = construct_counts[construct_counts >= 10].index.tolist()

print(f"Constructs with 10+ ligands: {len(full_constructs)}")

df_full = df[df['construct'].isin(full_constructs)]

# Create pivot table
pivot = df_full.pivot_table(
    index='construct',
    columns='ligand',
    values='vina_score',
    aggfunc='first'
)

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

# --- FILTERED VERSION: Replace positive scores with NaN for visualization ---
pivot_filtered = pivot.copy()
n_positive = (pivot_filtered > 0).sum().sum()
print(f"Positive scores (clashes) being filtered: {n_positive}")
pivot_filtered[pivot_filtered > 0] = np.nan

# Create figure - filtered heatmap
fig, ax = plt.subplots(figsize=(14, 6))

# Use clipped values for colormap
pivot_clipped = pivot_filtered.clip(lower=-6, upper=0)

im = ax.imshow(pivot_clipped.values, cmap='Blues_r', aspect='auto', vmin=-6, vmax=0)

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
                                fill=False, edgecolor='red', linewidth=2.5)
            ax.add_patch(rect)

# Add cell values
for i in range(len(pivot.index)):
    for j in range(len(pivot.columns)):
        val = pivot_filtered.iloc[i, j]
        if pd.notna(val):
            color = 'white' if val < -3 else 'black'
            ax.text(j, i, f'{val:.1f}', ha='center', va='center',
                   fontsize=6, color=color)
        elif pd.notna(pivot.iloc[i, j]):
            # Mark clashes with X
            ax.text(j, i, '✗', ha='center', va='center',
                   fontsize=8, color='red')

ax.set_xlabel('Amino Acid Ligand', fontsize=11)
ax.set_ylabel('Protein Construct', fontsize=11)
ax.set_title('AutoDock Vina Binding Scores (Filtered: clashes marked ✗)\n(Red boxes = Cognate; Darker blue = stronger binding)',
             fontsize=12, fontweight='bold')

plt.tight_layout()
outpath = OUTPUT_DIR / "FIG_vina_heatmap_filtered.png"
plt.savefig(outpath, dpi=150, bbox_inches='tight')
plt.savefig(outpath.with_suffix('.pdf'), bbox_inches='tight')
print(f"Saved: {outpath}")
plt.close()

# ============================================================
# 2. CORRELATION WITH AF3 ipTM
# ============================================================
print("\n" + "=" * 60)
print("2. VINA vs AF3 ipTM CORRELATION")
print("=" * 60)

# Filter to valid scores (negative Vina = reasonable binding)
df_valid = df[df['vina_score'].notna() & (df['vina_score'] < 0)].copy()
df_valid = df_valid[df_valid['iptm'].notna()]

print(f"Valid data points (negative Vina, has ipTM): {len(df_valid)}")

# Correlation
r, p = stats.pearsonr(df_valid['vina_score'], df_valid['iptm'])
print(f"Pearson r = {r:.3f}, p = {p:.2e}")

rho, p_rho = stats.spearmanr(df_valid['vina_score'], df_valid['iptm'])
print(f"Spearman rho = {rho:.3f}, p = {p_rho:.2e}")

# Scatter plot
fig, ax = plt.subplots(figsize=(8, 6))

# Color by enzyme type
colors = []
for _, row in df_valid.iterrows():
    enzyme = get_enzyme(row['construct']) if pd.notna(row.get('construct')) else get_enzyme(row['job_name'])
    if enzyme == 'ProRS':
        colors.append('blue')
    elif enzyme == 'ThrRS':
        colors.append('green')
    else:
        colors.append('gray')

ax.scatter(df_valid['vina_score'], df_valid['iptm'], c=colors, alpha=0.6, s=50)

# Regression line
z = np.polyfit(df_valid['vina_score'], df_valid['iptm'], 1)
p_line = np.poly1d(z)
x_range = np.linspace(df_valid['vina_score'].min(), df_valid['vina_score'].max(), 100)
ax.plot(x_range, p_line(x_range), 'r--', linewidth=2, label=f'r={r:.3f}')

ax.set_xlabel('Vina Score (kcal/mol)', fontsize=12)
ax.set_ylabel('AF3 ipTM', fontsize=12)
ax.set_title('Vina Score vs AF3 ipTM Correlation\n(Excluding clashes)', fontsize=12, fontweight='bold')
ax.legend()

# Add legend for enzyme types
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label='ProRS'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=10, label='ThrRS'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=10, label='Other')
]
ax.legend(handles=legend_elements, loc='upper left')

# Add correlation text
ax.text(0.95, 0.05, f'Pearson r = {r:.3f}\nSpearman ρ = {rho:.3f}',
        transform=ax.transAxes, ha='right', va='bottom',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
outpath = OUTPUT_DIR / "FIG_vina_vs_iptm.png"
plt.savefig(outpath, dpi=150, bbox_inches='tight')
plt.savefig(outpath.with_suffix('.pdf'), bbox_inches='tight')
print(f"Saved: {outpath}")
plt.close()

# ============================================================
# 3. ANCESTRAL vs MODERN SELECTIVITY
# ============================================================
print("\n" + "=" * 60)
print("3. ANCESTRAL vs MODERN SELECTIVITY")
print("=" * 60)

# Calculate selectivity for each construct
selectivity_data = []

for construct in pivot.index:
    enzyme = get_enzyme(construct)
    if not enzyme or enzyme not in COGNATES:
        continue

    cognate = COGNATES[enzyme]
    if cognate not in pivot.columns:
        continue

    cog_score = pivot.loc[construct, cognate]
    if pd.isna(cog_score) or cog_score > 0:
        continue

    non_cog = pivot.loc[construct, [c for c in pivot.columns if c != cognate]]
    # Only use negative (valid) non-cognate scores
    non_cog_valid = non_cog[non_cog < 0].dropna()

    if len(non_cog_valid) < 5:
        continue

    non_cog_mean = non_cog_valid.mean()
    delta = non_cog_mean - cog_score  # Positive = cognate binds better

    # Rank (lower score = better)
    row = pivot.loc[construct].dropna()
    row_valid = row[row < 0]
    rank = (row_valid < cog_score).sum() + 1

    era = 'ancestral' if 'anc' in construct.lower() else 'modern'

    selectivity_data.append({
        'construct': construct,
        'enzyme': enzyme,
        'era': era,
        'cognate': cognate,
        'cog_score': cog_score,
        'non_cog_mean': non_cog_mean,
        'selectivity': delta,
        'rank': rank,
        'n_valid': len(row_valid)
    })

sel_df = pd.DataFrame(selectivity_data)
print(sel_df.to_string(index=False))

# Compare ancestral vs modern
print("\n--- Summary by Era ---")
for enzyme in ['ProRS', 'ThrRS']:
    print(f"\n{enzyme}:")
    for era in ['ancestral', 'modern']:
        subset = sel_df[(sel_df['enzyme'] == enzyme) & (sel_df['era'] == era)]
        if len(subset) > 0:
            print(f"  {era}: selectivity = {subset['selectivity'].mean():.2f} ± {subset['selectivity'].std():.2f}")

# Bar plot
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for idx, enzyme in enumerate(['ProRS', 'ThrRS']):
    ax = axes[idx]
    subset = sel_df[sel_df['enzyme'] == enzyme].copy()

    if len(subset) == 0:
        continue

    # Sort by era then construct
    subset = subset.sort_values(['era', 'construct'])

    colors = ['#1f77b4' if e == 'ancestral' else '#ff7f0e' for e in subset['era']]

    bars = ax.bar(range(len(subset)), subset['selectivity'], color=colors)

    ax.set_xticks(range(len(subset)))
    ax.set_xticklabels(subset['construct'], rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Selectivity (ΔVina kcal/mol)')
    ax.set_title(f'{enzyme} Selectivity\n(cognate: {COGNATES[enzyme]})')
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#1f77b4', label='Ancestral'),
        Patch(facecolor='#ff7f0e', label='Modern')
    ]
    ax.legend(handles=legend_elements)

plt.tight_layout()
outpath = OUTPUT_DIR / "FIG_selectivity_anc_vs_modern.png"
plt.savefig(outpath, dpi=150, bbox_inches='tight')
plt.savefig(outpath.with_suffix('.pdf'), bbox_inches='tight')
print(f"\nSaved: {outpath}")
plt.close()

# ============================================================
# 4. SELECTIVITY AGREEMENT: Vina vs ipTM
# ============================================================
print("\n" + "=" * 60)
print("4. SELECTIVITY AGREEMENT: Vina vs ipTM")
print("=" * 60)

# Load comprehensive data with ipTM
comp_df = pd.read_csv("/storage/kiran-stuff/aaRS/phase2/figures/data/comprehensive_ligand_analysis.csv")
comp_df['job_base'] = comp_df['job_name'].str.replace(r'_\d{8}_\d{6}$', '', regex=True)

# Merge with Vina
merged = df.merge(comp_df[['job_base', 'ligand', 'AA_iptm']],
                  left_on=['job_name', 'ligand'],
                  right_on=['job_base', 'ligand'],
                  how='left')

# For each construct, calculate selectivity from both methods
selectivity_comparison = []

for construct in pivot.index:
    enzyme = get_enzyme(construct)
    if not enzyme or enzyme not in COGNATES:
        continue

    cognate = COGNATES[enzyme]

    # Vina selectivity
    if cognate in pivot.columns:
        cog_vina = pivot.loc[construct, cognate]
        non_cog_vina = pivot.loc[construct, [c for c in pivot.columns if c != cognate]]
        non_cog_vina_valid = non_cog_vina[non_cog_vina < 0].dropna()

        if pd.notna(cog_vina) and cog_vina < 0 and len(non_cog_vina_valid) >= 5:
            vina_selectivity = non_cog_vina_valid.mean() - cog_vina
        else:
            vina_selectivity = np.nan
    else:
        vina_selectivity = np.nan

    # ipTM selectivity
    construct_data = merged[merged['job_name'] == construct]
    cog_iptm = construct_data[construct_data['ligand'] == cognate]['AA_iptm'].values
    non_cog_iptm = construct_data[construct_data['ligand'] != cognate]['AA_iptm'].dropna()

    if len(cog_iptm) > 0 and len(non_cog_iptm) >= 5:
        iptm_selectivity = cog_iptm[0] - non_cog_iptm.mean()  # Higher ipTM = better
    else:
        iptm_selectivity = np.nan

    era = 'ancestral' if 'anc' in construct.lower() else 'modern'

    selectivity_comparison.append({
        'construct': construct,
        'enzyme': enzyme,
        'era': era,
        'vina_selectivity': vina_selectivity,
        'iptm_selectivity': iptm_selectivity
    })

selcomp_df = pd.DataFrame(selectivity_comparison)
selcomp_df = selcomp_df.dropna()

if len(selcomp_df) > 2:
    r_sel, p_sel = stats.pearsonr(selcomp_df['vina_selectivity'], selcomp_df['iptm_selectivity'])
    print(f"Selectivity correlation (Vina vs ipTM): r = {r_sel:.3f}, p = {p_sel:.3f}")

    # Scatter plot
    fig, ax = plt.subplots(figsize=(7, 6))

    colors = ['#1f77b4' if e == 'ancestral' else '#ff7f0e' for e in selcomp_df['era']]
    markers = ['o' if e == 'ProRS' else 's' for e in selcomp_df['enzyme']]

    for i, row in selcomp_df.iterrows():
        color = '#1f77b4' if row['era'] == 'ancestral' else '#ff7f0e'
        marker = 'o' if row['enzyme'] == 'ProRS' else 's'
        ax.scatter(row['vina_selectivity'], row['iptm_selectivity'],
                   c=color, marker=marker, s=100, alpha=0.7)

    ax.set_xlabel('Vina Selectivity (ΔVina for cognate)', fontsize=12)
    ax.set_ylabel('ipTM Selectivity (ΔipTM for cognate)', fontsize=12)
    ax.set_title(f'Selectivity Agreement: Vina vs AF3 ipTM\nr = {r_sel:.3f}', fontsize=12, fontweight='bold')
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#1f77b4', markersize=10, label='Ancestral ProRS'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#1f77b4', markersize=10, label='Ancestral ThrRS'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#ff7f0e', markersize=10, label='Modern ProRS'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#ff7f0e', markersize=10, label='Modern ThrRS'),
    ]
    ax.legend(handles=legend_elements, loc='best')

    plt.tight_layout()
    outpath = OUTPUT_DIR / "FIG_selectivity_vina_vs_iptm.png"
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.savefig(outpath.with_suffix('.pdf'), bbox_inches='tight')
    print(f"Saved: {outpath}")
    plt.close()

# Save summary CSV
sel_df.to_csv(OUTPUT_DIR / "selectivity_summary.csv", index=False)
selcomp_df.to_csv(OUTPUT_DIR / "selectivity_comparison.csv", index=False)
print(f"\nSaved: {OUTPUT_DIR / 'selectivity_summary.csv'}")
print(f"Saved: {OUTPUT_DIR / 'selectivity_comparison.csv'}")

print("\n" + "=" * 60)
print("ANALYSIS COMPLETE")
print("=" * 60)
