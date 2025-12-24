#!/usr/bin/env python3
"""
Create MM/GBSA heatmap with ΔΔG (relative binding energies).
Matches ipTM heatmap format with 6 conditions.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re

# Load MM/GBSA scores
df = pd.read_csv("energy_scoring/mmgbsa_scores.csv")

# Filter successful calculations
df = df[df['status'] == 'OK'].copy()

# Filter out extreme energies (numerical issues)
df = df[df['BE_mmgbsa'].abs() < 1e6].copy()

print(f"MM/GBSA data: {len(df)} structures with reasonable energies")

# Extract job name from file path
def extract_job_name(filepath):
    filename = filepath.split('/')[-1].replace('_model.cif', '')
    # Remove seed-X_sample-Y pattern if present
    job_name = re.sub(r'_seed-\d+_sample-\d+$', '', filename)
    return job_name

df['job_name'] = df['file'].apply(extract_job_name)

# Load ipTM data for merging
iptm_df = pd.read_csv("AF3_RESULTS_CORRECTED.csv")

# Merge
merged = df.merge(
    iptm_df[['job_name', 'AA_iptm', 'ptm', 'protein_len', 'ligands']],
    on='job_name',
    how='left'
)

print(f"Merged with ipTM: {len(merged)} structures")
print(f"With ipTM data: {merged['AA_iptm'].notna().sum()}")

# Filter valid data (same criteria as ipTM heatmap)
editing = merged[(merged['protein_len'] == 300) & (merged['AA_iptm'] >= 0.60)].copy()
other = merged[(merged['protein_len'] != 300) & (merged['ptm'] >= 0.65)].copy()
df_good = pd.concat([editing, other], ignore_index=True)

print(f"Filtered data: {len(df_good)} structures with good quality")

# Define amino acids
AMINO_ACIDS = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
               'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

# Define conditions
conditions = {
    'Anc ProRS\nCatalytic': (df_good['protein_len'] == 500),
    'Anc ProRS\nEditing': (df_good['protein_len'] == 300),
    'Modern ProRS\nCatalytic': (df_good['protein_len'] == 572),
    'Anc ThrRS\nno Zn': ((df_good['protein_len'] == 278) &
                        ~df_good['job_name'].str.contains('zn', case=False, na=False)),
    'Anc ThrRS\n+ Zn': ((df_good['protein_len'] == 278) &
                       df_good['job_name'].str.contains('zn', case=False, na=False)),
    'Modern ThrRS\n+ Zn': (df_good['protein_len'] == 401),
}

# Cognate amino acids
cognates = {
    'Anc ProRS\nCatalytic': 'PRO',
    'Anc ProRS\nEditing': 'THR',
    'Modern ProRS\nCatalytic': 'PRO',
    'Anc ThrRS\nno Zn': 'THR',
    'Anc ThrRS\n+ Zn': 'THR',
    'Modern ThrRS\n+ Zn': 'THR',
}

# Build heatmap data
abs_BE = pd.DataFrame(index=AMINO_ACIDS)
ddG = pd.DataFrame(index=AMINO_ACIDS)

for cond_name, mask in conditions.items():
    subset = df_good[mask].copy()

    # Get mean BE per ligand
    mean_BE = subset.groupby('lig_resname')['BE_mmgbsa'].mean()

    # Get cognate BE
    cognate = cognates[cond_name]
    cognate_BE = mean_BE.get(cognate, np.nan)

    # Fill heatmap
    for aa in AMINO_ACIDS:
        BE = mean_BE.get(aa, np.nan)
        abs_BE.loc[aa, cond_name] = BE

        # ΔΔG = BE(AA) - BE(cognate)
        # More positive = weaker binding (less favorable)
        if not np.isnan(BE) and not np.isnan(cognate_BE):
            ddG.loc[aa, cond_name] = BE - cognate_BE
        else:
            ddG.loc[aa, cond_name] = np.nan

    print(f"\n{cond_name}:")
    print(f"  Cognate: {cognate}, BE = {cognate_BE:.1f} kcal/mol")
    print(f"  Ligands with data: {subset['lig_resname'].nunique()}")

# Create figure
fig, axes = plt.subplots(1, 2, figsize=(20, 10))

# Panel A: Absolute Binding Energies
ax1 = axes[0]
sns.heatmap(abs_BE, annot=True, fmt='.0f', cmap='RdYlGn', center=-30000,
            vmin=-80000, vmax=20000, cbar_kws={'label': 'BE (kcal/mol)'},
            linewidths=0.5, linecolor='gray', ax=ax1)
ax1.set_title('A. Absolute Binding Energies (MM/GBSA)', fontweight='bold', fontsize=14)
ax1.set_xlabel('')
ax1.set_ylabel('Amino Acid', fontweight='bold')
ax1.axvline(x=2, color='black', linewidth=2, linestyle='-')

# Panel B: Relative Binding Energies (ΔΔG)
ax2 = axes[1]
sns.heatmap(ddG, annot=True, fmt='.0f', cmap='RdYlGn_r', center=0,
            vmin=-20000, vmax=20000, cbar_kws={'label': 'ΔΔG from cognate (kcal/mol)'},
            linewidths=0.5, linecolor='gray', ax=ax2)
ax2.set_title('B. Relative Binding Energies (ΔΔG from Cognate)', fontweight='bold', fontsize=14)
ax2.set_xlabel('')
ax2.set_ylabel('')
ax2.axvline(x=2, color='black', linewidth=2, linestyle='-')

# Rotate labels
for ax in axes:
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

plt.tight_layout()
plt.savefig('figures/outputs/15_mmgbsa_heatmap.png', dpi=300, bbox_inches='tight')
print("\n✓ Saved: figures/outputs/15_mmgbsa_heatmap.png")

# Save data
abs_BE.to_csv('energy_scoring/mmgbsa_heatmap_absolute.csv')
ddG.to_csv('energy_scoring/mmgbsa_heatmap_ddG.csv')
print("✓ Saved: energy_scoring/mmgbsa_heatmap_absolute.csv")
print("✓ Saved: energy_scoring/mmgbsa_heatmap_ddG.csv")

# Summary statistics
print("\n" + "="*80)
print("ΔΔG STATISTICS (Modern ThrRS + Zn)")
print("="*80)

modern_thrrs = ddG['Modern ThrRS\n+ Zn'].dropna()
if len(modern_thrrs) > 0:
    print(f"\nTHR (cognate): ΔΔG = 0.0 (reference)")

    # Key comparisons
    for aa in ['SER', 'ILE', 'VAL']:
        if aa in modern_thrrs.index:
            ddg_val = modern_thrrs[aa]
            print(f"{aa}: ΔΔG = {ddg_val:+.1f} kcal/mol")

    # Hydroxyl vs non-hydroxyl
    hydroxyl = ['THR', 'SER', 'TYR']
    non_hydroxyl = ['ILE', 'VAL', 'LEU', 'ALA']

    h_vals = [modern_thrrs.get(aa, np.nan) for aa in hydroxyl if aa in modern_thrrs.index]
    nh_vals = [modern_thrrs.get(aa, np.nan) for aa in non_hydroxyl if aa in modern_thrrs.index]

    if len(h_vals) > 0 and len(nh_vals) > 0:
        print(f"\nHydroxyl AAs: mean ΔΔG = {np.nanmean(h_vals):+.1f} kcal/mol")
        print(f"Non-hydroxyl AAs: mean ΔΔG = {np.nanmean(nh_vals):+.1f} kcal/mol")
        print(f"Difference: {np.nanmean(h_vals) - np.nanmean(nh_vals):+.1f} kcal/mol")

print("\n" + "="*80)
print("COMPLETE!")
print("="*80)
