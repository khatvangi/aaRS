#!/usr/bin/env python3
"""
MM/GBSA heatmap with corrected energies (4 Å cutoff, ε=20).
Only AF3 jobs from AF3_RESULTS_CORRECTED.csv.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re

# Load MM/GBSA scores (AF3 jobs only)
df = pd.read_csv("energy_scoring/mmgbsa_af3only.csv")

# Filter reasonable energies
df = df[(df['BE_mmgbsa'] > -50) & (df['BE_mmgbsa'] < 500)].copy()
print(f"MM/GBSA data: {len(df)} structures with reasonable energies")

# Extract job name
def extract_job_name(filepath):
    filename = filepath.split('/')[-1].replace('_model.cif', '')
    job_name = re.sub(r'_seed-\d+_sample-\d+$', '', filename)
    return job_name

df['job_name'] = df['file'].apply(extract_job_name)

# Load ipTM data
iptm_df = pd.read_csv("AF3_RESULTS_CORRECTED.csv")

# Merge
merged = df.merge(
    iptm_df[['job_name', 'AA_iptm', 'ptm', 'protein_len', 'ligands']],
    on='job_name',
    how='left'
)

print(f"Merged: {len(merged)} structures")

# Define amino acids
AMINO_ACIDS = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
               'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

# Define conditions
conditions = {
    'Anc ProRS\nCatalytic': (merged['protein_len'] == 500),
    'Anc ProRS\nEditing': (merged['protein_len'] == 300),
    'Modern ProRS\nCatalytic': (merged['protein_len'] == 572),
    'Anc ThrRS\nno Zn': ((merged['protein_len'] == 278) &
                        ~merged['job_name'].str.contains('zn', case=False, na=False)),
    'Anc ThrRS\n+ Zn': ((merged['protein_len'] == 278) &
                       merged['job_name'].str.contains('zn', case=False, na=False)),
    'Modern ThrRS\n+ Zn': (merged['protein_len'] == 401),
}

# Cognates
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
    subset = merged[mask].copy()

    # Mean BE per ligand
    mean_BE = subset.groupby('lig_resname')['BE_mmgbsa'].mean()

    # Cognate BE
    cognate = cognates[cond_name]
    cognate_BE = mean_BE.get(cognate, np.nan)

    for aa in AMINO_ACIDS:
        BE = mean_BE.get(aa, np.nan)
        abs_BE.loc[aa, cond_name] = BE

        # ΔΔG
        if not np.isnan(BE) and not np.isnan(cognate_BE):
            ddG.loc[aa, cond_name] = BE - cognate_BE
        else:
            ddG.loc[aa, cond_name] = np.nan

    print(f"\n{cond_name}:")
    print(f"  Cognate: {cognate}, BE = {cognate_BE:.1f} kcal/mol")
    print(f"  Ligands: {subset['lig_resname'].nunique()}")

# Create figure
fig, axes = plt.subplots(1, 2, figsize=(20, 10))

# Panel A: Absolute BE
ax1 = axes[0]
sns.heatmap(abs_BE, annot=True, fmt='.0f', cmap='RdYlGn', center=100,
            vmin=0, vmax=300, cbar_kws={'label': 'BE (kcal/mol)'},
            linewidths=0.5, linecolor='gray', ax=ax1)
ax1.set_title('A. Absolute Binding Energies (MM/GBSA, ε=20, 4Å cutoff)', fontweight='bold', fontsize=14)
ax1.set_xlabel('')
ax1.set_ylabel('Amino Acid', fontweight='bold')
ax1.axvline(x=2, color='black', linewidth=2, linestyle='-')

# Panel B: ΔΔG
ax2 = axes[1]
sns.heatmap(ddG, annot=True, fmt='.0f', cmap='RdYlGn_r', center=0,
            vmin=-100, vmax=100, cbar_kws={'label': 'ΔΔG from cognate (kcal/mol)'},
            linewidths=0.5, linecolor='gray', ax=ax2)
ax2.set_title('B. Relative Binding Energies (ΔΔG from Cognate)', fontweight='bold', fontsize=14)
ax2.set_xlabel('')
ax2.set_ylabel('')
ax2.axvline(x=2, color='black', linewidth=2, linestyle='-')

for ax in axes:
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

plt.tight_layout()
plt.savefig('figures/outputs/16_mmgbsa_corrected_heatmap.png', dpi=300, bbox_inches='tight')
print("\n✓ Saved: figures/outputs/16_mmgbsa_corrected_heatmap.png")

# Save data
abs_BE.to_csv('energy_scoring/mmgbsa_corrected_absolute.csv')
ddG.to_csv('energy_scoring/mmgbsa_corrected_ddG.csv')
print("✓ Saved CSV files")

# Modern ThrRS statistics
print("\n" + "="*80)
print("MODERN ThrRS + Zn RESULTS")
print("="*80)

col = 'Modern ThrRS\n+ Zn'
if col in ddG.columns:
    print(f"\nTHR (cognate): ΔΔG = 0.0 (reference)")
    for aa in ['SER', 'ILE', 'VAL']:
        if aa in ddG.index:
            val = ddG.loc[aa, col]
            if not np.isnan(val):
                print(f"{aa}: ΔΔG = {val:+.1f} kcal/mol")

    # Hydroxyl vs non-hydroxyl
    hydroxyl = ['THR', 'SER', 'TYR']
    non_hydroxyl = ['ILE', 'VAL', 'LEU', 'ALA']

    h_vals = [ddG.loc[aa, col] for aa in hydroxyl if aa in ddG.index and not np.isnan(ddG.loc[aa, col])]
    nh_vals = [ddG.loc[aa, col] for aa in non_hydroxyl if aa in ddG.index and not np.isnan(ddG.loc[aa, col])]

    if len(h_vals) > 0 and len(nh_vals) > 0:
        print(f"\nHydroxyl AAs: mean ΔΔG = {np.mean(h_vals):+.1f} kcal/mol")
        print(f"Non-hydroxyl AAs: mean ΔΔG = {np.mean(nh_vals):+.1f} kcal/mol")
        print(f"Difference: {np.mean(h_vals) - np.mean(nh_vals):+.1f} kcal/mol")

print("\n" + "="*80)
print("COMPLETE!")
print("="*80)
