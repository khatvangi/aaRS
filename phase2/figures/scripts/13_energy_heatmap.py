#!/usr/bin/env python3
"""
Create energy heatmap matching the ipTM heatmap format.
Shows relative energies (ΔE from cognate) for each condition.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load merged energy + ipTM data
df = pd.read_csv("energy_scoring/merged_scores_iptm.csv")

# Filter out extreme energies (unphysical, likely bad geometries)
# Energies > 10,000 kcal/mol are clearly wrong
df = df[df['Eint_kcal'].abs() < 10000].copy()
print(f"After filtering extreme energies: {len(df)} structures")

# Filter valid data (same criteria as ipTM heatmap)
# For editing domain (protein_len == 300), use AA_iptm filter
editing = df[(df['protein_len'] == 300) & (df['AA_iptm'] >= 0.60)].copy()

# For other enzymes, use pTM filter
other = df[(df['protein_len'] != 300) & (df['ptm'] >= 0.65)].copy()

# Combine
df_good = pd.concat([editing, other], ignore_index=True)

print(f"Filtered data: {len(df_good)} structures with good quality")

# Define amino acids
AMINO_ACIDS = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
               'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

# Define same conditions as ipTM heatmap
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

# Define cognate amino acids for each condition
cognates = {
    'Anc ProRS\nCatalytic': 'PRO',
    'Anc ProRS\nEditing': 'THR',  # Editing binds errors (THR), not cognate
    'Modern ProRS\nCatalytic': 'PRO',
    'Anc ThrRS\nno Zn': 'THR',
    'Anc ThrRS\n+ Zn': 'THR',
    'Modern ThrRS\n+ Zn': 'THR',
}

# Build heatmap data (absolute energies first)
abs_energy = pd.DataFrame(index=AMINO_ACIDS)
rel_energy = pd.DataFrame(index=AMINO_ACIDS)

for cond_name, mask in conditions.items():
    subset = df_good[mask].copy()

    # Get mean energy per ligand (average across samples)
    mean_energies = subset.groupby('lig_resname')['Eint_kcal'].mean()

    # Get cognate energy for this condition
    cognate = cognates[cond_name]
    cognate_E = mean_energies.get(cognate, np.nan)

    # Calculate absolute and relative energies
    for aa in AMINO_ACIDS:
        E = mean_energies.get(aa, np.nan)
        abs_energy.loc[aa, cond_name] = E

        # Relative energy (ΔE from cognate)
        if not np.isnan(E) and not np.isnan(cognate_E):
            rel_energy.loc[aa, cond_name] = E - cognate_E
        else:
            rel_energy.loc[aa, cond_name] = np.nan

    print(f"\n{cond_name}:")
    print(f"  Cognate: {cognate}, E = {cognate_E:.1f} kcal/mol")
    print(f"  Ligands with data: {subset['lig_resname'].nunique()}")

# Create figure with both absolute and relative energy heatmaps
fig, axes = plt.subplots(1, 2, figsize=(20, 10))

# Panel A: Absolute Energies
ax1 = axes[0]
sns.heatmap(abs_energy, annot=True, fmt='.0f', cmap='coolwarm', center=1500,
            vmin=0, vmax=3000, cbar_kws={'label': 'Eint (kcal/mol)'},
            linewidths=0.5, linecolor='gray', ax=ax1)
ax1.set_title('A. Absolute Interaction Energies', fontweight='bold', fontsize=14)
ax1.set_xlabel('')
ax1.set_ylabel('Amino Acid', fontweight='bold')
ax1.axvline(x=2, color='black', linewidth=2, linestyle='-')  # Separator

# Panel B: Relative Energies (ΔE from cognate)
ax2 = axes[1]
sns.heatmap(rel_energy, annot=True, fmt='.0f', cmap='RdYlGn_r', center=0,
            vmin=-1000, vmax=1000, cbar_kws={'label': 'ΔE from cognate (kcal/mol)'},
            linewidths=0.5, linecolor='gray', ax=ax2)
ax2.set_title('B. Relative Energies (ΔE from Cognate)', fontweight='bold', fontsize=14)
ax2.set_xlabel('')
ax2.set_ylabel('')
ax2.axvline(x=2, color='black', linewidth=2, linestyle='-')  # Separator

# Rotate x-axis labels
for ax in axes:
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

plt.tight_layout()
plt.savefig('figures/outputs/13_energy_heatmap.png', dpi=300, bbox_inches='tight')
print("\n✓ Saved: figures/outputs/13_energy_heatmap.png")

# Save data as CSV for reference
abs_energy.to_csv('energy_scoring/heatmap_absolute_energies.csv')
rel_energy.to_csv('energy_scoring/heatmap_relative_energies.csv')
print("✓ Saved: energy_scoring/heatmap_absolute_energies.csv")
print("✓ Saved: energy_scoring/heatmap_relative_energies.csv")

# Print summary statistics
print("\n" + "="*80)
print("RELATIVE ENERGY STATISTICS")
print("="*80)

for col in rel_energy.columns:
    valid = rel_energy[col].dropna()
    if len(valid) > 0:
        print(f"\n{col}:")
        print(f"  Mean ΔE: {valid.mean():.1f} ± {valid.std():.1f} kcal/mol")
        print(f"  Range: {valid.min():.1f} to {valid.max():.1f} kcal/mol")

        # Most/least favorable
        most_fav = valid.idxmin()
        least_fav = valid.idxmax()
        print(f"  Most favorable: {most_fav} (ΔE = {valid[most_fav]:.1f})")
        print(f"  Least favorable: {least_fav} (ΔE = {valid[least_fav]:.1f})")

print("\n" + "="*80)
print("COMPLETE!")
print("="*80)
