#!/usr/bin/env python3
"""
Compare energy vs ipTM correlation across all conditions.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

# Load merged data
df = pd.read_csv("energy_scoring/merged_scores_iptm.csv")

# Filter out extreme energies
df = df[df['Eint_kcal'].abs() < 10000].copy()

# Filter valid data
editing = df[(df['protein_len'] == 300) & (df['AA_iptm'] >= 0.60)].copy()
other = df[(df['protein_len'] != 300) & (df['ptm'] >= 0.65)].copy()
df_good = pd.concat([editing, other], ignore_index=True)

print("="*80)
print("ENERGY vs ipTM CORRELATION ANALYSIS")
print("="*80)

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

cognates = {
    'Anc ProRS\nCatalytic': 'PRO',
    'Anc ProRS\nEditing': 'THR',
    'Modern ProRS\nCatalytic': 'PRO',
    'Anc ThrRS\nno Zn': 'THR',
    'Anc ThrRS\n+ Zn': 'THR',
    'Modern ThrRS\n+ Zn': 'THR',
}

# Create 2x3 panel figure
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
axes = axes.flatten()

correlations = []

for idx, (cond_name, mask) in enumerate(conditions.items()):
    subset = df_good[mask].copy()
    cognate = cognates[cond_name]

    # Get per-ligand averages
    ligand_summary = subset.groupby('lig_resname').agg({
        'AA_iptm': 'mean',
        'Eint_kcal': 'mean'
    }).reset_index()

    # Calculate relative energy
    cognate_E = ligand_summary[ligand_summary['lig_resname'] == cognate]['Eint_kcal'].values
    if len(cognate_E) > 0:
        cognate_E = cognate_E[0]
        ligand_summary['dE_from_cognate'] = ligand_summary['Eint_kcal'] - cognate_E
    else:
        ligand_summary['dE_from_cognate'] = np.nan

    # Calculate correlations
    valid = ligand_summary[ligand_summary['AA_iptm'].notna() & ligand_summary['dE_from_cognate'].notna()]

    if len(valid) >= 3:
        r_abs, p_abs = pearsonr(valid['AA_iptm'], valid['Eint_kcal'])
        r_rel, p_rel = pearsonr(valid['AA_iptm'], valid['dE_from_cognate'])
        rho_rel, p_rho = spearmanr(valid['AA_iptm'], valid['dE_from_cognate'])
    else:
        r_abs, p_abs = np.nan, np.nan
        r_rel, p_rel = np.nan, np.nan
        rho_rel, p_rho = np.nan, np.nan

    correlations.append({
        'Condition': cond_name.replace('\n', ' '),
        'N': len(valid),
        'r_abs': r_abs,
        'p_abs': p_abs,
        'r_rel': r_rel,
        'p_rel': p_rel,
        'rho_rel': rho_rel,
        'p_rho': p_rho
    })

    # Plot
    ax = axes[idx]

    # Scatter plot with relative energy
    ax.scatter(valid['dE_from_cognate'], valid['AA_iptm'], s=100, alpha=0.7, c='steelblue', edgecolors='black')

    # Label points
    for _, row in valid.iterrows():
        if row['lig_resname'] in [cognate, 'SER', 'THR', 'PRO', 'ILE', 'VAL']:
            ax.annotate(row['lig_resname'], (row['dE_from_cognate'], row['AA_iptm']),
                       fontsize=8, fontweight='bold', ha='center', va='bottom')

    # Reference lines
    ax.axvline(x=0, color='red', linestyle='--', alpha=0.5, linewidth=1)
    ax.axhline(y=0.85, color='orange', linestyle='--', alpha=0.5, linewidth=1)

    # Add correlation text
    if not np.isnan(r_rel):
        ax.text(0.05, 0.95, f'r = {r_rel:.3f}\np = {p_rel:.3f}\nρ = {rho_rel:.3f}',
               transform=ax.transAxes, fontsize=9, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.set_xlabel('ΔE from cognate (kcal/mol)', fontweight='bold')
    ax.set_ylabel('ipTM Score', fontweight='bold')
    ax.set_title(cond_name, fontweight='bold', fontsize=12)
    ax.grid(True, alpha=0.3)

    print(f"\n{cond_name.replace(chr(10), ' ')}:")
    print(f"  N = {len(valid)} ligands")
    print(f"  Pearson r (abs E): {r_abs:.3f} (p = {p_abs:.4f})")
    print(f"  Pearson r (ΔE): {r_rel:.3f} (p = {p_rel:.4f})")
    print(f"  Spearman ρ (ΔE): {rho_rel:.3f} (p = {p_rho:.4f})")

plt.tight_layout()
plt.savefig('figures/outputs/14_energy_iptm_correlation.png', dpi=300, bbox_inches='tight')
print("\n✓ Saved: figures/outputs/14_energy_iptm_correlation.png")

# Save correlation summary
corr_df = pd.DataFrame(correlations)
corr_df.to_csv('energy_scoring/correlation_summary.csv', index=False)
print("✓ Saved: energy_scoring/correlation_summary.csv")

# Overall summary
print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print("\nCorrelation strength interpretation:")
print("  |r| < 0.3: Weak")
print("  0.3 ≤ |r| < 0.7: Moderate")
print("  |r| ≥ 0.7: Strong")
print("\nSignificance: p < 0.05 is statistically significant")
print("\n" + corr_df.to_string(index=False))

# Identify conditions with significant correlations
sig = corr_df[corr_df['p_rel'] < 0.05]
print(f"\n\nConditions with significant correlation (p < 0.05): {len(sig)}/{len(corr_df)}")
if len(sig) > 0:
    print(sig[['Condition', 'r_rel', 'p_rel']].to_string(index=False))

print("\n" + "="*80)
print("INTERPRETATION")
print("="*80)
print("""
The correlation between energy and ipTM varies by condition:

- Strong correlation: Energy and ipTM predict similar binding rankings
- Weak correlation: Energy captures different aspects than ipTM (e.g., electrostatics vs geometry)
- Negative correlation: Higher energy → lower ipTM (expected if energy is repulsive)

For ThrRS + Zn conditions, we expect:
  - Moderate correlation if hydroxyl mechanism dominates (both metrics capture it)
  - Weak correlation if energy captures electrostatics but ipTM captures geometry
""")

print("\n✓ COMPLETE!")
