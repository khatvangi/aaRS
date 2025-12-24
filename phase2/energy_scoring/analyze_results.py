#!/usr/bin/env python3
"""
Analyze energy scoring results and merge with ipTM data.
Create accept/reject classifications.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load data
print("="*80)
print("ENERGY + ipTM ANALYSIS")
print("="*80)

# Energy scores
scores = pd.read_csv("energy_scoring/scores_simple.csv")
scores_ok = scores[scores['status'] == 'OK'].copy()

# Filter numerical outliers (> 1e6 kcal/mol is clearly wrong)
scores_ok = scores_ok[scores_ok['Eint_kcal'] < 1e6]

print(f"\nEnergy scores: {len(scores_ok)} valid structures")

# ipTM data
iptm_df = pd.read_csv("AF3_RESULTS_CORRECTED.csv")

# Merge - extract job name from file path
# For paths like: .../modern_ecoli_full_pro/seed-1_sample-3/modern_ecoli_full_pro_seed-1_sample-3_model.cif
# Extract the job name (component -3, or from filename by removing seed-X_sample-Y_model.cif)
import re
def extract_job_name(filepath):
    # Get filename without extension
    filename = filepath.split('/')[-1].replace('_model.cif', '')
    # Remove seed-X_sample-Y pattern if present
    job_name = re.sub(r'_seed-\d+_sample-\d+$', '', filename)
    return job_name

scores_ok['job_name'] = scores_ok['file'].apply(extract_job_name)

# Merge
merged = scores_ok.merge(
    iptm_df[['job_name', 'AA_iptm', 'ptm', 'protein_len', 'ligands']],
    on='job_name',
    how='left'
)

print(f"Merged with ipTM: {len(merged)} structures")
print(f"With valid ipTM: {merged['AA_iptm'].notna().sum()}")

# Save merged dataset
merged.to_csv("energy_scoring/merged_scores_iptm.csv", index=False)
print(f"\n✓ Saved: energy_scoring/merged_scores_iptm.csv")

# Focus on key conditions
print("\n" + "="*80)
print("MODERN ThrRS + Zn ANALYSIS")
print("="*80)

modern_thrrs = merged[
    merged['job_name'].str.contains('modern_thrrs_ecoli_zn', case=False, na=False)
].copy()

if len(modern_thrrs) > 0:
    # Sort by ipTM
    modern_thrrs = modern_thrrs.sort_values('AA_iptm', ascending=False)

    print(f"\nStructures: {len(modern_thrrs)}")
    print("\nTop 10 by ipTM:")
    print(modern_thrrs[['lig_resname', 'AA_iptm', 'Eint_kcal', 'Evdw_kcal', 'Ecoul_kcal']].head(10).to_string(index=False))

    # Calculate relative energy (ΔE from THR)
    thr_energy = modern_thrrs[modern_thrrs['lig_resname'] == 'THR']['Eint_kcal'].values
    if len(thr_energy) > 0:
        thr_E = thr_energy[0]
        modern_thrrs['dE_from_THR'] = modern_thrrs['Eint_kcal'] - thr_E

        print(f"\n\nRelative energies (ΔE from THR = {thr_E:.1f} kcal/mol):")
        print(modern_thrrs[['lig_resname', 'AA_iptm', 'Eint_kcal', 'dE_from_THR']].head(10).to_string(index=False))

        # Correlation
        valid = modern_thrrs[modern_thrrs['AA_iptm'].notna()]
        if len(valid) > 5:
            corr_iptm_E = valid[['AA_iptm', 'Eint_kcal']].corr().iloc[0, 1]
            corr_iptm_dE = valid[['AA_iptm', 'dE_from_THR']].corr().iloc[0, 1]
            print(f"\nCorrelations:")
            print(f"  ipTM vs Eint: {corr_iptm_E:.3f}")
            print(f"  ipTM vs ΔE: {corr_iptm_dE:.3f}")

# Analyze editing domain
print("\n" + "="*80)
print("EDITING DOMAIN ANALYSIS")
print("="*80)

editing = merged[merged['protein_len'] == 300].copy()
if len(editing) > 0:
    # Get best per ligand
    editing_best = editing.loc[editing.groupby('lig_resname')['AA_iptm'].idxmax()]
    editing_best = editing_best.sort_values('AA_iptm', ascending=False)

    print(f"\nStructures: {len(editing_best)}")
    print("\nTop ligands:")
    print(editing_best[['lig_resname', 'AA_iptm', 'Eint_kcal', 'Evdw_kcal']].head(10).to_string(index=False))

    # Relative to PRO
    pro_energy = editing_best[editing_best['lig_resname'] == 'PRO']['Eint_kcal'].values
    if len(pro_energy) > 0:
        pro_E = pro_energy[0]
        editing_best['dE_from_PRO'] = editing_best['Eint_kcal'] - pro_E

        print(f"\n\nRelative energies (ΔE from PRO = {pro_E:.1f} kcal/mol):")
        print(editing_best[['lig_resname', 'AA_iptm', 'Eint_kcal', 'dE_from_PRO']].head(10).to_string(index=False))

# Create visualization
print("\n" + "="*80)
print("CREATING VISUALIZATION")
print("="*80)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel A: Modern ThrRS - ipTM vs Energy
ax1 = axes[0, 0]
if len(modern_thrrs) > 0:
    valid = modern_thrrs[modern_thrrs['AA_iptm'].notna()]
    ax1.scatter(valid['Eint_kcal'], valid['AA_iptm'], alpha=0.7, s=80)

    # Label key ligands
    for _, row in valid.iterrows():
        if row['lig_resname'] in ['THR', 'SER', 'ILE', 'VAL']:
            ax1.annotate(row['lig_resname'], (row['Eint_kcal'], row['AA_iptm']),
                        fontsize=9, fontweight='bold')

    ax1.set_xlabel('Interaction Energy (kcal/mol)', fontweight='bold')
    ax1.set_ylabel('ipTM Score', fontweight='bold')
    ax1.set_title('A. Modern ThrRS + Zn: Energy vs ipTM', fontweight='bold')
    ax1.grid(True, alpha=0.3)

# Panel B: Modern ThrRS - ΔE vs ipTM
ax2 = axes[0, 1]
if len(modern_thrrs) > 0 and 'dE_from_THR' in modern_thrrs.columns:
    valid = modern_thrrs[modern_thrrs['AA_iptm'].notna()]
    ax2.scatter(valid['dE_from_THR'], valid['AA_iptm'], alpha=0.7, s=80)

    for _, row in valid.iterrows():
        if row['lig_resname'] in ['THR', 'SER', 'ILE', 'VAL']:
            ax2.annotate(row['lig_resname'], (row['dE_from_THR'], row['AA_iptm']),
                        fontsize=9, fontweight='bold')

    ax2.axvline(x=0, color='red', linestyle='--', alpha=0.5)
    ax2.set_xlabel('ΔE from THR (kcal/mol)', fontweight='bold')
    ax2.set_ylabel('ipTM Score', fontweight='bold')
    ax2.set_title('B. Modern ThrRS + Zn: Relative Energy vs ipTM', fontweight='bold')
    ax2.grid(True, alpha=0.3)

# Panel C: Editing Domain
ax3 = axes[1, 0]
if len(editing_best) > 0:
    valid = editing_best[editing_best['AA_iptm'].notna()]
    ax3.scatter(valid['Eint_kcal'], valid['AA_iptm'], alpha=0.7, s=80, color='orange')

    for _, row in valid.iterrows():
        if row['lig_resname'] in ['THR', 'PRO', 'ALA', 'VAL']:
            ax3.annotate(row['lig_resname'], (row['Eint_kcal'], row['AA_iptm']),
                        fontsize=9, fontweight='bold')

    ax3.set_xlabel('Interaction Energy (kcal/mol)', fontweight='bold')
    ax3.set_ylabel('ipTM Score', fontweight='bold')
    ax3.set_title('C. Editing Domain: Energy vs ipTM', fontweight='bold')
    ax3.grid(True, alpha=0.3)

# Panel D: VdW vs Coulomb contributions
ax4 = axes[1, 1]
if len(modern_thrrs) > 0:
    ax4.scatter(modern_thrrs['Evdw_kcal'], modern_thrrs['Ecoul_kcal'], alpha=0.7, s=80)

    for _, row in modern_thrrs.iterrows():
        if row['lig_resname'] in ['THR', 'SER', 'ILE', 'VAL']:
            ax4.annotate(row['lig_resname'], (row['Evdw_kcal'], row['Ecoul_kcal']),
                        fontsize=9, fontweight='bold')

    ax4.set_xlabel('VdW Energy (kcal/mol)', fontweight='bold')
    ax4.set_ylabel('Coulomb Energy (kcal/mol)', fontweight='bold')
    ax4.set_title('D. Modern ThrRS: VdW vs Coulomb', fontweight='bold')
    ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('energy_scoring/energy_iptm_analysis.png', dpi=300, bbox_inches='tight')
print("\n✓ Saved: energy_scoring/energy_iptm_analysis.png")

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)
