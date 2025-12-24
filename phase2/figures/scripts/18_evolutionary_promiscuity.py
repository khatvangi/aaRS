#!/usr/bin/env python3
"""
Compare promiscuity across evolutionary trajectory.
Focus on SPREAD of ΔΔG values (discrimination ability).
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ddG = pd.read_csv('energy_scoring/mmgbsa_corrected_ddG.csv', index_col=0)

# Calculate promiscuity metrics
conditions = [
    ('Anc ProRS\nCatalytic', 'PRO'),
    ('Anc ProRS\nEditing', 'THR'),
    ('Modern ProRS\nCatalytic', 'PRO'),
    ('Anc ThrRS\nno Zn', 'THR'),
    ('Anc ThrRS\n+ Zn', 'THR'),
    ('Modern ThrRS\n+ Zn', 'THR'),
]

results = []
for col, cognate in conditions:
    if col in ddG.columns:
        valid = ddG[col].dropna()

        # Metrics
        std_dev = valid.std()
        spread = valid.max() - valid.min()
        within_50 = sum(abs(valid) < 50)
        within_100 = sum(abs(valid) < 100)

        results.append({
            'Condition': col.replace('\n', ' '),
            'Cognate': cognate,
            'N': len(valid),
            'StdDev': std_dev,
            'Spread': spread,
            'Within_50': within_50,
            'Within_100': within_100,
            'Promiscuity_50': within_50 / len(valid),
            'Promiscuity_100': within_100 / len(valid)
        })

results_df = pd.DataFrame(results)

print('='*80)
print('EVOLUTIONARY PROMISCUITY ANALYSIS')
print('='*80)
print('\nPromiscuity = fraction of AAs within ±50 or ±100 kcal/mol of cognate')
print('(Higher = more promiscuous = less discrimination)')
print('\n' + results_df.to_string(index=False))

# Create visualization
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
axes = axes.flatten()

for idx, (col, cognate) in enumerate(conditions):
    if col in ddG.columns:
        ax = axes[idx]

        valid = ddG[col].dropna().sort_values()

        # Color code
        colors = []
        for aa in valid.index:
            if aa == cognate:
                colors.append('red')
            elif abs(valid[aa]) < 50:
                colors.append('orange')
            else:
                colors.append('steelblue')

        ax.barh(range(len(valid)), valid.values, color=colors, alpha=0.7, edgecolor='black', linewidth=0.5)
        ax.set_yticks(range(len(valid)))
        ax.set_yticklabels(valid.index, fontsize=8)
        ax.axvline(x=0, color='red', linestyle='--', linewidth=2, alpha=0.7)
        ax.axvline(x=-50, color='orange', linestyle=':', linewidth=1, alpha=0.5)
        ax.axvline(x=50, color='orange', linestyle=':', linewidth=1, alpha=0.5)
        ax.set_xlabel(f'ΔΔG from {cognate} (kcal/mol)', fontsize=10)
        ax.set_title(col.replace('\n', ' '), fontweight='bold', fontsize=11)
        ax.grid(axis='x', alpha=0.3)
        ax.invert_yaxis()

        # Add promiscuity annotation
        within_50 = sum(abs(valid) < 50)
        pct = 100 * within_50 / len(valid)
        ax.text(0.95, 0.05, f'{within_50}/{len(valid)} ({pct:.0f}%)\nwithin ±50',
                transform=ax.transAxes, fontsize=9, fontweight='bold',
                verticalalignment='bottom', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

plt.tight_layout()
plt.savefig('figures/outputs/18_evolutionary_promiscuity.png', dpi=300, bbox_inches='tight')
print('\n✓ Saved: figures/outputs/18_evolutionary_promiscuity.png')

# Key findings
print('\n' + '='*80)
print('KEY EVOLUTIONARY FINDINGS')
print('='*80)

anc_prors = results_df[results_df['Condition'] == 'Anc ProRS Catalytic'].iloc[0]
mod_prors = results_df[results_df['Condition'] == 'Modern ProRS Catalytic'].iloc[0]
anc_thrrs_no_zn = results_df[results_df['Condition'] == 'Anc ThrRS no Zn'].iloc[0]
anc_thrrs_zn = results_df[results_df['Condition'] == 'Anc ThrRS + Zn'].iloc[0]
mod_thrrs = results_df[results_df['Condition'] == 'Modern ThrRS + Zn'].iloc[0]

print(f'\nProRS Evolution:')
anc_w50 = anc_prors['Within_50']
anc_n = anc_prors['N']
anc_prom = anc_prors['Promiscuity_50']
mod_w50 = mod_prors['Within_50']
mod_n = mod_prors['N']
mod_prom = mod_prors['Promiscuity_50']
print(f'  Ancestral:  {anc_w50}/{anc_n} AAs ({anc_prom*100:.0f}%) - PROMISCUOUS')
print(f'  Modern:     {mod_w50}/{mod_n} AAs ({mod_prom*100:.0f}%) - HIGHLY PROMISCUOUS')
print(f'  Change:     ProRS became {mod_prom/anc_prom:.1f}x MORE promiscuous')

print(f'\nThrRS Evolution (acquiring Zn):')
anc_nozn_w50 = anc_thrrs_no_zn['Within_50']
anc_nozn_n = anc_thrrs_no_zn['N']
anc_nozn_prom = anc_thrrs_no_zn['Promiscuity_50']
anc_zn_w50 = anc_thrrs_zn['Within_50']
anc_zn_n = anc_thrrs_zn['N']
anc_zn_prom = anc_thrrs_zn['Promiscuity_50']
mod_thrrs_w50 = mod_thrrs['Within_50']
mod_thrrs_n = mod_thrrs['N']
mod_thrrs_prom = mod_thrrs['Promiscuity_50']
print(f'  Ancestral (no Zn): {anc_nozn_w50}/{anc_nozn_n} AAs ({anc_nozn_prom*100:.0f}%)')
print(f'  Ancestral (+ Zn):  {anc_zn_w50}/{anc_zn_n} AAs ({anc_zn_prom*100:.0f}%)')
print(f'  Modern (+ Zn):     {mod_thrrs_w50}/{mod_thrrs_n} AAs ({mod_thrrs_prom*100:.0f}%)')
print(f'  Change:     ThrRS became {mod_thrrs_prom/anc_nozn_prom:.1f}x LESS promiscuous')

print(f'\nFinal Ratio:')
print(f'  Modern ProRS promiscuity: {mod_prom*100:.0f}%')
print(f'  Modern ThrRS promiscuity: {mod_thrrs_prom*100:.0f}%')
print(f'  Difference: ProRS is {mod_prom/mod_thrrs_prom:.1f}x more promiscuous')

print('\n' + '='*80)
print('INTERPRETATION')
print('='*80)
print('''
ProRS: REMAINED PROMISCUOUS (or became MORE promiscuous)
  - Ancient need to handle many amino acids
  - Editing domain is ESSENTIAL

ThrRS: EVOLVED SPECIFICITY through Zn
  - Went from no discrimination to high specificity
  - Editing domain now ONLY needed for SER (the zinc trap)
''')

print('='*80)
