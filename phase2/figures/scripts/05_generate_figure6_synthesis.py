#!/usr/bin/env python3
"""
Generate Figure 6: Evolutionary Synthesis

A comprehensive view comparing ThrRS and ProRS evolutionary strategies.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import numpy as np

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 9

# Load all data
df_anc_thrrs = pd.read_csv('figures/data/fig1c_anc_thrrs_no_zn.csv')
df_anc_prors = pd.read_csv('figures/data/fig1d_anc_prors.csv')
df_mod_thrrs = pd.read_csv('figures/data/fig4b_mod_thrrs_zn_all.csv')
df_mod_prors = pd.read_csv('figures/data/fig2b_mod_prors_catalytic.csv')

# Create comprehensive figure
fig = plt.figure(figsize=(16, 12))
gs = GridSpec(3, 2, figure=fig, hspace=0.35, wspace=0.3,
              height_ratios=[1, 1, 1.2])

# ==============================================================================
# ROW 1: ANCESTRAL STATE
# ==============================================================================

# Panel A: Ancestral ThrRS
ax1 = fig.add_subplot(gs[0, 0])
df_anc_thrrs_top10 = df_anc_thrrs.head(10)
x1 = np.arange(len(df_anc_thrrs_top10))
colors1 = ['#2ecc71' if aa == 'THR' else '#e74c3c' if rank < 8 else '#95a5a6'
           for aa, rank in zip(df_anc_thrrs_top10['primary_ligand'], df_anc_thrrs_top10['rank'])]

ax1.barh(x1, df_anc_thrrs_top10['AA_iptm'], color=colors1, edgecolor='black', linewidth=1.5)
ax1.set_yticks(x1)
ax1.set_yticklabels([f"#{int(r)}. {aa}" for aa, r in
                     zip(df_anc_thrrs_top10['primary_ligand'], df_anc_thrrs_top10['rank'])])
ax1.invert_yaxis()
ax1.set_xlabel('ipTM Score', fontweight='bold')
ax1.set_title('Ancestral ThrRS:\nTHR ranks #8 (Promiscuous)', fontweight='bold', fontsize=11)
ax1.set_xlim([0.75, 0.90])
ax1.grid(axis='x', alpha=0.3)

# Panel B: Ancestral ProRS
ax2 = fig.add_subplot(gs[0, 1])
df_anc_prors_top10 = df_anc_prors.head(10)
x2 = np.arange(len(df_anc_prors_top10))
colors2 = ['#2ecc71' if aa == 'PRO' else '#e74c3c' if rank < 3 else '#95a5a6'
           for aa, rank in zip(df_anc_prors_top10['primary_ligand'], df_anc_prors_top10['rank'])]

ax2.barh(x2, df_anc_prors_top10['AA_iptm'], color=colors2, edgecolor='black', linewidth=1.5)
ax2.set_yticks(x2)
ax2.set_yticklabels([f"#{int(r)}. {aa}" for aa, r in
                     zip(df_anc_prors_top10['primary_ligand'], df_anc_prors_top10['rank'])])
ax2.invert_yaxis()
ax2.set_xlabel('ipTM Score', fontweight='bold')
ax2.set_title('Ancestral ProRS:\nPRO ranks #3, GLU is #1!', fontweight='bold', fontsize=11)
ax2.set_xlim([0.75, 0.90])
ax2.grid(axis='x', alpha=0.3)

# ==============================================================================
# ROW 2: MODERN STATE
# ==============================================================================

# Panel C: Modern ThrRS
ax3 = fig.add_subplot(gs[1, 0])
df_mod_thrrs_key = df_mod_thrrs[df_mod_thrrs['primary_ligand'].isin(['THR', 'SER', 'ILE', 'VAL', 'ALA'])]
x3 = np.arange(len(df_mod_thrrs_key))
colors3 = ['#2ecc71', '#f39c12', '#3498db', '#3498db', '#3498db']

bars3 = ax3.bar(x3, df_mod_thrrs_key['AA_iptm'], color=colors3, edgecolor='black', linewidth=2)
ax3.set_xticks(x3)
ax3.set_xticklabels(df_mod_thrrs_key['primary_ligand'])
ax3.set_ylabel('ipTM Score', fontweight='bold')
ax3.set_title('Modern ThrRS + Zn:\nTHR #1, ILE rejected (85.6%)', fontweight='bold', fontsize=11)
ax3.set_ylim([0.75, 1.0])

for i, (idx, row) in enumerate(df_mod_thrrs_key.iterrows()):
    ax3.text(i, row['AA_iptm'] + 0.01, f"{row['pct_of_cognate']:.0f}%",
            ha='center', va='bottom', fontsize=9, fontweight='bold')

ax3.axhline(0.95, color='orange', linestyle=':', linewidth=2, alpha=0.7, label='95% threshold')
ax3.grid(axis='y', alpha=0.3)
ax3.legend()

# Panel D: Modern ProRS
ax4 = fig.add_subplot(gs[1, 1])
df_mod_prors_key = df_mod_prors[df_mod_prors['primary_ligand'].isin(['PRO', 'ALA', 'VAL', 'LEU', 'ILE'])]
x4 = np.arange(len(df_mod_prors_key))
colors4 = ['#2ecc71', '#f39c12', '#f39c12', '#f39c12', '#3498db']

bars4 = ax4.bar(x4, df_mod_prors_key['AA_iptm'], color=colors4, edgecolor='black', linewidth=2)
ax4.set_xticks(x4)
ax4.set_xticklabels(df_mod_prors_key['primary_ligand'])
ax4.set_ylabel('ipTM Score', fontweight='bold')
ax4.set_title('Modern ProRS:\nALA, VAL, LEU at 94-98%!', fontweight='bold', fontsize=11)
ax4.set_ylim([0.75, 1.0])

for i, (idx, row) in enumerate(df_mod_prors_key.iterrows()):
    ax4.text(i, row['AA_iptm'] + 0.01, f"{row['pct_of_cognate']:.0f}%",
            ha='center', va='bottom', fontsize=9, fontweight='bold')

ax4.axhline(0.95, color='orange', linestyle=':', linewidth=2, alpha=0.7, label='95% threshold')
ax4.grid(axis='y', alpha=0.3)
ax4.legend()

# ==============================================================================
# ROW 3: SUMMARY TABLE
# ==============================================================================

ax5 = fig.add_subplot(gs[2, :])
ax5.axis('off')

# Create summary table
table_data = [
    ['', 'ThrRS Pathway', 'ProRS Pathway'],
    ['Problem', 'Hydrophobic confusion\n(THR vs ILE/VAL)', 'Charge + size confusion\n(PRO vs GLU, ALA, VAL)'],
    ['Ancestral State', 'THR ranks #8/20\n(7 AAs bind better)', 'PRO ranks #3/20\n(GLU binds 5% better)'],
    ['Solution', 'STRUCTURAL\nZn-mediated filter', 'KINETIC\nEditing domain'],
    ['Modern Catalytic', 'THR: 0.97 (#1)\nILE: 0.83 (#15) ✓', 'PRO: 0.95 (#1)\nALA: 0.93 (98%) ✗\nVAL: 0.92 (97%) ✗'],
    ['Discrimination', 'THR/ILE: 1.17x\nTHR/SER: 1.02x (LEAK!)', 'PRO/ALA: 1.02x\nPRO/VAL: 1.03x'],
    ['Editing Role', 'SECONDARY\n(only for SER)', 'PRIMARY\n(essential for ALL)'],
    ['Chemical\nConstraint', 'Hydrophobics can\'t\ncoordinate Zn', 'No simple structural\nsolution for charge'],
]

# Create table
table = ax5.table(cellText=table_data, cellLoc='left',
                 loc='center', bbox=[0, 0, 1, 1])

table.auto_set_font_size(False)
table.set_fontsize(9)

# Style cells
for i, row in enumerate(table_data):
    for j, cell in enumerate(row):
        cell_obj = table[(i, j)]
        cell_obj.set_edgecolor('black')
        cell_obj.set_linewidth(1.5)

        if i == 0:  # Header
            cell_obj.set_facecolor('#34495e')
            cell_obj.set_text_props(weight='bold', color='white', fontsize=11)
            cell_obj.set_height(0.08)
        elif j == 0:  # Row labels
            cell_obj.set_facecolor('#ecf0f1')
            cell_obj.set_text_props(weight='bold', fontsize=9)
            cell_obj.set_width(0.15)
        elif j == 1:  # ThrRS column
            cell_obj.set_facecolor('#e8f8f5')
            cell_obj.set_height(0.12)
        elif j == 2:  # ProRS column
            cell_obj.set_facecolor('#fef5e7')
            cell_obj.set_height(0.12)

# Overall title
fig.suptitle('Figure 6: Evolutionary Synthesis - Two Solutions to Ancestral Promiscuity',
            fontsize=18, fontweight='bold', y=0.98)

# Save
output_path = 'figures/figure6/comprehensive_synthesis.png'
plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight', facecolor='white')
print(f"Saved: {output_path}")
plt.close()

print("\n" + "="*80)
print("FIGURE 6 COMPLETE")
print("="*80)
print("\nGenerated:")
print("  - figures/figure6/comprehensive_synthesis.png (+ PDF)")
print("\nThis figure shows:")
print("  - Ancestral promiscuity (both enzymes)")
print("  - Modern ThrRS: structural solution (Zn filter)")
print("  - Modern ProRS: kinetic solution (editing)")
print("  - Side-by-side comparison table")
