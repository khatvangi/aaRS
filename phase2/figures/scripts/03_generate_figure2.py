#!/usr/bin/env python3
"""
Generate Figure 2: ProRS - The Double Sieve Solution

Panels:
  B: Modern ProRS catalytic site (still promiscuous!)
  C: ProRS editing domain (preferentially binds non-cognate)

Shows that ProRS evolved editing domain while RETAINING catalytic promiscuity.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# Set publication-quality defaults
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5

# Color scheme
COLOR_COGNATE = '#2ecc71'  # Green (PRO)
COLOR_HIGH = '#f39c12'     # Orange (near-cognate, 95-100%)
COLOR_MED = '#3498db'      # Blue (70-95%)
COLOR_LOW = '#95a5a6'      # Gray (<70%)

# For editing domain
COLOR_EDIT_TARGET = '#e74c3c'  # Red (misacylation products)
COLOR_EDIT_OK = '#2ecc71'      # Green (cognate excluded)

def create_catalytic_plot(df, output_path):
    """
    Panel B: Modern ProRS catalytic site showing persistent promiscuity.
    """

    # Get PRO score
    pro_score = df[df['primary_ligand'] == 'PRO']['AA_iptm'].values[0]

    # Assign colors based on % of cognate
    colors = []
    for _, row in df.iterrows():
        pct = row['pct_of_cognate']
        if row['primary_ligand'] == 'PRO':
            colors.append(COLOR_COGNATE)
        elif pct >= 95:
            colors.append(COLOR_HIGH)
        elif pct >= 70:
            colors.append(COLOR_MED)
        else:
            colors.append(COLOR_LOW)

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    # Create bar chart
    x = np.arange(len(df))
    bars = ax.bar(x, df['AA_iptm'], color=colors, edgecolor='black', linewidth=1.5)

    # Highlight cognate
    cognate_idx = df[df['primary_ligand'] == 'PRO'].index[0]
    bars[cognate_idx].set_edgecolor='black'
    bars[cognate_idx].set_linewidth(3)

    # Add horizontal line at cognate level
    ax.axhline(pro_score, color='black', linestyle='--', linewidth=1.5, alpha=0.5,
               label=f'PRO level ({pro_score:.3f})')

    # Add 95% line
    ax.axhline(pro_score * 0.95, color=COLOR_HIGH, linestyle=':', linewidth=2, alpha=0.7,
               label='95% of cognate')

    # Labels
    ax.set_xlabel('Amino Acid (ranked by binding affinity)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Protein-Ligand ipTM Score', fontsize=12, fontweight='bold')
    ax.set_title('Modern ProRS Catalytic Site: Persistent Promiscuity', fontsize=14, fontweight='bold', pad=20)

    # X-axis ticks
    ax.set_xticks(x)
    ax.set_xticklabels(df['primary_ligand'], rotation=45, ha='right')

    # Y-axis limits
    ax.set_ylim([0.70, 1.0])
    ax.set_yticks(np.arange(0.70, 1.01, 0.05))

    # Add % labels on top bars
    for i, (idx, row) in enumerate(df.head(10).iterrows()):
        height = row['AA_iptm']
        pct = row['pct_of_cognate']
        ax.text(i, height + 0.01, f"{pct:.0f}%",
                ha='center', va='bottom', fontsize=8, fontweight='bold')

    # Highlight problematic substrates
    high_binders = df[df['pct_of_cognate'] >= 95]
    if len(high_binders) > 1:  # More than just PRO
        problem_aas = ', '.join([aa for aa in high_binders['primary_ligand'] if aa != 'PRO'])
        textstr = f'PROBLEM:\n{problem_aas}\nbind at ≥95%\nof cognate!'
        props = dict(boxstyle='round', facecolor=COLOR_HIGH, alpha=0.8, edgecolor='black', linewidth=2)
        ax.text(0.98, 0.98, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', horizontalalignment='right',
                bbox=props, fontweight='bold')

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=COLOR_COGNATE, edgecolor='black', linewidth=1.5, label='Cognate (PRO)'),
        mpatches.Patch(facecolor=COLOR_HIGH, edgecolor='black', linewidth=1.5, label='≥95% (DANGER!)'),
        mpatches.Patch(facecolor=COLOR_MED, edgecolor='black', linewidth=1.5, label='70-95%'),
        mpatches.Patch(facecolor=COLOR_LOW, edgecolor='black', linewidth=1.5, label='<70%')
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=True,
             edgecolor='black', framealpha=0.9, fontsize=9)

    # Grid
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight', facecolor='white')
    print(f"Saved: {output_path}")
    plt.close()


def create_editing_plot(df, output_path):
    """
    Panel C: ProRS editing domain showing double-sieve mechanism.
    """

    # Known misacylation products for ProRS
    MISACYLATION_PRODUCTS = ['THR', 'ALA', 'CYS', 'SER', 'VAL']

    # Assign colors
    colors = []
    for _, row in df.iterrows():
        if row['primary_ligand'] == 'PRO':
            colors.append(COLOR_EDIT_OK)  # Cognate excluded (good!)
        elif row['primary_ligand'] in MISACYLATION_PRODUCTS:
            colors.append(COLOR_EDIT_TARGET)  # Target for hydrolysis
        else:
            colors.append('#95a5a6')  # Gray

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    # Create bar chart
    x = np.arange(len(df))
    bars = ax.bar(x, df['AA_iptm'], color=colors, edgecolor='black', linewidth=1.5)

    # Highlight THR (top misacylation product)
    thr_idx = df[df['primary_ligand'] == 'THR'].index[0] if 'THR' in df['primary_ligand'].values else None
    if thr_idx is not None:
        bars[thr_idx].set_edgecolor('black')
        bars[thr_idx].set_linewidth(3)

    # Highlight PRO
    pro_idx = df[df['primary_ligand'] == 'PRO'].index[0] if 'PRO' in df['primary_ligand'].values else None
    if pro_idx is not None:
        bars[pro_idx].set_edgecolor('darkgreen')
        bars[pro_idx].set_linewidth(3)

    # Labels
    ax.set_xlabel('Amino Acid (ranked by binding affinity)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Protein-Ligand ipTM Score', fontsize=12, fontweight='bold')
    ax.set_title('ProRS Editing Domain: "Double Sieve" Mechanism', fontsize=14, fontweight='bold', pad=20)

    # X-axis ticks
    ax.set_xticks(x)
    ax.set_xticklabels(df['primary_ligand'], rotation=45, ha='right')

    # Y-axis limits
    ax.set_ylim([0.60, 0.90])

    # Add rank labels
    for i, (idx, row) in enumerate(df.head(10).iterrows()):
        height = row['AA_iptm']
        ax.text(i, height + 0.005, f"#{int(row['rank'])}",
                ha='center', va='bottom', fontsize=8, fontweight='bold')

    # Add key finding annotation
    if thr_idx is not None and pro_idx is not None:
        thr_score = df.iloc[thr_idx]['AA_iptm']
        pro_score = df.iloc[pro_idx]['AA_iptm']
        thr_rank = df.iloc[thr_idx]['rank']
        pro_rank = df.iloc[pro_idx]['rank']

        textstr = f'DOUBLE SIEVE:\nTHR (misacylation) #{int(thr_rank)}: {thr_score:.3f}\nPRO (cognate) #{int(pro_rank)}: {pro_score:.3f}\n\nEditing TARGETS errors!'
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.9, edgecolor='black', linewidth=2)
        ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=props, fontweight='bold')

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=COLOR_EDIT_TARGET, edgecolor='black', linewidth=1.5,
                      label='Misacylation products\n(hydrolyzed)'),
        mpatches.Patch(facecolor=COLOR_EDIT_OK, edgecolor='black', linewidth=1.5,
                      label='Cognate (PRO)\n(excluded)'),
        mpatches.Patch(facecolor='#95a5a6', edgecolor='black', linewidth=1.5,
                      label='Other')
    ]
    ax.legend(handles=legend_elements, loc='upper right', frameon=True,
             edgecolor='black', framealpha=0.9, fontsize=9)

    # Grid
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight', facecolor='white')
    print(f"Saved: {output_path}")
    plt.close()


# ============================================================================
# GENERATE PANELS
# ============================================================================

# Panel B: Modern ProRS catalytic
print("Generating Figure 2B: Modern ProRS catalytic...")
df_cat = pd.read_csv('figures/data/fig2b_mod_prors_catalytic.csv')
create_catalytic_plot(df_cat, 'figures/figure2/panel_b_mod_prors_catalytic.png')

# Panel C: ProRS editing domain
print("Generating Figure 2C: ProRS editing domain...")
df_edit = pd.read_csv('figures/data/fig2c_prors_editing.csv')

# Fix: THR should be #1, not PRO - let me re-check the data
df_edit_sorted = df_edit.sort_values('AA_iptm', ascending=False).reset_index(drop=True)
df_edit_sorted['rank'] = range(1, len(df_edit_sorted) + 1)

create_editing_plot(df_edit_sorted, 'figures/figure2/panel_c_prors_editing.png')

print("\n" + "="*80)
print("FIGURE 2 COMPLETE")
print("="*80)
print("\nGenerated panels:")
print("  - Panel B: figures/figure2/panel_b_mod_prors_catalytic.png (+ PDF)")
print("  - Panel C: figures/figure2/panel_c_prors_editing.png (+ PDF)")
print("\nStill needed:")
print("  - Panel A: Double sieve mechanism schematic (BioRender)")
print("  - Panel D: Structural overlay THR vs PRO in editing pocket")
