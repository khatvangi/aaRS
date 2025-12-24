#!/usr/bin/env python3
"""
Generate Figure 1: The Ancestral "Bucket" Problem

Panels:
  C: Ancestral ThrRS binding landscape (all 20 AAs)
  D: Ancestral ProRS binding landscape (all 20 AAs)

Shows that ancestral enzymes could not discriminate cognate from non-cognate.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# Set publication-quality defaults
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5

# Color scheme
COLOR_COGNATE = '#2ecc71'  # Green
COLOR_BETTER = '#e74c3c'   # Red (substrates better than cognate)
COLOR_WORSE = '#95a5a6'    # Gray

def create_bucket_plot(df, cognate_aa, title, output_path):
    """
    Create bar chart showing all 20 AAs ranked by binding affinity.

    Args:
        df: DataFrame with columns ['primary_ligand', 'AA_iptm', 'rank']
        cognate_aa: The cognate amino acid (e.g., 'THR' or 'PRO')
        title: Plot title
        output_path: Where to save the figure
    """

    # Get cognate score and rank
    cognate_row = df[df['primary_ligand'] == cognate_aa].iloc[0]
    cognate_score = cognate_row['AA_iptm']
    cognate_rank = cognate_row['rank']

    # Assign colors
    colors = []
    for _, row in df.iterrows():
        if row['primary_ligand'] == cognate_aa:
            colors.append(COLOR_COGNATE)
        elif row['AA_iptm'] > cognate_score:
            colors.append(COLOR_BETTER)
        else:
            colors.append(COLOR_WORSE)

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))

    # Create bar chart
    x = np.arange(len(df))
    bars = ax.bar(x, df['AA_iptm'], color=colors, edgecolor='black', linewidth=1.5)

    # Highlight cognate bar
    cognate_idx = df[df['primary_ligand'] == cognate_aa].index[0]
    bars[cognate_idx].set_edgecolor('black')
    bars[cognate_idx].set_linewidth(3)

    # Add horizontal line at cognate level
    ax.axhline(cognate_score, color='black', linestyle='--', linewidth=1.5, alpha=0.5)

    # Labels
    ax.set_xlabel('Amino Acid (ranked by binding affinity)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Protein-Ligand ipTM Score', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)

    # X-axis ticks
    ax.set_xticks(x)
    ax.set_xticklabels(df['primary_ligand'], rotation=45, ha='right')

    # Y-axis limits
    ax.set_ylim([0.70, 0.95])
    ax.set_yticks(np.arange(0.70, 0.96, 0.05))

    # Add rank labels on bars
    for i, (idx, row) in enumerate(df.iterrows()):
        height = row['AA_iptm']
        ax.text(i, height + 0.005, f"{int(row['rank'])}",
                ha='center', va='bottom', fontsize=8, fontweight='bold')

    # Add cognate annotation
    cognate_x = df[df['primary_ligand'] == cognate_aa].index[0]
    ax.annotate(f'{cognate_aa}\n(cognate)\nRank {int(cognate_rank)}/{len(df)}',
                xy=(cognate_x, cognate_score),
                xytext=(cognate_x, 0.75),
                fontsize=10, fontweight='bold',
                ha='center',
                bbox=dict(boxstyle='round,pad=0.5', facecolor=COLOR_COGNATE,
                         edgecolor='black', linewidth=2, alpha=0.8),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))

    # Add text box with key finding
    n_better = len(df[df['AA_iptm'] > cognate_score])
    textstr = f'KEY FINDING:\n{n_better} amino acids bind\nBETTER than cognate!'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8, edgecolor='black', linewidth=2)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=props, fontweight='bold')

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=COLOR_COGNATE, edgecolor='black', linewidth=1.5, label='Cognate'),
        mpatches.Patch(facecolor=COLOR_BETTER, edgecolor='black', linewidth=1.5, label='Binds better'),
        mpatches.Patch(facecolor=COLOR_WORSE, edgecolor='black', linewidth=1.5, label='Binds worse')
    ]
    ax.legend(handles=legend_elements, loc='upper right', frameon=True,
             edgecolor='black', framealpha=0.9, fontsize=9)

    # Grid
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

    # Tight layout
    plt.tight_layout()

    # Save
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight', facecolor='white')
    print(f"Saved: {output_path}")

    plt.close()


# ============================================================================
# GENERATE PANELS
# ============================================================================

# Panel C: Ancestral ThrRS
print("Generating Figure 1C: Ancestral ThrRS...")
df_thrrs = pd.read_csv('figures/data/fig1c_anc_thrrs_no_zn.csv')
create_bucket_plot(
    df_thrrs,
    'THR',
    'Ancestral ThrRS: Promiscuous "Bucket" Enzyme',
    'figures/figure1/panel_c_anc_thrrs.png'
)

# Panel D: Ancestral ProRS
print("Generating Figure 1D: Ancestral ProRS...")
df_prors = pd.read_csv('figures/data/fig1d_anc_prors.csv')
create_bucket_plot(
    df_prors,
    'PRO',
    'Ancestral ProRS: Inverted Selectivity',
    'figures/figure1/panel_d_anc_prors.png'
)

print("\n" + "="*80)
print("FIGURE 1 COMPLETE")
print("="*80)
print("\nGenerated panels:")
print("  - Panel C: figures/figure1/panel_c_anc_thrrs.png (+ PDF)")
print("  - Panel D: figures/figure1/panel_d_anc_prors.png (+ PDF)")
print("\nStill needed:")
print("  - Panel A: Phylogenetic tree (from IQ-TREE)")
print("  - Panel B: Domain architecture (BioRender)")
