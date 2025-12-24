#!/usr/bin/env python3
"""
Generate Zinc Evolution Figures (3, 4, 5)

Figure 3: Zinc Disconnect - Zn binding ≠ Zn discrimination
Figure 4: Zinc Filter - How modern ThrRS uses Zn
Figure 5: Zinc Trap - SER escapes the filter
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.gridspec import GridSpec

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.5

# ==========================================================================
# FIGURE 3: ZINC DISCONNECT
# ==========================================================================

def create_competition_plot(df, output_path):
    """
    Show competition experiments: ancestral vs modern ThrRS with Zn.
    """

    # Parse competition data
    competitions = []
    for _, row in df.iterrows():
        job = row['job_name']
        ligands = row['ligands'].split(',')

        # Parse chain_pair_iptm matrix
        import ast
        cp_matrix = ast.literal_eval(row['chain_pair_iptm'])

        # For competition runs: [Protein, Lig1, Lig2, Zn]
        # cp[0][1] = Protein-Lig1, cp[0][2] = Protein-Lig2

        if 'anc_thrrs_THR_vs_ILE' in job:
            competitions.append({
                'name': 'Ancestral ThrRS\n+ Zn',
                'thr': cp_matrix[0][1],
                'ile': cp_matrix[0][2],
                'era': 'ancestral'
            })
        elif 'modern_thrrs_THR_vs_ILE' in job:
            competitions.append({
                'name': 'Modern ThrRS\n+ Zn',
                'thr': cp_matrix[0][1],
                'ile': cp_matrix[0][2],
                'era': 'modern'
            })

    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    for idx, comp in enumerate(competitions):
        ax = axes[idx]

        # Data
        ligands = ['THR\n(cognate)', 'ILE\n(non-cognate)']
        scores = [comp['thr'], comp['ile']]
        colors = ['#2ecc71', '#e74c3c'] if comp['thr'] > comp['ile'] * 1.05 else ['#f39c12', '#f39c12']

        # Bars
        x = np.arange(len(ligands))
        bars = ax.bar(x, scores, color=colors, edgecolor='black', linewidth=2, width=0.6)

        # Add scores on bars
        for i, (lig, score) in enumerate(zip(ligands, scores)):
            ax.text(i, score + 0.02, f'{score:.3f}',
                   ha='center', va='bottom', fontsize=12, fontweight='bold')

        # Discrimination ratio
        ratio = comp['thr'] / comp['ile']
        ax.axhline(comp['ile'], color='gray', linestyle='--', linewidth=1, alpha=0.5)

        # Title and labels
        ax.set_title(comp['name'], fontsize=14, fontweight='bold', pad=15)
        ax.set_ylabel('Protein-Ligand ipTM Score', fontsize=11, fontweight='bold')
        ax.set_ylim([0.75, 1.0])
        ax.set_xticks(x)
        ax.set_xticklabels(ligands, fontsize=11)

        # Result annotation
        if ratio > 1.05:
            result = f"DISCRIMINATES\nRatio: {ratio:.2f}x"
            bbox_color = '#2ecc71'
        else:
            result = f"DEAD HEAT\nRatio: {ratio:.2f}x"
            bbox_color = '#e74c3c'

        props = dict(boxstyle='round', facecolor=bbox_color, alpha=0.8,
                    edgecolor='black', linewidth=2)
        ax.text(0.5, 0.95, result, transform=ax.transAxes,
               fontsize=11, fontweight='bold', ha='center', va='top',
               bbox=props)

        # Grid
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)

    # Overall title
    fig.suptitle('Figure 3: The "Zinc Disconnect" - Zn Binding ≠ Zn Discrimination',
                fontsize=16, fontweight='bold', y=1.02)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight', facecolor='white')
    print(f"Saved: {output_path}")
    plt.close()


# ==============================================================================
# FIGURE 4: ZINC FILTER HEATMAP
# ==============================================================================

def create_zinc_filter_heatmap(df, output_path):
    """
    Heatmap showing all 20 AAs with modern ThrRS + Zn.
    """

    # Amino acid properties for grouping
    HYDROPHOBIC = ['ALA', 'VAL', 'LEU', 'ILE', 'PRO', 'PHE', 'TRP', 'MET']
    POLAR = ['SER', 'THR', 'CYS', 'ASN', 'GLN', 'TYR', 'GLY']
    CHARGED = ['ARG', 'LYS', 'ASP', 'GLU', 'HIS']

    # Group AAs
    df['group'] = df['primary_ligand'].apply(lambda x:
        'Hydrophobic' if x in HYDROPHOBIC else
        'Polar' if x in POLAR else
        'Charged')

    # Sort by group then score
    df_sorted = df.sort_values(['group', 'AA_iptm'], ascending=[True, False])

    # Create figure with grid
    fig = plt.figure(figsize=(14, 8))
    gs = GridSpec(2, 2, height_ratios=[3, 1], width_ratios=[1, 0.05],
                 hspace=0.3, wspace=0.05)

    # Main heatmap
    ax_heat = fig.add_subplot(gs[0, 0])

    # Prepare data for heatmap
    aa_scores = df_sorted[['primary_ligand', 'AA_iptm', 'Zn_iptm']].values

    # Create 2-row heatmap (AA binding, Zn binding)
    heat_data = np.array([
        [row[1] for row in aa_scores],  # AA ipTM
        [row[2] if not np.isnan(row[2]) else 0 for row in aa_scores]  # Zn ipTM
    ])

    im = ax_heat.imshow(heat_data, cmap='RdYlGn', aspect='auto', vmin=0.75, vmax=1.0)

    # Axes
    ax_heat.set_xticks(np.arange(len(aa_scores)))
    ax_heat.set_xticklabels([row[0] for row in aa_scores], rotation=45, ha='right')
    ax_heat.set_yticks([0, 1])
    ax_heat.set_yticklabels(['AA Binding', 'Zn Binding'], fontsize=11, fontweight='bold')

    # Add scores as text
    for i in range(heat_data.shape[0]):
        for j in range(heat_data.shape[1]):
            text = ax_heat.text(j, i, f'{heat_data[i, j]:.2f}',
                               ha="center", va="center", color="black", fontsize=8, fontweight='bold')

    # Highlight cognate
    thr_idx = list(df_sorted['primary_ligand']).index('THR')
    rect = mpatches.Rectangle((thr_idx - 0.5, -0.5), 1, 2, linewidth=4,
                              edgecolor='blue', facecolor='none')
    ax_heat.add_patch(rect)

    # Colorbar
    ax_cbar = fig.add_subplot(gs[0, 1])
    cbar = plt.colorbar(im, cax=ax_cbar)
    cbar.set_label('ipTM Score', fontsize=11, fontweight='bold')

    # Bar chart below
    ax_bar = fig.add_subplot(gs[1, 0])

    x_bar = np.arange(len(df_sorted))
    colors_bar = ['#2ecc71' if aa == 'THR' else
                  '#f39c12' if pct >= 95 else
                  '#3498db' if pct >= 85 else
                  '#95a5a6'
                  for aa, pct in zip(df_sorted['primary_ligand'], df_sorted['pct_of_cognate'])]

    ax_bar.bar(x_bar, df_sorted['pct_of_cognate'], color=colors_bar, edgecolor='black', linewidth=1)
    ax_bar.set_xticks(x_bar)
    ax_bar.set_xticklabels(df_sorted['primary_ligand'], rotation=45, ha='right')
    ax_bar.set_ylabel('% of Cognate', fontsize=10, fontweight='bold')
    ax_bar.axhline(100, color='black', linestyle='--', linewidth=1.5)
    ax_bar.axhline(95, color='orange', linestyle=':', linewidth=1.5)
    ax_bar.set_ylim([75, 105])
    ax_bar.grid(axis='y', alpha=0.3)

    # Overall title
    fig.suptitle('Figure 4: The Zinc Filter - Structural Discrimination',
                fontsize=16, fontweight='bold')

    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight', facecolor='white')
    print(f"Saved: {output_path}")
    plt.close()


# ==============================================================================
# FIGURE 5: ZINC TRAP
# ==============================================================================

def create_zinc_trap_plot(df, output_path):
    """
    Focus on THR, SER, ILE to show the zinc trap.
    """

    # Filter to relevant AAs
    df_trap = df[df['primary_ligand'].isin(['THR', 'SER', 'ILE'])].copy()
    df_trap = df_trap.sort_values('AA_iptm', ascending=False)

    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Panel A: AA binding scores
    x = np.arange(len(df_trap))
    colors = ['#2ecc71', '#f39c12', '#3498db']  # THR=green, SER=orange, ILE=blue

    bars1 = ax1.bar(x, df_trap['AA_iptm'], color=colors, edgecolor='black', linewidth=2, width=0.6)

    # Add scores
    for i, (idx, row) in enumerate(df_trap.iterrows()):
        ax1.text(i, row['AA_iptm'] + 0.01, f"{row['AA_iptm']:.3f}",
                ha='center', va='bottom', fontsize=12, fontweight='bold')

    ax1.set_xticks(x)
    ax1.set_xticklabels(df_trap['primary_ligand'], fontsize=12)
    ax1.set_ylabel('Protein-AA ipTM Score', fontsize=11, fontweight='bold')
    ax1.set_title('AA Binding Affinity', fontsize=13, fontweight='bold')
    ax1.set_ylim([0.75, 1.0])
    ax1.grid(axis='y', alpha=0.3)

    # Panel B: Zn coordination comparison
    bars2 = ax2.bar(x, df_trap['Zn_iptm'], color=colors, edgecolor='black', linewidth=2, width=0.6)

    # Add scores
    for i, (idx, row) in enumerate(df_trap.iterrows()):
        ax2.text(i, row['Zn_iptm'] + 0.005, f"{row['Zn_iptm']:.3f}",
                ha='center', va='bottom', fontsize=12, fontweight='bold')

    ax2.set_xticks(x)
    ax2.set_xticklabels(df_trap['primary_ligand'], fontsize=12)
    ax2.set_ylabel('Protein-Zn ipTM Score', fontsize=11, fontweight='bold')
    ax2.set_title('Zn Coordination', fontsize=13, fontweight='bold')
    ax2.set_ylim([0.75, 1.0])
    ax2.grid(axis='y', alpha=0.3)

    # Add key findings
    thr_score = df_trap[df_trap['primary_ligand'] == 'THR']['AA_iptm'].values[0]
    ser_score = df_trap[df_trap['primary_ligand'] == 'SER']['AA_iptm'].values[0]
    ile_score = df_trap[df_trap['primary_ligand'] == 'ILE']['AA_iptm'].values[0]

    thr_ser_ratio = thr_score / ser_score
    thr_ile_ratio = thr_score / ile_score

    textstr = f"THE ZINC TRAP:\n\nTHR/SER: {thr_ser_ratio:.2f}x\n(TOO CLOSE!)\n\nTHR/ILE: {thr_ile_ratio:.2f}x\n(GOOD)\n\nSER coordinates Zn\nvia hydroxyl → TRAPPED!"
    props = dict(boxstyle='round', facecolor='#f39c12', alpha=0.9,
                edgecolor='black', linewidth=2)
    ax1.text(0.98, 0.02, textstr, transform=ax1.transAxes,
            fontsize=10, fontweight='bold', va='bottom', ha='right', bbox=props)

    # Overall title
    fig.suptitle('Figure 5: The "Zinc Trap" - Why Editing is Mandatory',
                fontsize=16, fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight', facecolor='white')
    print(f"Saved: {output_path}")
    plt.close()


# ==============================================================================
# MAIN
# ==============================================================================

print("="*80)
print("GENERATING ZINC EVOLUTION FIGURES")
print("="*80)

# Figure 3: Competition experiments
print("\nGenerating Figure 3: Zinc Disconnect...")
df_comp = pd.read_csv('figures/data/fig3_competitions.csv')
create_competition_plot(df_comp, 'figures/figure3/panel_b_competitions.png')

# Figure 4: Zinc filter heatmap
print("\nGenerating Figure 4: Zinc Filter...")
df_zn_all = pd.read_csv('figures/data/fig4b_mod_thrrs_zn_all.csv')
create_zinc_filter_heatmap(df_zn_all, 'figures/figure4/panel_b_zinc_filter_heatmap.png')

# Figure 5: Zinc trap
print("\nGenerating Figure 5: Zinc Trap...")
df_trap = pd.read_csv('figures/data/fig5b_zinc_trap.csv')
create_zinc_trap_plot(df_trap, 'figures/figure5/panel_ab_zinc_trap.png')

print("\n" + "="*80)
print("ZINC FIGURES COMPLETE")
print("="*80)
print("\nGenerated:")
print("  - Figure 3: figures/figure3/panel_b_competitions.png (+ PDF)")
print("  - Figure 4: figures/figure4/panel_b_zinc_filter_heatmap.png (+ PDF)")
print("  - Figure 5: figures/figure5/panel_ab_zinc_trap.png (+ PDF)")
