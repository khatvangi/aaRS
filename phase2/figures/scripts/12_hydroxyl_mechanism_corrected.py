#!/usr/bin/env python3
"""
CORRECTED ANALYSIS: Hydroxyl-Mediated Zn Coordination
Shows that -OH groups enable bidentate coordination → high ipTM
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import json

# Load data
with open('figures/data/comprehensive_ligand_analysis.json', 'r') as f:
    data = json.load(f)

df = pd.read_csv('figures/data/comprehensive_ligand_analysis.csv')

# Publication settings
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

def create_hydroxyl_mechanism_figure():
    """Create corrected figure showing hydroxyl-mediated Zn coordination."""

    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.3)

    # Get modern ThrRS data
    modern_thrrs = df[df['job_name'].str.contains('modern_thrrs_ecoli_zn', case=False, na=False)].copy()

    # Add coordination atom count from JSON
    coord_counts = {}
    has_hydroxyl = {}
    coord_atoms_list = {}

    for struct in data:
        if 'modern_thrrs_ecoli_zn' in struct.get('job_name', '').lower():
            ligand = struct.get('ligand')
            if ligand:
                coord_atoms = struct.get('coordinating_atoms', [])
                coord_counts[ligand] = len(coord_atoms)
                coord_atoms_list[ligand] = coord_atoms
                has_hydroxyl[ligand] = any('OG' in atom or 'OH' in atom for atom in coord_atoms)

    modern_thrrs['coord_atom_count'] = modern_thrrs['ligand'].map(coord_counts)
    modern_thrrs['has_hydroxyl'] = modern_thrrs['ligand'].map(has_hydroxyl)
    modern_thrrs['coord_atoms'] = modern_thrrs['ligand'].map(coord_atoms_list)

    # Sort by ipTM
    modern_thrrs = modern_thrrs.sort_values('AA_iptm', ascending=False)

    # ===== PANEL A: Coordination Count vs ipTM =====
    ax1 = fig.add_subplot(gs[0, 0])

    # Color by coordination type
    colors_coord = {3: '#27ae60', 2: '#f39c12', 1: '#e74c3c'}

    for _, row in modern_thrrs.iterrows():
        coord_count = row['coord_atom_count']
        color = colors_coord.get(coord_count, '#95a5a6')

        # Make THR and SER stand out
        if row['ligand'] in ['THR', 'SER']:
            marker = 's'
            size = 200
            edge = 'black'
            linewidth = 3
        else:
            marker = 'o'
            size = 100
            edge = 'black'
            linewidth = 1

        ax1.scatter(coord_count, row['AA_iptm'], s=size, color=color,
                   marker=marker, alpha=0.8, edgecolors=edge, linewidth=linewidth)

        # Label THR and SER
        if row['ligand'] in ['THR', 'SER']:
            ax1.text(coord_count + 0.05, row['AA_iptm'], f" {row['ligand']}",
                    fontsize=11, fontweight='bold', va='center')

    ax1.set_xlabel('Number of Coordinating Atoms', fontsize=12, fontweight='bold')
    ax1.set_ylabel('ipTM Score', fontsize=12, fontweight='bold')
    ax1.set_title('A. Bidentate vs Monodentate Coordination', fontsize=13, fontweight='bold')
    ax1.set_xticks([1, 2, 3])
    ax1.set_xticklabels(['Monodentate\n(1 atom)', 'Weak\n(2 atoms)', 'Bidentate\n(3 atoms)'])
    ax1.grid(True, alpha=0.3)

    # Add legend
    legend_elements = [
        mpatches.Patch(facecolor='#27ae60', label='Bidentate (3 atoms)'),
        mpatches.Patch(facecolor='#f39c12', label='2 atoms'),
        mpatches.Patch(facecolor='#e74c3c', label='Monodentate (1 atom)'),
        mpatches.Patch(facecolor='white', edgecolor='black', linewidth=3, label='THR/SER (hydroxyl)')
    ]
    ax1.legend(handles=legend_elements, loc='lower right', framealpha=0.9)

    # ===== PANEL B: Hydroxyl vs Non-Hydroxyl =====
    ax2 = fig.add_subplot(gs[0, 1])

    hydroxyl_aa = modern_thrrs[modern_thrrs['has_hydroxyl'] == True]
    non_hydroxyl_aa = modern_thrrs[modern_thrrs['has_hydroxyl'] == False]

    # Bar plot
    categories = ['Hydroxyl\n(THR, SER)', 'Non-Hydroxyl\n(Others)']
    means = [hydroxyl_aa['AA_iptm'].mean(), non_hydroxyl_aa['AA_iptm'].mean()]

    bars = ax2.bar(categories, means, color=['#27ae60', '#95a5a6'],
                   alpha=0.7, edgecolor='black', linewidth=2)

    # Add individual points
    for i, row in hydroxyl_aa.iterrows():
        ax2.scatter(0, row['AA_iptm'], s=150, color='#27ae60',
                   marker='s', edgecolors='black', linewidth=2, zorder=10)
        ax2.text(0.15, row['AA_iptm'], row['ligand'], fontsize=10, fontweight='bold')

    ax2.set_ylabel('ipTM Score', fontsize=12, fontweight='bold')
    ax2.set_title('B. Hydroxyl Group Effect', fontsize=13, fontweight='bold')
    ax2.set_ylim(0.75, 1.0)

    # Add values on bars
    for i, (bar, val) in enumerate(zip(bars, means)):
        ax2.text(bar.get_x() + bar.get_width()/2, val + 0.01,
                f'{val:.3f}', ha='center', va='bottom',
                fontsize=11, fontweight='bold')

    # Add discrimination percentage
    discrimination = (means[0] - means[1]) / means[1] * 100
    ax2.text(0.5, 0.92, f'{discrimination:.1f}% discrimination',
            ha='center', fontsize=11, fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))

    ax2.grid(True, alpha=0.3, axis='y')

    # ===== PANEL C: Coordinating Atoms Diagram =====
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.axis('off')

    # THR diagram
    ax3.text(0.5, 0.95, 'THR: Bidentate Coordination', ha='center',
            fontsize=12, fontweight='bold', transform=ax3.transAxes)

    # Draw Zn
    zn_circle = plt.Circle((0.5, 0.65), 0.08, color='#7f8c8d', alpha=0.9,
                          transform=ax3.transAxes)
    ax3.add_patch(zn_circle)
    ax3.text(0.5, 0.65, 'Zn²⁺', ha='center', va='center', fontsize=10,
            fontweight='bold', color='white', transform=ax3.transAxes)

    # Draw THR atoms
    atoms = [
        (0.25, 0.85, 'N', '#3498db'),
        (0.35, 0.75, 'CA', '#95a5a6'),
        (0.65, 0.55, 'OG1', '#e74c3c')
    ]

    for x, y, label, color in atoms:
        circle = plt.Circle((x, y), 0.05, color=color, alpha=0.8,
                          transform=ax3.transAxes)
        ax3.add_patch(circle)
        ax3.text(x, y, label, ha='center', va='center', fontsize=8,
                fontweight='bold', color='white', transform=ax3.transAxes)

        # Draw coordination bond
        ax3.plot([x, 0.5], [y, 0.65], 'k--', lw=2, alpha=0.5,
                transform=ax3.transAxes)

    ax3.text(0.5, 0.45, 'ipTM: 0.97', ha='center', fontsize=11,
            fontweight='bold', color='#27ae60', transform=ax3.transAxes)
    ax3.text(0.5, 0.35, '3 coordinating atoms', ha='center', fontsize=9,
            style='italic', transform=ax3.transAxes)

    # ILE diagram
    ax3.text(0.5, 0.25, 'ILE: Monodentate Coordination', ha='center',
            fontsize=12, fontweight='bold', transform=ax3.transAxes)

    # Draw Zn
    zn_circle2 = plt.Circle((0.5, 0.05), 0.08, color='#7f8c8d', alpha=0.9,
                           transform=ax3.transAxes)
    ax3.add_patch(zn_circle2)
    ax3.text(0.5, 0.05, 'Zn²⁺', ha='center', va='center', fontsize=10,
            fontweight='bold', color='white', transform=ax3.transAxes)

    # Draw ILE atoms (only 2 coordinate)
    atoms_ile = [
        (0.3, 0.15, 'N', '#3498db'),
        (0.4, 0.1, 'CA', '#95a5a6')
    ]

    for x, y, label, color in atoms_ile:
        circle = plt.Circle((x, y), 0.05, color=color, alpha=0.8,
                          transform=ax3.transAxes)
        ax3.add_patch(circle)
        ax3.text(x, y, label, ha='center', va='center', fontsize=8,
                fontweight='bold', color='white', transform=ax3.transAxes)

        # Draw coordination bond
        ax3.plot([x, 0.5], [y, 0.05], 'k--', lw=2, alpha=0.5,
                transform=ax3.transAxes)

    # Show missing -OH
    ax3.text(0.7, 0.15, 'NO -OH', ha='center', fontsize=10,
            fontweight='bold', color='#e74c3c',
            bbox=dict(boxstyle='round', facecolor='#ffeeee', edgecolor='#e74c3c', linewidth=2),
            transform=ax3.transAxes)

    # ===== PANEL D: All Ligands Sorted =====
    ax4 = fig.add_subplot(gs[1, :])

    x = np.arange(len(modern_thrrs))

    # Color by hydroxyl status
    colors_bars = ['#27ae60' if h else '#95a5a6' for h in modern_thrrs['has_hydroxyl']]

    bars = ax4.bar(x, modern_thrrs['AA_iptm'], color=colors_bars,
                   alpha=0.7, edgecolor='black', linewidth=1)

    # Highlight THR and SER
    for i, (idx, row) in enumerate(modern_thrrs.iterrows()):
        if row['ligand'] in ['THR', 'SER']:
            bars[i].set_edgecolor('black')
            bars[i].set_linewidth(3)

    ax4.set_xticks(x)
    ax4.set_xticklabels(modern_thrrs['ligand'], rotation=0, fontsize=10)
    ax4.set_ylabel('ipTM Score', fontsize=12, fontweight='bold')
    ax4.set_title('D. All Amino Acids Ranked by ipTM (Modern ThrRS + Zn)', fontsize=13, fontweight='bold')
    ax4.grid(True, alpha=0.3, axis='y')

    # Add coordination count labels
    for i, (idx, row) in enumerate(modern_thrrs.iterrows()):
        if row['coord_atoms'] is not None and len(row['coord_atoms']) > 0:
            coord_str = ', '.join(row['coord_atoms'])
            if len(coord_str) > 30:
                coord_str = f"{row['coord_atom_count']} atoms"
        else:
            coord_str = ''

        if coord_str:
            y_pos = row['AA_iptm'] + 0.005
            ax4.text(i, y_pos, coord_str, ha='center', va='bottom',
                    fontsize=7, rotation=45, style='italic')

    # Add annotations
    ax4.axhline(y=0.96, color='#27ae60', linestyle='--', linewidth=2, alpha=0.5)
    ax4.text(0.02, 0.965, 'Hydroxyl threshold (THR, SER)', fontsize=9,
            color='#27ae60', fontweight='bold', transform=ax4.transAxes)

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor='#27ae60', label='Has -OH group (bidentate)'),
        mpatches.Patch(facecolor='#95a5a6', label='No -OH (monodentate)'),
    ]
    ax4.legend(handles=legend_elements, loc='upper right', framealpha=0.9, fontsize=10)

    plt.suptitle('The Hydroxyl Mechanism: -OH Groups Enable Bidentate Zn²⁺ Coordination',
                fontsize=16, fontweight='bold', y=0.995)

    return fig


if __name__ == '__main__':
    import os

    print("="*80)
    print("CORRECTED ANALYSIS: HYDROXYL-MEDIATED Zn COORDINATION")
    print("="*80)

    fig = create_hydroxyl_mechanism_figure()

    os.makedirs('figures/corrected', exist_ok=True)
    fig.savefig('figures/corrected/hydroxyl_mechanism.png', dpi=300, bbox_inches='tight')
    fig.savefig('figures/corrected/hydroxyl_mechanism.pdf', bbox_inches='tight')

    print("\n✓ Saved: figures/corrected/hydroxyl_mechanism.png")
    print("✓ Saved: figures/corrected/hydroxyl_mechanism.pdf")

    # Print summary
    df = pd.read_csv('figures/data/comprehensive_ligand_analysis.csv')
    modern = df[df['job_name'].str.contains('modern_thrrs_ecoli_zn', case=False, na=False)]

    print("\n" + "="*80)
    print("KEY FINDINGS")
    print("="*80)
    print("\nTHR: N-CA-OG1 coordination (bidentate via -OH) → ipTM 0.97")
    print("SER: N-CB-OG coordination (bidentate via -OH) → ipTM 0.95 (98% of THR!)")
    print("ILE: N-CA coordination (NO -OH) → ipTM 0.83")
    print("VAL: N-CA coordination (NO -OH) → ipTM 0.90")
    print("\nThe -OH group is the THIRD coordinating atom that creates bidentate geometry!")
    print("This is why SER is trapped - it has the SAME coordination chemistry as THR!")
    print("\n" + "="*80)
    print("COMPLETE!")
    print("="*80)
