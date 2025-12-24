#!/usr/bin/env python3
"""
Physical Mechanism Analysis Figure
Shows the molecular details behind ipTM scores:
- Zn coordination evolution
- Contact and H-bond networks
- Editing domain selectivity mechanism
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load data
df = pd.read_csv('figures/data/comprehensive_ligand_analysis.csv')

# Publication settings
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 9

def create_mechanism_figure():
    """Create comprehensive physical mechanism figure."""

    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.4, wspace=0.3)

    # ===== PANEL A: Zn Coordination Evolution =====
    ax1 = fig.add_subplot(gs[0, :])

    # Get ThrRS structures
    thrrs = df[df['job_name'].str.contains('thrrs', case=False, na=False)].copy()
    thrrs = thrrs[thrrs['ligand'].isin(['THR', 'SER', 'ILE', 'VAL', 'PRO'])]

    # Categorize
    thrrs['condition'] = 'Other'
    thrrs.loc[thrrs['job_name'].str.contains('modern', case=False) &
              thrrs['has_zn'], 'condition'] = 'Modern + Zn'
    thrrs.loc[~thrrs['job_name'].str.contains('modern', case=False) &
              thrrs['has_zn'], 'condition'] = 'Ancestral + Zn'
    thrrs.loc[~thrrs['job_name'].str.contains('modern', case=False) &
              ~thrrs['has_zn'], 'condition'] = 'Ancestral no Zn'

    thrrs = thrrs[thrrs['condition'] != 'Other']

    # Plot
    conditions = ['Ancestral no Zn', 'Ancestral + Zn', 'Modern + Zn']
    x_pos = {c: i for i, c in enumerate(conditions)}

    colors = {'THR': '#27ae60', 'SER': '#f39c12', 'ILE': '#e74c3c',
              'VAL': '#9b59b6', 'PRO': '#3498db'}

    for ligand in ['THR', 'SER', 'ILE', 'VAL', 'PRO']:
        subset = thrrs[thrrs['ligand'] == ligand]
        if len(subset) == 0:
            continue

        x_vals = []
        y_vals = []

        for cond in conditions:
            cond_data = subset[subset['condition'] == cond]
            if len(cond_data) > 0:
                x_vals.append(x_pos[cond])
                if cond == 'Ancestral no Zn':
                    # No Zn distance for this condition
                    y_vals.append(35)  # Off scale
                else:
                    y_vals.append(cond_data.iloc[0]['zn_min_distance'])

        if x_vals:
            ax1.plot(x_vals, y_vals, 'o-', color=colors[ligand],
                    label=ligand, linewidth=2, markersize=8, alpha=0.8)

    ax1.axhline(y=3.0, color='red', linestyle='--', linewidth=2, alpha=0.7)
    ax1.text(2.5, 3.5, 'Coordination cutoff (3Å)', fontsize=9, color='red')

    ax1.set_xticks(range(len(conditions)))
    ax1.set_xticklabels(conditions, fontsize=11)
    ax1.set_ylabel('Distance to Zn²⁺ (Å)', fontsize=12, fontweight='bold')
    ax1.set_title('A. Evolution of Zn Coordination', fontsize=14, fontweight='bold', pad=10)
    ax1.set_ylim(0, 40)
    ax1.legend(loc='upper right', ncol=5, framealpha=0.9)
    ax1.grid(True, alpha=0.3)

    # Add annotations
    ax1.text(0.5, 32, 'No Zn coordination\n(no Zn in structure)', ha='center',
            fontsize=9, style='italic', color='gray')
    ax1.text(1, 28, 'Zn present but DISTANT\n(not coordinating)', ha='center',
            fontsize=9, style='italic', color='gray')
    ax1.text(2, 8, 'Zn integrated into\nactive site', ha='center',
            fontsize=9, style='italic', color='green', fontweight='bold')

    # ===== PANEL B: Contact Atoms vs ipTM (Modern ThrRS) =====
    ax2 = fig.add_subplot(gs[1, 0])

    modern_thrrs = df[df['job_name'].str.contains('modern_thrrs_ecoli_zn', case=False, na=False)].copy()

    # Get best per ligand
    modern_best = modern_thrrs.loc[modern_thrrs.groupby('ligand')['AA_iptm'].idxmax()]

    for _, row in modern_best.iterrows():
        ligand = row['ligand']
        color = colors.get(ligand, '#95a5a6')
        ax2.scatter(row['contact_atoms_4A'], row['AA_iptm'], s=100,
                   color=color, alpha=0.7, edgecolors='black', linewidth=1)
        ax2.text(row['contact_atoms_4A'], row['AA_iptm'], f" {ligand}",
                fontsize=8, va='center')

    ax2.set_xlabel('Contact Atoms (4Å cutoff)', fontsize=10, fontweight='bold')
    ax2.set_ylabel('ipTM Score', fontsize=10, fontweight='bold')
    ax2.set_title('B. Modern ThrRS: Contacts vs Affinity', fontsize=11, fontweight='bold')
    ax2.grid(True, alpha=0.3)

    # Correlation
    corr = modern_best[['contact_atoms_4A', 'AA_iptm']].corr().iloc[0, 1]
    ax2.text(0.05, 0.95, f'R = {corr:.3f}', transform=ax2.transAxes,
            fontsize=9, va='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # ===== PANEL C: H-bonds vs ipTM (Modern ThrRS) =====
    ax3 = fig.add_subplot(gs[1, 1])

    for _, row in modern_best.iterrows():
        ligand = row['ligand']
        color = colors.get(ligand, '#95a5a6')
        ax3.scatter(row['hbond_count'], row['AA_iptm'], s=100,
                   color=color, alpha=0.7, edgecolors='black', linewidth=1)
        ax3.text(row['hbond_count'], row['AA_iptm'], f" {ligand}",
                fontsize=8, va='center')

    ax3.set_xlabel('H-bond Count', fontsize=10, fontweight='bold')
    ax3.set_ylabel('ipTM Score', fontsize=10, fontweight='bold')
    ax3.set_title('C. Modern ThrRS: H-bonds vs Affinity', fontsize=11, fontweight='bold')
    ax3.grid(True, alpha=0.3)

    # Correlation
    corr = modern_best[['hbond_count', 'AA_iptm']].corr().iloc[0, 1]
    ax3.text(0.05, 0.95, f'R = {corr:.3f}', transform=ax3.transAxes,
            fontsize=9, va='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # ===== PANEL D: Zn Distance vs ipTM (Modern ThrRS) =====
    ax4 = fig.add_subplot(gs[1, 2])

    for _, row in modern_best.iterrows():
        ligand = row['ligand']
        color = colors.get(ligand, '#95a5a6')
        ax4.scatter(row['zn_min_distance'], row['AA_iptm'], s=100,
                   color=color, alpha=0.7, edgecolors='black', linewidth=1)
        ax4.text(row['zn_min_distance'], row['AA_iptm'], f" {ligand}",
                fontsize=8, va='center')

    ax4.axvline(x=3.0, color='red', linestyle='--', linewidth=2, alpha=0.7)

    ax4.set_xlabel('Zn Distance (Å)', fontsize=10, fontweight='bold')
    ax4.set_ylabel('ipTM Score', fontsize=10, fontweight='bold')
    ax4.set_title('D. Modern ThrRS: Zn Coordination', fontsize=11, fontweight='bold')
    ax4.grid(True, alpha=0.3)

    # Correlation
    corr = modern_best[['zn_min_distance', 'AA_iptm']].corr().iloc[0, 1]
    ax4.text(0.05, 0.95, f'R = {corr:.3f}', transform=ax4.transAxes,
            fontsize=9, va='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # ===== PANEL E: Editing Domain Mechanism =====
    ax5 = fig.add_subplot(gs[2, :])

    editing = df[df['protein_len'] == 300].copy()
    editing_best = editing.loc[editing.groupby('ligand')['AA_iptm'].idxmax()]
    editing_best = editing_best.sort_values('AA_iptm', ascending=False).head(12)

    x = np.arange(len(editing_best))
    width = 0.35

    bars1 = ax5.bar(x - width/2, editing_best['AA_iptm'], width,
                    label='ipTM Score', color='#3498db', alpha=0.7, edgecolor='black')
    ax5_twin = ax5.twinx()
    bars2 = ax5_twin.bar(x + width/2, editing_best['contact_atoms_4A'], width,
                         label='Contact Atoms', color='#e74c3c', alpha=0.7, edgecolor='black')

    ax5.set_xlabel('Amino Acid', fontsize=12, fontweight='bold')
    ax5.set_ylabel('ipTM Score', fontsize=12, fontweight='bold', color='#3498db')
    ax5_twin.set_ylabel('Contact Atoms (4Å)', fontsize=12, fontweight='bold', color='#e74c3c')
    ax5.set_title('E. Editing Domain: THR Binds Better Than Cognate PRO', fontsize=14, fontweight='bold', pad=10)

    ax5.set_xticks(x)
    ax5.set_xticklabels(editing_best['ligand'], fontsize=10)
    ax5.tick_params(axis='y', labelcolor='#3498db')
    ax5_twin.tick_params(axis='y', labelcolor='#e74c3c')

    # Highlight THR and PRO
    thr_idx = list(editing_best['ligand']).index('THR')
    pro_idx = list(editing_best['ligand']).index('PRO')

    ax5.axvline(x=thr_idx, color='#27ae60', linestyle='--', linewidth=2, alpha=0.5)
    ax5.text(thr_idx, 0.92, 'THR\n(error)', ha='center', fontsize=9,
            color='#27ae60', fontweight='bold')

    ax5.axvline(x=pro_idx, color='#9b59b6', linestyle='--', linewidth=2, alpha=0.5)
    ax5.text(pro_idx, 0.92, 'PRO\n(cognate)', ha='center', fontsize=9,
            color='#9b59b6', fontweight='bold')

    # Add annotations
    ax5.text(0.5, 0.05, 'INVERTED SELECTIVITY: Editing domain binds THR (error) better than PRO (cognate)',
            transform=ax5.transAxes, ha='center', fontsize=11, fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='#fef9e7', edgecolor='#f39c12', linewidth=2))

    ax5.legend(loc='upper right', framealpha=0.9)
    ax5_twin.legend(loc='upper left', framealpha=0.9)
    ax5.grid(True, alpha=0.3, axis='y')

    plt.suptitle('Physical Mechanisms Behind ipTM Scores', fontsize=16, fontweight='bold', y=0.995)

    return fig


if __name__ == '__main__':
    import os
    os.makedirs('figures/mechanisms', exist_ok=True)

    print("="*80)
    print("GENERATING PHYSICAL MECHANISM FIGURE")
    print("="*80)

    fig = create_mechanism_figure()

    fig.savefig('figures/mechanisms/physical_mechanisms.png', dpi=300, bbox_inches='tight')
    fig.savefig('figures/mechanisms/physical_mechanisms.pdf', bbox_inches='tight')

    print("\n✓ Saved: figures/mechanisms/physical_mechanisms.png")
    print("✓ Saved: figures/mechanisms/physical_mechanisms.pdf")

    print("\n" + "="*80)
    print("COMPLETE!")
    print("="*80)
