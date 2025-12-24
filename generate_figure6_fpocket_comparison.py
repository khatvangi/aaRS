#!/usr/bin/env python3
"""
generate_figure6_fpocket_comparison.py
======================================
Figure 6: fpocket Binding Pocket Analysis - LUCA vs Modern

Compares binding pocket properties from fpocket analysis:
- Pocket volume
- Surface area
- Druggability score
- Polarity

Shows how binding pocket evolved from LUCA to Modern.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np
import seaborn as sns

# Set style
plt.style.use('default')
sns.set_palette("husl")

# Color scheme
COLORS = {
    'luca': '#E91E63',  # Pink for LUCA
    'modern': '#2196F3',  # Blue for Modern
    'volume': '#4CAF50',  # Green
    'sasa': '#FF9800',  # Orange
    'druggability': '#9C27B0',  # Purple
}

# Data from fpocket analysis (Pocket 1 - the main binding pocket)
data = {
    'LUCA ProRS': {
        'Volume': 749.6,  # Ų
        'Total_SASA': 148.9,  # ų
        'Polar_SASA': 93.3,  # ų
        'Apolar_SASA': 55.6,  # ų
        'Score': 0.321,
        'Druggability': 0.033,
        'N_spheres': 47,
        'Polarity': 66.67,  # % polar atoms
        'Hydrophobicity': -13.636,
    },
    'Modern E. coli ProRS': {
        'Volume': 1847.0,  # Ų
        'Total_SASA': 380.6,  # ų
        'Polar_SASA': 305.8,  # ų
        'Apolar_SASA': 74.9,  # ų
        'Score': 1.479,
        'Druggability': 0.000,
        'N_spheres': 126,
        'Polarity': 68.48,  # % polar atoms
        'Hydrophobicity': -10.955,
    }
}

def create_figure6():
    """Create fpocket comparison figure"""

    fig = plt.figure(figsize=(16, 10))

    # Create grid layout
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.25,
                          left=0.08, right=0.96, top=0.92, bottom=0.06)

    # =========================================================================
    # Panel A: Pocket Volume Comparison
    # =========================================================================
    ax1 = fig.add_subplot(gs[0, 0])

    enzymes = ['LUCA\nProRS', 'Modern\nE. coli ProRS']
    volumes = [data['LUCA ProRS']['Volume'], data['Modern E. coli ProRS']['Volume']]
    colors_list = [COLORS['luca'], COLORS['modern']]

    bars = ax1.bar(enzymes, volumes, color=colors_list, alpha=0.85,
                   edgecolor='black', linewidth=2, width=0.6)

    # Add value labels
    for bar, vol in zip(bars, volumes):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 50,
                f'{vol:.0f} Ų', ha='center', va='bottom',
                fontsize=13, fontweight='bold')

    # Show fold-change
    fold_change = volumes[1] / volumes[0]
    ax1.text(0.5, max(volumes) * 0.7, f'{fold_change:.1f}× larger',
            ha='center', va='center', fontsize=14, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow',
                     alpha=0.7, edgecolor='black', linewidth=2))

    ax1.set_ylabel('Pocket Volume (Ų)', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Enzyme', fontsize=14, fontweight='bold')
    ax1.set_ylim(0, max(volumes) * 1.15)
    ax1.grid(axis='y', alpha=0.3, linestyle='--')
    ax1.set_title('A. Binding Pocket Volume:\nModern is 2.5× Larger than LUCA',
                 fontsize=14, fontweight='bold', pad=15)

    # =========================================================================
    # Panel B: Surface Area Comparison (Polar vs Apolar)
    # =========================================================================
    ax2 = fig.add_subplot(gs[0, 1])

    labels = ['LUCA ProRS', 'Modern E. coli ProRS']
    polar_sasa = [data['LUCA ProRS']['Polar_SASA'],
                  data['Modern E. coli ProRS']['Polar_SASA']]
    apolar_sasa = [data['LUCA ProRS']['Apolar_SASA'],
                   data['Modern E. coli ProRS']['Apolar_SASA']]

    x = np.arange(len(labels))
    width = 0.35

    bars1 = ax2.bar(x - width/2, polar_sasa, width, label='Polar SASA',
                    color='#2196F3', alpha=0.85, edgecolor='black', linewidth=1.5)
    bars2 = ax2.bar(x + width/2, apolar_sasa, width, label='Apolar SASA',
                    color='#FF9800', alpha=0.85, edgecolor='black', linewidth=1.5)

    # Add value labels
    for bar1, bar2 in zip(bars1, bars2):
        height1 = bar1.get_height()
        height2 = bar2.get_height()
        ax2.text(bar1.get_x() + bar1.get_width()/2., height1 + 5,
                f'{height1:.0f}', ha='center', va='bottom',
                fontsize=11, fontweight='bold')
        ax2.text(bar2.get_x() + bar2.get_width()/2., height2 + 5,
                f'{height2:.0f}', ha='center', va='bottom',
                fontsize=11, fontweight='bold')

    ax2.set_ylabel('Surface Area (ų)', fontsize=14, fontweight='bold')
    ax2.set_xlabel('Enzyme', fontsize=14, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, fontsize=12, fontweight='bold')
    ax2.set_ylim(0, max(polar_sasa + apolar_sasa) * 1.2)
    ax2.legend(loc='upper left', fontsize=12, frameon=True, shadow=True)
    ax2.grid(axis='y', alpha=0.3, linestyle='--')
    ax2.set_title('B. Solvent Accessible Surface Area:\nBoth Primarily Polar (Hydrophilic)',
                 fontsize=14, fontweight='bold', pad=15)

    # =========================================================================
    # Panel C: Multiple Pocket Properties Radar Chart
    # =========================================================================
    ax3 = fig.add_subplot(gs[1, 0], projection='polar')

    # Properties to compare (normalized to 0-1 scale)
    categories = ['Volume\n(normalized)', 'Surface Area\n(normalized)',
                  'Polarity\n(%)', 'Score\n(normalized)', 'Alpha Spheres\n(normalized)']
    N = len(categories)

    # Normalize data to 0-1 scale
    luca_values = [
        data['LUCA ProRS']['Volume'] / 2000,  # Normalize by max expected
        data['LUCA ProRS']['Total_SASA'] / 400,
        data['LUCA ProRS']['Polarity'] / 100,
        data['LUCA ProRS']['Score'] / 1.5,
        data['LUCA ProRS']['N_spheres'] / 130,
    ]

    modern_values = [
        data['Modern E. coli ProRS']['Volume'] / 2000,
        data['Modern E. coli ProRS']['Total_SASA'] / 400,
        data['Modern E. coli ProRS']['Polarity'] / 100,
        data['Modern E. coli ProRS']['Score'] / 1.5,
        data['Modern E. coli ProRS']['N_spheres'] / 130,
    ]

    # Angles for each property
    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    luca_values += luca_values[:1]
    modern_values += modern_values[:1]
    angles += angles[:1]

    # Plot
    ax3.plot(angles, luca_values, 'o-', linewidth=2, label='LUCA ProRS',
            color=COLORS['luca'], markersize=8)
    ax3.fill(angles, luca_values, alpha=0.25, color=COLORS['luca'])

    ax3.plot(angles, modern_values, 'o-', linewidth=2, label='Modern E. coli ProRS',
            color=COLORS['modern'], markersize=8)
    ax3.fill(angles, modern_values, alpha=0.25, color=COLORS['modern'])

    # Labels
    ax3.set_xticks(angles[:-1])
    ax3.set_xticklabels(categories, fontsize=11, fontweight='bold')
    ax3.set_ylim(0, 1)
    ax3.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax3.set_yticklabels(['0.25', '0.5', '0.75', '1.0'], fontsize=9)
    ax3.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize=11)
    ax3.grid(True, alpha=0.3)
    ax3.set_title('C. Multi-Property Comparison\n(Normalized Radar Chart)',
                 fontsize=14, fontweight='bold', pad=20)

    # =========================================================================
    # Panel D: Key Findings Summary Table
    # =========================================================================
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis('off')

    summary_data = [
        ['Property', 'LUCA ProRS', 'Modern E. coli', 'Fold Change'],
        ['─────────', '──────────', '──────────────', '────────────'],
        ['Volume (Ų)', f"{data['LUCA ProRS']['Volume']:.0f}",
         f"{data['Modern E. coli ProRS']['Volume']:.0f}",
         f"{data['Modern E. coli ProRS']['Volume']/data['LUCA ProRS']['Volume']:.1f}×"],
        ['Total SASA (ų)', f"{data['LUCA ProRS']['Total_SASA']:.0f}",
         f"{data['Modern E. coli ProRS']['Total_SASA']:.0f}",
         f"{data['Modern E. coli ProRS']['Total_SASA']/data['LUCA ProRS']['Total_SASA']:.1f}×"],
        ['Polar SASA (ų)', f"{data['LUCA ProRS']['Polar_SASA']:.0f}",
         f"{data['Modern E. coli ProRS']['Polar_SASA']:.0f}",
         f"{data['Modern E. coli ProRS']['Polar_SASA']/data['LUCA ProRS']['Polar_SASA']:.1f}×"],
        ['Polarity (%)', f"{data['LUCA ProRS']['Polarity']:.1f}",
         f"{data['Modern E. coli ProRS']['Polarity']:.1f}", '~Equal'],
        ['fpocket Score', f"{data['LUCA ProRS']['Score']:.2f}",
         f"{data['Modern E. coli ProRS']['Score']:.2f}",
         f"{data['Modern E. coli ProRS']['Score']/data['LUCA ProRS']['Score']:.1f}×"],
        ['Alpha Spheres', f"{data['LUCA ProRS']['N_spheres']}",
         f"{data['Modern E. coli ProRS']['N_spheres']}",
         f"{data['Modern E. coli ProRS']['N_spheres']/data['LUCA ProRS']['N_spheres']:.1f}×"],
    ]

    # Create table
    table = ax4.table(cellText=summary_data, cellLoc='left',
                     loc='center', bbox=[0, 0, 1, 1])

    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 2.5)

    # Style header row
    for i in range(4):
        cell = table[(0, i)]
        cell.set_facecolor('#E0E0E0')
        cell.set_text_props(weight='bold', fontsize=12)

    # Color code fold changes
    for i in range(2, 8):
        cell = table[(i, 3)]
        if '×' in cell.get_text().get_text():
            try:
                fold = float(cell.get_text().get_text().replace('×', ''))
                if fold > 2.0:
                    cell.set_facecolor('#FFEBEE')  # Light red - big increase
                elif fold > 1.5:
                    cell.set_facecolor('#FFF9C4')  # Light yellow - moderate increase
            except:
                pass

    ax4.set_title('D. Summary Table: Modern Pocket is Larger & Higher Quality',
                 fontsize=14, fontweight='bold', pad=15)

    # Add interpretation box
    interpretation = """
KEY FINDINGS:

1. VOLUME EXPANSION: Modern pocket is 2.5× larger (1847 vs 750 Ų)
   → More spacious substrate binding site
   → Potentially accommodates diverse conformations

2. SURFACE AREA: Modern has 2.6× more surface area (381 vs 149 ų)
   → Increased protein-ligand interactions
   → More contact points for substrate recognition

3. POLARITY MAINTAINED: Both ~67-68% polar
   → Hydrophilic character conserved across evolution
   → Consistent with binding charged amino acids (Pro, Thr)

4. QUALITY SCORE: Modern scores 4.6× higher (1.48 vs 0.32)
   → Better-defined binding pocket geometry
   → More optimal for substrate binding

5. COMPLEXITY: Modern has 2.7× more alpha spheres (126 vs 47)
   → More sophisticated pocket architecture
   → Refined geometric organization

BIOLOGICAL INTERPRETATION:
Evolution EXPANDED the binding pocket while maintaining its chemical
properties (polarity). Larger volume may allow promiscuous binding of
both Pro and Thr, consistent with ipTM data (Figure 3).
"""

    ax4.text(0.5, -0.35, interpretation, transform=ax4.transAxes,
            fontsize=9, verticalalignment='top', horizontalalignment='center',
            family='monospace',
            bbox=dict(boxstyle='round,pad=1', facecolor='lightyellow',
                     edgecolor='orange', linewidth=2, alpha=0.9))

    # Overall title
    fig.suptitle('Figure 6: fpocket Binding Pocket Analysis - Evolution from LUCA to Modern',
                fontsize=17, fontweight='bold', y=0.98)

    return fig

# =============================================================================
# GENERATE AND SAVE
# =============================================================================

if __name__ == '__main__':
    fig = create_figure6()

    # Save in multiple formats and locations
    output_dirs = [
        '/storage/kiran-stuff/aaRS/manuscript_figures',
        '/storage/kiran-stuff/aaRS/final_figures',
        '/storage/kiran-stuff/aaRS/figures'
    ]

    for output_dir in output_dirs:
        fig.savefig(f'{output_dir}/Figure6_fpocket_Pocket_Comparison.png',
                    dpi=300, bbox_inches='tight', facecolor='white')
        fig.savefig(f'{output_dir}/Figure6_fpocket_Pocket_Comparison.pdf',
                    bbox_inches='tight', facecolor='white')
        fig.savefig(f'{output_dir}/Figure6_fpocket_Pocket_Comparison.svg',
                    bbox_inches='tight', facecolor='white')
        print(f"✓ Saved to {output_dir}/")

    plt.close()
    print("\n✅ Figure 6 generation complete!")
    print("\nKEY FINDING: Modern pocket is 2.5× LARGER than LUCA!")
    print("This expanded volume may explain maintained promiscuity (Figure 3)")
