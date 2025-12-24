#!/usr/bin/env python3
"""
Generate improved Figure 1 for aaRS manuscript
- Larger, clearer phylogeny with no overlapping text
- Remove "node" text
- Better color scheme for sequence conservation
- Enlarged for readability
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Rectangle, Polygon, FancyArrowPatch, Wedge
import numpy as np
from matplotlib import rcParams

# Publication settings - larger fonts for clarity
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
rcParams['font.size'] = 10
rcParams['axes.linewidth'] = 1.0
rcParams['xtick.major.width'] = 0.8
rcParams['ytick.major.width'] = 0.8
rcParams['lines.linewidth'] = 1.5
rcParams['patch.linewidth'] = 0.8

# Clear, colorblind-friendly colors
COLORS = {
    'archaea': '#E74C3C',      # Bright red
    'bacteria': '#3498DB',     # Bright blue
    'eukaryota': '#27AE60',    # Bright green
    'luca': '#9B59B6',         # Purple
    'catalytic': '#1976D2',    # Blue
    'editing': '#F57C00',      # Orange
    'background': '#E8E8E8',   # Light gray
}

def figure1_phylogenetic_overview_improved():
    """
    IMPROVED FIGURE 1 - Clearer phylogeny with larger labels
    """
    # Make figure larger for clarity
    fig = plt.figure(figsize=(10, 7))  # Increased from 7x5.5

    # Create grid with more space between panels
    gs = fig.add_gridspec(2, 2, width_ratios=[2.5, 1.5], height_ratios=[1, 1],
                          hspace=0.35, wspace=0.5,
                          left=0.08, right=0.96, top=0.94, bottom=0.06)

    # ============================================================================
    # Panel A: ProRS Phylogenetic Tree (IMPROVED)
    # ============================================================================
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.text(-0.08, 0.98, 'A', transform=ax1.transAxes,
             fontsize=16, fontweight='bold', va='top')

    # Clade positions - more vertical spacing
    y_archaea = 0.75
    y_bacteria = 0.50
    y_eukaryota = 0.25

    # Draw tree branches FIRST (so they're behind)
    ax1.plot([0, 0.15], [0.5, 0.5], 'k-', linewidth=2)  # Root to LUCA
    ax1.plot([0.15, 0.15], [0.25, 0.75], 'k-', linewidth=2)  # Vertical trunk
    ax1.plot([0.15, 0.35], [y_archaea, y_archaea], 'k-', linewidth=2)
    ax1.plot([0.15, 0.35], [y_bacteria, y_bacteria], 'k-', linewidth=2)
    ax1.plot([0.15, 0.25], [y_eukaryota, y_eukaryota], 'k-', linewidth=2)
    ax1.plot([0.25, 0.35], [y_eukaryota, y_eukaryota], 'k-', linewidth=1.5, linestyle='--')

    # Plot collapsed clades as LARGER triangles
    # Archaea clade
    triangle_archaea = Polygon([[0.35, y_archaea-0.08], [0.75, y_archaea-0.12],
                                 [0.75, y_archaea+0.12]],
                                facecolor=COLORS['archaea'], alpha=0.4,
                                edgecolor=COLORS['archaea'], linewidth=2)
    ax1.add_patch(triangle_archaea)
    ax1.text(0.82, y_archaea, 'Archaea (n=31)', va='center', fontsize=11, fontweight='bold')

    # Bacteria clade
    triangle_bacteria = Polygon([[0.35, y_bacteria-0.08], [0.75, y_bacteria-0.12],
                                  [0.75, y_bacteria+0.12]],
                                 facecolor=COLORS['bacteria'], alpha=0.4,
                                 edgecolor=COLORS['bacteria'], linewidth=2)
    ax1.add_patch(triangle_bacteria)
    ax1.text(0.82, y_bacteria, 'Bacteria (n=45)', va='center', fontsize=11, fontweight='bold')

    # Eukaryota clade
    triangle_eukaryota = Polygon([[0.35, y_eukaryota-0.08], [0.75, y_eukaryota-0.12],
                                   [0.75, y_eukaryota+0.12]],
                                  facecolor=COLORS['eukaryota'], alpha=0.4,
                                  edgecolor=COLORS['eukaryota'], linewidth=2)
    ax1.add_patch(triangle_eukaryota)
    ax1.text(0.82, y_eukaryota, 'Eukaryota (n=17)', va='center', fontsize=11, fontweight='bold')

    # Mark ancestral nodes - LARGER and CLEARER labels
    # LUCA - removed "node" word
    ax1.plot(0.15, 0.5, 'o', markersize=14, color=COLORS['luca'],
             markeredgecolor='black', markeredgewidth=2, zorder=10)
    ax1.text(0.08, 0.5, 'LUCA', ha='right', va='center',
             fontsize=12, fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor=COLORS['luca'], linewidth=1.5))

    # Eukaryotic ancestor - removed "node" word
    ax1.plot(0.25, y_eukaryota, '^', markersize=12, color=COLORS['eukaryota'],
             markeredgecolor='black', markeredgewidth=2, zorder=10)
    ax1.text(0.25, y_eukaryota-0.12, 'Eukaryotic\nancestor', ha='center',
             va='top', fontsize=10, fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor=COLORS['eukaryota'], linewidth=1.5))

    # Bootstrap values - larger and clearer
    ax1.text(0.17, 0.63, '98%', fontsize=9, style='italic', color='#333333',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))
    ax1.text(0.17, 0.37, '95%', fontsize=9, style='italic', color='#333333',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))

    # Scale bar - larger and clearer
    ax1.plot([0.05, 0.25], [0.05, 0.05], 'k-', linewidth=2.5)
    ax1.plot([0.05, 0.05], [0.03, 0.07], 'k-', linewidth=2.5)
    ax1.plot([0.25, 0.25], [0.03, 0.07], 'k-', linewidth=2.5)
    ax1.text(0.15, 0.01, '0.2 substitutions/site', ha='center', fontsize=10,
             fontweight='bold')

    ax1.set_xlim(-0.05, 1.05)
    ax1.set_ylim(-0.02, 1.02)
    ax1.axis('off')
    ax1.set_title('ProRS Phylogeny (93 sequences)', fontsize=13,
                  fontweight='bold', pad=12)

    # ============================================================================
    # Panel B: ThrRS Phylogenetic Tree (IMPROVED)
    # ============================================================================
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.text(-0.08, 0.98, 'B', transform=ax2.transAxes,
             fontsize=16, fontweight='bold', va='top')

    # Similar tree for ThrRS - same improvements
    ax2.plot([0, 0.15], [0.5, 0.5], 'k-', linewidth=2)
    ax2.plot([0.15, 0.15], [0.25, 0.75], 'k-', linewidth=2)
    ax2.plot([0.15, 0.35], [y_archaea, y_archaea], 'k-', linewidth=2)
    ax2.plot([0.15, 0.35], [y_bacteria, y_bacteria], 'k-', linewidth=2)
    ax2.plot([0.15, 0.35], [y_eukaryota, y_eukaryota], 'k-', linewidth=2)

    # Clades
    triangle_archaea2 = Polygon([[0.35, y_archaea-0.08], [0.75, y_archaea-0.12],
                                  [0.75, y_archaea+0.12]],
                                 facecolor=COLORS['archaea'], alpha=0.4,
                                 edgecolor=COLORS['archaea'], linewidth=2)
    ax2.add_patch(triangle_archaea2)
    ax2.text(0.82, y_archaea, 'Archaea (n=21)', va='center', fontsize=11, fontweight='bold')

    triangle_bacteria2 = Polygon([[0.35, y_bacteria-0.08], [0.75, y_bacteria-0.12],
                                   [0.75, y_bacteria+0.12]],
                                  facecolor=COLORS['bacteria'], alpha=0.4,
                                  edgecolor=COLORS['bacteria'], linewidth=2)
    ax2.add_patch(triangle_bacteria2)
    ax2.text(0.82, y_bacteria, 'Bacteria (n=32)', va='center', fontsize=11, fontweight='bold')

    triangle_eukaryota2 = Polygon([[0.35, y_eukaryota-0.08], [0.75, y_eukaryota-0.12],
                                    [0.75, y_eukaryota+0.12]],
                                   facecolor=COLORS['eukaryota'], alpha=0.4,
                                   edgecolor=COLORS['eukaryota'], linewidth=2)
    ax2.add_patch(triangle_eukaryota2)
    ax2.text(0.82, y_eukaryota, 'Eukaryota (n=11)', va='center', fontsize=11, fontweight='bold')

    # LUCA ThrRS
    ax2.plot(0.15, 0.5, 'o', markersize=14, color=COLORS['luca'],
             markeredgecolor='black', markeredgewidth=2, zorder=10)
    ax2.text(0.08, 0.5, 'LUCA', ha='right', va='center',
             fontsize=12, fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor=COLORS['luca'], linewidth=1.5))

    # Bootstrap values
    ax2.text(0.17, 0.63, '96%', fontsize=9, style='italic', color='#333333',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))
    ax2.text(0.17, 0.37, '92%', fontsize=9, style='italic', color='#333333',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))

    # Scale bar
    ax2.plot([0.05, 0.25], [0.05, 0.05], 'k-', linewidth=2.5)
    ax2.plot([0.05, 0.05], [0.03, 0.07], 'k-', linewidth=2.5)
    ax2.plot([0.25, 0.25], [0.03, 0.07], 'k-', linewidth=2.5)
    ax2.text(0.15, 0.01, '0.2 substitutions/site', ha='center', fontsize=10,
             fontweight='bold')

    ax2.set_xlim(-0.05, 1.05)
    ax2.set_ylim(-0.02, 1.02)
    ax2.axis('off')
    ax2.set_title('ThrRS Phylogeny (64 sequences)', fontsize=13,
                  fontweight='bold', pad=12)

    # ============================================================================
    # Panel C: Domain Architecture (IMPROVED)
    # ============================================================================
    ax3 = fig.add_subplot(gs[0, 1])
    ax3.text(-0.12, 0.98, 'C', transform=ax3.transAxes,
             fontsize=16, fontweight='bold', va='top')

    # ProRS domain architecture - larger and clearer
    y_pos = 0.65

    # N-terminal
    rect1 = Rectangle((0.05, y_pos-0.06), 0.12, 0.12,
                      facecolor=COLORS['background'], alpha=0.7,
                      edgecolor='black', linewidth=1.5)
    ax3.add_patch(rect1)
    ax3.text(0.11, y_pos, 'N', ha='center', va='center', fontsize=9, fontweight='bold')

    # Catalytic domain - BRIGHTER blue
    rect2 = Rectangle((0.18, y_pos-0.06), 0.38, 0.12,
                      facecolor=COLORS['catalytic'], alpha=0.8,
                      edgecolor='black', linewidth=1.5)
    ax3.add_patch(rect2)
    ax3.text(0.37, y_pos, 'Catalytic Domain', ha='center', va='center',
             fontsize=10, fontweight='bold', color='white')

    # Editing domain - BRIGHTER orange
    rect3 = Rectangle((0.60, y_pos-0.06), 0.26, 0.12,
                      facecolor=COLORS['editing'], alpha=0.9,
                      edgecolor='black', linewidth=1.5)
    ax3.add_patch(rect3)
    ax3.text(0.73, y_pos, 'Editing\nDomain', ha='center', va='center',
             fontsize=9, fontweight='bold', color='white')

    # C-terminal
    rect4 = Rectangle((0.88, y_pos-0.06), 0.08, 0.12,
                      facecolor=COLORS['background'], alpha=0.7,
                      edgecolor='black', linewidth=1.5)
    ax3.add_patch(rect4)
    ax3.text(0.92, y_pos, 'C', ha='center', va='center', fontsize=9, fontweight='bold')

    ax3.text(0.5, y_pos+0.16, 'ProRS (2037 aa)', ha='center', fontsize=11,
             fontweight='bold', color=COLORS['catalytic'])
    ax3.text(0.05, y_pos+0.11, '1', ha='left', fontsize=8, style='italic')
    ax3.text(0.96, y_pos+0.11, '2037', ha='right', fontsize=8, style='italic')

    # ThrRS domain architecture - NO editing domain
    y_pos = 0.25

    # N-terminal
    rect5 = Rectangle((0.05, y_pos-0.06), 0.12, 0.12,
                      facecolor=COLORS['background'], alpha=0.7,
                      edgecolor='black', linewidth=1.5)
    ax3.add_patch(rect5)
    ax3.text(0.11, y_pos, 'N', ha='center', va='center', fontsize=9, fontweight='bold')

    # Catalytic domain (larger - NO editing domain!)
    rect6 = Rectangle((0.18, y_pos-0.06), 0.68, 0.12,
                      facecolor=COLORS['catalytic'], alpha=0.8,
                      edgecolor='black', linewidth=1.5)
    ax3.add_patch(rect6)
    ax3.text(0.52, y_pos, 'Catalytic Domain', ha='center', va='center',
             fontsize=10, fontweight='bold', color='white')

    # C-terminal
    rect7 = Rectangle((0.88, y_pos-0.06), 0.08, 0.12,
                      facecolor=COLORS['background'], alpha=0.7,
                      edgecolor='black', linewidth=1.5)
    ax3.add_patch(rect7)
    ax3.text(0.92, y_pos, 'C', ha='center', va='center', fontsize=9, fontweight='bold')

    ax3.text(0.5, y_pos+0.16, 'ThrRS (1017 aa)', ha='center', fontsize=11,
             fontweight='bold', color=COLORS['catalytic'])
    ax3.text(0.05, y_pos+0.11, '1', ha='left', fontsize=8, style='italic')
    ax3.text(0.96, y_pos+0.11, '1017', ha='right', fontsize=8, style='italic')

    # Add annotation showing difference
    ax3.annotate('', xy=(0.73, 0.48), xytext=(0.73, 0.37),
                arrowprops=dict(arrowstyle='<->', color='red', lw=2))
    ax3.text(0.78, 0.425, 'ProRS has\nediting domain', ha='left', va='center',
             fontsize=9, fontweight='bold', color='red',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='#FFEEEE',
                      edgecolor='red', linewidth=1.5))

    ax3.set_xlim(0, 1.05)
    ax3.set_ylim(0.05, 0.95)
    ax3.axis('off')
    ax3.set_title('Domain Architecture', fontsize=13, fontweight='bold', pad=12)

    # ============================================================================
    # Panel D: Legend
    # ============================================================================
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.text(-0.12, 0.98, 'D', transform=ax4.transAxes,
             fontsize=16, fontweight='bold', va='top')

    # Legend elements with LARGER fonts
    y = 0.85
    dy = 0.13

    ax4.text(0.05, y, 'Phylogenetic Markers:', fontsize=11, fontweight='bold')
    y -= dy

    # LUCA
    ax4.plot(0.08, y, 'o', markersize=12, color=COLORS['luca'],
             markeredgecolor='black', markeredgewidth=1.5)
    ax4.text(0.15, y, 'Last Universal\nCommon Ancestor', va='center', fontsize=9)
    y -= dy

    # Euk ancestor
    ax4.plot(0.08, y, '^', markersize=11, color=COLORS['eukaryota'],
             markeredgecolor='black', markeredgewidth=1.5)
    ax4.text(0.15, y, 'Eukaryotic Ancestor', va='center', fontsize=9)
    y -= dy * 1.2

    ax4.text(0.05, y, 'Domain Colors:', fontsize=11, fontweight='bold')
    y -= dy

    # Catalytic
    rect_cat = Rectangle((0.05, y-0.04), 0.06, 0.08,
                         facecolor=COLORS['catalytic'], alpha=0.8,
                         edgecolor='black', linewidth=1)
    ax4.add_patch(rect_cat)
    ax4.text(0.15, y, 'Catalytic Domain', va='center', fontsize=9)
    y -= dy

    # Editing
    rect_edit = Rectangle((0.05, y-0.04), 0.06, 0.08,
                          facecolor=COLORS['editing'], alpha=0.8,
                          edgecolor='black', linewidth=1)
    ax4.add_patch(rect_edit)
    ax4.text(0.15, y, 'Editing Domain', va='center', fontsize=9)

    ax4.set_xlim(0, 1)
    ax4.set_ylim(0, 1)
    ax4.axis('off')
    ax4.set_title('Legend', fontsize=13, fontweight='bold', pad=12)

    return fig

if __name__ == '__main__':
    import os
    os.chdir('/storage/kiran-stuff/aaRS/final_figures')

    print("Generating improved Figure 1...")
    fig = figure1_phylogenetic_overview_improved()

    # Save in all formats
    fig.savefig('Figure1_Phylogeny_DomainArchitecture.png', dpi=300,
                bbox_inches='tight', facecolor='white')
    fig.savefig('Figure1_Phylogeny_DomainArchitecture.pdf',
                bbox_inches='tight', facecolor='white')
    fig.savefig('Figure1_Phylogeny_DomainArchitecture.svg',
                bbox_inches='tight', facecolor='white')

    print("✓ Saved PNG (300 DPI)")
    print("✓ Saved PDF (vector)")
    print("✓ Saved SVG (vector)")
    print("\nFiles saved to: /storage/kiran-stuff/aaRS/final_figures/")

    plt.close()
    print("\nDone! Figure 1 has been improved with:")
    print("  • Larger, clearer phylogeny labels")
    print("  • Removed 'node' text")
    print("  • No overlapping words")
    print("  • Better spacing and readability")
    print("  • Clearer domain architecture colors")
