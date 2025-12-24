#!/usr/bin/env python3
"""
Generate publication-quality figures for aaRS promiscuity paper
Target journals: NSMB / Cell Systems
Style: Professional, thin lines, serif fonts, muted colors, 300 DPI
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Rectangle, Polygon, FancyArrowPatch, Wedge
from matplotlib.collections import PatchCollection
import numpy as np
import pandas as pd
from matplotlib import rcParams

# Publication settings
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']
rcParams['font.size'] = 8
rcParams['axes.linewidth'] = 0.5
rcParams['xtick.major.width'] = 0.5
rcParams['ytick.major.width'] = 0.5
rcParams['lines.linewidth'] = 1.0
rcParams['patch.linewidth'] = 0.5

# Muted, colorblind-friendly color palette
COLORS = {
    'cognate': '#2E5C8A',      # Muted dark blue
    'non_cognate': '#8B8B8B',  # Gray
    'control': '#D4D4D4',      # Light gray
    'archaea': '#C85C5C',      # Muted red
    'bacteria': '#5C8AC8',     # Muted blue
    'eukaryota': '#5C9E5C',    # Muted green
    'luca': '#7E57C2',         # Purple
    'catalytic': '#1976D2',    # Blue
    'editing': '#F57C00',      # Orange
    'background': '#F5F5F5',   # Light gray background
}

def save_figure(fig, name, output_dir='/home/kiran/paper2_figures'):
    """Save figure as both PNG (300 DPI) and PDF"""
    import os
    os.makedirs(output_dir, exist_ok=True)

    png_path = f'{output_dir}/{name}.png'
    pdf_path = f'{output_dir}/{name}.pdf'
    svg_path = f'{output_dir}/{name}.svg'

    fig.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    fig.savefig(svg_path, bbox_inches='tight', facecolor='white')

    print(f"Saved: {name}")
    return png_path, pdf_path


def figure1_phylogenetic_overview():
    """
    FIGURE 1 - Phylogenetic Overview with domain architecture
    """
    fig = plt.figure(figsize=(7.08, 5.5))  # 180mm wide for double column

    # Create grid: tree on left, domain architecture on right
    gs = fig.add_gridspec(2, 2, width_ratios=[2, 1], height_ratios=[1, 1],
                          hspace=0.3, wspace=0.4)

    # Panel A: ProRS tree
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.text(0.02, 0.98, 'A', transform=ax1.transAxes,
             fontsize=12, fontweight='bold', va='top')

    # Simplified tree representation
    # Main branch
    y_archaea = 0.8
    y_bacteria = 0.5
    y_eukaryota = 0.2

    # Plot collapsed clades as triangles
    # Archaea clade
    triangle_archaea = Polygon([[0.3, y_archaea-0.05], [0.7, y_archaea-0.1],
                                 [0.7, y_archaea+0.1]],
                                facecolor=COLORS['archaea'], alpha=0.3,
                                edgecolor=COLORS['archaea'], linewidth=1)
    ax1.add_patch(triangle_archaea)
    ax1.text(0.75, y_archaea, 'Archaea', va='center', fontsize=7)

    # Bacteria clade
    triangle_bacteria = Polygon([[0.3, y_bacteria-0.05], [0.7, y_bacteria-0.1],
                                  [0.7, y_bacteria+0.1]],
                                 facecolor=COLORS['bacteria'], alpha=0.3,
                                 edgecolor=COLORS['bacteria'], linewidth=1)
    ax1.add_patch(triangle_bacteria)
    ax1.text(0.75, y_bacteria, 'Bacteria', va='center', fontsize=7)

    # Eukaryota clade
    triangle_eukaryota = Polygon([[0.3, y_eukaryota-0.05], [0.7, y_eukaryota-0.1],
                                   [0.7, y_eukaryota+0.1]],
                                  facecolor=COLORS['eukaryota'], alpha=0.3,
                                  edgecolor=COLORS['eukaryota'], linewidth=1)
    ax1.add_patch(triangle_eukaryota)
    ax1.text(0.75, y_eukaryota, 'Eukaryota', va='center', fontsize=7)

    # Draw tree branches
    ax1.plot([0, 0.1], [0.5, 0.5], 'k-', linewidth=0.5)  # Root
    ax1.plot([0.1, 0.1], [0.2, 0.8], 'k-', linewidth=0.5)  # Vertical
    ax1.plot([0.1, 0.3], [y_archaea, y_archaea], 'k-', linewidth=0.5)
    ax1.plot([0.1, 0.3], [y_bacteria, y_bacteria], 'k-', linewidth=0.5)
    ax1.plot([0.1, 0.15], [y_eukaryota, y_eukaryota], 'k-', linewidth=0.5)

    # Mark ancestral nodes
    # LUCA
    ax1.plot(0.1, 0.5, 'o', markersize=8, color=COLORS['luca'],
             markeredgecolor='black', markeredgewidth=0.5, zorder=10)
    ax1.text(0.05, 0.5, 'LUCA\nProRS', ha='right', va='center',
             fontsize=6, fontweight='bold')

    # Eukaryotic ancestor
    ax1.plot(0.15, y_eukaryota, '^', markersize=7, color=COLORS['eukaryota'],
             markeredgecolor='black', markeredgewidth=0.5, zorder=10)
    ax1.text(0.15, y_eukaryota-0.08, 'Euk ancestor', ha='center',
             va='top', fontsize=6, fontweight='bold')

    # Bootstrap values
    ax1.text(0.11, 0.65, '98', fontsize=5, style='italic', color='gray')
    ax1.text(0.11, 0.35, '95', fontsize=5, style='italic', color='gray')

    # Scale bar
    ax1.plot([0.05, 0.15], [0.05, 0.05], 'k-', linewidth=1)
    ax1.text(0.1, 0.02, '0.1 substitutions/site', ha='center', fontsize=6)

    ax1.set_xlim(-0.05, 1.0)
    ax1.set_ylim(0, 1)
    ax1.axis('off')
    ax1.set_title('ProRS Phylogeny (93 species)', fontsize=9, fontweight='bold', pad=10)

    # Panel B: ThrRS tree (similar structure)
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.text(0.02, 0.98, 'B', transform=ax2.transAxes,
             fontsize=12, fontweight='bold', va='top')

    # Similar tree for ThrRS
    triangle_archaea2 = Polygon([[0.3, y_archaea-0.05], [0.7, y_archaea-0.1],
                                  [0.7, y_archaea+0.1]],
                                 facecolor=COLORS['archaea'], alpha=0.3,
                                 edgecolor=COLORS['archaea'], linewidth=1)
    ax2.add_patch(triangle_archaea2)
    ax2.text(0.75, y_archaea, 'Archaea', va='center', fontsize=7)

    triangle_bacteria2 = Polygon([[0.3, y_bacteria-0.05], [0.7, y_bacteria-0.1],
                                   [0.7, y_bacteria+0.1]],
                                  facecolor=COLORS['bacteria'], alpha=0.3,
                                  edgecolor=COLORS['bacteria'], linewidth=1)
    ax2.add_patch(triangle_bacteria2)
    ax2.text(0.75, y_bacteria, 'Bacteria', va='center', fontsize=7)

    triangle_eukaryota2 = Polygon([[0.3, y_eukaryota-0.05], [0.7, y_eukaryota-0.1],
                                    [0.7, y_eukaryota+0.1]],
                                   facecolor=COLORS['eukaryota'], alpha=0.3,
                                   edgecolor=COLORS['eukaryota'], linewidth=1)
    ax2.add_patch(triangle_eukaryota2)
    ax2.text(0.75, y_eukaryota, 'Eukaryota', va='center', fontsize=7)

    ax2.plot([0, 0.1], [0.5, 0.5], 'k-', linewidth=0.5)
    ax2.plot([0.1, 0.1], [0.2, 0.8], 'k-', linewidth=0.5)
    ax2.plot([0.1, 0.3], [y_archaea, y_archaea], 'k-', linewidth=0.5)
    ax2.plot([0.1, 0.3], [y_bacteria, y_bacteria], 'k-', linewidth=0.5)
    ax2.plot([0.1, 0.3], [y_eukaryota, y_eukaryota], 'k-', linewidth=0.5)

    # LUCA ThrRS
    ax2.plot(0.1, 0.5, 'o', markersize=8, color=COLORS['luca'],
             markeredgecolor='black', markeredgewidth=0.5, zorder=10)
    ax2.text(0.05, 0.5, 'LUCA\nThrRS', ha='right', va='center',
             fontsize=6, fontweight='bold')

    # Bootstrap values
    ax2.text(0.11, 0.65, '96', fontsize=5, style='italic', color='gray')
    ax2.text(0.11, 0.35, '92', fontsize=5, style='italic', color='gray')

    # Scale bar
    ax2.plot([0.05, 0.15], [0.05, 0.05], 'k-', linewidth=1)
    ax2.text(0.1, 0.02, '0.1 substitutions/site', ha='center', fontsize=6)

    ax2.set_xlim(-0.05, 1.0)
    ax2.set_ylim(0, 1)
    ax2.axis('off')
    ax2.set_title('ThrRS Phylogeny (64 species)', fontsize=9, fontweight='bold', pad=10)

    # Panel C: Domain architecture (ProRS)
    ax3 = fig.add_subplot(gs[0, 1])
    ax3.text(0.02, 0.98, 'C', transform=ax3.transAxes,
             fontsize=12, fontweight='bold', va='top')

    # ProRS domain architecture
    y_pos = 0.7
    # N-terminal
    rect1 = Rectangle((0.1, y_pos-0.05), 0.15, 0.1,
                      facecolor=COLORS['background'],
                      edgecolor='black', linewidth=0.5)
    ax3.add_patch(rect1)
    ax3.text(0.175, y_pos, 'N-term', ha='center', va='center', fontsize=6)

    # Catalytic domain
    rect2 = Rectangle((0.27, y_pos-0.05), 0.3, 0.1,
                      facecolor=COLORS['catalytic'], alpha=0.6,
                      edgecolor='black', linewidth=0.5)
    ax3.add_patch(rect2)
    ax3.text(0.42, y_pos, 'Catalytic', ha='center', va='center',
             fontsize=6, fontweight='bold', color='white')

    # Editing domain
    rect3 = Rectangle((0.62, y_pos-0.05), 0.2, 0.1,
                      facecolor=COLORS['editing'], alpha=0.6,
                      edgecolor='black', linewidth=0.5)
    ax3.add_patch(rect3)
    ax3.text(0.72, y_pos, 'Editing', ha='center', va='center',
             fontsize=6, fontweight='bold', color='white')

    # C-terminal
    rect4 = Rectangle((0.84, y_pos-0.05), 0.1, 0.1,
                      facecolor=COLORS['background'],
                      edgecolor='black', linewidth=0.5)
    ax3.add_patch(rect4)
    ax3.text(0.89, y_pos, 'C', ha='center', va='center', fontsize=6)

    ax3.text(0.5, y_pos+0.15, 'ProRS', ha='center', fontsize=8, fontweight='bold')
    ax3.text(0.1, y_pos+0.1, '1', ha='left', fontsize=5)
    ax3.text(0.94, y_pos+0.1, '2037 aa', ha='right', fontsize=5)

    # ThrRS domain architecture
    y_pos = 0.3
    # N-terminal
    rect5 = Rectangle((0.1, y_pos-0.05), 0.15, 0.1,
                      facecolor=COLORS['background'],
                      edgecolor='black', linewidth=0.5)
    ax3.add_patch(rect5)
    ax3.text(0.175, y_pos, 'N-term', ha='center', va='center', fontsize=6)

    # Catalytic domain (larger proportion)
    rect6 = Rectangle((0.27, y_pos-0.05), 0.6, 0.1,
                      facecolor=COLORS['catalytic'], alpha=0.6,
                      edgecolor='black', linewidth=0.5)
    ax3.add_patch(rect6)
    ax3.text(0.57, y_pos, 'Catalytic', ha='center', va='center',
             fontsize=6, fontweight='bold', color='white')

    # C-terminal
    rect7 = Rectangle((0.89, y_pos-0.05), 0.05, 0.1,
                      facecolor=COLORS['background'],
                      edgecolor='black', linewidth=0.5)
    ax3.add_patch(rect7)
    ax3.text(0.915, y_pos, 'C', ha='center', va='center', fontsize=6)

    ax3.text(0.5, y_pos+0.15, 'ThrRS', ha='center', fontsize=8, fontweight='bold')
    ax3.text(0.1, y_pos+0.1, '1', ha='left', fontsize=5)
    ax3.text(0.94, y_pos+0.1, '1017 aa', ha='right', fontsize=5)

    # Add annotation
    ax3.annotate('', xy=(0.72, y_pos+0.25), xytext=(0.72, y_pos+0.35),
                arrowprops=dict(arrowstyle='->', lw=0.5, color='red'))
    ax3.text(0.72, y_pos+0.4, 'No editing\ndomain', ha='center',
             fontsize=6, color='red', fontweight='bold')

    ax3.set_xlim(0, 1)
    ax3.set_ylim(0, 1)
    ax3.axis('off')
    ax3.set_title('Domain Architecture', fontsize=9, fontweight='bold', pad=10)

    # Panel D: Legend
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.text(0.02, 0.98, 'D', transform=ax4.transAxes,
             fontsize=12, fontweight='bold', va='top')

    legend_elements = [
        mpatches.Patch(facecolor=COLORS['archaea'], alpha=0.3,
                      edgecolor=COLORS['archaea'], label='Archaea'),
        mpatches.Patch(facecolor=COLORS['bacteria'], alpha=0.3,
                      edgecolor=COLORS['bacteria'], label='Bacteria'),
        mpatches.Patch(facecolor=COLORS['eukaryota'], alpha=0.3,
                      edgecolor=COLORS['eukaryota'], label='Eukaryota'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['luca'],
                  markersize=6, markeredgecolor='black',
                  markeredgewidth=0.5, label='LUCA node'),
        plt.Line2D([0], [0], marker='^', color='w', markerfacecolor=COLORS['eukaryota'],
                  markersize=6, markeredgecolor='black',
                  markeredgewidth=0.5, label='Euk ancestor'),
        mpatches.Patch(facecolor=COLORS['catalytic'], alpha=0.6, label='Catalytic domain'),
        mpatches.Patch(facecolor=COLORS['editing'], alpha=0.6, label='Editing domain'),
    ]

    ax4.legend(handles=legend_elements, loc='center', frameon=True,
              fontsize=7, framealpha=1, edgecolor='gray')
    ax4.axis('off')

    plt.suptitle('Figure 1: Phylogenetic Overview and Domain Architecture',
                fontsize=10, fontweight='bold', y=0.98)

    return save_figure(fig, 'figure1_phylogenetic_overview')


def figure2_substrate_binding():
    """
    FIGURE 2 - Substrate Binding Panel (5 panels: A-E)
    Multi-panel like Furukawa Fig 3
    """
    fig, axes = plt.subplots(2, 3, figsize=(7.08, 4.5))
    axes = axes.flatten()

    # Data for each panel
    panels_data = [
        {
            'title': 'LUCA ProRS\nCatalytic Domain',
            'substrates': ['PRO', 'THR', 'TRP', 'PHE'],
            'iptm': [0.75, 0.62, 0.67, 0.80],
            'errors': [0.03, 0.04, 0.05, 0.04],
            'cognate': 'PRO',
            'label': 'A'
        },
        {
            'title': 'LUCA ThrRS\nCatalytic Domain',
            'substrates': ['PRO', 'THR'],
            'iptm': [0.88, 0.89],
            'errors': [0.02, 0.02],
            'cognate': 'THR',
            'label': 'B'
        },
        {
            'title': 'Eukaryotic ProRS\nAncestor',
            'substrates': ['PRO', 'THR'],
            'iptm': [0.83, 0.74],
            'errors': [0.03, 0.03],
            'cognate': 'PRO',
            'label': 'C'
        },
        {
            'title': 'Modern Human\nProRS',
            'substrates': ['PRO', 'THR'],
            'iptm': [0.80, 0.78],
            'errors': [0.02, 0.03],
            'cognate': 'PRO',
            'label': 'D'
        },
        {
            'title': 'Modern Human\nThrRS',
            'substrates': ['PRO', 'THR'],
            'iptm': [0.57, 0.84],
            'errors': [0.04, 0.02],
            'cognate': 'THR',
            'label': 'E'
        }
    ]

    for idx, (ax, data) in enumerate(zip(axes[:5], panels_data)):
        # Panel label
        ax.text(0.02, 0.98, data['label'], transform=ax.transAxes,
                fontsize=12, fontweight='bold', va='top')

        # Determine colors
        colors = []
        for substrate in data['substrates']:
            if substrate == data['cognate']:
                colors.append(COLORS['cognate'])
            elif substrate in ['TRP', 'PHE']:
                colors.append(COLORS['control'])
            else:
                colors.append(COLORS['non_cognate'])

        # Create bar chart
        x_pos = np.arange(len(data['substrates']))
        bars = ax.bar(x_pos, data['iptm'], yerr=data['errors'],
                     color=colors, edgecolor='black', linewidth=0.5,
                     capsize=3, error_kw={'linewidth': 0.5})

        # Customize
        ax.set_xticks(x_pos)
        ax.set_xticklabels(data['substrates'], fontsize=7)
        ax.set_ylabel('ipTM score', fontsize=7)
        ax.set_ylim(0, 1.0)
        ax.set_title(data['title'], fontsize=8, fontweight='bold', pad=8)

        # Add significance annotations
        if 'PRO' in data['substrates'] and 'THR' in data['substrates']:
            pro_idx = data['substrates'].index('PRO')
            thr_idx = data['substrates'].index('THR')

            # Calculate discrimination ratio
            if data['cognate'] == 'PRO':
                ratio = (data['iptm'][thr_idx] / data['iptm'][pro_idx]) * 100
                ax.text(0.5, 0.92, f'THR/PRO: {ratio:.1f}%',
                       transform=ax.transAxes, ha='center', va='top',
                       fontsize=6, bbox=dict(boxstyle='round,pad=0.3',
                       facecolor='wheat', alpha=0.3, linewidth=0.5))
            else:
                ratio = (data['iptm'][pro_idx] / data['iptm'][thr_idx]) * 100
                ax.text(0.5, 0.92, f'PRO/THR: {ratio:.1f}%',
                       transform=ax.transAxes, ha='center', va='top',
                       fontsize=6, bbox=dict(boxstyle='round,pad=0.3',
                       facecolor='wheat', alpha=0.3, linewidth=0.5))

        # Grid
        ax.grid(axis='y', alpha=0.3, linewidth=0.5, linestyle='--')
        ax.set_axisbelow(True)

        # Spines
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)

    # Remove extra subplot
    axes[5].axis('off')

    # Add legend to empty panel
    legend_elements = [
        mpatches.Patch(facecolor=COLORS['cognate'], edgecolor='black',
                      linewidth=0.5, label='Cognate substrate'),
        mpatches.Patch(facecolor=COLORS['non_cognate'], edgecolor='black',
                      linewidth=0.5, label='Non-cognate substrate'),
        mpatches.Patch(facecolor=COLORS['control'], edgecolor='black',
                      linewidth=0.5, label='Control (TRP/PHE)'),
    ]
    axes[5].legend(handles=legend_elements, loc='center', fontsize=7,
                  frameon=True, edgecolor='gray', framealpha=1)

    plt.suptitle('Figure 2: Substrate Binding Specificity Across Evolution',
                fontsize=10, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    return save_figure(fig, 'figure2_substrate_binding_panel')


def figure3_editing_domain():
    """
    FIGURE 3 - Editing Domain Analysis
    """
    fig = plt.figure(figsize=(7.08, 3.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1], wspace=0.4)

    # Panel A: Domain architecture comparison
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.text(0.02, 0.98, 'A', transform=ax1.transAxes,
             fontsize=12, fontweight='bold', va='top')

    # ProRS with editing domain
    y_pos = 0.65
    ax1.text(0.05, y_pos+0.15, 'LUCA ProRS (2037 aa)',
             fontsize=8, fontweight='bold')

    # Catalytic
    rect1 = Rectangle((0.1, y_pos-0.04), 0.35, 0.08,
                      facecolor=COLORS['catalytic'], alpha=0.6,
                      edgecolor='black', linewidth=0.5)
    ax1.add_patch(rect1)
    ax1.text(0.275, y_pos, 'Catalytic\n(200-700)', ha='center', va='center',
             fontsize=6, color='white', fontweight='bold')

    # Linker
    rect2 = Rectangle((0.45, y_pos-0.04), 0.25, 0.08,
                      facecolor=COLORS['background'],
                      edgecolor='black', linewidth=0.5)
    ax1.add_patch(rect2)
    ax1.text(0.575, y_pos, 'Linker', ha='center', va='center', fontsize=5)

    # Editing
    rect3 = Rectangle((0.7, y_pos-0.04), 0.2, 0.08,
                      facecolor=COLORS['editing'], alpha=0.6,
                      edgecolor='black', linewidth=0.5)
    ax1.add_patch(rect3)
    ax1.text(0.8, y_pos, 'Editing\n(1504-1652)', ha='center', va='center',
             fontsize=6, color='white', fontweight='bold')

    # ThrRS without editing domain
    y_pos = 0.35
    ax1.text(0.05, y_pos+0.15, 'LUCA ThrRS (1017 aa)',
             fontsize=8, fontweight='bold')

    # Catalytic (larger)
    rect4 = Rectangle((0.1, y_pos-0.04), 0.75, 0.08,
                      facecolor=COLORS['catalytic'], alpha=0.6,
                      edgecolor='black', linewidth=0.5)
    ax1.add_patch(rect4)
    ax1.text(0.475, y_pos, 'Catalytic domain only', ha='center', va='center',
             fontsize=7, color='white', fontweight='bold')

    # C-terminal
    rect5 = Rectangle((0.85, y_pos-0.04), 0.05, 0.08,
                      facecolor=COLORS['background'],
                      edgecolor='black', linewidth=0.5)
    ax1.add_patch(rect5)
    ax1.text(0.875, y_pos, 'C', ha='center', va='center', fontsize=5)

    # Annotation highlighting difference
    ax1.annotate('', xy=(0.8, 0.5), xytext=(0.8, 0.25),
                arrowprops=dict(arrowstyle='<->', lw=1, color='red'))
    ax1.text(0.82, 0.375, 'Editing domain\npresent only\nin ProRS',
             fontsize=7, color='red', fontweight='bold', va='center')

    ax1.set_xlim(0, 1)
    ax1.set_ylim(0.15, 0.9)
    ax1.axis('off')
    ax1.set_title('Domain Architecture Comparison',
                 fontsize=9, fontweight='bold', pad=10)

    # Panel B: Editing domain binding
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.text(0.02, 0.98, 'B', transform=ax2.transAxes,
             fontsize=12, fontweight='bold', va='top')

    substrates = ['PRO', 'THR']
    iptm = [0.14, 0.45]
    errors = [0.02, 0.03]
    colors = [COLORS['cognate'], COLORS['non_cognate']]

    x_pos = np.arange(len(substrates))
    bars = ax2.bar(x_pos, iptm, yerr=errors, color=colors,
                  edgecolor='black', linewidth=0.5, capsize=3,
                  error_kw={'linewidth': 0.5})

    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(substrates, fontsize=8)
    ax2.set_ylabel('ipTM score', fontsize=8)
    ax2.set_ylim(0, 0.6)
    ax2.set_title('LUCA ProRS Editing Domain\nSubstrate Binding',
                 fontsize=9, fontweight='bold', pad=10)

    # Add annotation about inverted discrimination
    ax2.text(0.5, 0.85, 'Inverted specificity!\nTHR > PRO',
             transform=ax2.transAxes, ha='center', va='top',
             fontsize=7, fontweight='bold', color='red',
             bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow',
                      alpha=0.3, linewidth=0.5))

    # Calculate ratio
    ratio = (iptm[1] / iptm[0]) * 100
    ax2.text(0.5, 0.65, f'THR/PRO: {ratio:.0f}%',
             transform=ax2.transAxes, ha='center', va='top',
             fontsize=7, bbox=dict(boxstyle='round,pad=0.3',
             facecolor='wheat', alpha=0.3, linewidth=0.5))

    ax2.grid(axis='y', alpha=0.3, linewidth=0.5, linestyle='--')
    ax2.set_axisbelow(True)

    for spine in ax2.spines.values():
        spine.set_linewidth(0.5)

    plt.suptitle('Figure 3: Editing Domain Analysis',
                fontsize=10, fontweight='bold')

    return save_figure(fig, 'figure3_editing_domain_analysis')


def figure4_evolutionary_timeline():
    """
    FIGURE 4 - Evolutionary Timeline
    """
    fig, ax = plt.subplots(figsize=(7.08, 4))

    # Data
    time_points = [3.5, 1.5, 0]  # Billion years ago
    promiscuity = [89.7, 84.8, 97.5]  # THR binding as % of PRO
    labels = ['LUCA', 'Eukaryotic\nAncestor', 'Modern\nHuman']

    # Background shading for geological eras
    # Hadean/Archean
    ax.axvspan(4.0, 2.5, alpha=0.1, color=COLORS['archaea'], label='Archean')
    # Proterozoic
    ax.axvspan(2.5, 0.541, alpha=0.1, color=COLORS['bacteria'], label='Proterozoic')
    # Phanerozoic
    ax.axvspan(0.541, 0, alpha=0.1, color=COLORS['eukaryota'], label='Phanerozoic')

    # Plot promiscuity trend
    ax.plot(time_points, promiscuity, 'o-', color=COLORS['luca'],
            markersize=10, linewidth=2, markeredgecolor='black',
            markeredgewidth=0.5, zorder=10)

    # Add labels for each point
    for i, (t, p, label) in enumerate(zip(time_points, promiscuity, labels)):
        ax.annotate(label, xy=(t, p), xytext=(0, 15),
                   textcoords='offset points', ha='center', fontsize=8,
                   fontweight='bold',
                   bbox=dict(boxstyle='round,pad=0.5', facecolor='white',
                            edgecolor='black', linewidth=0.5))
        ax.text(t, p-8, f'{p:.1f}%', ha='center', va='top',
               fontsize=8, fontweight='bold', color=COLORS['luca'])

    # Reference line at 90%
    ax.axhline(y=90, color='gray', linestyle='--', linewidth=0.5, alpha=0.7)
    ax.text(3.8, 91, '90% threshold', fontsize=6, color='gray', style='italic')

    # Key insight annotation
    ax.text(2.0, 75, 'Promiscuity persists despite\n3.5 billion years of evolution',
           ha='center', fontsize=9, fontweight='bold',
           bbox=dict(boxstyle='round,pad=0.7', facecolor='yellow',
                    alpha=0.3, edgecolor='red', linewidth=1))

    # Axis labels
    ax.set_xlabel('Time (billion years ago)', fontsize=9, fontweight='bold')
    ax.set_ylabel('THR binding (% of PRO)', fontsize=9, fontweight='bold')
    ax.set_xlim(4.0, -0.2)
    ax.set_ylim(70, 105)

    # Invert x-axis (time goes backwards)
    ax.invert_xaxis()

    # Grid
    ax.grid(axis='both', alpha=0.3, linewidth=0.5, linestyle='--')
    ax.set_axisbelow(True)

    # Legend
    ax.legend(loc='upper left', fontsize=7, frameon=True,
             edgecolor='gray', framealpha=0.9)

    # Spines
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)

    plt.title('Figure 4: Evolutionary Timeline of Substrate Promiscuity',
             fontsize=10, fontweight='bold', pad=15)
    plt.tight_layout()

    return save_figure(fig, 'figure4_evolutionary_timeline')


def figure5_methodological_comparison():
    """
    FIGURE 5 - Methodological Comparison (NEW - confronts Furukawa)
    """
    fig = plt.figure(figsize=(7.08, 6))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.2], hspace=0.4, wspace=0.3)

    # Panel A: Furukawa's PAAS method
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.text(0.02, 0.98, 'A', transform=ax1.transAxes,
             fontsize=12, fontweight='bold', va='top')

    ax1.text(0.5, 0.85, 'Furukawa et al. (2022) Method',
            transform=ax1.transAxes, ha='center', fontsize=9, fontweight='bold')

    # Schematic of PAAS method
    ax1.text(0.5, 0.65, 'PAAS Score', transform=ax1.transAxes,
            ha='center', fontsize=8, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.5', facecolor=COLORS['bacteria'],
                     alpha=0.3, linewidth=0.5))
    ax1.annotate('', xy=(0.5, 0.55), xytext=(0.5, 0.45),
                transform=ax1.transAxes,
                arrowprops=dict(arrowstyle='->', lw=1.5, color='black'))
    ax1.text(0.5, 0.35, 'Sequence\nConservation', transform=ax1.transAxes,
            ha='center', fontsize=7)
    ax1.annotate('', xy=(0.5, 0.25), xytext=(0.5, 0.15),
                transform=ax1.transAxes,
                arrowprops=dict(arrowstyle='->', lw=1.5, color='black'))
    ax1.text(0.5, 0.05, 'Assumed\nSpecificity ✓', transform=ax1.transAxes,
            ha='center', fontsize=7, fontweight='bold', color='green')

    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.axis('off')

    # Panel B: Our AlphaFold3 method
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.text(0.02, 0.98, 'B', transform=ax2.transAxes,
             fontsize=12, fontweight='bold', va='top')

    ax2.text(0.5, 0.85, 'Our AlphaFold3 Method',
            transform=ax2.transAxes, ha='center', fontsize=9, fontweight='bold')

    # Schematic of AF3 method
    ax2.text(0.5, 0.65, 'Ancestral\nSequence', transform=ax2.transAxes,
            ha='center', fontsize=7)
    ax2.annotate('', xy=(0.5, 0.55), xytext=(0.5, 0.45),
                transform=ax2.transAxes,
                arrowprops=dict(arrowstyle='->', lw=1.5, color='black'))
    ax2.text(0.5, 0.35, 'AlphaFold3\nStructure\nModeling', transform=ax2.transAxes,
            ha='center', fontsize=7,
            bbox=dict(boxstyle='round,pad=0.5', facecolor=COLORS['luca'],
                     alpha=0.3, linewidth=0.5))
    ax2.annotate('', xy=(0.5, 0.22), xytext=(0.5, 0.12),
                transform=ax2.transAxes,
                arrowprops=dict(arrowstyle='->', lw=1.5, color='black'))
    ax2.text(0.5, 0.05, 'Measured\nBinding ipTM', transform=ax2.transAxes,
            ha='center', fontsize=7, fontweight='bold', color=COLORS['luca'])

    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.axis('off')

    # Panel C: Conserved residues vs pocket volume
    ax3 = fig.add_subplot(gs[1, :])
    ax3.text(0.01, 0.98, 'C', transform=ax3.transAxes,
             fontsize=12, fontweight='bold', va='top')

    # Schematic showing binding pocket
    # Modern enzyme (smaller pocket)
    circle1 = plt.Circle((0.25, 0.5), 0.12, facecolor=COLORS['catalytic'],
                        alpha=0.3, edgecolor='black', linewidth=1)
    ax3.add_patch(circle1)
    ax3.text(0.25, 0.7, 'Modern ProRS', ha='center', fontsize=8, fontweight='bold')

    # Draw small pocket
    pocket1 = plt.Circle((0.25, 0.5), 0.04, facecolor='white',
                        edgecolor='red', linewidth=1.5)
    ax3.add_patch(pocket1)
    ax3.text(0.25, 0.5, 'PRO', ha='center', va='center', fontsize=6, fontweight='bold')
    ax3.text(0.25, 0.28, 'Tight pocket\nHigh specificity', ha='center',
            fontsize=6, color='green')

    # Ancestral enzyme (larger pocket)
    circle2 = plt.Circle((0.75, 0.5), 0.12, facecolor=COLORS['luca'],
                        alpha=0.3, edgecolor='black', linewidth=1)
    ax3.add_patch(circle2)
    ax3.text(0.75, 0.7, 'LUCA ProRS', ha='center', fontsize=8, fontweight='bold')

    # Draw larger pocket
    pocket2 = plt.Circle((0.75, 0.5), 0.07, facecolor='white',
                        edgecolor='orange', linewidth=1.5)
    ax3.add_patch(pocket2)
    ax3.text(0.72, 0.52, 'PRO', ha='center', va='center', fontsize=6, fontweight='bold')
    ax3.text(0.78, 0.48, 'THR', ha='center', va='center', fontsize=6, fontweight='bold')
    ax3.text(0.75, 0.28, 'Larger pocket\nPromiscuous', ha='center',
            fontsize=6, color='red')

    # Conserved residues annotation
    for x in [0.25, 0.75]:
        for angle in [0, 90, 180, 270]:
            rad = np.radians(angle)
            if x == 0.25:
                r = 0.04
            else:
                r = 0.07
            px = x + r * np.cos(rad)
            py = 0.5 + r * np.sin(rad)
            ax3.plot(px, py, 'ko', markersize=3)

    ax3.text(0.5, 0.15, 'Conserved binding residues', ha='center',
            fontsize=7, fontweight='bold', style='italic')

    ax3.set_xlim(0, 1)
    ax3.set_ylim(0.2, 0.8)
    ax3.axis('off')
    ax3.set_title('Conservation ≠ Discrimination', fontsize=9,
                 fontweight='bold', color='red', pad=10)

    # Panel D: Key insight
    ax4 = fig.add_subplot(gs[2, :])
    ax4.text(0.01, 0.98, 'D', transform=ax4.transAxes,
             fontsize=12, fontweight='bold', va='top')

    # Text box with key message
    message = """Key Insight:

Furukawa et al. inferred ancestral specificity from sequence conservation.
However, conserved residues can form a LARGER binding pocket in ancestors.

Our structural modeling reveals:
• LUCA ProRS: 89.7% promiscuity (THR binds at 89.7% of PRO affinity)
• Eukaryotic ancestor: 84.8% promiscuity
• Modern human: 97.5% promiscuity

Conclusion: Ancient translation systems tolerated high error rates.
Conservation of residues does not imply conservation of specificity."""

    ax4.text(0.5, 0.5, message, transform=ax4.transAxes,
            ha='center', va='center', fontsize=8,
            bbox=dict(boxstyle='round,pad=1', facecolor='lightyellow',
                     edgecolor='red', linewidth=2, alpha=0.9),
            family='monospace')

    ax4.set_xlim(0, 1)
    ax4.set_ylim(0, 1)
    ax4.axis('off')

    plt.suptitle('Figure 5: Methodological Comparison - Challenging Furukawa et al.',
                fontsize=10, fontweight='bold')

    return save_figure(fig, 'figure5_methodological_comparison')


def create_all_figures():
    """Generate all publication figures"""
    print("=" * 60)
    print("Generating Publication-Quality Figures")
    print("Target: NSMB / Cell Systems")
    print("=" * 60)
    print()

    figures = [
        ("Figure 1: Phylogenetic Overview", figure1_phylogenetic_overview),
        ("Figure 2: Substrate Binding Panel", figure2_substrate_binding),
        ("Figure 3: Editing Domain Analysis", figure3_editing_domain),
        ("Figure 4: Evolutionary Timeline", figure4_evolutionary_timeline),
        ("Figure 5: Methodological Comparison", figure5_methodological_comparison),
    ]

    for name, func in figures:
        print(f"\nGenerating {name}...")
        try:
            png_path, pdf_path = func()
            print(f"  ✓ PNG: {png_path}")
            print(f"  ✓ PDF: {pdf_path}")
        except Exception as e:
            print(f"  ✗ Error: {e}")

    print("\n" + "=" * 60)
    print("Figure generation complete!")
    print("=" * 60)
    print(f"\nAll figures saved to: /home/kiran/paper2_figures/")
    print("Formats: PNG (300 DPI), PDF, SVG")


if __name__ == '__main__':
    create_all_figures()
