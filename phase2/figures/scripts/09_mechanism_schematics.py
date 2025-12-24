#!/usr/bin/env python3
"""
Figure 2A: Double Sieve Mechanism Schematic
Creates a publication-quality diagram of the ProRS two-stage quality control system.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
import numpy as np

def create_double_sieve_figure():
    """Create the Double Sieve mechanism schematic."""

    fig, ax = plt.subplots(1, 1, figsize=(10, 12))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 12)
    ax.axis('off')

    # Colors
    purple = '#9b59b6'  # Catalytic domain
    orange = '#f39c12'  # Editing domain
    green = '#27ae60'   # Correct/Pass
    red = '#e74c3c'     # Error
    gray = '#95a5a6'    # Neutral
    dark = '#2c3e50'    # Text

    # ===== TITLE =====
    ax.text(5, 11.5, 'DOUBLE SIEVE MECHANISM', fontsize=22, fontweight='bold',
            ha='center', color=dark)
    ax.text(5, 11.0, 'ProRS Two-Stage Quality Control', fontsize=14,
            ha='center', color=gray, style='italic')

    # ===== AMINO ACID POOL (TOP) =====
    pool_box = FancyBboxPatch((2.5, 9.8), 5, 0.8, boxstyle="round,pad=0.05",
                               facecolor='#ecf0f1', edgecolor=dark, linewidth=2)
    ax.add_patch(pool_box)
    ax.text(5, 10.2, 'AMINO ACID POOL', fontsize=12, fontweight='bold',
            ha='center', color=dark)
    ax.text(5, 10.0, 'PRO, THR, ALA, VAL, LEU, ...', fontsize=10,
            ha='center', color=gray)

    # Arrow down
    ax.annotate('', xy=(5, 9.3), xytext=(5, 9.7),
                arrowprops=dict(arrowstyle='->', color=dark, lw=2))

    # ===== SIEVE 1: CATALYTIC DOMAIN =====
    sieve1_box = FancyBboxPatch((1.5, 7.3), 7, 1.8, boxstyle="round,pad=0.05",
                                 facecolor='#f5eef8', edgecolor=purple, linewidth=3)
    ax.add_patch(sieve1_box)

    ax.text(5, 8.8, 'SIEVE 1: CATALYTIC DOMAIN', fontsize=14, fontweight='bold',
            ha='center', color=purple)
    ax.text(5, 8.4, 'COARSE FILTER', fontsize=12, ha='center', color=purple)
    ax.text(5, 7.9, 'Binds most substrates with similar affinity', fontsize=10,
            ha='center', color=gray, style='italic')
    ax.text(5, 7.5, 'ALA = 98% | VAL = 97% | THR = 92% of PRO', fontsize=9,
            ha='center', color=purple)

    # Sieve holes (decorative)
    for x in [2.5, 3.5, 4.5, 5.5, 6.5, 7.5]:
        circle = plt.Circle((x, 7.35), 0.08, color=purple, alpha=0.5)
        ax.add_patch(circle)

    # ===== SUBSTRATES PASSING THROUGH =====
    # Arrow down
    ax.annotate('', xy=(5, 6.8), xytext=(5, 7.2),
                arrowprops=dict(arrowstyle='->', color=dark, lw=2))

    # Three pathways
    y_sub = 6.3

    # PRO (correct)
    pro_circle = plt.Circle((2.5, y_sub), 0.35, color=green, alpha=0.8)
    ax.add_patch(pro_circle)
    ax.text(2.5, y_sub, 'PRO', fontsize=10, fontweight='bold',
            ha='center', va='center', color='white')
    ax.text(2.5, y_sub - 0.55, '✓ Cognate', fontsize=9, ha='center', color=green)

    # THR (error)
    thr_circle = plt.Circle((5, y_sub), 0.35, color=red, alpha=0.8)
    ax.add_patch(thr_circle)
    ax.text(5, y_sub, 'THR', fontsize=10, fontweight='bold',
            ha='center', va='center', color='white')
    ax.text(5, y_sub - 0.55, '✗ Error', fontsize=9, ha='center', color=red)

    # ALA (error)
    ala_circle = plt.Circle((7.5, y_sub), 0.35, color=red, alpha=0.8)
    ax.add_patch(ala_circle)
    ax.text(7.5, y_sub, 'ALA', fontsize=10, fontweight='bold',
            ha='center', va='center', color='white')
    ax.text(7.5, y_sub - 0.55, '✗ Error', fontsize=9, ha='center', color=red)

    # Arrows to Sieve 2
    for x in [2.5, 5, 7.5]:
        ax.annotate('', xy=(x, 5.0), xytext=(x, 5.6),
                    arrowprops=dict(arrowstyle='->', color=dark, lw=1.5))

    # ===== SIEVE 2: EDITING DOMAIN =====
    sieve2_box = FancyBboxPatch((1.5, 3.0), 7, 1.8, boxstyle="round,pad=0.05",
                                 facecolor='#fef9e7', edgecolor=orange, linewidth=3)
    ax.add_patch(sieve2_box)

    ax.text(5, 4.5, 'SIEVE 2: EDITING DOMAIN', fontsize=14, fontweight='bold',
            ha='center', color=orange)
    ax.text(5, 4.1, 'FINE FILTER', fontsize=12, ha='center', color=orange)
    ax.text(5, 3.6, 'Preferentially binds error substrates', fontsize=10,
            ha='center', color=gray, style='italic')
    ax.text(5, 3.2, 'THR ranks #1 (0.87) | PRO ranks #2 (0.82)', fontsize=9,
            ha='center', color=orange)

    # Sieve holes (decorative)
    for x in [2.5, 3.5, 4.5, 5.5, 6.5, 7.5]:
        circle = plt.Circle((x, 3.05), 0.08, color=orange, alpha=0.5)
        ax.add_patch(circle)

    # ===== FINAL OUTPUTS =====
    # Arrows down
    for x in [2.5, 5, 7.5]:
        ax.annotate('', xy=(x, 2.0), xytext=(x, 2.9),
                    arrowprops=dict(arrowstyle='->', color=dark, lw=1.5))

    y_out = 1.5

    # PRO output (correct product)
    pro_out = FancyBboxPatch((1.7, y_out - 0.4), 1.6, 0.8, boxstyle="round,pad=0.02",
                              facecolor=green, edgecolor='white', linewidth=2)
    ax.add_patch(pro_out)
    ax.text(2.5, y_out, 'PRO-tRNA', fontsize=11, fontweight='bold',
            ha='center', va='center', color='white')
    ax.text(2.5, y_out - 0.65, '✓ PRODUCT', fontsize=9, ha='center',
            color=green, fontweight='bold')

    # THR output (hydrolyzed)
    thr_out = FancyBboxPatch((4.2, y_out - 0.4), 1.6, 0.8, boxstyle="round,pad=0.02",
                              facecolor=gray, edgecolor='white', linewidth=2, alpha=0.6)
    ax.add_patch(thr_out)
    ax.text(5, y_out, 'THR', fontsize=11, ha='center', va='center', color='white')
    ax.text(5, y_out - 0.65, '✗ HYDROLYZED', fontsize=9, ha='center',
            color=gray, fontweight='bold')

    # ALA output (hydrolyzed)
    ala_out = FancyBboxPatch((6.7, y_out - 0.4), 1.6, 0.8, boxstyle="round,pad=0.02",
                              facecolor=gray, edgecolor='white', linewidth=2, alpha=0.6)
    ax.add_patch(ala_out)
    ax.text(7.5, y_out, 'ALA', fontsize=11, ha='center', va='center', color='white')
    ax.text(7.5, y_out - 0.65, '✗ HYDROLYZED', fontsize=9, ha='center',
            color=gray, fontweight='bold')

    # ===== RESULT BOX =====
    result_box = FancyBboxPatch((2, 0.1), 6, 0.5, boxstyle="round,pad=0.02",
                                 facecolor='#d5f5e3', edgecolor=green, linewidth=2)
    ax.add_patch(result_box)
    ax.text(5, 0.35, 'RESULT: Only correct PRO-tRNA produced', fontsize=11,
            fontweight='bold', ha='center', color=dark)

    # ===== KEY INSIGHT BOX =====
    ax.text(0.3, 5.5, 'KEY:', fontsize=10, fontweight='bold', color=dark)
    ax.text(0.3, 5.1, 'Catalytic site is', fontsize=9, color=purple)
    ax.text(0.3, 4.8, 'PROMISCUOUS', fontsize=9, fontweight='bold', color=purple)
    ax.text(0.3, 4.4, 'Editing site is', fontsize=9, color=orange)
    ax.text(0.3, 4.1, 'SELECTIVE', fontsize=9, fontweight='bold', color=orange)

    plt.tight_layout()
    return fig

def create_zinc_filter_figure():
    """Create the Zinc Filter mechanism schematic."""

    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 10)
    ax.axis('off')

    # Colors
    blue = '#3498db'    # ThrRS theme
    green = '#27ae60'   # Accepted
    red = '#e74c3c'     # Rejected
    orange = '#f39c12'  # Trapped (SER)
    gray = '#7f8c8d'    # Zinc
    dark = '#2c3e50'    # Text

    # ===== TITLE =====
    ax.text(6, 9.5, 'ZINC FILTER MECHANISM', fontsize=22, fontweight='bold',
            ha='center', color=dark)
    ax.text(6, 9.0, 'ThrRS Structural Discrimination via Metal Coordination', fontsize=14,
            ha='center', color=gray, style='italic')

    # ===== ACTIVE SITE POCKET =====
    pocket = FancyBboxPatch((3, 3), 6, 5, boxstyle="round,pad=0.1",
                             facecolor='#ebf5fb', edgecolor=blue, linewidth=3)
    ax.add_patch(pocket)
    ax.text(6, 7.7, 'ACTIVE SITE POCKET', fontsize=12, fontweight='bold',
            ha='center', color=blue)

    # Zinc ion in center
    zn_circle = plt.Circle((6, 5.5), 0.5, color=gray, alpha=0.9)
    ax.add_patch(zn_circle)
    ax.text(6, 5.5, 'Zn²⁺', fontsize=14, fontweight='bold',
            ha='center', va='center', color='white')
    ax.text(6, 4.8, 'FILTER', fontsize=10, ha='center', color=gray, fontweight='bold')

    # Coordination bonds (dashed lines)
    ax.plot([6, 6], [6.0, 6.8], 'k--', lw=2, alpha=0.5)
    ax.plot([5.3, 6], [5.5, 5.5], 'k--', lw=2, alpha=0.5)
    ax.plot([6.7, 6], [5.5, 5.5], 'k--', lw=2, alpha=0.5)
    ax.text(6.1, 6.5, 'coordination', fontsize=8, color=gray, rotation=90)

    # ===== LEFT SIDE: ACCEPTED =====
    ax.text(1.5, 7.5, 'ACCEPTED ✓', fontsize=14, fontweight='bold', color=green)
    ax.text(1.5, 7.0, '(hydroxyl groups coordinate Zn)', fontsize=9, color=gray)

    # THR
    thr_box = FancyBboxPatch((0.5, 5.8), 2, 1, boxstyle="round,pad=0.02",
                              facecolor=green, edgecolor='white', linewidth=2)
    ax.add_patch(thr_box)
    ax.text(1.5, 6.5, 'THR', fontsize=12, fontweight='bold',
            ha='center', va='center', color='white')
    ax.text(1.5, 6.1, '-OH bidentate', fontsize=9, ha='center', color='white')
    ax.text(1.5, 5.5, 'ipTM: 0.97', fontsize=10, ha='center', color=green, fontweight='bold')

    # SER (the trap!)
    ser_box = FancyBboxPatch((0.5, 3.8), 2, 1, boxstyle="round,pad=0.02",
                              facecolor=orange, edgecolor='white', linewidth=2)
    ax.add_patch(ser_box)
    ax.text(1.5, 4.5, 'SER', fontsize=12, fontweight='bold',
            ha='center', va='center', color='white')
    ax.text(1.5, 4.1, '-OH bidentate', fontsize=9, ha='center', color='white')
    ax.text(1.5, 3.5, 'ipTM: 0.95 ⚠️', fontsize=10, ha='center', color=orange, fontweight='bold')
    ax.text(1.5, 3.1, 'TRAPPED!', fontsize=9, ha='center', color=orange, fontweight='bold')

    # ===== RIGHT SIDE: REJECTED =====
    ax.text(10.5, 7.5, 'REJECTED ✗', fontsize=14, fontweight='bold', color=red)
    ax.text(10.5, 7.0, '(hydrophobics cannot coordinate)', fontsize=9, color=gray)

    # ILE
    ile_box = FancyBboxPatch((9.5, 5.8), 2, 1, boxstyle="round,pad=0.02",
                              facecolor=red, edgecolor='white', linewidth=2, alpha=0.7)
    ax.add_patch(ile_box)
    ax.text(10.5, 6.5, 'ILE', fontsize=12, fontweight='bold',
            ha='center', va='center', color='white')
    ax.text(10.5, 6.1, '-CH₃ no coord', fontsize=9, ha='center', color='white')
    ax.text(10.5, 5.5, 'ipTM: 0.83', fontsize=10, ha='center', color=red, fontweight='bold')

    # VAL
    val_box = FancyBboxPatch((9.5, 4.3), 2, 1, boxstyle="round,pad=0.02",
                              facecolor=red, edgecolor='white', linewidth=2, alpha=0.7)
    ax.add_patch(val_box)
    ax.text(10.5, 5.0, 'VAL', fontsize=12, fontweight='bold',
            ha='center', va='center', color='white')
    ax.text(10.5, 4.6, '-CH₃ no coord', fontsize=9, ha='center', color='white')
    ax.text(10.5, 4.0, 'ipTM: 0.90', fontsize=10, ha='center', color=red, fontweight='bold')

    # ===== ARROWS =====
    # THR to Zn
    ax.annotate('', xy=(5.5, 5.7), xytext=(2.5, 6.3),
                arrowprops=dict(arrowstyle='->', color=green, lw=2))
    # SER to Zn
    ax.annotate('', xy=(5.5, 5.3), xytext=(2.5, 4.3),
                arrowprops=dict(arrowstyle='->', color=orange, lw=2))
    # ILE blocked
    ax.annotate('', xy=(9.5, 6.3), xytext=(6.7, 5.7),
                arrowprops=dict(arrowstyle='-|>', color=red, lw=2,
                               connectionstyle='arc3,rad=-0.2'))
    ax.text(8.2, 6.3, '✗', fontsize=20, color=red)

    # ===== KEY MESSAGE BOX =====
    msg_box = FancyBboxPatch((2.5, 0.5), 7, 1.5, boxstyle="round,pad=0.05",
                              facecolor='#fef9e7', edgecolor=orange, linewidth=2)
    ax.add_patch(msg_box)
    ax.text(6, 1.7, 'THE ZINC TRAP', fontsize=14, fontweight='bold',
            ha='center', color=orange)
    ax.text(6, 1.2, 'Zn filter rejects hydrophobics (ILE, VAL) ✓', fontsize=11,
            ha='center', color=dark)
    ax.text(6, 0.8, 'But SER also coordinates Zn → EDITING REQUIRED', fontsize=11,
            ha='center', color=orange, fontweight='bold')

    plt.tight_layout()
    return fig


if __name__ == '__main__':
    # Create Double Sieve figure
    fig1 = create_double_sieve_figure()
    fig1.savefig('fig2a_double_sieve_mechanism.png', dpi=300, bbox_inches='tight',
                 facecolor='white', edgecolor='none')
    fig1.savefig('fig2a_double_sieve_mechanism.pdf', bbox_inches='tight',
                 facecolor='white', edgecolor='none')
    print("✓ Created: fig2a_double_sieve_mechanism.png/pdf")

    # Create Zinc Filter figure
    fig2 = create_zinc_filter_figure()
    fig2.savefig('fig4a_zinc_filter_mechanism.png', dpi=300, bbox_inches='tight',
                 facecolor='white', edgecolor='none')
    fig2.savefig('fig4a_zinc_filter_mechanism.pdf', bbox_inches='tight',
                 facecolor='white', edgecolor='none')
    print("✓ Created: fig4a_zinc_filter_mechanism.png/pdf")

    plt.show()
