#!/usr/bin/env python3
"""
Figure 3: The ancestral editing domain displays inverted substrate preference
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np

# Data from editing domain AF3 predictions
data = {
    'Pro-tRNA^Pro (cognate)': 0.14,
    'Thr-tRNA^Pro (non-cognate)': 0.45,
}

fig = plt.figure(figsize=(14, 6))
gs = fig.add_gridspec(1, 2, wspace=0.3, left=0.08, right=0.96, top=0.88, bottom=0.12)

# =====================================================================
# Panel A: Domain Architecture
# =====================================================================
ax1 = fig.add_subplot(gs[0, 0])
ax1.set_xlim(0, 10)
ax1.set_ylim(0, 5)
ax1.axis('off')

# Draw ProRS protein
# Catalytic domain (left, larger)
catalytic = Rectangle((0.5, 1.5), 4, 2, facecolor='#1976D2', alpha=0.7,
                       edgecolor='black', linewidth=2, label='Catalytic domain')
ax1.add_patch(catalytic)
ax1.text(2.5, 2.5, 'Catalytic\nDomain', ha='center', va='center',
         fontsize=12, fontweight='bold', color='white')

# Editing domain (right, smaller, highlighted)
editing = Rectangle((5, 1.8), 2.5, 1.4, facecolor='#F57C00', alpha=0.9,
                     edgecolor='red', linewidth=3, label='Editing domain')
ax1.add_patch(editing)
ax1.text(6.25, 2.5, 'Editing\nDomain', ha='center', va='center',
         fontsize=11, fontweight='bold', color='white')

# Annotation arrow pointing to editing domain
ax1.annotate('3.2× prefers\nThr over Pro', xy=(6.25, 3.5), xytext=(7.5, 4.2),
             fontsize=11, fontweight='bold', color='red',
             arrowprops=dict(arrowstyle='->', color='red', lw=2),
             bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow',
                      edgecolor='red', linewidth=2, alpha=0.9))

# Residue numbers
ax1.text(0.5, 1.2, '1', ha='center', fontsize=10, fontweight='bold')
ax1.text(4.5, 1.2, '~1500', ha='center', fontsize=10, fontweight='bold')
ax1.text(7.5, 1.2, '2037', ha='center', fontsize=10, fontweight='bold')

# Title
ax1.text(4, 4.5, 'LUCA ProRS Domain Architecture', ha='center',
         fontsize=14, fontweight='bold')

ax1.set_title('A. Domain Architecture', fontsize=13, fontweight='bold', pad=15)

# =====================================================================
# Panel B: Editing Domain ipTM - Inverted Preference
# =====================================================================
ax2 = fig.add_subplot(gs[0, 1])

substrates = list(data.keys())
iptm_values = list(data.values())
colors = ['#4CAF50', '#FF9800']  # Green for Pro, Orange for Thr

bars = ax2.bar(substrates, iptm_values, color=colors, alpha=0.85,
               edgecolor='black', linewidth=2, width=0.6)

# Add value labels
for bar, val in zip(bars, iptm_values):
    ax2.text(bar.get_x() + bar.get_width()/2., val + 0.02,
            f'{val:.2f}', ha='center', va='bottom',
            fontsize=13, fontweight='bold')

# Show fold preference
fold_pref = data['Thr-tRNA^Pro (non-cognate)'] / data['Pro-tRNA^Pro (cognate)']
ax2.text(0.5, 0.35, f'{fold_pref:.1f}× higher', ha='center', va='center',
         fontsize=14, fontweight='bold', color='red',
         bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow',
                  edgecolor='red', linewidth=2))

# Arrow showing preference
ax2.annotate('', xy=(1, data['Thr-tRNA^Pro (non-cognate)']),
             xytext=(0, data['Pro-tRNA^Pro (cognate)']),
             arrowprops=dict(arrowstyle='->', color='red', lw=3))

# Add interpretation box
interpretation = """
INVERTED PREFERENCE:

Editing domain PREFERS the
WRONG substrate (Thr-tRNA^Pro)
over the correct one (Pro-tRNA^Pro).

This is OPPOSITE to the catalytic
domain, which binds both equally
(promiscuity).

BIOLOGICAL FUNCTION:
The editing domain's job is to
HYDROLYZE mis-acylated tRNA
(Thr-tRNA^Pro), so it SHOULD
preferentially bind the error
product, not the correct product.

This inverted preference validates
the editing domain's proofreading
function in ancestral ProRS.
"""

ax2.text(1.5, 0.25, interpretation, transform=ax2.transData,
         fontsize=9, verticalalignment='top', family='monospace',
         bbox=dict(boxstyle='round,pad=0.8', facecolor='lightblue',
                  edgecolor='blue', linewidth=2, alpha=0.8))

ax2.set_ylabel('Pocket ipTM Score', fontsize=14, fontweight='bold')
ax2.set_xlabel('Substrate-tRNA Complex', fontsize=14, fontweight='bold')
ax2.set_ylim(0, 0.6)
ax2.set_xticklabels(substrates, fontsize=11, fontweight='bold', rotation=15, ha='right')
ax2.grid(axis='y', alpha=0.3, linestyle='--')
ax2.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, linewidth=1)
ax2.set_title('B. Editing Domain Shows 3.2× Preference\nfor Non-Cognate Substrate',
             fontsize=13, fontweight='bold', pad=15)

# Overall title
fig.suptitle('Figure 3: Ancestral Editing Domain Displays Inverted Substrate Preference',
            fontsize=16, fontweight='bold', y=0.96)

plt.tight_layout()

# Save
for output_dir in ['manuscript_figures', 'final_figures', 'figures']:
    for fmt in ['png', 'pdf', 'svg']:
        plt.savefig(f'{output_dir}/Figure3_Editing_Domain_Inverted.{fmt}',
                    dpi=300, bbox_inches='tight', facecolor='white')

print("✓ Figure 3 created: manuscript_figures/Figure3_Editing_Domain_Inverted.pdf")
print(f"\nKEY FINDING: Editing domain shows 3.2× preference for Thr-tRNA (0.45)")
print(f"             over Pro-tRNA (0.14) - INVERTED specificity!")
