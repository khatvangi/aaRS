#!/usr/bin/env python3
"""
Figure 4: H-bond evolution - Steric permissiveness and evolutionary pruning
"""
import matplotlib.pyplot as plt
import numpy as np

# Data from H-bond analysis
data = {
    'LUCA (deep)': {'Pro': 227, 'Thr': 218},
    'Modern E. coli': {'Pro': 42, 'Thr': 39},
}

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel A: Evolution of interaction networks
ax = axes[0]
enzymes = list(data.keys())
pro_bonds = [data[e]['Pro'] for e in enzymes]
thr_bonds = [data[e]['Thr'] for e in enzymes]

x = np.arange(len(enzymes))
width = 0.35
bars1 = ax.bar(x - width/2, pro_bonds, width, label='Proline', color='#4CAF50', alpha=0.85, edgecolor='black', linewidth=2)
bars2 = ax.bar(x + width/2, thr_bonds, width, label='Threonine', color='#FF9800', alpha=0.85, edgecolor='black', linewidth=2)

for bar in bars1 + bars2:
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + 5, f'{int(height)}',
            ha='center', va='bottom', fontsize=12, fontweight='bold')

ax.set_ylabel('Total H-bonds', fontsize=14, fontweight='bold')
ax.set_xlabel('Enzyme', fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(enzymes, fontsize=12, fontweight='bold')
ax.set_ylim(0, 250)
ax.legend(fontsize=12, frameon=True, shadow=True)
ax.grid(axis='y', alpha=0.3)
ax.set_title('A. Evolution of Interaction Networks\nReduction from LUCA (~223) to Modern (~44)',
             fontsize=13, fontweight='bold', pad=10)

# Panel B: No compensation
ax = axes[1]
luca_pro = data['LUCA (deep)']['Pro']
luca_thr = data['LUCA (deep)']['Thr']
delta = luca_thr - luca_pro

categories = ['LUCA ProRS\n+ Proline', 'LUCA ProRS\n+ Threonine']
values = [luca_pro, luca_thr]
colors = ['#4CAF50', '#FF9800']

bars = ax.bar(categories, values, color=colors, alpha=0.85, edgecolor='black', linewidth=2)
for bar, val in zip(bars, values):
    ax.text(bar.get_x() + bar.get_width()/2., val + 5, f'{int(val)}',
            ha='center', va='bottom', fontsize=12, fontweight='bold')

# Annotate delta
ax.annotate('', xy=(1, luca_thr), xytext=(1, luca_pro),
            arrowprops=dict(arrowstyle='<->', color='red', lw=2))
ax.text(1.15, (luca_pro + luca_thr)/2, f'Δ = {delta}\nbonds',
        fontsize=11, fontweight='bold', color='red')

ax.set_ylabel('H-bonds', fontsize=14, fontweight='bold')
ax.set_ylim(0, 250)
ax.grid(axis='y', alpha=0.3)
ax.set_title('B. No Compensation\nNet Loss (Δ = -9 bonds) for Threonine',
             fontsize=13, fontweight='bold', pad=10)

# Panel C: Interpretation text
ax = axes[2]
ax.axis('off')

text = """
STERIC PERMISSIVENESS:

LUCA ProRS (Ancestral):
  • Proline: 227 H-bonds
  • Threonine: 218 H-bonds (-9)
  • Total ~223 bonds average

Modern E. coli ProRS:
  • Proline: 42 H-bonds
  • Threonine: 39 H-bonds (-3)
  • Total ~41 bonds average

EVOLUTIONARY PRUNING:
  → 81% reduction (223 → 41)
  → Despite fewer contacts,
    promiscuity MAINTAINED

INTERPRETATION:
LUCA had promiscuous binding
through STERIC PERMISSIVENESS
(large pocket fits both Pro/Thr)
rather than specific chemical
engagement.

Modern enzymes PRUNED contacts
while maintaining promiscuity
via POCKET EXPANSION (2.5×).

Thr β-hydroxyl forms NO specific
compensating interactions (0 bonds).
"""

ax.text(0.1, 0.95, text, transform=ax.transAxes,
        fontsize=10, verticalalignment='top', family='monospace',
        bbox=dict(boxstyle='round,pad=1', facecolor='lightyellow',
                 edgecolor='orange', linewidth=2, alpha=0.9))

ax.set_title('C. Structural Basis of Permissiveness',
             fontsize=13, fontweight='bold', pad=10)

fig.suptitle('Figure 4: Steric Permissiveness and Evolutionary Pruning',
             fontsize=16, fontweight='bold')
plt.tight_layout()

# Save
for fmt in ['png', 'pdf', 'svg']:
    plt.savefig(f'manuscript_figures/Figure4_Hbond_Evolution.{fmt}',
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(f'final_figures/Figure4_Hbond_Evolution.{fmt}',
                dpi=300, bbox_inches='tight', facecolor='white')

print("✓ Figure 4 created: manuscript_figures/Figure4_Hbond_Evolution.pdf")
