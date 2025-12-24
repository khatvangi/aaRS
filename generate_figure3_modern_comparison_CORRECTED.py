#!/usr/bin/env python3
"""
Figure 3: Modern vs Ancestral Comparison (CORRECTED DATA)

IMPORTANT: This replaces the previous version with CORRECT ipTM values
extracted directly from raw JSON files.

KEY FINDING: ThrRS evolved strong specificity (47% preference) while
ProRS maintained ancestral promiscuity (9% preference).
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# CORRECTED DATA - directly from JSON files
data = {
    'LUCA ProRS': {
        'Pro': 0.78,    # cognate
        'Thr': 0.70,    # non-cognate (11% lower)
        'pLDDT': 63.1,
        'preference': 11  # percent
    },
    'LUCA ThrRS': {
        'Pro': 0.88,    # non-cognate
        'Thr': 0.89,    # cognate (1% higher - truly promiscuous)
        'pLDDT': 62.8,
        'preference': 1   # percent
    },
    'Modern E. coli ProRS': {
        'Pro': 0.95,    # cognate
        'Thr': 0.87,    # non-cognate (9% lower)
        'pLDDT': 94.1,
        'preference': 9   # percent
    },
    'Modern ThrRS': {
        'Pro': 0.57,    # non-cognate
        'Thr': 0.84,    # cognate (47% higher - STRONG specificity)
        'pLDDT': 92.3,
        'preference': 47  # percent
    },
}

fig = plt.figure(figsize=(18, 10))
gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.3,
                      left=0.08, right=0.96, top=0.92, bottom=0.08)

# Colors
color_pro = '#4CAF50'  # Green for Pro
color_thr = '#FF9800'  # Orange for Thr

# =====================================================================
# Panel A: ipTM Comparison - SHOWS DIFFERENTIAL SPECIFICITY
# =====================================================================
ax1 = fig.add_subplot(gs[0, 0])

enzymes = list(data.keys())
pro_scores = [data[e]['Pro'] for e in enzymes]
thr_scores = [data[e]['Thr'] for e in enzymes]

x = np.arange(len(enzymes))
width = 0.35

bars1 = ax1.bar(x - width/2, pro_scores, width, label='Proline',
                color=color_pro, alpha=0.85, edgecolor='black', linewidth=2)
bars2 = ax1.bar(x + width/2, thr_scores, width, label='Threonine',
                color=color_thr, alpha=0.85, edgecolor='black', linewidth=2)

# Add value labels
for bar in bars1 + bars2:
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height + 0.02,
            f'{height:.2f}', ha='center', va='bottom',
            fontsize=11, fontweight='bold')

# Highlight key findings with arrows
# ThrRS evolved strong specificity
ax1.annotate('47% preference\nevolved!',
             xy=(3 + width/2, data['Modern ThrRS']['Thr']),
             xytext=(3.5, 0.75),
             fontsize=10, fontweight='bold', color='red',
             arrowprops=dict(arrowstyle='->', color='red', lw=2),
             bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow',
                      edgecolor='red', linewidth=2, alpha=0.9))

# ProRS maintained promiscuity
ax1.annotate('Maintained\npromiscuity',
             xy=(2, (data['Modern E. coli ProRS']['Pro'] +
                     data['Modern E. coli ProRS']['Thr'])/2),
             xytext=(1.5, 0.50),
             fontsize=10, fontweight='bold', color='blue',
             arrowprops=dict(arrowstyle='->', color='blue', lw=2),
             bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue',
                      edgecolor='blue', linewidth=2, alpha=0.9))

ax1.set_ylabel('Pocket ipTM Score', fontsize=14, fontweight='bold')
ax1.set_xlabel('Enzyme', fontsize=14, fontweight='bold')
ax1.set_xticks(x)
ax1.set_xticklabels(enzymes, fontsize=11, fontweight='bold', rotation=20, ha='right')
ax1.set_ylim(0, 1.05)
ax1.legend(fontsize=12, loc='upper left', frameon=True, shadow=True)
ax1.grid(axis='y', alpha=0.3, linestyle='--')
ax1.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, linewidth=1,
            label='Confidence threshold')
ax1.set_title('A. Pocket ipTM: Asymmetric Evolution of Specificity',
             fontsize=13, fontweight='bold', pad=10)

# =====================================================================
# Panel B: Substrate Preference (% difference)
# =====================================================================
ax2 = fig.add_subplot(gs[0, 1])

preferences = [data[e]['preference'] for e in enzymes]
colors_pref = ['#FFB74D', '#81C784', '#FFB74D', '#F44336']  # Yellow/green/yellow/red

bars = ax2.bar(enzymes, preferences, color=colors_pref, alpha=0.85,
               edgecolor='black', linewidth=2)

for bar, pref in zip(bars, preferences):
    ax2.text(bar.get_x() + bar.get_width()/2., pref + 2,
            f'{pref}%', ha='center', va='bottom',
            fontsize=13, fontweight='bold')

# Highlight zones
ax2.axhspan(0, 10, alpha=0.1, color='green', label='Promiscuous (<10%)')
ax2.axhspan(10, 50, alpha=0.1, color='red', label='Specific (>10%)')

ax2.set_ylabel('Cognate Preference (%)', fontsize=14, fontweight='bold')
ax2.set_xlabel('Enzyme', fontsize=14, fontweight='bold')
ax2.set_xticklabels(enzymes, fontsize=11, fontweight='bold', rotation=20, ha='right')
ax2.set_ylim(0, 55)
ax2.grid(axis='y', alpha=0.3, linestyle='--')
ax2.legend(fontsize=10, loc='upper left')
ax2.set_title('B. Substrate Discrimination: ThrRS Evolved 47% Preference',
             fontsize=13, fontweight='bold', pad=10)

# =====================================================================
# Panel C: pLDDT (Overall Structure Quality)
# =====================================================================
ax3 = fig.add_subplot(gs[0, 2])

plddt_scores = [data[e]['pLDDT'] for e in enzymes]
colors_plddt = ['#FF9800', '#FF9800', '#4CAF50', '#4CAF50']  # Orange for LUCA, Green for Modern

bars = ax3.bar(enzymes, plddt_scores, color=colors_plddt, alpha=0.85,
               edgecolor='black', linewidth=2)

for bar, score in zip(bars, plddt_scores):
    ax3.text(bar.get_x() + bar.get_width()/2., score + 2,
            f'{score:.1f}', ha='center', va='bottom',
            fontsize=12, fontweight='bold')

ax3.axhline(y=70, color='red', linestyle='--', linewidth=2,
            label='High confidence threshold')
ax3.set_ylabel('Average pLDDT', fontsize=14, fontweight='bold')
ax3.set_xlabel('Enzyme', fontsize=14, fontweight='bold')
ax3.set_xticklabels(enzymes, fontsize=11, fontweight='bold', rotation=20, ha='right')
ax3.set_ylim(0, 105)
ax3.grid(axis='y', alpha=0.3, linestyle='--')
ax3.legend(fontsize=10)
ax3.set_title('C. Structure Quality: Modern >> LUCA',
             fontsize=13, fontweight='bold', pad=10)

# =====================================================================
# Panel D: Evolutionary Trajectories (Bottom Left - Large)
# =====================================================================
ax4 = fig.add_subplot(gs[1, :2])

# Plot trajectories
# ProRS trajectory
pro_pref = [data['LUCA ProRS']['preference'], data['Modern E. coli ProRS']['preference']]
ax4.plot([0, 1], pro_pref, 'o-', color='#4CAF50', linewidth=3,
         markersize=15, label='ProRS (maintained promiscuity)',
         markeredgecolor='black', markeredgewidth=2)

# ThrRS trajectory
thr_pref = [data['LUCA ThrRS']['preference'], data['Modern ThrRS']['preference']]
ax4.plot([0, 1], thr_pref, 'o-', color='#F44336', linewidth=3,
         markersize=15, label='ThrRS (evolved specificity)',
         markeredgecolor='black', markeredgewidth=2)

# Annotate values
ax4.text(-0.05, pro_pref[0], f'{pro_pref[0]}%', fontsize=13, fontweight='bold', ha='right', va='center')
ax4.text(1.05, pro_pref[1], f'{pro_pref[1]}%', fontsize=13, fontweight='bold', ha='left', va='center')
ax4.text(-0.05, thr_pref[0], f'{thr_pref[0]}%', fontsize=13, fontweight='bold', ha='right', va='center')
ax4.text(1.05, thr_pref[1], f'{thr_pref[1]}%', fontsize=13, fontweight='bold', ha='left', va='center')

# Shade promiscuous vs specific zones
ax4.axhspan(0, 10, alpha=0.15, color='green')
ax4.axhspan(10, 50, alpha=0.15, color='red')
ax4.text(0.5, 5, 'PROMISCUOUS', ha='center', fontsize=12,
         fontweight='bold', color='darkgreen', alpha=0.7)
ax4.text(0.5, 30, 'SPECIFIC', ha='center', fontsize=12,
         fontweight='bold', color='darkred', alpha=0.7)

# Add arrow showing ThrRS evolution
ax4.annotate('', xy=(1, thr_pref[1]), xytext=(0, thr_pref[0]),
            arrowprops=dict(arrowstyle='->', color='red', lw=4, alpha=0.5))

ax4.set_xticks([0, 1])
ax4.set_xticklabels(['LUCA\n(~3.5 Ga)', 'Modern\n(Present)'],
                    fontsize=13, fontweight='bold')
ax4.set_ylabel('Cognate Preference (%)', fontsize=14, fontweight='bold')
ax4.set_ylim(0, 52)
ax4.set_xlim(-0.15, 1.15)
ax4.legend(fontsize=12, loc='upper left', frameon=True, shadow=True)
ax4.grid(axis='y', alpha=0.3, linestyle='--')
ax4.set_title('D. Evolutionary Trajectories: Asymmetric Refinement',
             fontsize=13, fontweight='bold', pad=10)

# =====================================================================
# Panel E: Interpretation Summary (Bottom Right)
# =====================================================================
ax5 = fig.add_subplot(gs[1, 2])
ax5.axis('off')

interpretation = """
ASYMMETRIC EVOLUTION:

LUCA (Ancestral):
  • ProRS: 11% preference (mild)
  • ThrRS: 1% preference (promiscuous)
  → Both relatively promiscuous

Modern:
  • ProRS: 9% preference (MAINTAINED)
  • ThrRS: 47% preference (EVOLVED!)
  → Divergent fidelity strategies

KEY FINDINGS:

1. ProRS MAINTAINED promiscuity
   → Likely due to editing domain
     compensation (3.2× preference
     for Thr-tRNA^Pro)

2. ThrRS EVOLVED specificity
   → 47-fold increase in preference
   → May lack effective editing

3. Differential selection pressure
   → Thr→Pro errors more deleterious?
   → ProRS can tolerate catalytic
     promiscuity with editing backup

BIOLOGICAL SIGNIFICANCE:
ProRS uses TWO-STAGE fidelity:
  1) Catalytic promiscuity (9%)
  2) Editing specificity (320%)
  = Overall high fidelity

ThrRS evolved SINGLE-STAGE:
  - High catalytic specificity (47%)
  - Less editing dependency
"""

ax5.text(0.05, 0.95, interpretation, transform=ax5.transAxes,
        fontsize=10, verticalalignment='top', family='monospace',
        bbox=dict(boxstyle='round,pad=1', facecolor='lightyellow',
                 edgecolor='orange', linewidth=2, alpha=0.9))

ax5.set_title('E. Biological Interpretation',
             fontsize=13, fontweight='bold', pad=10)

# Overall title
fig.suptitle('Figure 3: Asymmetric Evolution of Substrate Specificity (CORRECTED)',
            fontsize=16, fontweight='bold', y=0.96)

plt.tight_layout()

# Save
for output_dir in ['manuscript_figures', 'final_figures', 'figures']:
    for fmt in ['png', 'pdf', 'svg']:
        plt.savefig(f'{output_dir}/Figure3_Modern_Ancestral_CORRECTED.{fmt}',
                    dpi=300, bbox_inches='tight', facecolor='white')

print("✓ Figure 3 (CORRECTED) created")
print("\nKEY FINDING: Asymmetric evolution of specificity")
print("  - ProRS maintained promiscuity: LUCA 11% → Modern 9%")
print("  - ThrRS evolved specificity: LUCA 1% → Modern 47%")
print("\nThis explains why ProRS can tolerate catalytic promiscuity:")
print("  → Editing domain provides backup quality control")
