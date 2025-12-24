#!/usr/bin/env python3
"""
Generate ipTM bar charts for ProRS manuscript
Figure 5: Binding Affinity Comparison
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import sys

# Read the master ipTM data
data_file = '/storage/kiran-stuff/aaRS/figures/table_master_iptm_data.csv'
df = pd.read_csv(data_file)

print("Loaded ipTM data:")
print(df.head())
print(f"\nTotal models: {len(df)}")

# Extract relevant data
# LUCA ProRS catalytic: deep_domain_pro (PRO), deep_domain_thr (THR)
# Modern ProRS catalytic: shallow_domain_pro (PRO), shallow_domain_thr (THR)
# LUCA Editing: deep_editing_pro (PRO), deep_editing_thr (THR)

data = {
    'LUCA ProRS': {
        'Pro': df[df['model'] == 'deep_domain_pro']['iptm'].values[0],
        'Thr': df[df['model'] == 'deep_domain_thr']['iptm'].values[0]
    },
    'Modern ProRS': {
        'Pro': df[df['model'] == 'shallow_domain_pro']['iptm'].values[0],
        'Thr': df[df['model'] == 'shallow_domain_thr']['iptm'].values[0]
    },
    'LUCA Editing': {
        'Pro': df[df['model'] == 'deep_editing_pro']['iptm'].values[0],
        'Thr': df[df['model'] == 'deep_editing_thr']['iptm'].values[0]
    }
}

print("\nData to plot:")
for key, vals in data.items():
    print(f"{key}: PRO={vals['Pro']:.2f}, THR={vals['Thr']:.2f}")
    cross_react = vals['Thr'] / vals['Pro'] * 100
    print(f"  Cross-reactivity: {cross_react:.1f}%")

# Create figure with 3 subplots
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Color scheme
color_pro = '#2E86AB'  # Blue for proline
color_thr = '#A23B72'  # Purple/magenta for threonine

for ax, (title, values) in zip(axes, data.items()):
    x = np.arange(2)
    bars = ax.bar(x, [values['Pro'], values['Thr']],
                  color=[color_pro, color_thr], width=0.6,
                  edgecolor='black', linewidth=1.5)

    # Axis settings
    ax.set_xticks(x)
    ax.set_xticklabels(['Proline\n(Cognate)', 'Threonine\n(Non-cognate)'],
                       fontsize=12, fontweight='bold')
    ax.set_ylabel('ipTM Score', fontsize=14, fontweight='bold')
    ax.set_title(title, fontsize=16, fontweight='bold', pad=15)
    ax.set_ylim(0, 1.0)

    # Add horizontal reference line at 0.5
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, linewidth=1.5)

    # Add value labels on top of bars
    for bar, val in zip(bars, [values['Pro'], values['Thr']]):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, height + 0.03,
                f'{val:.2f}', ha='center', va='bottom',
                fontsize=13, fontweight='bold')

    # Calculate and display cross-reactivity percentage
    cross_react = values['Thr'] / values['Pro'] * 100

    # Position the annotation box
    if title == 'LUCA Editing':
        # For editing domain, THR > PRO, so show it differently
        annotation = f'THR binds {cross_react:.0f}%\nbetter than PRO'
        y_pos = 0.85
        bgcolor = 'lightcoral'
    else:
        annotation = f'{cross_react:.1f}%\ncross-reactivity'
        y_pos = 0.90
        bgcolor = 'wheat'

    ax.text(0.5, y_pos, annotation,
            transform=ax.transAxes, ha='center', fontsize=11,
            bbox=dict(boxstyle='round,pad=0.5', facecolor=bgcolor,
                     alpha=0.7, edgecolor='black', linewidth=1.5),
            fontweight='bold')

    # Grid for readability
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

    # Spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)

# Add overall title
fig.suptitle('Substrate Binding Specificity in Ancestral vs Modern ProRS',
             fontsize=18, fontweight='bold', y=1.02)

# Adjust layout
plt.tight_layout()

# Save in multiple formats
output_base = '/storage/kiran-stuff/aaRS/figure5_iptm_bars'
plt.savefig(f'{output_base}.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(f'{output_base}.pdf', bbox_inches='tight', facecolor='white')
plt.savefig(f'{output_base}.svg', bbox_inches='tight', facecolor='white')

print(f"\nFigures saved:")
print(f"  {output_base}.png")
print(f"  {output_base}.pdf")
print(f"  {output_base}.svg")

# Create a second figure with more comprehensive data
# Include LUCA ThrRS and Modern ThrRS
fig2, axes2 = plt.subplots(2, 3, figsize=(15, 10))
axes2 = axes2.flatten()

# Extended dataset
extended_data = {
    'LUCA ProRS\nCatalytic': {
        'Pro': df[df['model'] == 'deep_domain_pro']['iptm'].values[0],
        'Thr': df[df['model'] == 'deep_domain_thr']['iptm'].values[0]
    },
    'LUCA ThrRS\nCatalytic': {
        'Pro': df[df['model'] == 'deep_thrrs_pro']['iptm'].values[0],
        'Thr': df[df['model'] == 'deep_thrrs_thr']['iptm'].values[0]
    },
    'LUCA ProRS\nEditing': {
        'Pro': df[df['model'] == 'deep_editing_pro']['iptm'].values[0],
        'Thr': df[df['model'] == 'deep_editing_thr']['iptm'].values[0]
    },
    'Modern ProRS\nCatalytic': {
        'Pro': df[df['model'] == 'shallow_domain_pro']['iptm'].values[0],
        'Thr': df[df['model'] == 'shallow_domain_thr']['iptm'].values[0]
    },
    'Modern ProRS\n(Human)': {
        'Pro': df[df['model'] == 'modern_prours_pro']['iptm'].values[0],
        'Thr': df[df['model'] == 'modern_prours_thr']['iptm'].values[0]
    },
    'Modern ThrRS\n(Human)': {
        'Pro': df[df['model'] == 'modern_thrrs_pro']['iptm'].values[0],
        'Thr': df[df['model'] == 'modern_thrrs_thr']['iptm'].values[0]
    }
}

for i, (ax, (title, values)) in enumerate(zip(axes2, extended_data.items())):
    x = np.arange(2)
    bars = ax.bar(x, [values['Pro'], values['Thr']],
                  color=[color_pro, color_thr], width=0.6,
                  edgecolor='black', linewidth=1.5)

    ax.set_xticks(x)
    ax.set_xticklabels(['PRO', 'THR'], fontsize=11, fontweight='bold')
    ax.set_ylabel('ipTM Score', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=13, fontweight='bold', pad=10)
    ax.set_ylim(0, 1.0)
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, linewidth=1)

    # Value labels
    for bar, val in zip(bars, [values['Pro'], values['Thr']]):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, height + 0.02,
                f'{val:.2f}', ha='center', va='bottom',
                fontsize=11, fontweight='bold')

    # Grid and spines
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)

fig2.suptitle('Comprehensive Substrate Binding Analysis',
              fontsize=18, fontweight='bold', y=0.995)

plt.tight_layout()

output_base2 = '/storage/kiran-stuff/aaRS/figure5_iptm_comprehensive'
plt.savefig(f'{output_base2}.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(f'{output_base2}.pdf', bbox_inches='tight', facecolor='white')

print(f"\nComprehensive figure saved:")
print(f"  {output_base2}.png")
print(f"  {output_base2}.pdf")

print("\nFigure 5 generation complete!")
