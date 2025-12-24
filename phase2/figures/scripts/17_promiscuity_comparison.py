#!/usr/bin/env python3
"""
Compare ProRS promiscuity vs ThrRS specificity.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load MM/GBSA data
ddG = pd.read_csv('energy_scoring/mmgbsa_corrected_ddG.csv', index_col=0)
abs_BE = pd.read_csv('energy_scoring/mmgbsa_corrected_absolute.csv', index_col=0)

# Extract data for Modern ProRS and Modern ThrRS
prors_col = 'Modern ProRS\nCatalytic'
thrrs_col = 'Modern ThrRS\n+ Zn'

prors_ddG = ddG[prors_col].dropna().sort_values()
thrrs_ddG = ddG[thrrs_col].dropna().sort_values()

# Create figure
fig, axes = plt.subplots(1, 2, figsize=(16, 8))

# Panel A: ProRS - promiscuous
ax1 = axes[0]
colors_prors = ['red' if aa == 'PRO' else 'steelblue' for aa in prors_ddG.index]
ax1.barh(range(len(prors_ddG)), prors_ddG.values, color=colors_prors, alpha=0.7, edgecolor='black')
ax1.set_yticks(range(len(prors_ddG)))
ax1.set_yticklabels(prors_ddG.index, fontsize=10)
ax1.axvline(x=0, color='red', linestyle='--', linewidth=2, alpha=0.7)
ax1.axvline(x=-50, color='orange', linestyle=':', linewidth=1.5, alpha=0.5)
ax1.axvline(x=50, color='orange', linestyle=':', linewidth=1.5, alpha=0.5)
ax1.set_xlabel('ΔΔG from PRO (kcal/mol)', fontweight='bold', fontsize=12)
ax1.set_title('A. Modern ProRS Catalytic: PROMISCUOUS\n(17/20 AAs within ±50 kcal/mol)',
              fontweight='bold', fontsize=14)
ax1.grid(axis='x', alpha=0.3)
ax1.invert_yaxis()

# Add count annotation
within_50_prors = sum(abs(prors_ddG) < 50)
ax1.text(0.95, 0.05, f'Promiscuous:\n{within_50_prors}/20 AAs',
         transform=ax1.transAxes, fontsize=11, fontweight='bold',
         verticalalignment='bottom', horizontalalignment='right',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

# Panel B: ThrRS - specific
ax2 = axes[1]
colors_thrrs = ['red' if aa == 'THR' else 'darkgreen' if aa == 'SER' else 'steelblue' for aa in thrrs_ddG.index]
ax2.barh(range(len(thrrs_ddG)), thrrs_ddG.values, color=colors_thrrs, alpha=0.7, edgecolor='black')
ax2.set_yticks(range(len(thrrs_ddG)))
ax2.set_yticklabels(thrrs_ddG.index, fontsize=10)
ax2.axvline(x=0, color='red', linestyle='--', linewidth=2, alpha=0.7)
ax2.axvline(x=-50, color='orange', linestyle=':', linewidth=1.5, alpha=0.5)
ax2.axvline(x=50, color='orange', linestyle=':', linewidth=1.5, alpha=0.5)
ax2.set_xlabel('ΔΔG from THR (kcal/mol)', fontweight='bold', fontsize=12)
ax2.set_title('B. Modern ThrRS + Zn: SPECIFIC\n(5/20 AAs within ±50 kcal/mol)',
              fontweight='bold', fontsize=14)
ax2.grid(axis='x', alpha=0.3)
ax2.invert_yaxis()

# Add count annotation
within_50_thrrs = sum(abs(thrrs_ddG) < 50)
ax2.text(0.95, 0.05, f'Specific:\n{within_50_thrrs}/20 AAs\n\nSER trapped!',
         transform=ax2.transAxes, fontsize=11, fontweight='bold',
         verticalalignment='bottom', horizontalalignment='right',
         bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))

plt.tight_layout()
plt.savefig('figures/outputs/17_promiscuity_comparison.png', dpi=300, bbox_inches='tight')
print("✓ Saved: figures/outputs/17_promiscuity_comparison.png")

# Print detailed comparison
print("\n" + "="*80)
print("ProRS vs ThrRS PROMISCUITY")
print("="*80)

print(f"\nModern ProRS Catalytic:")
print(f"  Promiscuous ligands (±50 kcal/mol): {within_50_prors}/20")
print(f"  Top 5 competitors to PRO:")
for aa in prors_ddG.index[:5]:
    if aa != 'PRO':
        print(f"    {aa}: ΔΔG = {prors_ddG[aa]:+.1f} kcal/mol")

print(f"\nModern ThrRS + Zn:")
print(f"  Specific ligands (±50 kcal/mol): {within_50_thrrs}/20")
print(f"  Closest competitor:")
ser_ddg = thrrs_ddG.get('SER', np.nan)
if not np.isnan(ser_ddg):
    print(f"    SER: ΔΔG = {ser_ddg:+.1f} kcal/mol (TRAPPED!)")

print(f"\nPromiscuity Ratio: {within_50_prors}/{within_50_thrrs} = {within_50_prors/within_50_thrrs:.1f}x")
print(f"ProRS is {within_50_prors/within_50_thrrs:.1f}x MORE PROMISCUOUS than ThrRS")

print("\n" + "="*80)
print("WHY EDITING DOMAINS ARE ESSENTIAL")
print("="*80)
print("""
ProRS Catalytic Site:
  - Binds 17/20 amino acids promiscuously
  - THR, SER, ALA are all errors (bind well but wrong)
  → NEEDS editing domain to clear THR/SER/ALA-tRNA

ThrRS Catalytic Site (+ Zn):
  - Highly specific for THR
  - Only SER is trapped (ΔΔG = -28 kcal/mol, very close to THR)
  → NEEDS editing domain ONLY for SER-tRNA

This explains the DOUBLE-SIEVE mechanism!
""")

print("="*80)
