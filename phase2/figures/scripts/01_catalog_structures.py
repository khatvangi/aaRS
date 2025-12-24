#!/usr/bin/env python3
"""
Catalog all AF3 structures and organize by figure requirements.
Creates a comprehensive inventory for figure generation.
"""

import pandas as pd
import json
from pathlib import Path
from collections import defaultdict

# Load master CSV
df = pd.read_csv('AF3_RESULTS_CORRECTED.csv')

# Filter high-quality predictions
df_clean = df[(df['ptm'] >= 0.40) & (df['has_rna'] == False)].copy()

# Extract primary ligand
df_clean['primary_ligand'] = df_clean['ligands'].str.split(',').str[0]

# Categorize by enzyme and era
def categorize_run(row):
    """Categorize each run by enzyme, era, domain, and has_zinc."""
    name = row['job_name'].lower()

    # Enzyme
    if 'thrrs' in name or 'thr' in name:
        enzyme = 'ThrRS'
    elif 'prors' in name or 'pro' in name:
        enzyme = 'ProRS'
    else:
        enzyme = 'Unknown'

    # Era
    if 'anc' in name or 'ancestral' in name or 'deep' in name or 'luca' in name:
        era = 'ancestral'
    elif 'modern' in name or 'ecoli' in name:
        era = 'modern'
    else:
        era = 'unknown'

    # Domain
    if 'edit' in name:
        domain = 'editing'
    elif 'cat' in name or 'domain' in name:
        domain = 'catalytic'
    elif 'full' in name:
        domain = 'full'
    else:
        domain = 'catalytic'  # default

    # Zinc
    has_zinc = '_zn' in name or 'zinc' in name or row['Zn_iptm'] > 0

    # Competition
    is_competition = 'competition' in name

    return pd.Series({
        'enzyme': enzyme,
        'era': era,
        'domain': domain,
        'has_zinc': has_zinc,
        'is_competition': is_competition
    })

# Apply categorization
df_categorized = df_clean.join(df_clean.apply(categorize_run, axis=1))

# Save categorized data
df_categorized.to_csv('figures/data/categorized_predictions.csv', index=False)

print("="*80)
print("STRUCTURE CATALOG FOR FIGURE GENERATION")
print("="*80)

# ============================================================================
# FIGURE 1: ANCESTRAL BUCKET PROBLEM
# ============================================================================
print("\n" + "="*80)
print("FIGURE 1: ANCESTRAL BUCKET PROBLEM")
print("="*80)

# Fig 1C: Ancestral ThrRS (no Zn)
anc_thrrs_no_zn = df_categorized[
    (df_categorized['enzyme'] == 'ThrRS') &
    (df_categorized['era'] == 'ancestral') &
    (df_categorized['domain'] == 'catalytic') &
    (df_categorized['has_zinc'] == False) &
    (df_categorized['protein_len'] == 278)
].copy()

print(f"\nFig 1C - Ancestral ThrRS (no Zn): {len(anc_thrrs_no_zn)} predictions")
if not anc_thrrs_no_zn.empty:
    anc_thrrs_no_zn = anc_thrrs_no_zn.sort_values('AA_iptm', ascending=False)
    anc_thrrs_no_zn['rank'] = range(1, len(anc_thrrs_no_zn) + 1)

    print("\nTop 10:")
    for _, row in anc_thrrs_no_zn.head(10).iterrows():
        print(f"  {row['rank']:2d}. {row['primary_ligand']:3s}: {row['AA_iptm']:.3f}")

    thr_rank = anc_thrrs_no_zn[anc_thrrs_no_zn['primary_ligand'] == 'THR']['rank'].values
    if len(thr_rank) > 0:
        print(f"\n  → THR ranks {thr_rank[0]} out of {len(anc_thrrs_no_zn)}")

    anc_thrrs_no_zn.to_csv('figures/data/fig1c_anc_thrrs_no_zn.csv', index=False)

# Fig 1D: Ancestral ProRS
anc_prors = df_categorized[
    (df_categorized['enzyme'] == 'ProRS') &
    (df_categorized['era'] == 'ancestral') &
    (df_categorized['domain'] == 'catalytic') &
    (df_categorized['protein_len'] == 500)
].copy()

print(f"\nFig 1D - Ancestral ProRS: {len(anc_prors)} predictions")
if not anc_prors.empty:
    anc_prors = anc_prors.sort_values('AA_iptm', ascending=False)
    anc_prors['rank'] = range(1, len(anc_prors) + 1)

    print("\nTop 10:")
    for _, row in anc_prors.head(10).iterrows():
        print(f"  {row['rank']:2d}. {row['primary_ligand']:3s}: {row['AA_iptm']:.3f}")

    pro_rank = anc_prors[anc_prors['primary_ligand'] == 'PRO']['rank'].values
    if len(pro_rank) > 0:
        print(f"\n  → PRO ranks {pro_rank[0]} out of {len(anc_prors)}")

    anc_prors.to_csv('figures/data/fig1d_anc_prors.csv', index=False)

# ============================================================================
# FIGURE 2: ProRS DOUBLE SIEVE
# ============================================================================
print("\n" + "="*80)
print("FIGURE 2: ProRS DOUBLE SIEVE")
print("="*80)

# Fig 2B: Modern ProRS catalytic
mod_prors = df_categorized[
    (df_categorized['enzyme'] == 'ProRS') &
    (df_categorized['era'] == 'modern') &
    (df_categorized['protein_len'] == 572)
].copy()

print(f"\nFig 2B - Modern ProRS catalytic: {len(mod_prors)} predictions")
if not mod_prors.empty:
    mod_prors = mod_prors.sort_values('AA_iptm', ascending=False)
    mod_prors['rank'] = range(1, len(mod_prors) + 1)

    pro_score = mod_prors[mod_prors['primary_ligand'] == 'PRO']['AA_iptm'].values[0]
    mod_prors['pct_of_cognate'] = (mod_prors['AA_iptm'] / pro_score) * 100

    print("\nTop 10:")
    for _, row in mod_prors.head(10).iterrows():
        print(f"  {row['rank']:2d}. {row['primary_ligand']:3s}: {row['AA_iptm']:.3f} ({row['pct_of_cognate']:.1f}%)")

    mod_prors.to_csv('figures/data/fig2b_mod_prors_catalytic.csv', index=False)

# Fig 2C: ProRS editing domain
prors_edit = df_categorized[
    (df_categorized['enzyme'] == 'ProRS') &
    (df_categorized['domain'] == 'editing') &
    (df_categorized['protein_len'] == 300)
].copy()

print(f"\nFig 2C - ProRS editing domain: {len(prors_edit)} predictions")
if not prors_edit.empty:
    prors_edit = prors_edit.sort_values('AA_iptm', ascending=False)
    prors_edit['rank'] = range(1, len(prors_edit) + 1)

    print("\nTop 10:")
    for _, row in prors_edit.head(10).iterrows():
        print(f"  {row['rank']:2d}. {row['primary_ligand']:3s}: {row['AA_iptm']:.3f}")

    prors_edit.to_csv('figures/data/fig2c_prors_editing.csv', index=False)

# ============================================================================
# FIGURE 3: ZINC DISCONNECT
# ============================================================================
print("\n" + "="*80)
print("FIGURE 3: ZINC DISCONNECT")
print("="*80)

# Competition experiments
competitions = df_categorized[df_categorized['is_competition'] == True].copy()

print(f"\nCompetition experiments: {len(competitions)}")
for _, row in competitions.iterrows():
    print(f"  {row['job_name']}")
    print(f"    Ligands: {row['ligands']}, AA_iptm: {row['AA_iptm']:.3f}")

competitions.to_csv('figures/data/fig3_competitions.csv', index=False)

# ============================================================================
# FIGURE 4: ZINC FILTER
# ============================================================================
print("\n" + "="*80)
print("FIGURE 4: ZINC FILTER")
print("="*80)

# Modern ThrRS with Zn (all 20 AAs)
mod_thrrs_zn = df_categorized[
    (df_categorized['enzyme'] == 'ThrRS') &
    (df_categorized['era'] == 'modern') &
    (df_categorized['has_zinc'] == True) &
    (df_categorized['is_competition'] == False) &
    (df_categorized['protein_len'] == 401)
].copy()

print(f"\nFig 4B - Modern ThrRS + Zn: {len(mod_thrrs_zn)} predictions")
if not mod_thrrs_zn.empty:
    mod_thrrs_zn = mod_thrrs_zn.sort_values('AA_iptm', ascending=False)
    mod_thrrs_zn['rank'] = range(1, len(mod_thrrs_zn) + 1)

    thr_score = mod_thrrs_zn[mod_thrrs_zn['primary_ligand'] == 'THR']['AA_iptm'].values[0]
    mod_thrrs_zn['pct_of_cognate'] = (mod_thrrs_zn['AA_iptm'] / thr_score) * 100

    print("\nAll substrates:")
    for _, row in mod_thrrs_zn.iterrows():
        print(f"  {row['rank']:2d}. {row['primary_ligand']:3s}: {row['AA_iptm']:.3f} ({row['pct_of_cognate']:.1f}%)")

    mod_thrrs_zn.to_csv('figures/data/fig4b_mod_thrrs_zn_all.csv', index=False)

# ============================================================================
# FIGURE 5: ZINC TRAP
# ============================================================================
print("\n" + "="*80)
print("FIGURE 5: ZINC TRAP")
print("="*80)

# Focus on THR, SER, ILE
zinc_trap = mod_thrrs_zn[mod_thrrs_zn['primary_ligand'].isin(['THR', 'SER', 'ILE'])].copy()

print(f"\nFig 5B - Zinc trap comparison:")
for _, row in zinc_trap.iterrows():
    print(f"  {row['primary_ligand']:3s}: AA={row['AA_iptm']:.3f}, Zn={row['Zn_iptm']:.3f}")

zinc_trap.to_csv('figures/data/fig5b_zinc_trap.csv', index=False)

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================
print("\n" + "="*80)
print("SUMMARY STATISTICS")
print("="*80)

summary = {
    'total_predictions': len(df_categorized),
    'ancestral_thrrs': len(df_categorized[(df_categorized['enzyme'] == 'ThrRS') & (df_categorized['era'] == 'ancestral')]),
    'modern_thrrs': len(df_categorized[(df_categorized['enzyme'] == 'ThrRS') & (df_categorized['era'] == 'modern')]),
    'ancestral_prors': len(df_categorized[(df_categorized['enzyme'] == 'ProRS') & (df_categorized['era'] == 'ancestral')]),
    'modern_prors': len(df_categorized[(df_categorized['enzyme'] == 'ProRS') & (df_categorized['era'] == 'modern')]),
    'with_zinc': len(df_categorized[df_categorized['has_zinc'] == True]),
    'competitions': len(competitions)
}

print(json.dumps(summary, indent=2))

with open('figures/data/catalog_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)

print("\n" + "="*80)
print("CATALOG COMPLETE")
print("="*80)
print(f"\nData files saved to: figures/data/")
print(f"  - categorized_predictions.csv (full dataset)")
print(f"  - fig1c_anc_thrrs_no_zn.csv")
print(f"  - fig1d_anc_prors.csv")
print(f"  - fig2b_mod_prors_catalytic.csv")
print(f"  - fig2c_prors_editing.csv")
print(f"  - fig3_competitions.csv")
print(f"  - fig4b_mod_thrrs_zn_all.csv")
print(f"  - fig5b_zinc_trap.csv")
print(f"  - catalog_summary.json")
