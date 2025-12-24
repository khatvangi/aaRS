#!/usr/bin/env python3
"""
Analyze geometry-based metrics for aaRS-ligand binding.
Focus on promiscuity, specificity, and physical interactions.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set style
sns.set_style('whitegrid')
plt.rcParams['figure.dpi'] = 150

def load_data():
    """Load geometry metrics."""
    df = pd.read_csv('geometry_metrics.csv')
    print(f"Loaded {len(df)} structures")
    return df

def analyze_by_enzyme_ligand(df):
    """Create enzyme x ligand heatmaps for key metrics."""

    # Filter to main enzymes
    df_main = df[df['enzyme'].isin(['ProRS', 'ThrRS'])].copy()

    # Define amino acids
    amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                   'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                   'THR', 'TRP', 'TYR', 'VAL']

    metrics = {
        'pocket_iptm': 'Pocket ipTM',
        'contacts_4A': 'Contacts (4Å)',
        'polar_close_contacts_3p5A': 'Polar close contacts (≤3.5Å)',
        'polar_fraction': 'Polar fraction',
        'clashes_2.2A': 'Clashes (<2.2Å)'
    }

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()

    for idx, (metric, label) in enumerate(metrics.items()):
        ax = axes[idx]

        # Create pivot table
        pivot = df_main.pivot_table(
            values=metric,
            index='enzyme',
            columns='ligand',
            aggfunc='mean'
        )

        # Reorder columns to standard AA order
        cols = [aa for aa in amino_acids if aa in pivot.columns]
        pivot = pivot[cols]

        # Plot heatmap
        if metric == 'clashes_2.2A':
            cmap = 'Reds'  # More clashes = worse
        else:
            cmap = 'viridis'  # More = better (generally)

        sns.heatmap(pivot, annot=True, fmt='.2f', cmap=cmap, ax=ax,
                   cbar_kws={'label': label})
        ax.set_title(label, fontsize=14, fontweight='bold')
        ax.set_xlabel('Ligand', fontsize=12)
        ax.set_ylabel('Enzyme', fontsize=12)

    # Remove extra subplot
    fig.delaxes(axes[5])

    plt.tight_layout()
    plt.savefig('heatmap_geometry_metrics.png', dpi=300, bbox_inches='tight')
    print("Saved heatmap_geometry_metrics.png")
    plt.close()

def analyze_specificity(df):
    """Analyze enzyme specificity using geometry metrics."""

    print("\n" + "="*70)
    print("SPECIFICITY ANALYSIS (Geometry-Based)")
    print("="*70)

    # ProRS specificity
    prors = df[df['enzyme'] == 'ProRS'].copy()
    prors_pro = prors[prors['ligand'] == 'PRO']
    prors_other = prors[prors['ligand'] != 'PRO']

    print("\nProRS Specificity (PRO vs other):")
    print(f"  PRO (n={len(prors_pro)}):")
    print(f"    Pocket ipTM: {prors_pro['pocket_iptm'].mean():.3f} ± {prors_pro['pocket_iptm'].std():.3f}")
    print(f"    Contacts:    {prors_pro['contacts_4A'].mean():.1f} ± {prors_pro['contacts_4A'].std():.1f}")
    print(f"    Polar contacts (≤3.5Å): {prors_pro['polar_close_contacts_3p5A'].mean():.1f} ± {prors_pro['polar_close_contacts_3p5A'].std():.1f}")

    print(f"  Other (n={len(prors_other)}):")
    print(f"    Pocket ipTM: {prors_other['pocket_iptm'].mean():.3f} ± {prors_other['pocket_iptm'].std():.3f}")
    print(f"    Contacts:    {prors_other['contacts_4A'].mean():.1f} ± {prors_other['contacts_4A'].std():.1f}")
    print(f"    Polar contacts (≤3.5Å): {prors_other['polar_close_contacts_3p5A'].mean():.1f} ± {prors_other['polar_close_contacts_3p5A'].std():.1f}")

    print(f"  Δ (PRO - Other):")
    print(f"    Pocket ipTM: {prors_pro['pocket_iptm'].mean() - prors_other['pocket_iptm'].mean():+.3f}")
    print(f"    Contacts:    {prors_pro['contacts_4A'].mean() - prors_other['contacts_4A'].mean():+.1f}")
    print(f"    Polar contacts: {prors_pro['polar_close_contacts_3p5A'].mean() - prors_other['polar_close_contacts_3p5A'].mean():+.1f}")

    # ThrRS specificity
    thrrs = df[df['enzyme'] == 'ThrRS'].copy()
    thrrs_thr = thrrs[thrrs['ligand'] == 'THR']
    thrrs_other = thrrs[thrrs['ligand'] != 'THR']

    print("\nThrRS Specificity (THR vs other):")
    print(f"  THR (n={len(thrrs_thr)}):")
    print(f"    Pocket ipTM: {thrrs_thr['pocket_iptm'].mean():.3f} ± {thrrs_thr['pocket_iptm'].std():.3f}")
    print(f"    Contacts:    {thrrs_thr['contacts_4A'].mean():.1f} ± {thrrs_thr['contacts_4A'].std():.1f}")
    print(f"    Polar contacts (≤3.5Å): {thrrs_thr['polar_close_contacts_3p5A'].mean():.1f} ± {thrrs_thr['polar_close_contacts_3p5A'].std():.1f}")

    print(f"  Other (n={len(thrrs_other)}):")
    print(f"    Pocket ipTM: {thrrs_other['pocket_iptm'].mean():.3f} ± {thrrs_other['pocket_iptm'].std():.3f}")
    print(f"    Contacts:    {thrrs_other['contacts_4A'].mean():.1f} ± {thrrs_other['contacts_4A'].std():.1f}")
    print(f"    Polar contacts (≤3.5Å): {thrrs_other['polar_close_contacts_3p5A'].mean():.1f} ± {thrrs_other['polar_close_contacts_3p5A'].std():.1f}")

    print(f"  Δ (THR - Other):")
    print(f"    Pocket ipTM: {thrrs_thr['pocket_iptm'].mean() - thrrs_other['pocket_iptm'].mean():+.3f}")
    print(f"    Contacts:    {thrrs_thr['contacts_4A'].mean() - thrrs_other['contacts_4A'].mean():+.1f}")
    print(f"    Polar contacts: {thrrs_thr['polar_close_contacts_3p5A'].mean() - thrrs_other['polar_close_contacts_3p5A'].mean():+.1f}")

def analyze_zn_effect(df):
    """Analyze effect of Zn on binding metrics."""

    print("\n" + "="*70)
    print("ZINC EFFECT ANALYSIS")
    print("="*70)

    with_zn = df[df['zn_present'] == True]
    without_zn = df[df['zn_present'] == False]

    print(f"\nWithout Zn (n={len(without_zn)}):")
    print(f"  Pocket ipTM: {without_zn['pocket_iptm'].mean():.3f} ± {without_zn['pocket_iptm'].std():.3f}")
    print(f"  Contacts:    {without_zn['contacts_4A'].mean():.1f} ± {without_zn['contacts_4A'].std():.1f}")
    print(f"  Polar contacts (≤3.5Å): {without_zn['polar_close_contacts_3p5A'].mean():.1f} ± {without_zn['polar_close_contacts_3p5A'].std():.1f}")
    print(f"  Clashes:     {without_zn['clashes_2.2A'].mean():.1f} ± {without_zn['clashes_2.2A'].std():.1f}")

    print(f"\nWith Zn (n={len(with_zn)}):")
    print(f"  Pocket ipTM: {with_zn['pocket_iptm'].mean():.3f} ± {with_zn['pocket_iptm'].std():.3f}")
    print(f"  Contacts:    {with_zn['contacts_4A'].mean():.1f} ± {with_zn['contacts_4A'].std():.1f}")
    print(f"  Polar contacts (≤3.5Å): {with_zn['polar_close_contacts_3p5A'].mean():.1f} ± {with_zn['polar_close_contacts_3p5A'].std():.1f}")
    print(f"  Clashes:     {with_zn['clashes_2.2A'].mean():.1f} ± {with_zn['clashes_2.2A'].std():.1f}")

    # Zn coordination analysis
    zn_df = with_zn[with_zn['zn_min_dist'].notna()].copy()
    print(f"\nZn Coordination (n={len(zn_df)}):")
    print(f"  Min Zn-ligand dist: {zn_df['zn_min_dist'].mean():.2f} ± {zn_df['zn_min_dist'].std():.2f} Å")
    print(f"  Coordination count: {zn_df['zn_coordination_count'].mean():.1f} ± {zn_df['zn_coordination_count'].std():.1f}")
    print(f"  Ligand contacts:    {zn_df['zn_ligand_contacts'].mean():.1f} ± {zn_df['zn_ligand_contacts'].std():.1f}")

def analyze_promiscuity(df):
    """Compute promiscuity index for each enzyme variant."""

    print("\n" + "="*70)
    print("PROMISCUITY INDEX (based on pocket ipTM)")
    print("="*70)

    # Group by job_name (unique structure)
    df_grouped = df.groupby('job_name').agg({
        'enzyme': 'first',
        'epoch': 'first',
        'domain': 'first',
        'pocket_iptm': 'mean',
        'contacts_4A': 'mean',
        'ligand': 'count'  # number of ligands tested
    }).reset_index()

    # Filter to structures with multiple ligands tested
    df_multi = df_grouped[df_grouped['ligand'] >= 3].copy()

    print(f"\nStructures with ≥3 ligands tested: {len(df_multi)}")

    # Compute promiscuity: mean pocket ipTM across ligands
    # Higher mean = binds many ligands well = more promiscuous
    top_promiscuous = df_multi.nlargest(10, 'pocket_iptm')

    print("\nTop 10 most promiscuous (high pocket ipTM across ligands):")
    print(top_promiscuous[['job_name', 'enzyme', 'epoch', 'domain', 'pocket_iptm', 'ligand']].to_string(index=False))

    # Most specific (low mean pocket ipTM)
    top_specific = df_multi.nsmallest(10, 'pocket_iptm')
    print("\nTop 10 most specific (low pocket ipTM across ligands):")
    print(top_specific[['job_name', 'enzyme', 'epoch', 'domain', 'pocket_iptm', 'ligand']].to_string(index=False))

def create_comparison_plots(df):
    """Create comparison plots for key metrics."""

    # Filter to main enzymes
    df_main = df[df['enzyme'].isin(['ProRS', 'ThrRS'])].copy()

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1. Pocket ipTM by enzyme and epoch
    ax = axes[0, 0]
    df_plot = df_main[df_main['pocket_iptm'].notna()]
    sns.boxplot(data=df_plot, x='enzyme', y='pocket_iptm', hue='epoch', ax=ax)
    ax.set_title('Pocket ipTM by Enzyme and Epoch', fontweight='bold')
    ax.set_ylabel('Pocket ipTM')
    ax.legend(title='Epoch')

    # 2. Contacts vs pocket ipTM
    ax = axes[0, 1]
    for enzyme in ['ProRS', 'ThrRS']:
        df_enzyme = df_main[df_main['enzyme'] == enzyme]
        ax.scatter(df_enzyme['contacts_4A'], df_enzyme['pocket_iptm'],
                  label=enzyme, alpha=0.6, s=50)
    ax.set_xlabel('Contacts (4Å)')
    ax.set_ylabel('Pocket ipTM')
    ax.set_title('Contacts vs Pocket ipTM', fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)

    # 3. Polar close contacts by Zn presence
    ax = axes[1, 0]
    sns.violinplot(data=df_main, x='zn_present', y='polar_close_contacts_3p5A', ax=ax)
    ax.set_xlabel('Zn Present')
    ax.set_ylabel('Polar close contacts (≤3.5Å)')
    ax.set_title('Polar close contacts by Zn Presence', fontweight='bold')
    ax.set_xticklabels(['No', 'Yes'])

    # 4. Clashes vs pocket ipTM
    ax = axes[1, 1]
    df_plot = df_main[df_main['clashes_2.2A'] < 20]  # exclude extreme outliers
    for enzyme in ['ProRS', 'ThrRS']:
        df_enzyme = df_plot[df_plot['enzyme'] == enzyme]
        ax.scatter(df_enzyme['clashes_2.2A'], df_enzyme['pocket_iptm'],
                  label=enzyme, alpha=0.6, s=50)
    ax.set_xlabel('Clashes (<2.2Å)')
    ax.set_ylabel('Pocket ipTM')
    ax.set_title('Clashes vs Pocket ipTM', fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig('comparison_plots_geometry.png', dpi=300, bbox_inches='tight')
    print("\nSaved comparison_plots_geometry.png")
    plt.close()

def main():
    df = load_data()

    print("\n" + "="*70)
    print("GEOMETRY-BASED ANALYSIS")
    print("="*70)

    # Create heatmaps
    analyze_by_enzyme_ligand(df)

    # Analyze specificity
    analyze_specificity(df)

    # Analyze Zn effect
    analyze_zn_effect(df)

    # Analyze promiscuity
    analyze_promiscuity(df)

    # Create comparison plots
    create_comparison_plots(df)

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print("\nGenerated files:")
    print("  - heatmap_geometry_metrics.png")
    print("  - comparison_plots_geometry.png")
    print("\nData file:")
    print("  - geometry_metrics.csv")

if __name__ == '__main__':
    main()
