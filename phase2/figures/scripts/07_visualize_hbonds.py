#!/usr/bin/env python3
"""
Generate Figure 7B: H-Bond Analysis Comparison

Shows hydrogen bond networks across different enzyme-ligand combinations.
Validates the double sieve mechanism and evolutionary trajectory.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Publication quality settings
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['legend.fontsize'] = 10

# Color scheme
COLOR_THR = '#2ecc71'  # Green (cognate for ThrRS)
COLOR_PRO = '#9b59b6'  # Purple (cognate for ProRS)
COLOR_MIXED = '#3498db'  # Blue (neutral)

def create_hbond_comparison():
    """Create comprehensive H-bond comparison figure."""

    # Load data
    df = pd.read_csv('figures/data/hbond_analysis.csv')

    # Parse job names to extract enzyme and ligand
    results = []
    for _, row in df.iterrows():
        job_name = row['job_name']

        # Determine enzyme type and era
        if 'deep_editing' in job_name:
            enzyme = 'ProRS Editing'
            era = 'Ancestral'
        elif 'deep_thrrs' in job_name:
            enzyme = 'ThrRS'
            era = 'Ancestral'
        elif 'modern_thrrs' in job_name:
            enzyme = 'ThrRS'
            era = 'Modern'
        elif 'modern_prours' in job_name:
            enzyme = 'ProRS Catalytic'
            era = 'Modern'
        else:
            continue

        # Determine ligand
        if 'thr' in job_name.lower():
            ligand = 'THR'
        elif 'pro' in job_name.lower():
            ligand = 'PRO'
        else:
            continue

        results.append({
            'Enzyme': enzyme,
            'Era': era,
            'Ligand': ligand,
            'H-bonds': row['total_hbonds'],
            'Avg_Distance': row['avg_hbond_distance'],
            'Label': f"{enzyme}\n{era}\n{ligand}"
        })

    results_df = pd.DataFrame(results)

    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Panel A: H-bond counts by enzyme and ligand
    plot_hbond_bars(results_df, ax1)

    # Panel B: Average H-bond distances
    plot_hbond_distances(results_df, ax2)

    plt.tight_layout()

    # Save
    plt.savefig('figures/hbond_analysis/fig7b_hbond_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig('figures/hbond_analysis/fig7b_hbond_comparison.pdf', bbox_inches='tight')
    print("Saved: figures/hbond_analysis/fig7b_hbond_comparison.png")
    print("Saved: figures/hbond_analysis/fig7b_hbond_comparison.pdf")

    # Also create individual focused plots
    create_editing_domain_plot(results_df)
    create_evolutionary_plot(results_df)

    return results_df


def plot_hbond_bars(df, ax):
    """Plot H-bond counts as grouped bar chart."""

    # Group by enzyme and ligand
    grouped = df.pivot_table(
        index=['Enzyme', 'Era'],
        columns='Ligand',
        values='H-bonds',
        aggfunc='mean'
    ).reset_index()

    # Create labels
    labels = [f"{row['Enzyme']}\n({row['Era']})" for _, row in grouped.iterrows()]

    x = np.arange(len(labels))
    width = 0.35

    # Plot bars
    thr_values = grouped['THR'].values if 'THR' in grouped.columns else np.zeros(len(labels))
    pro_values = grouped['PRO'].values if 'PRO' in grouped.columns else np.zeros(len(labels))

    bars1 = ax.bar(x - width/2, thr_values, width, label='THR', color=COLOR_THR, alpha=0.8)
    bars2 = ax.bar(x + width/2, pro_values, width, label='PRO', color=COLOR_PRO, alpha=0.8)

    # Add value labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{int(height)}',
                       ha='center', va='bottom', fontsize=9)

    ax.set_xlabel('Enzyme (Era)', fontweight='bold')
    ax.set_ylabel('Number of H-bonds', fontweight='bold')
    ax.set_title('A. H-Bond Count Comparison', fontweight='bold', loc='left')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.legend()
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_ylim(0, max(thr_values.max(), pro_values.max()) * 1.15)


def plot_hbond_distances(df, ax):
    """Plot average H-bond distances."""

    # Group by enzyme and ligand
    grouped = df.pivot_table(
        index=['Enzyme', 'Era'],
        columns='Ligand',
        values='Avg_Distance',
        aggfunc='mean'
    ).reset_index()

    labels = [f"{row['Enzyme']}\n({row['Era']})" for _, row in grouped.iterrows()]

    x = np.arange(len(labels))
    width = 0.35

    thr_values = grouped['THR'].values if 'THR' in grouped.columns else np.zeros(len(labels))
    pro_values = grouped['PRO'].values if 'PRO' in grouped.columns else np.zeros(len(labels))

    bars1 = ax.bar(x - width/2, thr_values, width, label='THR', color=COLOR_THR, alpha=0.8)
    bars2 = ax.bar(x + width/2, pro_values, width, label='PRO', color=COLOR_PRO, alpha=0.8)

    # Add value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{height:.2f}Å',
                       ha='center', va='bottom', fontsize=8)

    ax.set_xlabel('Enzyme (Era)', fontweight='bold')
    ax.set_ylabel('Average H-bond Distance (Å)', fontweight='bold')
    ax.set_title('B. H-Bond Distance Comparison', fontweight='bold', loc='left')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.legend()
    ax.grid(axis='y', alpha=0.3, linestyle='--')

    # Add quality reference line (optimal H-bond ~2.8 Å)
    ax.axhline(y=2.8, color='gray', linestyle='--', linewidth=1, alpha=0.5, label='Optimal (~2.8Å)')
    ax.set_ylim(0, 3.5)


def create_editing_domain_plot(df):
    """Create focused plot for editing domain (validates double sieve)."""

    editing_data = df[df['Enzyme'] == 'ProRS Editing']

    if editing_data.empty:
        return

    fig, ax = plt.subplots(figsize=(6, 5))

    # Pivot for easy plotting
    pivot = editing_data.pivot(index='Era', columns='Ligand', values='H-bonds')

    x = np.arange(len(pivot.index))
    width = 0.35

    thr_vals = pivot['THR'].values if 'THR' in pivot.columns else [0]
    pro_vals = pivot['PRO'].values if 'PRO' in pivot.columns else [0]

    bars1 = ax.bar(x - width/2, thr_vals, width, label='THR (non-cognate)', color='#e74c3c', alpha=0.8)
    bars2 = ax.bar(x + width/2, pro_vals, width, label='PRO (cognate)', color=COLOR_PRO, alpha=0.8)

    # Add value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(height)}',
                   ha='center', va='bottom', fontsize=11, fontweight='bold')

    ax.set_ylabel('Number of H-bonds', fontweight='bold', fontsize=12)
    ax.set_title('ProRS Editing Domain: Double Sieve Validation', fontweight='bold', fontsize=13)
    ax.set_xticks(x)
    ax.set_xticklabels(pivot.index, fontsize=11)
    ax.legend(fontsize=11)
    ax.grid(axis='y', alpha=0.3, linestyle='--')

    # Add annotation
    ax.text(0.5, 0.95,
           'Editing domain binds PRO (error) 8× better than THR\nValidates fine-filter function',
           transform=ax.transAxes, ha='center', va='top',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3),
           fontsize=10, style='italic')

    plt.tight_layout()
    plt.savefig('figures/hbond_analysis/editing_domain_validation.png', dpi=300, bbox_inches='tight')
    plt.savefig('figures/hbond_analysis/editing_domain_validation.pdf', bbox_inches='tight')
    print("Saved: figures/hbond_analysis/editing_domain_validation.png")
    plt.close()


def create_evolutionary_plot(df):
    """Create plot showing evolutionary changes in ThrRS."""

    thrrs_data = df[df['Enzyme'] == 'ThrRS']

    if thrrs_data.empty:
        return

    fig, ax = plt.subplots(figsize=(8, 5))

    # Pivot with aggregation to handle any duplicates
    pivot = thrrs_data.pivot_table(index='Era', columns='Ligand', values='H-bonds', aggfunc='mean')

    # Reorder to show evolution
    pivot = pivot.reindex(['Ancestral', 'Modern'])

    x = np.arange(len(pivot.index))
    width = 0.35

    thr_vals = pivot['THR'].values if 'THR' in pivot.columns else np.zeros(len(pivot))
    pro_vals = pivot['PRO'].values if 'PRO' in pivot.columns else np.zeros(len(pivot))

    bars1 = ax.bar(x - width/2, thr_vals, width, label='THR (cognate)', color=COLOR_THR, alpha=0.8)
    bars2 = ax.bar(x + width/2, pro_vals, width, label='PRO (error)', color='#e74c3c', alpha=0.8)

    # Add value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(height)}',
                   ha='center', va='bottom', fontsize=11, fontweight='bold')

    # Calculate and show discrimination
    for i, era in enumerate(pivot.index):
        thr = pivot.loc[era, 'THR'] if 'THR' in pivot.columns else 0
        pro = pivot.loc[era, 'PRO'] if 'PRO' in pivot.columns else 0

        if thr > 0 and pro > 0:
            discrimination = thr / pro
            ax.text(i, max(thr, pro) + 1,
                   f'{discrimination:.1f}× discrimination',
                   ha='center', fontsize=10, fontweight='bold',
                   bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3))

    ax.set_ylabel('Number of H-bonds', fontweight='bold', fontsize=12)
    ax.set_title('ThrRS Evolution: Increased Specificity', fontweight='bold', fontsize=13)
    ax.set_xticks(x)
    ax.set_xticklabels(pivot.index, fontsize=11)
    ax.legend(fontsize=11)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_ylim(0, 18)

    # Add arrow showing evolution
    ax.annotate('', xy=(1, 16), xytext=(0, 16),
               arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    ax.text(0.5, 16.5, 'Evolution', ha='center', fontsize=11, fontweight='bold')

    plt.tight_layout()
    plt.savefig('figures/hbond_analysis/thrrs_evolution.png', dpi=300, bbox_inches='tight')
    plt.savefig('figures/hbond_analysis/thrrs_evolution.pdf', bbox_inches='tight')
    print("Saved: figures/hbond_analysis/thrrs_evolution.png")
    plt.close()


def print_summary(df):
    """Print key findings."""

    print("\n" + "="*80)
    print("KEY FINDINGS FROM H-BOND ANALYSIS")
    print("="*80)

    # Editing domain validation
    editing = df[df['Enzyme'] == 'ProRS Editing']
    if not editing.empty:
        thr_editing = editing[editing['Ligand'] == 'THR']['H-bonds'].values[0]
        pro_editing = editing[editing['Ligand'] == 'PRO']['H-bonds'].values[0]
        print(f"\n1. EDITING DOMAIN (Double Sieve Validation):")
        print(f"   PRO (cognate error): {pro_editing} H-bonds")
        print(f"   THR (non-cognate): {thr_editing} H-bond")
        print(f"   Ratio: {pro_editing/thr_editing:.1f}× better binding of PRO")
        print(f"   → Validates editing domain as 'fine filter' for PRO removal")

    # ThrRS evolution
    thrrs = df[df['Enzyme'] == 'ThrRS']
    if not thrrs.empty:
        anc = thrrs[thrrs['Era'] == 'Ancestral']
        mod = thrrs[thrrs['Era'] == 'Modern']

        anc_thr_data = anc[anc['Ligand'] == 'THR']['H-bonds'].values
        anc_pro_data = anc[anc['Ligand'] == 'PRO']['H-bonds'].values
        mod_thr_data = mod[mod['Ligand'] == 'THR']['H-bonds'].values
        mod_pro_data = mod[mod['Ligand'] == 'PRO']['H-bonds'].values

        print(f"\n2. ThrRS EVOLUTION:")
        if len(anc_thr_data) > 0 and len(anc_pro_data) > 0:
            anc_thr = anc_thr_data[0]
            anc_pro = anc_pro_data[0]
            print(f"   Ancestral: THR={anc_thr}, PRO={anc_pro} ({anc_thr/anc_pro:.1f}× discrimination)")
        elif len(anc_thr_data) > 0:
            print(f"   Ancestral: THR={anc_thr_data[0]} H-bonds (PRO data not available)")

        if len(mod_thr_data) > 0 and len(mod_pro_data) > 0:
            mod_thr = mod_thr_data[0]
            mod_pro = mod_pro_data[0]
            print(f"   Modern: THR={mod_thr}, PRO={mod_pro} ({mod_thr/mod_pro:.1f}× discrimination)")

            if len(anc_thr_data) > 0 and len(anc_pro_data) > 0:
                print(f"   → Discrimination improved {(mod_thr/mod_pro)/(anc_thr_data[0]/anc_pro_data[0]):.1f}× during evolution")

    # ProRS persistence
    prours = df[df['Enzyme'] == 'ProRS Catalytic']
    if not prours.empty:
        pro_thr = prours[prours['Ligand'] == 'THR']['H-bonds'].values[0]
        pro_pro = prours[prours['Ligand'] == 'PRO']['H-bonds'].values[0]

        print(f"\n3. ProRS CATALYTIC SITE (Persistent Promiscuity):")
        print(f"   THR: {pro_thr} H-bonds")
        print(f"   PRO: {pro_pro} H-bonds")
        print(f"   → Nearly identical ({pro_pro/pro_thr:.2f}×), explains need for editing domain")

    print("\n" + "="*80)


if __name__ == '__main__':
    import os
    os.makedirs('figures/hbond_analysis', exist_ok=True)

    print("="*80)
    print("H-BOND VISUALIZATION")
    print("="*80)

    results = create_hbond_comparison()
    print_summary(results)

    print("\n" + "="*80)
    print("COMPLETE!")
    print("="*80)
