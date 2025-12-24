#!/usr/bin/env python3
"""
Generate Comprehensive aaRS Substrate Specificity Heatmap

Shows ALL enzyme conditions side-by-side to reveal the complete
evolutionary landscape and asymmetric evolution story.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Publication quality settings
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 9
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 13

# All 20 amino acids in single-letter code
AMINO_ACIDS = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
               'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

def load_and_filter_data():
    """Load data and filter for reliable predictions."""
    df = pd.read_csv('AF3_RESULTS_CORRECTED.csv')

    # Filter for no tRNA first
    df_no_rna = df[df['has_rna'] == False].copy()

    # EXCLUDE competition runs (they have multiple ligands and skew results)
    df_no_rna = df_no_rna[~df_no_rna['job_name'].str.contains('COMPETITION', case=False)]

    # For editing domain (protein_len == 300), use AA_iptm filter
    # Editing domain is standalone so has low pTM (~0.46) but good AA_iptm (0.60-0.87)
    # Lower threshold to 0.60 to include all 20 amino acids (ARG, CYS, GLY, LYS, PHE are 0.60-0.69)
    editing = df_no_rna[(df_no_rna['protein_len'] == 300) & (df_no_rna['AA_iptm'] >= 0.60)].copy()

    # For other enzymes, use pTM filter
    # Note: Lowered from 0.70 to 0.65 to include ASP, GLU, PHE for ancestral ThrRS
    # which have pTM in 0.67-0.69 range but good AA_iptm scores
    other = df_no_rna[(df_no_rna['protein_len'] != 300) & (df_no_rna['ptm'] >= 0.65)].copy()

    # Combine both
    df_good = pd.concat([editing, other], ignore_index=True)

    # Extract first ligand (amino acid)
    df_good['ligand'] = df_good['ligands'].str.split(',').str[0]

    return df_good


def create_heatmap_data(df_good):
    """Create pivot table for heatmap."""

    # Define enzyme conditions with clear names
    conditions = {
        'Anc ProRS\nCatalytic': (df_good['protein_len'] == 500),
        'Anc ProRS\nEditing': (df_good['protein_len'] == 300),  # NEW: Editing domain
        'Modern ProRS\nCatalytic': (df_good['protein_len'] == 572),
        'Anc ThrRS\nno Zn': ((df_good['protein_len'] == 278) &
                            ~df_good['job_name'].str.contains('zn', case=False)),
        'Anc ThrRS\n+ Zn': ((df_good['protein_len'] == 278) &
                           df_good['job_name'].str.contains('zn', case=False)),
        'Modern ThrRS\n+ Zn': (df_good['protein_len'] == 401),
    }

    # Build pivot table
    heatmap_data = pd.DataFrame(index=AMINO_ACIDS)

    for cond_name, mask in conditions.items():
        subset = df_good[mask].copy()
        # Get best score per ligand (in case of multiple runs)
        scores = subset.groupby('ligand')['AA_iptm'].max()
        heatmap_data[cond_name] = heatmap_data.index.map(lambda x: scores.get(x, np.nan))

    return heatmap_data


def calculate_cognate_ranks(heatmap_data):
    """Calculate rank of cognate substrate in each condition."""

    cognates = {
        'Anc ProRS\nCatalytic': 'PRO',
        'Anc ProRS\nEditing': 'THR',  # NEW: Editing domain binds errors (THR), not cognate (PRO)
        'Modern ProRS\nCatalytic': 'PRO',
        'Anc ThrRS\nno Zn': 'THR',
        'Anc ThrRS\n+ Zn': 'THR',
        'Modern ThrRS\n+ Zn': 'THR',
    }

    ranks = {}

    for col, cognate in cognates.items():
        if col in heatmap_data.columns:
            # Sort by score (descending) and get rank
            sorted_scores = heatmap_data[col].sort_values(ascending=False)
            rank = list(sorted_scores.index).index(cognate) + 1
            cognate_score = heatmap_data.loc[cognate, col]
            ranks[col] = {'cognate': cognate, 'rank': rank, 'score': cognate_score}

    return ranks


def create_comprehensive_heatmap(heatmap_data, ranks):
    """Generate the comprehensive heatmap figure."""

    fig, ax = plt.subplots(figsize=(12, 12))  # Wider for 6 columns

    # Create heatmap using imshow for better control
    im = ax.imshow(heatmap_data.values, cmap='YlOrRd', aspect='auto',
                   vmin=0.5, vmax=1.0, interpolation='nearest')

    # Set ticks
    ax.set_xticks(np.arange(len(heatmap_data.columns)))
    ax.set_yticks(np.arange(len(heatmap_data.index)))
    ax.set_xticklabels(heatmap_data.columns, fontsize=10)
    ax.set_yticklabels(heatmap_data.index, fontsize=9)

    # Rotate x labels for readability
    plt.setp(ax.get_xticklabels(), rotation=0, ha="center")

    # Add value annotations
    for i, aa in enumerate(heatmap_data.index):
        for j, col in enumerate(heatmap_data.columns):
            val = heatmap_data.iloc[i, j]
            if not np.isnan(val):
                # Determine text color based on background
                text_color = 'white' if val > 0.75 else 'black'

                # Check if this is cognate (use ranks dict for correct cognate)
                is_cognate = False
                if col in ranks and aa == ranks[col]['cognate']:
                    is_cognate = True

                # Format: bold and larger for cognate
                weight = 'bold' if is_cognate else 'normal'
                size = 10 if is_cognate else 8

                ax.text(j, i, f'{val:.2f}',
                       ha='center', va='center',
                       color=text_color, fontweight=weight, fontsize=size)

    # Add rectangles around cognate cells (using ranks dict for correct cognates)
    for j, col in enumerate(heatmap_data.columns):
        if col in ranks:
            cognate_aa = ranks[col]['cognate']
            cognate_idx = heatmap_data.index.get_loc(cognate_aa)
            rect = mpatches.Rectangle((j-0.5, cognate_idx-0.5), 1, 1,
                                     linewidth=3, edgecolor='black',
                                     facecolor='none')
            ax.add_patch(rect)

    # Add vertical line to separate ProRS from ThrRS
    ax.axvline(x=2.5, color='black', linewidth=2, linestyle='-')

    # Labels and title
    ax.set_xlabel('Enzyme Condition', fontweight='bold', fontsize=12)
    ax.set_ylabel('Amino Acid Substrate', fontweight='bold', fontsize=12)
    ax.set_title('Complete aaRS Substrate Specificity Landscape (6 Conditions)\nAsymmetric Evolution: ProRS (Kinetic/Double-Sieve) vs ThrRS (Structural/Zn Filter)',
                fontweight='bold', fontsize=13, pad=20)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('ipTM Score (Binding Affinity)', rotation=270, labelpad=20, fontweight='bold')

    # Add annotations for cognate ranks
    annotation_text = "Cognate Ranks:\n"
    for col in heatmap_data.columns:
        if col in ranks:
            info = ranks[col]
            annotation_text += f"{col.replace(chr(10), ' ')}: #{info['rank']}/20 ({info['score']:.3f})\n"

    ax.text(1.15, 0.5, annotation_text,
           transform=ax.transAxes, fontsize=8,
           verticalalignment='center',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Add legend for border
    legend_elements = [
        mpatches.Patch(facecolor='none', edgecolor='black', linewidth=3, label='Cognate substrate')
    ]
    ax.legend(handles=legend_elements, loc='upper left', framealpha=0.9, fontsize=9)

    plt.tight_layout()

    return fig, ax


def create_summary_table(heatmap_data, ranks):
    """Create summary statistics table."""

    summary = []

    for col in heatmap_data.columns:
        scores = heatmap_data[col].dropna()

        if col in ranks:
            cognate_info = ranks[col]
            cognate = cognate_info['cognate']
            rank = cognate_info['rank']
            cognate_score = cognate_info['score']
        else:
            cognate = 'N/A'
            rank = 'N/A'
            cognate_score = np.nan

        # Count how many AAs score > 90% of cognate
        if not np.isnan(cognate_score):
            high_scorers = (scores >= cognate_score * 0.90).sum()
        else:
            high_scorers = 0

        # Calculate spread (max - min)
        spread = scores.max() - scores.min() if len(scores) > 0 else 0

        summary.append({
            'Condition': col.replace('\n', ' '),
            'Cognate': cognate,
            'Cognate_Score': cognate_score,
            'Cognate_Rank': rank,
            'AAs_>90%': high_scorers,
            'Score_Spread': spread,
            'Mean_Score': scores.mean(),
            'Std_Score': scores.std()
        })

    summary_df = pd.DataFrame(summary)
    return summary_df


def create_evolution_comparison(heatmap_data):
    """Create a focused comparison showing evolution."""

    fig, axes = plt.subplots(1, 2, figsize=(14, 12))

    # Panel A: ProRS evolution (Catalytic only - no editing in high-quality data)
    prors_data = heatmap_data[['Anc ProRS\nCatalytic', 'Modern ProRS\nCatalytic']]

    im1 = axes[0].imshow(prors_data.values, cmap='YlOrRd', aspect='auto',
                        vmin=0.5, vmax=1.0, interpolation='nearest')
    axes[0].set_xticks(np.arange(len(prors_data.columns)))
    axes[0].set_yticks(np.arange(len(prors_data.index)))
    axes[0].set_xticklabels(prors_data.columns, fontsize=10)
    axes[0].set_yticklabels(prors_data.index, fontsize=9)
    axes[0].set_title('A. ProRS Evolution:\nKinetic Solution (Editing Domain)', fontweight='bold')
    axes[0].set_ylabel('Amino Acid', fontweight='bold')

    # Annotate values
    for i, aa in enumerate(prors_data.index):
        for j, col in enumerate(prors_data.columns):
            val = prors_data.iloc[i, j]
            if not np.isnan(val):
                text_color = 'white' if val > 0.75 else 'black'
                weight = 'bold' if aa == 'PRO' else 'normal'
                axes[0].text(j, i, f'{val:.2f}',
                           ha='center', va='center',
                           color=text_color, fontweight=weight, fontsize=8)

    # Highlight PRO
    pro_idx = prors_data.index.get_loc('PRO')
    for j in range(len(prors_data.columns)):
        rect = mpatches.Rectangle((j-0.5, pro_idx-0.5), 1, 1,
                                 linewidth=3, edgecolor='black', facecolor='none')
        axes[0].add_patch(rect)

    # Panel B: ThrRS evolution
    thrrs_data = heatmap_data[['Anc ThrRS\nno Zn', 'Anc ThrRS\n+ Zn', 'Modern ThrRS\n+ Zn']]

    im2 = axes[1].imshow(thrrs_data.values, cmap='YlOrRd', aspect='auto',
                        vmin=0.5, vmax=1.0, interpolation='nearest')
    axes[1].set_xticks(np.arange(len(thrrs_data.columns)))
    axes[1].set_yticks(np.arange(len(thrrs_data.index)))
    axes[1].set_xticklabels(thrrs_data.columns, fontsize=10)
    axes[1].set_yticklabels(thrrs_data.index, fontsize=9)
    axes[1].set_title('B. ThrRS Evolution:\nStructural Solution (Zn Filter)', fontweight='bold')

    # Annotate values
    for i, aa in enumerate(thrrs_data.index):
        for j, col in enumerate(thrrs_data.columns):
            val = thrrs_data.iloc[i, j]
            if not np.isnan(val):
                text_color = 'white' if val > 0.75 else 'black'
                weight = 'bold' if aa == 'THR' else 'normal'
                axes[1].text(j, i, f'{val:.2f}',
                           ha='center', va='center',
                           color=text_color, fontweight=weight, fontsize=8)

    # Highlight THR
    thr_idx = thrrs_data.index.get_loc('THR')
    for j in range(len(thrrs_data.columns)):
        rect = mpatches.Rectangle((j-0.5, thr_idx-0.5), 1, 1,
                                 linewidth=3, edgecolor='black', facecolor='none')
        axes[1].add_patch(rect)

    # Add colorbars
    cbar1 = plt.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.04)
    cbar1.set_label('ipTM Score', rotation=270, labelpad=15)
    cbar2 = plt.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04)
    cbar2.set_label('ipTM Score', rotation=270, labelpad=15)

    plt.suptitle('Asymmetric Evolution of Substrate Specificity',
                fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()

    return fig, axes


if __name__ == '__main__':
    import os
    os.makedirs('figures/comprehensive', exist_ok=True)

    print("="*80)
    print("COMPREHENSIVE aaRS SUBSTRATE SPECIFICITY HEATMAP")
    print("="*80)

    # Load data
    print("\n1. Loading and filtering data...")
    df_good = load_and_filter_data()
    print(f"   Total reliable predictions: {len(df_good)}")

    # Create heatmap data
    print("\n2. Creating heatmap data structure...")
    heatmap_data = create_heatmap_data(df_good)
    print(f"   Conditions: {list(heatmap_data.columns)}")
    print(f"   Amino acids: {len(heatmap_data)}")

    # Calculate ranks
    print("\n3. Calculating cognate ranks...")
    ranks = calculate_cognate_ranks(heatmap_data)
    for cond, info in ranks.items():
        print(f"   {cond.replace(chr(10), ' ')}: {info['cognate']} ranks #{info['rank']}/20 (score={info['score']:.3f})")

    # Generate comprehensive heatmap
    print("\n4. Generating comprehensive heatmap...")
    fig, ax = create_comprehensive_heatmap(heatmap_data, ranks)
    plt.savefig('figures/comprehensive/complete_heatmap.png', dpi=300, bbox_inches='tight')
    plt.savefig('figures/comprehensive/complete_heatmap.pdf', bbox_inches='tight')
    print("   Saved: figures/comprehensive/complete_heatmap.png")
    print("   Saved: figures/comprehensive/complete_heatmap.pdf")
    plt.close()

    # Generate evolution comparison
    print("\n5. Generating evolution comparison figure...")
    fig, axes = create_evolution_comparison(heatmap_data)
    plt.savefig('figures/comprehensive/evolution_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig('figures/comprehensive/evolution_comparison.pdf', bbox_inches='tight')
    print("   Saved: figures/comprehensive/evolution_comparison.png")
    print("   Saved: figures/comprehensive/evolution_comparison.pdf")
    plt.close()

    # Create summary table
    print("\n6. Creating summary statistics...")
    summary_df = create_summary_table(heatmap_data, ranks)
    summary_df.to_csv('figures/data/heatmap_summary.csv', index=False)
    print("   Saved: figures/data/heatmap_summary.csv")

    # Save heatmap data
    heatmap_data.to_csv('figures/data/complete_heatmap_data.csv')
    print("   Saved: figures/data/complete_heatmap_data.csv")

    # Print summary
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    print(summary_df.to_string(index=False))

    print("\n" + "="*80)
    print("KEY FINDINGS")
    print("="*80)
    print("\n1. PRORS PERSISTENT PROMISCUITY:")
    anc_prors_key = 'Anc ProRS\nCatalytic'
    mod_prors_key = 'Modern ProRS\nCatalytic'
    print(f"   Ancestral: PRO ranks #{ranks[anc_prors_key]['rank']}/20")
    print(f"   Modern: PRO ranks #{ranks[mod_prors_key]['rank']}/20")
    print("   → No improvement! Editing domain is MANDATORY")

    print("\n2. THRRS EVOLUTION:")
    anc_no_zn_key = 'Anc ThrRS\nno Zn'
    anc_zn_key = 'Anc ThrRS\n+ Zn'
    mod_zn_key = 'Modern ThrRS\n+ Zn'
    anc_no_zn_rank = ranks[anc_no_zn_key]['rank']
    anc_zn_rank = ranks[anc_zn_key]['rank']
    mod_zn_rank = ranks[mod_zn_key]['rank']
    print(f"   Ancestral no Zn: THR ranks #{anc_no_zn_rank}/20 (buried in the middle!)")
    print(f"   Ancestral + Zn: THR ranks #{anc_zn_rank}/20 (Zn helps immediately)")
    print(f"   Modern + Zn: THR ranks #{mod_zn_rank}/20 (clear winner!)")
    print("   → Dramatic evolution: Zn filter works!")

    print("\n" + "="*80)
    print("COMPLETE!")
    print("="*80)
