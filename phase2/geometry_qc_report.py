#!/usr/bin/env python3
"""
Rigorous geometry-based analysis with proper stratification and normalization.

Fixes:
1. Rename hbonds → polar_close_contacts_3p5A (no angle check)
2. Split Zn into engaged vs floating
3. Normalize contacts per heavy atom
4. Stratify by (enzyme, epoch, domain, zn_engaged)
5. QC gates for clashes
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set style
sns.set_style('whitegrid')
plt.rcParams['figure.dpi'] = 150

def add_derived_features(df):
    """Add derived features with proper normalization and Zn classification."""

    df = df.copy()

    # 1. Rename H-bonds to polar close contacts (no angle check)
    df.rename(columns={'hbonds': 'polar_close_contacts_3p5A'}, inplace=True)

    # 2. Normalize by ligand size
    df['contacts_per_atom'] = df['contacts_4A'] / df['ligand_heavy_atoms'].clip(lower=1)
    df['polar_contacts_per_atom'] = df['polar_contacts'] / df['ligand_heavy_atoms'].clip(lower=1)
    df['clash_rate'] = df['clashes_2.2A'] / df['ligand_heavy_atoms'].clip(lower=1)
    df['polar_close_contacts_per_atom'] = df['polar_close_contacts_3p5A'] / df['ligand_heavy_atoms'].clip(lower=1)

    # 3. Zn classification: engaged vs floating
    df['zn_engaged'] = (
        (df['zn_present'] == True) &
        (df['zn_min_dist'] <= 3.0) &
        (df['zn_coordination_count'] >= 2)
    )

    df['zn_floating'] = (
        (df['zn_present'] == True) &
        (~df['zn_engaged'])
    )

    # 4. Create condition label
    df['condition'] = df.apply(
        lambda row: f"{row['epoch']}_{row['enzyme']}_{row['domain']}" +
                    ("_Zn" if row['zn_engaged'] else ""),
        axis=1
    )

    return df

def create_clean_dataset(df, clash_threshold=2):
    """Filter to low-clash structures."""
    return df[df['clashes_2.2A'] <= clash_threshold].copy()

def permutation_test(x, y, n_perm=10000, seed=0):
    """
    Two-sample permutation test for difference in means.
    Returns p-value (two-sided).
    """
    rng = np.random.default_rng(seed)
    x = np.asarray(x)
    y = np.asarray(y)
    observed = np.abs(x.mean() - y.mean())
    pooled = np.concatenate([x, y])
    n_x = len(x)
    count = 0
    for _ in range(n_perm):
        rng.shuffle(pooled)
        perm_diff = np.abs(pooled[:n_x].mean() - pooled[n_x:].mean())
        if perm_diff >= observed:
            count += 1
    return (count + 1) / (n_perm + 1)


def benjamini_hochberg(pvals):
    """Apply Benjamini-Hochberg FDR correction. Returns q-values."""
    pvals = np.asarray(pvals)
    n = len(pvals)
    if n == 0:
        return np.array([])
    order = np.argsort(pvals)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, n + 1)
    qvals = pvals * n / ranks
    # Ensure monotonicity
    qvals_sorted = qvals[order]
    for i in range(n - 2, -1, -1):
        qvals_sorted[i] = min(qvals_sorted[i], qvals_sorted[i + 1])
    qvals[order] = qvals_sorted
    return np.clip(qvals, 0, 1)


def compute_cognate_vs_noncognate(df, n_bootstrap=10000, n_perm=10000, min_n_cognate=3):
    """
    Compute cognate vs non-cognate effect for each condition.
    Returns DataFrame with Δ values, 95% CI, permutation p-values, and BH-corrected q-values.

    IMPORTANT: Conditions with n_cognate < min_n_cognate are excluded from inferential
    statistics (p-values, q-values) because bootstrap/permutation tests are unreliable
    with tiny samples.
    """

    results = []

    # Define cognate pairs
    cognate_map = {
        'ProRS': 'PRO',
        'ThrRS': 'THR'
    }

    conditions = df['condition'].unique()

    for condition in conditions:
        df_cond = df[df['condition'] == condition]

        # Get enzyme type
        enzyme = df_cond['enzyme'].iloc[0]
        if enzyme not in cognate_map:
            continue

        cognate_aa = cognate_map[enzyme]

        # Split cognate vs non-cognate
        df_cognate = df_cond[df_cond['ligand'] == cognate_aa]
        df_noncognate = df_cond[df_cond['ligand'] != cognate_aa]

        if len(df_cognate) == 0 or len(df_noncognate) == 0:
            continue

        # Metrics to analyze
        metrics = ['pocket_iptm', 'contacts_per_atom', 'polar_contacts_per_atom',
                  'polar_close_contacts_per_atom', 'clash_rate']

        n_cog = len(df_cognate)
        n_noncog = len(df_noncognate)

        # Flag: is sample size sufficient for inferential stats?
        sufficient_n = (n_cog >= min_n_cognate)

        row = {
            'condition': condition,
            'enzyme': enzyme,
            'cognate_aa': cognate_aa,
            'n_cognate': n_cog,
            'n_noncognate': n_noncog,
            'sufficient_n': sufficient_n
        }

        for metric in metrics:
            cog_vals = df_cognate[metric].values
            noncog_vals = df_noncognate[metric].values

            # Mean values
            cognate_mean = cog_vals.mean()
            noncognate_mean = noncog_vals.mean()
            delta = cognate_mean - noncognate_mean

            # Bootstrap 95% CI for delta (correct method: resample both groups)
            rng = np.random.default_rng(0)
            bootstrap_deltas = []
            for _ in range(n_bootstrap):
                cog_sample = rng.choice(cog_vals, size=len(cog_vals), replace=True).mean()
                noncog_sample = rng.choice(noncog_vals, size=len(noncog_vals), replace=True).mean()
                bootstrap_deltas.append(cog_sample - noncog_sample)

            ci_lower = np.percentile(bootstrap_deltas, 2.5)
            ci_upper = np.percentile(bootstrap_deltas, 97.5)

            # Permutation p-value (only if sufficient sample size)
            if sufficient_n:
                pval = permutation_test(cog_vals, noncog_vals, n_perm=n_perm, seed=0)
            else:
                pval = np.nan  # Not computed for small samples

            row[f'{metric}_cognate'] = cognate_mean
            row[f'{metric}_noncognate'] = noncognate_mean
            row[f'{metric}_delta'] = delta
            row[f'{metric}_ci_lower'] = ci_lower
            row[f'{metric}_ci_upper'] = ci_upper
            row[f'{metric}_pval'] = pval

        results.append(row)

    df_results = pd.DataFrame(results)

    # Apply BH-FDR correction across all tests (per metric)
    metrics = ['pocket_iptm', 'contacts_per_atom', 'polar_contacts_per_atom',
              'polar_close_contacts_per_atom', 'clash_rate']

    for metric in metrics:
        pvals = df_results[f'{metric}_pval'].values
        # Only apply FDR to non-NaN p-values
        valid_mask = ~np.isnan(pvals)
        qvals = np.full_like(pvals, np.nan)
        if valid_mask.sum() > 0:
            qvals[valid_mask] = benjamini_hochberg(pvals[valid_mask])
        df_results[f'{metric}_qval'] = qvals
        # Significance based on FDR-corrected q < 0.05
        df_results[f'{metric}_significant'] = df_results[f'{metric}_qval'] < 0.05

    return df_results

def create_stratified_heatmaps(df, output_file='heatmap_stratified_conditions.png'):
    """Create heatmaps stratified by condition (not just enzyme)."""

    # Standard amino acid order
    amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                   'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                   'THR', 'TRP', 'TYR', 'VAL']

    # Metrics to plot (normalized)
    metrics = {
        'pocket_iptm': 'Pocket ipTM\n(AF3 local support)',
        'contacts_per_atom': 'Contacts per atom\n(shape complementarity)',
        'polar_contacts_per_atom': 'Polar contacts per atom',
        'clash_rate': 'Clash rate\n(steric problems)'
    }

    # Get top conditions by sample size
    condition_counts = df.groupby('condition').size().sort_values(ascending=False)
    top_conditions = condition_counts.head(12).index.tolist()

    fig, axes = plt.subplots(4, 3, figsize=(20, 16))
    axes = axes.flatten()

    for idx, (metric, label) in enumerate(metrics.items()):
        ax = axes[idx]

        # Create pivot table
        pivot = df[df['condition'].isin(top_conditions)].pivot_table(
            values=metric,
            index='condition',
            columns='ligand',
            aggfunc='mean'
        )

        # Reorder columns to standard AA order
        cols = [aa for aa in amino_acids if aa in pivot.columns]
        pivot = pivot[cols]

        # Sort rows by condition name
        pivot = pivot.sort_index()

        # Plot heatmap
        if metric == 'clash_rate':
            cmap = 'Reds'  # More clashes = worse
            vmax = pivot.values.max() if pivot.values.max() > 0 else 1
            vmin = 0
        else:
            cmap = 'viridis'
            vmax = None
            vmin = None

        sns.heatmap(pivot, annot=True, fmt='.2f', cmap=cmap, ax=ax,
                   cbar_kws={'label': label}, vmin=vmin, vmax=vmax)
        ax.set_title(label, fontsize=12, fontweight='bold')
        ax.set_xlabel('Ligand', fontsize=10)
        ax.set_ylabel('Condition', fontsize=10)
        ax.tick_params(axis='y', labelsize=8)
        ax.tick_params(axis='x', labelsize=8)

    # Remove extra subplots
    for idx in range(len(metrics), 12):
        fig.delaxes(axes[idx])

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved {output_file}")
    plt.close()

def plot_cognate_effects(df_effects, output_file='cognate_vs_noncognate_effects.png'):
    """
    Plot cognate vs non-cognate effects with confidence intervals.

    IMPORTANT: Significance coloring is based on FDR-corrected q < 0.05 (permutation test),
    NOT on whether CI excludes 0. Conditions with n_cognate < 3 are marked as
    "insufficient n" and do not receive inferential markers.
    """

    if len(df_effects) == 0:
        print("No cognate/non-cognate comparisons available")
        return

    metrics = ['pocket_iptm', 'contacts_per_atom', 'polar_contacts_per_atom']

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    for idx, metric in enumerate(metrics):
        ax = axes[idx]

        # Sort by delta
        df_plot = df_effects.sort_values(f'{metric}_delta', ascending=False)

        # Plot
        x = range(len(df_plot))
        y = df_plot[f'{metric}_delta'].values
        yerr_lower = y - df_plot[f'{metric}_ci_lower'].values
        yerr_upper = df_plot[f'{metric}_ci_upper'].values - y

        # Color scheme:
        # - Green: significant (q < 0.05 after FDR)
        # - Gray: not significant but sufficient n
        # - Light red: insufficient sample size (n_cognate < 3)
        colors = []
        for _, row in df_plot.iterrows():
            if not row['sufficient_n']:
                colors.append('lightcoral')  # Insufficient n
            elif row[f'{metric}_significant']:
                colors.append('green')  # Significant (q < 0.05)
            else:
                colors.append('steelblue')  # Not significant

        ax.barh(x, y, color=colors, alpha=0.7)
        ax.errorbar(y, x, xerr=[yerr_lower, yerr_upper], fmt='none',
                   color='black', capsize=3, linewidth=1)

        ax.axvline(0, color='black', linestyle='--', linewidth=1)
        ax.set_yticks(x)

        # Label with sample size
        labels = [f"{row['condition']} (n={row['n_cognate']})"
                 for _, row in df_plot.iterrows()]
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel(f'Δ {metric.replace("_", " ")}\n(cognate − non-cognate)', fontsize=10)
        ax.set_title(f'{metric.replace("_", " ").title()}', fontweight='bold')
        ax.grid(axis='x', alpha=0.3)

        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='green', alpha=0.7, label='Significant (q < 0.05, FDR)'),
            Patch(facecolor='steelblue', alpha=0.7, label='Not significant'),
            Patch(facecolor='lightcoral', alpha=0.7, label='Insufficient n (< 3)')
        ]
        ax.legend(handles=legend_elements, loc='best', fontsize=8)

    # Add note about methodology
    fig.text(0.5, 0.01,
             'Significance: permutation test (10k) + Benjamini-Hochberg FDR correction. '
             'CIs: bootstrap (10k) on difference.',
             ha='center', fontsize=9, style='italic')

    plt.tight_layout(rect=[0, 0.03, 1, 1])
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved {output_file}")
    plt.close()

def plot_zn_diagnostics(df, output_file='zn_diagnostics.png'):
    """Create Zn diagnostic plots."""

    df_zn = df[df['zn_present'] == True].copy()

    if len(df_zn) == 0:
        print("No Zn-containing structures")
        return

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # 1. Histogram of zn_min_dist
    ax = axes[0]
    df_zn_valid = df_zn[df_zn['zn_min_dist'].notna()]
    ax.hist(df_zn_valid['zn_min_dist'], bins=30, color='steelblue', alpha=0.7, edgecolor='black')
    ax.axvline(3.0, color='red', linestyle='--', linewidth=2, label='Engaged cutoff (3.0 Å)')
    ax.set_xlabel('Min Zn-ligand distance (Å)', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Zn-Ligand Distance Distribution', fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)

    # 2. Pocket ipTM vs zn_min_dist scatter
    ax = axes[1]
    df_plot = df_zn[df_zn['zn_min_dist'].notna()]

    engaged = df_plot['zn_engaged']
    ax.scatter(df_plot[engaged]['zn_min_dist'], df_plot[engaged]['pocket_iptm'],
              c='green', s=50, alpha=0.6, label=f'Engaged (n={engaged.sum()})')
    ax.scatter(df_plot[~engaged]['zn_min_dist'], df_plot[~engaged]['pocket_iptm'],
              c='red', s=50, alpha=0.6, label=f'Floating (n={(~engaged).sum()})')

    ax.axvline(3.0, color='red', linestyle='--', linewidth=2, alpha=0.5)
    ax.set_xlabel('Min Zn-ligand distance (Å)', fontsize=12)
    ax.set_ylabel('Pocket ipTM', fontsize=12)
    ax.set_title('Pocket ipTM vs Zn Distance', fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)

    # 3. Zn engagement stats
    ax = axes[2]

    zn_stats = pd.DataFrame({
        'Category': ['Engaged', 'Floating'],
        'Count': [df_zn['zn_engaged'].sum(), df_zn['zn_floating'].sum()],
        'Mean pocket ipTM': [
            df_zn[df_zn['zn_engaged']]['pocket_iptm'].mean(),
            df_zn[df_zn['zn_floating']]['pocket_iptm'].mean()
        ]
    })

    x = range(len(zn_stats))
    ax.bar(x, zn_stats['Count'], color=['green', 'red'], alpha=0.7)

    for i, row in zn_stats.iterrows():
        ax.text(i, row['Count'] + 1,
               f"n={int(row['Count'])}\nipTM={row['Mean pocket ipTM']:.3f}",
               ha='center', fontsize=10)

    ax.set_xticks(x)
    ax.set_xticklabels(zn_stats['Category'])
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Zn Engagement Classification', fontweight='bold')
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved {output_file}")
    plt.close()

def main():
    print("="*70)
    print("RIGOROUS GEOMETRY QC REPORT")
    print("="*70)

    # 1. Load data
    df = pd.read_csv('geometry_metrics.csv')
    print(f"\nLoaded {len(df)} structures from geometry_metrics.csv")

    # 2. Add derived features
    print("\nAdding derived features...")
    df = add_derived_features(df)

    # Print Zn classification
    print(f"\nZn classification:")
    print(f"  Zn engaged:  {df['zn_engaged'].sum()} structures")
    print(f"  Zn floating: {df['zn_floating'].sum()} structures")
    print(f"  No Zn:       {(~df['zn_present']).sum()} structures")

    # 3. Save derived dataset
    output_derived = 'geometry_metrics_derived.csv'
    df.to_csv(output_derived, index=False)
    print(f"\nSaved {output_derived}")

    # 4. Create clean dataset
    print("\nCreating clean dataset (clashes ≤ 2)...")
    df_clean = create_clean_dataset(df, clash_threshold=2)
    print(f"  Clean dataset: {len(df_clean)}/{len(df)} structures ({100*len(df_clean)/len(df):.1f}%)")

    output_clean = 'geometry_metrics_clean.csv'
    df_clean.to_csv(output_clean, index=False)
    print(f"Saved {output_clean}")

    # 5. Create stratified heatmaps (use clean dataset)
    print("\nCreating stratified heatmaps...")
    create_stratified_heatmaps(df_clean)

    # 6. Compute cognate vs non-cognate effects
    print("\nComputing cognate vs non-cognate effects (with bootstrap 95% CI)...")
    df_effects = compute_cognate_vs_noncognate(df_clean, n_bootstrap=1000)

    if len(df_effects) > 0:
        output_effects = 'cognate_vs_noncognate_effects.csv'
        df_effects.to_csv(output_effects, index=False)
        print(f"Saved {output_effects}")

        # Print summary
        print("\n" + "="*70)
        print("COGNATE SPECIFICITY SUMMARY")
        print("="*70)
        print("(Significance: q < 0.05 after BH-FDR correction)")
        print("(Conditions with n_cognate < 3 excluded from inferential tests)")

        for _, row in df_effects.iterrows():
            print(f"\n{row['condition']}:")
            print(f"  Cognate: {row['cognate_aa']} (n={row['n_cognate']})")
            print(f"  Non-cognate: n={row['n_noncognate']}")

            if not row['sufficient_n']:
                print(f"  ⚠️  Insufficient sample size (n < 3) - no inferential stats")

            # Pocket ipTM
            if row['sufficient_n']:
                sig = "✓ SIG" if row['pocket_iptm_significant'] else "✗ NS"
                qval = row['pocket_iptm_qval']
                print(f"  Δ pocket_iptm: {row['pocket_iptm_delta']:+.3f} "
                      f"[{row['pocket_iptm_ci_lower']:+.3f}, {row['pocket_iptm_ci_upper']:+.3f}] q={qval:.3f} {sig}")
            else:
                print(f"  Δ pocket_iptm: {row['pocket_iptm_delta']:+.3f} "
                      f"[{row['pocket_iptm_ci_lower']:+.3f}, {row['pocket_iptm_ci_upper']:+.3f}] (no test)")

            # Contacts per atom
            if row['sufficient_n']:
                sig = "✓ SIG" if row['contacts_per_atom_significant'] else "✗ NS"
                qval = row['contacts_per_atom_qval']
                print(f"  Δ contacts/atom: {row['contacts_per_atom_delta']:+.3f} "
                      f"[{row['contacts_per_atom_ci_lower']:+.3f}, {row['contacts_per_atom_ci_upper']:+.3f}] q={qval:.3f} {sig}")
            else:
                print(f"  Δ contacts/atom: {row['contacts_per_atom_delta']:+.3f} "
                      f"[{row['contacts_per_atom_ci_lower']:+.3f}, {row['contacts_per_atom_ci_upper']:+.3f}] (no test)")

        # Plot effects
        plot_cognate_effects(df_effects)

    # 7. Zn diagnostics
    print("\nCreating Zn diagnostics...")
    plot_zn_diagnostics(df)

    # 8. Summary statistics
    print("\n" + "="*70)
    print("DATASET STATISTICS")
    print("="*70)

    print(f"\nAll structures (n={len(df)}):")
    print(df[['pocket_iptm', 'contacts_per_atom', 'polar_contacts_per_atom',
             'clash_rate']].describe())

    print(f"\nClean structures (n={len(df_clean)}):")
    print(df_clean[['pocket_iptm', 'contacts_per_atom', 'polar_contacts_per_atom',
                    'clash_rate']].describe())

    # Condition breakdown
    print(f"\nBy condition (clean dataset):")
    print(df_clean.groupby('condition').size().sort_values(ascending=False))

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print("\nGenerated files:")
    print("  - geometry_metrics_derived.csv")
    print("  - geometry_metrics_clean.csv")
    print("  - cognate_vs_noncognate_effects.csv")
    print("  - heatmap_stratified_conditions.png")
    print("  - cognate_vs_noncognate_effects.png")
    print("  - zn_diagnostics.png")

    print("\n" + "="*70)
    print("INTERPRETATION RULES")
    print("="*70)
    print("✓ pocket_iptm = AF3 local support (NOT binding energy)")
    print("✓ contacts_per_atom = shape complementarity proxy")
    print("✓ polar_contacts_per_atom = polar engagement proxy")
    print("✓ Use Zn results ONLY in zn_engaged subset")
    print("✓ Treat zn_floating as AF3 placement artifact")
    print("="*70)

if __name__ == '__main__':
    main()
