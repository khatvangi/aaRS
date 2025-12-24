#!/usr/bin/env python3
"""
Summarize Zn engagement rules by condition.
"""
import pandas as pd
import numpy as np
from pathlib import Path

def main():
    print("="*70)
    print("ZN ENGAGEMENT RULES REPORT")
    print("="*70)

    # Load data - prefer recomputed if available
    if Path('qc/recomputed_run1.csv').exists():
        print("\nLoading qc/recomputed_run1.csv...")
        df = pd.read_csv('qc/recomputed_run1.csv')
    elif Path('geometry_metrics_recomputed.csv').exists():
        print("\nLoading geometry_metrics_recomputed.csv...")
        df = pd.read_csv('geometry_metrics_recomputed.csv')
    else:
        print("\nLoading geometry_metrics_clean.csv...")
        df = pd.read_csv('geometry_metrics_clean.csv')

    # Filter to Zn-containing structures
    df_zn = df[df['zn_present'] == True].copy()

    print(f"  Total structures: {len(df)}")
    print(f"  Zn-containing: {len(df_zn)}")

    if len(df_zn) == 0:
        print("\nNo Zn-containing structures found!")
        return

    # Summarize by condition
    print("\n" + "="*70)
    print("ZN ENGAGEMENT BY CONDITION")
    print("="*70)

    summary_rows = []

    for condition in sorted(df_zn['condition'].unique()):
        df_cond = df_zn[df_zn['condition'] == condition]

        n_total = len(df_cond)
        n_engaged = df_cond['zn_engaged'].sum()
        n_floating = df_cond['zn_floating'].sum()

        # Zn distance stats
        zn_dists = df_cond['zn_min_dist_hetero'].dropna()

        if len(zn_dists) > 0:
            mean_dist = zn_dists.mean()
            std_dist = zn_dists.std()
            min_dist = zn_dists.min()
            max_dist = zn_dists.max()
        else:
            mean_dist = std_dist = min_dist = max_dist = None

        summary_rows.append({
            'condition': condition,
            'n_total': n_total,
            'n_engaged': n_engaged,
            'n_floating': n_floating,
            'pct_engaged': 100 * n_engaged / n_total if n_total > 0 else 0,
            'mean_zn_dist': mean_dist,
            'std_zn_dist': std_dist,
            'min_zn_dist': min_dist,
            'max_zn_dist': max_dist
        })

    df_summary = pd.DataFrame(summary_rows)
    df_summary = df_summary.sort_values('pct_engaged', ascending=False)

    # Save
    df_summary.to_csv('zn/zn_engagement_by_condition.csv', index=False)
    print("\nSaved zn/zn_engagement_by_condition.csv")

    print("\nEngagement by condition:")
    print(df_summary.to_string(index=False))

    # Overall stats
    print("\n" + "="*70)
    print("OVERALL ZN STATISTICS")
    print("="*70)

    n_engaged_total = df_zn['zn_engaged'].sum()
    n_floating_total = df_zn['zn_floating'].sum()

    print(f"\nTotal Zn structures: {len(df_zn)}")
    print(f"  Engaged (≤3.0 Å): {n_engaged_total} ({100*n_engaged_total/len(df_zn):.1f}%)")
    print(f"  Floating (>3.0 Å): {n_floating_total} ({100*n_floating_total/len(df_zn):.1f}%)")

    # Distance distribution
    zn_dists_all = df_zn['zn_min_dist_hetero'].dropna()
    if len(zn_dists_all) > 0:
        print(f"\nZn-ligand distance (heteroatoms):")
        print(f"  Mean: {zn_dists_all.mean():.2f} Å")
        print(f"  Std:  {zn_dists_all.std():.2f} Å")
        print(f"  Min:  {zn_dists_all.min():.2f} Å")
        print(f"  Max:  {zn_dists_all.max():.2f} Å")
        print(f"  Median: {zn_dists_all.median():.2f} Å")

    # Engaged vs floating comparison
    if n_engaged_total > 0 and n_floating_total > 0:
        print("\n" + "="*70)
        print("ENGAGED vs FLOATING COMPARISON")
        print("="*70)

        df_engaged = df_zn[df_zn['zn_engaged'] == True]
        df_floating = df_zn[df_zn['zn_floating'] == True]

        metrics = ['pocket_iptm', 'contacts_per_atom', 'polar_contacts_per_atom', 'clash_rate']

        for metric in metrics:
            if metric not in df.columns:
                continue

            engaged_vals = df_engaged[metric].dropna()
            floating_vals = df_floating[metric].dropna()

            if len(engaged_vals) > 0 and len(floating_vals) > 0:
                print(f"\n{metric}:")
                print(f"  Engaged:  {engaged_vals.mean():.3f} ± {engaged_vals.std():.3f}")
                print(f"  Floating: {floating_vals.mean():.3f} ± {floating_vals.std():.3f}")
                print(f"  Δ (engaged - floating): {engaged_vals.mean() - floating_vals.mean():+.3f}")

    # Write markdown report
    print("\n" + "="*70)
    print("WRITING RULES DOCUMENT")
    print("="*70)

    md_content = f"""# Zn Engagement Rules

**Cutoff:** zn_min_dist_hetero ≤ 3.0 Å (Zn to ligand O/N/S atoms)

---

## Classification

| Category | Count | Percentage |
|----------|-------|------------|
| **Engaged** | {n_engaged_total} | {100*n_engaged_total/len(df_zn):.1f}% |
| **Floating** | {n_floating_total} | {100*n_floating_total/len(df_zn):.1f}% |

---

## Distance Distribution

- **Mean:** {zn_dists_all.mean():.2f} ± {zn_dists_all.std():.2f} Å
- **Median:** {zn_dists_all.median():.2f} Å
- **Range:** [{zn_dists_all.min():.2f}, {zn_dists_all.max():.2f}] Å

---

## Engagement by Condition

"""

    md_content += df_summary[['condition', 'n_engaged', 'n_floating', 'pct_engaged', 'mean_zn_dist']].to_markdown(index=False)

    md_content += """

---

## Interpretation

### Engaged Zn (≤3.0 Å)
- Zn is within coordination distance of ligand heteroatoms
- Likely represents true Zn-ligand interaction
- **USE for mechanistic interpretation**

### Floating Zn (>3.0 Å)
- Zn is distant from ligand (mean ~10-15 Å)
- AF3 modeling artifact - Zn placed for protein but not ligand
- **EXCLUDE from biological claims**

---

## Recommendation

**Only use structures with `zn_engaged=True` for:**
- Claims about Zn-mediated binding
- Comparisons of Zn effect on selectivity
- Mechanistic interpretation of metal coordination

**Flag `zn_floating` structures as:**
- AF3 structural predictions with incomplete Zn coordination
- Not suitable for Zn-ligand binding analysis
"""

    with open('zn/zn_engagement_rules.md', 'w') as f:
        f.write(md_content)

    print("Saved zn/zn_engagement_rules.md")

    print("\n" + "="*70)
    print("ZN RULES REPORT COMPLETE")
    print("="*70)

if __name__ == '__main__':
    main()
