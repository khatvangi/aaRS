#!/usr/bin/env python3
"""
Generate figures for Paper 2: Ancestral aaRS Promiscuity
Outputs: SVG, PDF, PNG formats + data tables
"""

import json
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from pathlib import Path

# Set publication-quality defaults
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 9
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.bbox'] = 'tight'

# Color palette (colorblind-friendly)
COLORS = {
    'cognate': '#2E7D32',      # Dark green
    'non_cognate': '#D32F2F',  # Dark red
    'control': '#757575',      # Gray
    'luca': '#1976D2',         # Blue
    'eukaryotic': '#388E3C',   # Green
    'modern': '#7B1FA2',       # Purple
    'domain': '#3498db',       # Light blue
    'fulllength': '#e74c3c',   # Red
    'archaea': '#e74c3c',      # Red
    'bacteria': '#3498db',     # Blue
    'eukaryota': '#27ae60',    # Green
}

# Output directory
OUTPUT_DIR = Path('/storage/kiran-stuff/aaRS/figures')
OUTPUT_DIR.mkdir(exist_ok=True)

def extract_iptm_from_json(json_path):
    """Extract ipTM and other metrics from AF3 summary JSON.

    IMPORTANT: Uses chain_pair_iptm[0][2] for protein-ligand binding,
    NOT the overall iptm which measures whole complex structure.
    """
    with open(json_path) as f:
        data = json.load(f)

    # Extract protein-ligand interface confidence from chain_pair_iptm matrix
    # Matrix is [protein, tRNA/other, ligand] where indices are [0, 1, 2]
    # chain_pair_iptm[0][2] = protein-ligand interaction score
    chain_pair_iptm = data.get('chain_pair_iptm', [[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    protein_ligand_iptm = chain_pair_iptm[0][2] if len(chain_pair_iptm) > 0 and len(chain_pair_iptm[0]) > 2 else 0

    return {
        'iptm': protein_ligand_iptm,  # Use protein-ligand binding score
        'iptm_overall': data.get('iptm', 0),  # Keep overall for reference
        'ptm': data.get('ptm', 0),
        'chain_ptm': data.get('chain_ptm', []),
        'fraction_disordered': data.get('fraction_disordered', 0),
        'has_clash': data.get('has_clash', False),
    }

def find_best_sample(model_dir, model_name):
    """Find the best sample based on AF3 ranking_score (overall model quality)."""
    best_ranking = -1
    best_path = None

    # Check all samples
    for sample_num in range(5):
        sample_path = model_dir / f'seed-1_sample-{sample_num}' / f'{model_name}_seed-1_sample-{sample_num}_summary_confidences.json'
        if sample_path.exists():
            with open(sample_path) as f:
                data = json.load(f)
            ranking = data.get('ranking_score', 0)
            if ranking > best_ranking:
                best_ranking = ranking
                best_path = sample_path

    return best_path

def collect_all_iptm_data():
    """Collect ipTM values from all models using BEST sample by ranking."""
    base_path = Path('/storage/kiran-stuff/aaRS/phase2')

    # Define all models with their categories
    models = {
        # Domain models - LUCA ProRS
        'deep_domain_pro': {'enzyme': 'LUCA ProRS', 'ligand': 'PRO', 'type': 'domain', 'category': 'catalytic'},
        'deep_domain_thr': {'enzyme': 'LUCA ProRS', 'ligand': 'THR', 'type': 'domain', 'category': 'catalytic'},
        'deep_catalytic_trp': {'enzyme': 'LUCA ProRS', 'ligand': 'TRP', 'type': 'domain', 'category': 'control'},
        'deep_catalytic_phe': {'enzyme': 'LUCA ProRS', 'ligand': 'PHE', 'type': 'domain', 'category': 'control'},
        'deep_cat_trp': {'enzyme': 'LUCA ProRS', 'ligand': 'TRP', 'type': 'domain', 'category': 'control_alt'},
        'deep_cat_phe': {'enzyme': 'LUCA ProRS', 'ligand': 'PHE', 'type': 'domain', 'category': 'control_alt'},

        # Editing domain
        'deep_editing_pro': {'enzyme': 'LUCA ProRS Editing', 'ligand': 'PRO', 'type': 'domain', 'category': 'editing'},
        'deep_editing_thr': {'enzyme': 'LUCA ProRS Editing', 'ligand': 'THR', 'type': 'domain', 'category': 'editing'},

        # LUCA ThrRS
        'deep_thrrs_pro': {'enzyme': 'LUCA ThrRS', 'ligand': 'PRO', 'type': 'domain', 'category': 'catalytic'},
        'deep_thrrs_thr': {'enzyme': 'LUCA ThrRS', 'ligand': 'THR', 'type': 'domain', 'category': 'catalytic'},

        # Eukaryotic ProRS
        'shallow_domain_pro': {'enzyme': 'Eukaryotic ProRS', 'ligand': 'PRO', 'type': 'domain', 'category': 'catalytic'},
        'shallow_domain_thr': {'enzyme': 'Eukaryotic ProRS', 'ligand': 'THR', 'type': 'domain', 'category': 'catalytic'},

        # Modern enzymes
        'modern_prours_pro': {'enzyme': 'Modern ProRS', 'ligand': 'PRO', 'type': 'domain', 'category': 'modern'},
        'modern_prours_thr': {'enzyme': 'Modern ProRS', 'ligand': 'THR', 'type': 'domain', 'category': 'modern'},
        'modern_thrrs_pro': {'enzyme': 'Modern ThrRS', 'ligand': 'PRO', 'type': 'domain', 'category': 'modern'},
        'modern_thrrs_thr': {'enzyme': 'Modern ThrRS', 'ligand': 'THR', 'type': 'domain', 'category': 'modern'},

        # Full-length models
        'fulllength_deep_pro': {'enzyme': 'LUCA ProRS', 'ligand': 'PRO', 'type': 'fulllength', 'category': 'catalytic'},
        'fulllength_deep_thr': {'enzyme': 'LUCA ProRS', 'ligand': 'THR', 'type': 'fulllength', 'category': 'catalytic'},
        'fulllength_shallow_pro': {'enzyme': 'Eukaryotic ProRS', 'ligand': 'PRO', 'type': 'fulllength', 'category': 'catalytic'},
        'fulllength_shallow_thr': {'enzyme': 'Eukaryotic ProRS', 'ligand': 'THR', 'type': 'fulllength', 'category': 'catalytic'},
    }

    results = []

    for model_name, info in models.items():
        # Try different directory patterns
        dirs_to_try = [
            base_path / 'outputs' / model_name,
            base_path / 'outputs' / model_name / model_name,
            base_path / 'af3_output_full' / model_name,
        ]

        found = False
        for model_dir in dirs_to_try:
            if model_dir.exists():
                # Find best sample
                best_path = find_best_sample(model_dir, model_name)
                if best_path:
                    metrics = extract_iptm_from_json(best_path)
                    results.append({
                        'model': model_name,
                        'enzyme': info['enzyme'],
                        'ligand': info['ligand'],
                        'type': info['type'],
                        'category': info['category'],
                        'iptm': metrics['iptm'],
                        'ptm': metrics['ptm'],
                        'fraction_disordered': metrics['fraction_disordered'],
                        'has_clash': metrics['has_clash'],
                        'sample_path': str(best_path),
                    })
                    found = True
                    break

        if not found:
            print(f"Warning: Could not find JSON for {model_name}")

    return pd.DataFrame(results)

def save_figure(fig, name, formats=['svg', 'pdf', 'png']):
    """Save figure in multiple formats."""
    for fmt in formats:
        filepath = OUTPUT_DIR / f'{name}.{fmt}'
        fig.savefig(filepath, format=fmt, bbox_inches='tight', dpi=300)
        print(f"Saved: {filepath}")

def create_figure_2b(df):
    """Figure 2B: ipTM Score Comparison (Domain-Level)"""

    # Filter for domain-level catalytic models
    domain_df = df[(df['type'] == 'domain') & (df['category'].isin(['catalytic', 'control', 'editing']))].copy()

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    # Define enzyme order and colors
    enzymes = ['LUCA ProRS', 'LUCA ThrRS', 'Eukaryotic ProRS', 'Modern ProRS', 'Modern ThrRS']
    ligand_colors = {'PRO': COLORS['cognate'], 'THR': COLORS['non_cognate'],
                     'TRP': COLORS['control'], 'PHE': '#9E9E9E'}

    # Plot grouped bars
    x = np.arange(len(enzymes))
    width = 0.2

    for i, ligand in enumerate(['PRO', 'THR', 'TRP', 'PHE']):
        values = []
        for enzyme in enzymes:
            mask = (domain_df['enzyme'] == enzyme) & (domain_df['ligand'] == ligand)
            if mask.any():
                values.append(domain_df[mask]['iptm'].values[0])
            else:
                values.append(0)

        offset = (i - 1.5) * width
        bars = ax.bar(x + offset, values, width, label=ligand, color=ligand_colors[ligand],
                      edgecolor='black', linewidth=0.5)

        # Add value labels on bars
        for bar, val in zip(bars, values):
            if val > 0:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                       f'{val:.2f}', ha='center', va='bottom', fontsize=7)

    # Customize plot
    ax.set_xlabel('Enzyme')
    ax.set_ylabel('ipTM Score')
    ax.set_title('Domain-Level Substrate Binding Affinity')
    ax.set_xticks(x)
    ax.set_xticklabels(enzymes, rotation=15, ha='right')
    ax.set_ylim(0, 1.1)
    ax.legend(title='Ligand', loc='upper right')
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.3)

    # Add grid
    ax.yaxis.grid(True, alpha=0.3)
    ax.set_axisbelow(True)

    plt.tight_layout()
    save_figure(fig, 'figure_2b_domain_iptm_comparison')
    plt.close()

    return domain_df

def create_figure_3b(df):
    """Figure 3B: Domain vs Full-Length Comparison (KEY FIGURE)"""

    # Calculate promiscuity ratios
    data = []

    for enzyme_base in ['LUCA ProRS', 'Eukaryotic ProRS']:
        for model_type in ['domain', 'fulllength']:
            pro_mask = (df['enzyme'] == enzyme_base) & (df['ligand'] == 'PRO') & (df['type'] == model_type)
            thr_mask = (df['enzyme'] == enzyme_base) & (df['ligand'] == 'THR') & (df['type'] == model_type)

            if pro_mask.any() and thr_mask.any():
                pro_iptm = df[pro_mask]['iptm'].values[0]
                thr_iptm = df[thr_mask]['iptm'].values[0]
                ratio = (thr_iptm / pro_iptm * 100) if pro_iptm > 0 else 0

                data.append({
                    'Enzyme': enzyme_base,
                    'Type': 'Domain' if model_type == 'domain' else 'Full-length',
                    'Ratio': ratio,
                    'PRO_iptm': pro_iptm,
                    'THR_iptm': thr_iptm,
                })

    ratio_df = pd.DataFrame(data)

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot grouped bars
    enzymes = ['LUCA ProRS', 'Eukaryotic ProRS']
    x = np.arange(len(enzymes))
    width = 0.35

    domain_vals = [ratio_df[(ratio_df['Enzyme'] == e) & (ratio_df['Type'] == 'Domain')]['Ratio'].values[0]
                   for e in enzymes]
    full_vals = [ratio_df[(ratio_df['Enzyme'] == e) & (ratio_df['Type'] == 'Full-length')]['Ratio'].values[0]
                 for e in enzymes]

    bars1 = ax.bar(x - width/2, domain_vals, width, label='Domain-only',
                   color=COLORS['domain'], edgecolor='black', linewidth=1)
    bars2 = ax.bar(x + width/2, full_vals, width, label='Full-length',
                   color=COLORS['fulllength'], edgecolor='black', linewidth=1)

    # Add value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2, height + 1,
                   f'{height:.1f}%', ha='center', va='bottom', fontweight='bold', fontsize=10)

    # Add arrows showing increase
    for i, (d, f) in enumerate(zip(domain_vals, full_vals)):
        diff = f - d
        mid_x = x[i]
        ax.annotate(f'+{diff:.1f}pp',
                   xy=(mid_x, max(d, f) + 8),
                   ha='center', fontsize=9, color='darkred', fontweight='bold')

    # Reference line at 90%
    ax.axhline(y=90, color='gray', linestyle='--', alpha=0.5, label='90% threshold')

    # Customize plot
    ax.set_xlabel('Enzyme')
    ax.set_ylabel('THR Binding (% of PRO)')
    ax.set_title('Full-Length Context Enhances Promiscuity', fontweight='bold', fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(enzymes)
    ax.set_ylim(0, 110)
    ax.legend(loc='lower right')

    # Add grid
    ax.yaxis.grid(True, alpha=0.3)
    ax.set_axisbelow(True)

    # Add text box with key finding
    textstr = 'Full-length context\nENHANCES promiscuity'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', bbox=props, fontweight='bold')

    plt.tight_layout()
    save_figure(fig, 'figure_3b_domain_vs_fulllength')
    plt.close()

    return ratio_df

def create_figure_5b(df):
    """Figure 5B: Binding Affinity Heatmap"""

    # Filter for catalytic domain models
    heatmap_df = df[df['category'].isin(['catalytic', 'control', 'modern'])].copy()

    # Create pivot table
    enzymes = ['LUCA ProRS', 'LUCA ThrRS', 'Eukaryotic ProRS', 'Modern ProRS', 'Modern ThrRS']
    ligands = ['PRO', 'THR', 'TRP', 'PHE']

    # Build matrix
    matrix = np.zeros((len(enzymes), len(ligands)))
    for i, enzyme in enumerate(enzymes):
        for j, ligand in enumerate(ligands):
            mask = (heatmap_df['enzyme'] == enzyme) & (heatmap_df['ligand'] == ligand) & (heatmap_df['type'] == 'domain')
            if mask.any():
                matrix[i, j] = heatmap_df[mask]['iptm'].values[0]
            else:
                matrix[i, j] = np.nan

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))

    # Create heatmap
    cmap = LinearSegmentedColormap.from_list('custom', ['white', '#FFCDD2', '#EF5350', '#C62828'])

    im = ax.imshow(matrix, cmap=cmap, aspect='auto', vmin=0, vmax=1)

    # Add colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('ipTM Score', rotation=-90, va="bottom")

    # Add text annotations
    for i in range(len(enzymes)):
        for j in range(len(ligands)):
            val = matrix[i, j]
            if not np.isnan(val):
                text_color = 'white' if val > 0.5 else 'black'
                ax.text(j, i, f'{val:.2f}', ha='center', va='center',
                       color=text_color, fontsize=10, fontweight='bold')

    # Customize
    ax.set_xticks(np.arange(len(ligands)))
    ax.set_yticks(np.arange(len(enzymes)))
    ax.set_xticklabels(ligands)
    ax.set_yticklabels(enzymes)
    ax.set_xlabel('Ligand')
    ax.set_ylabel('Enzyme')
    ax.set_title('Substrate Binding Affinity Matrix', fontweight='bold')

    # Add borders for cognate pairs
    # ProRS cognate = PRO (column 0), ThrRS cognate = THR (column 1)
    for i, enzyme in enumerate(enzymes):
        if 'ProRS' in enzyme:
            rect = plt.Rectangle((0-0.5, i-0.5), 1, 1, fill=False,
                                 edgecolor='green', linewidth=3)
            ax.add_patch(rect)
        elif 'ThrRS' in enzyme:
            rect = plt.Rectangle((1-0.5, i-0.5), 1, 1, fill=False,
                                 edgecolor='green', linewidth=3)
            ax.add_patch(rect)

    plt.tight_layout()
    save_figure(fig, 'figure_5b_binding_heatmap')
    plt.close()

    return matrix

def create_figure_4(df):
    """Figure 4: Evolutionary Timeline"""

    # Calculate promiscuity values at each time point
    timeline_data = []

    # LUCA (3.5 Gya)
    luca_pro = df[(df['enzyme'] == 'LUCA ProRS') & (df['ligand'] == 'PRO') & (df['type'] == 'fulllength')]['iptm']
    luca_thr = df[(df['enzyme'] == 'LUCA ProRS') & (df['ligand'] == 'THR') & (df['type'] == 'fulllength')]['iptm']
    if len(luca_pro) > 0 and len(luca_thr) > 0:
        luca_ratio = luca_thr.values[0] / luca_pro.values[0] * 100
        timeline_data.append({'time': 3.5, 'label': 'LUCA', 'ratio': luca_ratio})

    # Eukaryotic (1.5 Gya)
    euk_pro = df[(df['enzyme'] == 'Eukaryotic ProRS') & (df['ligand'] == 'PRO') & (df['type'] == 'fulllength')]['iptm']
    euk_thr = df[(df['enzyme'] == 'Eukaryotic ProRS') & (df['ligand'] == 'THR') & (df['type'] == 'fulllength')]['iptm']
    if len(euk_pro) > 0 and len(euk_thr) > 0:
        euk_ratio = euk_thr.values[0] / euk_pro.values[0] * 100
        timeline_data.append({'time': 1.5, 'label': 'Eukaryotic\nAncestor', 'ratio': euk_ratio})

    # Modern (0 Gya)
    mod_pro = df[(df['enzyme'] == 'Modern ProRS') & (df['ligand'] == 'PRO') & (df['type'] == 'domain')]['iptm']
    mod_thr = df[(df['enzyme'] == 'Modern ProRS') & (df['ligand'] == 'THR') & (df['type'] == 'domain')]['iptm']
    if len(mod_pro) > 0 and len(mod_thr) > 0:
        mod_ratio = mod_thr.values[0] / mod_pro.values[0] * 100
        timeline_data.append({'time': 0, 'label': 'Modern', 'ratio': mod_ratio})

    timeline_df = pd.DataFrame(timeline_data)

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 5))

    # Plot timeline
    times = timeline_df['time'].values
    ratios = timeline_df['ratio'].values
    labels = timeline_df['label'].values

    # Background gradient
    gradient = np.linspace(0, 1, 256).reshape(1, -1)
    ax.imshow(gradient, extent=[3.8, -0.3, 0, 120], aspect='auto',
             cmap='YlOrRd_r', alpha=0.2)

    # Plot line and points
    ax.plot(times, ratios, 'ko-', markersize=15, linewidth=3, markerfacecolor='white',
           markeredgewidth=2)

    # Add ratio labels
    for t, r, l in zip(times, ratios, labels):
        ax.annotate(f'{r:.1f}%', xy=(t, r), xytext=(0, 15),
                   textcoords='offset points', ha='center', fontsize=12, fontweight='bold')
        ax.annotate(l, xy=(t, r), xytext=(0, -30),
                   textcoords='offset points', ha='center', fontsize=10)

    # Customize
    ax.set_xlim(4, -0.5)
    ax.set_ylim(40, 110)
    ax.set_xlabel('Billions of Years Ago', fontsize=12)
    ax.set_ylabel('THR Binding (% of PRO)', fontsize=12)
    ax.set_title('Promiscuity Persists Across Evolutionary Time', fontweight='bold', fontsize=14)

    # Add key events
    ax.axvline(x=3.5, color='gray', linestyle=':', alpha=0.5)
    ax.text(3.5, 45, 'Origin\nof Life', ha='center', fontsize=8, alpha=0.7)

    ax.axvline(x=2.0, color='gray', linestyle=':', alpha=0.5)
    ax.text(2.0, 45, 'GOE', ha='center', fontsize=8, alpha=0.7)

    ax.axhline(y=90, color='red', linestyle='--', alpha=0.3)
    ax.text(-0.2, 91, '90%', fontsize=8, color='red', alpha=0.7)

    # Add text box
    textstr = 'Minimal improvement\nover 3.5 billion years'
    props = dict(boxstyle='round', facecolor='lightyellow', alpha=0.8)
    ax.text(0.75, 0.15, textstr, transform=ax.transAxes, fontsize=10,
           verticalalignment='bottom', bbox=props)

    plt.tight_layout()
    save_figure(fig, 'figure_4_evolutionary_timeline')
    plt.close()

    return timeline_df

def create_figure_1ab():
    """Figure 1A/B: Phylogenetic Trees"""
    try:
        from ete3 import Tree, TreeStyle, NodeStyle, TextFace, CircleFace
    except ImportError:
        print("ete3 not available. Creating simplified tree visualization.")
        return create_simplified_trees()

    tree_files = {
        'ProRS': '/storage/kiran-stuff/aaRS/phase1b/results/ProRS_deep.treefile',
        'ThrRS': '/storage/kiran-stuff/aaRS/phase1b/results/ThrRS_deep.treefile',
    }

    for name, tree_file in tree_files.items():
        if not os.path.exists(tree_file):
            print(f"Tree file not found: {tree_file}")
            continue

        # Load tree
        t = Tree(tree_file)

        # Color nodes by domain
        archaea_keywords = ['Methan', 'Pyro', 'Thermo', 'Halo', 'Sulfo', 'Archaeo']
        bacteria_keywords = ['Escherichia', 'Bacillus', 'Pseudo', 'Strepto', 'Salmonella', 'Clostrid']
        eukaryota_keywords = ['Homo', 'Mus', 'Saccharomyces', 'Drosophila', 'Arabidopsis', 'Xenopus']

        for node in t.traverse():
            nstyle = NodeStyle()
            if node.is_leaf():
                node_name = node.name.lower()
                if any(kw.lower() in node_name for kw in archaea_keywords):
                    nstyle['fgcolor'] = COLORS['archaea']
                elif any(kw.lower() in node_name for kw in bacteria_keywords):
                    nstyle['fgcolor'] = COLORS['bacteria']
                elif any(kw.lower() in node_name for kw in eukaryota_keywords):
                    nstyle['fgcolor'] = COLORS['eukaryota']
                else:
                    nstyle['fgcolor'] = 'gray'
            node.set_style(nstyle)

        # Tree style
        ts = TreeStyle()
        ts.show_leaf_name = False
        ts.mode = 'r'
        ts.scale = 100

        # Render
        output_path = OUTPUT_DIR / f'figure_1{"a" if name == "ProRS" else "b"}_{name.lower()}_tree'
        t.render(str(output_path) + '.png', tree_style=ts, w=800)
        t.render(str(output_path) + '.svg', tree_style=ts, w=800)
        print(f"Saved: {output_path}")

def create_simplified_trees():
    """Create simplified tree visualization using matplotlib when ete3 not available."""

    tree_files = {
        'ProRS': '/storage/kiran-stuff/aaRS/phase1b/results/ProRS_deep.treefile',
        'ThrRS': '/storage/kiran-stuff/aaRS/phase1b/results/ThrRS_deep.treefile',
    }

    for panel, (name, tree_file) in enumerate(tree_files.items()):
        if not os.path.exists(tree_file):
            print(f"Tree file not found: {tree_file}")
            continue

        # Read tree file to count species
        with open(tree_file) as f:
            tree_str = f.read()

        # Count tips (approximate by counting commas + 1)
        n_tips = tree_str.count(',') + 1

        # Create schematic tree figure
        fig, ax = plt.subplots(figsize=(8, 8))

        # Draw schematic tree
        # Main trunk
        ax.plot([0.5, 0.5], [0, 0.3], 'k-', linewidth=3)

        # Three main clades
        colors = [COLORS['archaea'], COLORS['bacteria'], COLORS['eukaryota']]
        labels = ['Archaea', 'Bacteria', 'Eukaryota']
        x_positions = [0.2, 0.5, 0.8]

        for i, (x, color, label) in enumerate(zip(x_positions, colors, labels)):
            # Branch to clade
            ax.plot([0.5, x], [0.3, 0.5], 'k-', linewidth=2)
            # Clade expansion
            for j in range(5):
                y_end = 0.6 + j * 0.08
                x_var = x + (j - 2) * 0.03
                ax.plot([x, x_var], [0.5, y_end], '-', color=color, linewidth=1, alpha=0.7)
            # Label
            ax.text(x, 0.98, label, ha='center', fontsize=10, color=color, fontweight='bold')

        # Mark LUCA node
        ax.plot(0.5, 0.3, 'k*', markersize=20)
        ax.text(0.52, 0.28, 'LUCA', fontsize=9, fontweight='bold')

        # Title and info
        ax.set_title(f'{name} Phylogenetic Tree ({n_tips} species)', fontweight='bold', fontsize=12)
        ax.text(0.5, 0.05, 'Schematic representation', ha='center', fontsize=8, style='italic', alpha=0.7)

        # Remove axes
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')

        # Save
        panel_letter = 'a' if name == 'ProRS' else 'b'
        save_figure(fig, f'figure_1{panel_letter}_{name.lower()}_tree_schematic')
        plt.close()

def create_editing_domain_figure(df):
    """Create supplementary figure showing editing domain results."""

    editing_df = df[df['category'] == 'editing'].copy()

    if editing_df.empty:
        print("No editing domain data found")
        return

    fig, ax = plt.subplots(figsize=(6, 4))

    ligands = ['PRO', 'THR']
    values = [editing_df[editing_df['ligand'] == l]['iptm'].values[0] if len(editing_df[editing_df['ligand'] == l]) > 0 else 0
              for l in ligands]

    colors = [COLORS['cognate'], COLORS['non_cognate']]
    bars = ax.bar(ligands, values, color=colors, edgecolor='black', linewidth=1)

    # Add value labels
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
               f'{val:.2f}', ha='center', va='bottom', fontweight='bold')

    # Calculate ratio
    if values[0] > 0:
        ratio = values[1] / values[0] * 100
        ax.text(0.5, max(values) + 0.08, f'THR = {ratio:.0f}% of PRO',
               ha='center', fontsize=10, fontweight='bold', color='darkred')

    ax.set_ylabel('ipTM Score')
    ax.set_title('LUCA ProRS Editing Domain\nFails to Discriminate', fontweight='bold')
    ax.set_ylim(0, max(values) * 1.3 if values else 1)

    plt.tight_layout()
    save_figure(fig, 'figure_2d_editing_domain')
    plt.close()

def create_model_quality_table(df):
    """Create supplementary table with all model quality metrics."""

    quality_df = df[['model', 'enzyme', 'ligand', 'type', 'iptm', 'ptm',
                     'fraction_disordered', 'has_clash']].copy()
    quality_df = quality_df.round(3)

    # Save as CSV
    csv_path = OUTPUT_DIR / 'table_s3_model_quality_metrics.csv'
    quality_df.to_csv(csv_path, index=False)
    print(f"Saved: {csv_path}")

    return quality_df

def main():
    """Main function to generate all figures."""

    print("="*60)
    print("Generating Paper 2 Figures: Ancestral aaRS Promiscuity")
    print("="*60)

    # Collect all data
    print("\n1. Collecting ipTM data from all models...")
    df = collect_all_iptm_data()

    # Save master data table
    master_table_path = OUTPUT_DIR / 'table_master_iptm_data.csv'
    df.to_csv(master_table_path, index=False)
    print(f"Saved master data: {master_table_path}")
    print(f"   Total models found: {len(df)}")

    # Generate figures
    print("\n2. Creating Figure 3B (Domain vs Full-Length)...")
    ratio_df = create_figure_3b(df)
    ratio_df.to_csv(OUTPUT_DIR / 'table_figure_3b_data.csv', index=False)

    print("\n3. Creating Figure 2B (Domain ipTM Comparison)...")
    domain_df = create_figure_2b(df)

    print("\n4. Creating Figure 5B (Binding Heatmap)...")
    matrix = create_figure_5b(df)

    print("\n5. Creating Figure 4 (Evolutionary Timeline)...")
    timeline_df = create_figure_4(df)
    timeline_df.to_csv(OUTPUT_DIR / 'table_figure_4_data.csv', index=False)

    print("\n6. Creating Figure 1A/B (Phylogenetic Trees)...")
    create_simplified_trees()

    print("\n7. Creating Figure 2D (Editing Domain)...")
    create_editing_domain_figure(df)

    print("\n8. Creating model quality table...")
    create_model_quality_table(df)

    print("\n" + "="*60)
    print("Figure generation complete!")
    print(f"Output directory: {OUTPUT_DIR}")
    print("="*60)

if __name__ == '__main__':
    main()
