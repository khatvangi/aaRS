#!/usr/bin/env python3
"""
plot_competition_heatmap.py

Generate a heatmap visualization of competition tournament results.
Rows = Protein constructs, Columns = Competitor ligands
Color = Cognate vs Competitor margin (positive = cognate wins)
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Load data
df = pd.read_csv("/storage/kiran-stuff/aaRS/phase3_balanced/competition_tournament.csv")

print("Competition Tournament Data:")
print(df[["protein_construct", "enzyme", "era", "cognate_ligand", "competitor_ligand",
          "cognate_iptm", "competitor_iptm", "margin", "winner"]].to_string(index=False))

# Create a pivot table for the heatmap
# Rows = protein construct, Columns = competitor ligand
# Value = margin (cognate_iptm - competitor_iptm)

# Create construct label
df["construct_label"] = df["era"].str.capitalize() + " " + df["enzyme"]

# Pivot
pivot = df.pivot_table(
    index="construct_label",
    columns="competitor_ligand",
    values="margin",
    aggfunc="first"  # Only one value per cell expected
)

print("\nPivot table (margins):")
print(pivot)

# Sort constructs by enzyme then era
construct_order = ["Ancestral ProRS", "Modern ProRS", "Ancestral ThrRS", "Modern ThrRS"]
construct_order = [c for c in construct_order if c in pivot.index]
pivot = pivot.reindex(construct_order)

# Sort competitors
competitor_order = ["ALA", "GLY", "SER", "THR", "PRO", "VAL", "ILE", "LEU", "MET", "CYS", "TYR", "GLU"]
competitor_order = [c for c in competitor_order if c in pivot.columns]
pivot = pivot[competitor_order]

print("\nFinal pivot (sorted):")
print(pivot)

# Create figure
fig, ax = plt.subplots(figsize=(10, 5))

# Create mask for NaN values
mask = pivot.isna()

# Plot heatmap
im = ax.imshow(
    np.where(mask, 0, pivot.values),
    cmap="RdBu",  # Red = competitor wins, Blue = cognate wins
    vmin=-0.5, vmax=0.5,
    aspect="auto"
)

# Set axis labels
ax.set_xticks(range(len(pivot.columns)))
ax.set_xticklabels(pivot.columns, fontsize=11)
ax.set_yticks(range(len(pivot.index)))
ax.set_yticklabels(pivot.index, fontsize=11)

# Add colorbar
cbar = plt.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label("Cognate - Competitor ipTM\n(Blue = Cognate Wins)", fontsize=10)

# Add text annotations
for i in range(len(pivot.index)):
    for j in range(len(pivot.columns)):
        val = pivot.iloc[i, j]
        if pd.notna(val):
            color = "white" if abs(val) > 0.2 else "black"
            symbol = "✓" if val > 0 else "✗"
            ax.text(j, i, f"{symbol}\n{val:.2f}", ha="center", va="center",
                    fontsize=10, fontweight="bold", color=color)
        else:
            ax.text(j, i, "—", ha="center", va="center",
                    fontsize=14, color="gray")

# Title
ax.set_title("AF3 Competition Tournament: Cognate vs Non-Cognate\n"
             "(✓ = Cognate wins pocket, ✗ = Competitor wins)",
             fontsize=12, fontweight="bold")

ax.set_xlabel("Competitor Ligand", fontsize=11)
ax.set_ylabel("Protein Construct", fontsize=11)

plt.tight_layout()

# Save
outpath = Path("/storage/kiran-stuff/aaRS/phase3_balanced/FIG_competition_tournament.png")
plt.savefig(outpath, dpi=150, bbox_inches="tight")
print(f"\nSaved figure to {outpath}")

# Also save PDF
plt.savefig(outpath.with_suffix(".pdf"), bbox_inches="tight")
print(f"Saved PDF to {outpath.with_suffix('.pdf')}")

plt.close()

# Print summary stats
print("\n" + "="*60)
print("SUMMARY STATISTICS")
print("="*60)
print(f"\nTotal matchups: {len(df)}")
print(f"Cognate wins: {df['cognate_wins'].sum()} ({df['cognate_wins'].mean():.0%})")
print(f"Mean margin when cognate wins: {df[df['winner']=='cognate']['margin'].mean():.3f}")

# Per-enzyme summary
print("\nPer-enzyme:")
for enzyme in df["enzyme"].unique():
    subset = df[df["enzyme"] == enzyme]
    print(f"  {enzyme}: {subset['cognate_wins'].sum()}/{len(subset)} wins "
          f"(mean Δ = {subset['margin'].mean():.3f})")
