#!/usr/bin/env python3
"""
Generate graphical abstract figure for ProRS evolution.

Outputs:
- figures/fig_graphical_abstract_proRS_evolution.png
- figures/fig_graphical_abstract_proRS_evolution.svg
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch
from pathlib import Path

def draw_protein(ax, x, y, label, domains, color_main, color_secondary, pocket_label=None):
    """
    Draw a simple domain cartoon.
    domains: list of (width, height, offset_y) relative to (x, y)
    """
    for i, (w, h, dy) in enumerate(domains):
        color = color_main if i == 0 else color_secondary
        rect = Rectangle(
            (x + sum(d[0] for d in domains[:i]), y + dy),
            w,
            h,
            linewidth=1.0,
            edgecolor="black",
            facecolor=color,
        )
        ax.add_patch(rect)

    total_width = sum(d[0] for d in domains)
    pocket = Circle(
        (x + total_width * 0.6, y + 0.6),
        radius=0.25,
        linewidth=1.0,
        edgecolor="black",
        facecolor="none",
        alpha=0.8,
    )
    ax.add_patch(pocket)

    ax.text(x + total_width / 2, y - 0.2, label, ha="center", va="top", fontsize=9)

    if pocket_label:
        ax.text(x + total_width * 0.6, y + 1.0, pocket_label, ha="center", va="bottom", fontsize=7)


def draw_pro_thr_icons(ax, x, y):
    c_pro = Circle((x, y), radius=0.15, edgecolor="black", facecolor="#aaaaaa")
    ax.add_patch(c_pro)
    ax.text(x, y, "Pro", ha="center", va="center", fontsize=7)

    c_thr = Circle((x + 0.22, y + 0.05), radius=0.15, edgecolor="black", facecolor="#88aadd")
    ax.add_patch(c_thr)
    ax.text(x + 0.22, y + 0.05, "Thr", ha="center", va="center", fontsize=7)

    ax.text(x + 0.11, y - 0.25, "Pro/Thr\nconfusion", ha="center", va="top", fontsize=7)


def draw_chem_limit(ax, x, y):
    pro = Circle((x, y), radius=0.2, edgecolor="black", facecolor="#aaaaaa")
    ax.add_patch(pro)
    ax.text(x, y, "Pro", ha="center", va="center", fontsize=8)

    thr = Circle((x + 0.8, y), radius=0.2, edgecolor="black", facecolor="#88aadd")
    ax.add_patch(thr)
    ax.text(x + 0.8, y, "Thr", ha="center", va="center", fontsize=8)

    oh = Circle((x + 1.1, y + 0.25), radius=0.08, edgecolor="black", facecolor="#ff7777")
    ax.add_patch(oh)
    ax.text(x + 1.1, y + 0.25, "OH", ha="center", va="center", fontsize=6)

    arr1 = FancyArrowPatch((x + 0.8, y + 0.2), (x + 1.0, y + 0.22), arrowstyle="-|>", mutation_scale=10, linewidth=1.0)
    ax.add_patch(arr1)
    arr2 = FancyArrowPatch((x + 0.8, y - 0.2), (x + 1.05, y + 0.18), arrowstyle="-|>", mutation_scale=10, linewidth=1.0)
    ax.add_patch(arr2)

    ax.text(
        x + 0.4,
        y - 0.6,
        "Thr OH forms\ncompensatory H-bonds\n→ intrinsic limit",
        ha="left",
        va="top",
        fontsize=7,
    )


def main():
    gold = "#D4A017"
    muted_red = "#C06060"
    steel_blue = "#4A6FA5"
    light_gray = "#dddddd"

    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_axes([0.03, 0.25, 0.94, 0.65])
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4)
    ax.axis("off")

    ax.text(0.5, 3.8, "LUCA ProRS", fontsize=11, fontweight="bold", ha="center")
    ax.text(4.0, 3.8, "Ancestral editing-lineage", fontsize=11, fontweight="bold", ha="center")
    ax.text(7.5, 3.8, "Modern ProRS", fontsize=11, fontweight="bold", ha="center")

    draw_protein(
        ax,
        x=0.2,
        y=2.0,
        label="Compact, promiscuous",
        domains=[(1.0, 1.0, 0.0), (0.8, 0.8, 0.1)],
        color_main=gold,
        color_secondary=light_gray,
        pocket_label="Small active site",
    )
    draw_pro_thr_icons(ax, x=0.9, y=1.4)

    draw_protein(
        ax,
        x=3.4,
        y=2.0,
        label="Domain accretion\n+ editing module",
        domains=[(1.0, 1.0, 0.0), (0.8, 0.9, 0.0), (0.6, 0.8, 0.1)],
        color_main=muted_red,
        color_secondary=light_gray,
        pocket_label="Elaborated pocket",
    )
    draw_pro_thr_icons(ax, x=4.3, y=1.4)

    draw_protein(
        ax,
        x=6.6,
        y=2.0,
        label="More domains,\nmore PUs,\nhigher ambiguity",
        domains=[(1.0, 1.0, 0.0), (0.8, 1.1, -0.05), (0.8, 0.9, 0.05)],
        color_main=steel_blue,
        color_secondary=light_gray,
        pocket_label="Crowded pocket,\nmore contacts",
    )
    draw_pro_thr_icons(ax, x=7.5, y=1.4)

    arrow = FancyArrowPatch((1.8, 3.3), (9.0, 3.3), arrowstyle="-|>", mutation_scale=20, linewidth=1.5)
    ax.add_patch(arrow)
    ax.text(5.4, 3.45, "Evolutionary time →", ha="center", va="bottom", fontsize=9)

    draw_chem_limit(ax, x=6.5, y=0.8)
    ax.text(6.5, 1.6, "Chemical limit", fontsize=9, fontweight="bold", ha="left")

    ax2 = fig.add_axes([0.03, 0.05, 0.94, 0.13])
    ax2.axis("off")
    rect = Rectangle((0, 0), 1, 1, transform=ax2.transAxes, facecolor="#f5f5f5", edgecolor="none")
    ax2.add_patch(rect)
    ax2.text(
        0.5,
        0.5,
        "ProRS evolution adds domains, Protein Units, and structural ambiguity,\n"
        "but fails to overcome the intrinsic Pro/Thr specificity limit.",
        ha="center",
        va="center",
        fontsize=10,
    )

    fig.suptitle(
        "Structural elaboration without improved Pro/Thr specificity in ProRS evolution",
        fontsize=12,
        fontweight="bold",
    )

    outdir = Path("figures")
    outdir.mkdir(parents=True, exist_ok=True)
    png = outdir / "fig_graphical_abstract_proRS_evolution.png"
    svg = outdir / "fig_graphical_abstract_proRS_evolution.svg"
    fig.savefig(png, dpi=400)
    fig.savefig(svg)
    print(f"Saved {png}")
    print(f"Saved {svg}")


if __name__ == "__main__":
    main()
