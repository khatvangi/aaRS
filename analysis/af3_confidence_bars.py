#!/usr/bin/env python3
"""
Plot AF3 confidence (ipTM / pTM / mean pLDDT) for ProRS complexes.

Usage:
PYTHONPATH=. python -m analysis.af3_confidence_bars \
    --config metadata/af3_promiscuity_metrics.json \
    --out-fig figures/af3_confidence_bars.png
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def load_metrics(label: str, dirpath: Path) -> dict:
    dirpath = Path(dirpath)
    if not dirpath.exists():
        raise FileNotFoundError(f"Metrics directory missing: {dirpath}")

    summary_candidates = sorted(dirpath.glob("*summary_confidences.json"))
    if not summary_candidates:
        summary_candidates = sorted(dirpath.glob("*.json"))
    if not summary_candidates:
        raise FileNotFoundError(f"No JSON metrics found in {dirpath}")
    summary_path = summary_candidates[0]

    with open(summary_path) as f:
        data = json.load(f)

    iptm = data.get("iptm", data.get("ipTM"))
    ptm = data.get("ptm", data.get("pTM"))

    mean_plddt: Optional[float] = None
    confid_path = dirpath / (dirpath.name + "_confidences.json")
    if confid_path.exists():
        try:
            with open(confid_path) as f:
                conf = json.load(f)
            if isinstance(conf, dict) and "plddt" in conf:
                vals = conf["plddt"]
                if vals:
                    mean_plddt = float(np.mean(vals)) / 100.0
        except Exception:
            mean_plddt = None

    return {
        "label": label,
        "iptm": iptm,
        "ptm": ptm,
        "mean_plddt": mean_plddt,
    }


def main():
    ap = argparse.ArgumentParser(description="Plot AF3 confidence metrics.")
    ap.add_argument("--config", required=True, help="JSON list of {label, metrics_dir}")
    ap.add_argument("--out-fig", required=True)
    args = ap.parse_args()

    cfg = pd.read_json(args.config)
    rows = [load_metrics(row["label"], Path(row["metrics_dir"])) for _, row in cfg.iterrows()]
    df = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=(6, 4))
    x = np.arange(len(df))
    width = 0.25

    ax.bar(x - width, df["iptm"], width, label="ipTM")
    if df["ptm"].notna().any():
        ax.bar(x, df["ptm"], width, label="pTM")
    if df["mean_plddt"].notna().any():
        ax.bar(x + width, df["mean_plddt"], width, label="mean pLDDT (scaled)")

    ax.set_xticks(x)
    ax.set_xticklabels(df["label"], rotation=45, ha="right")
    ax.set_ylabel("Score")
    ax.set_ylim(0, 1.05)
    ax.set_title("AF3 confidence for ProRS ligand complexes")
    ax.legend()
    fig.tight_layout()

    out = Path(args.out_fig)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
