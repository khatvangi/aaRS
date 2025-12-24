#!/usr/bin/env python3
"""
02_sensitivity_sweep.py

Sensitivity sweep over geometry metric recomputation parameters.
Computes cognate vs non-cognate effects with:
  - effect = mean(cognate) - mean(noncognate)
  - bootstrap CI on the *difference* (independent resampling of both groups)
  - permutation p-values on the same statistic
  - BH-FDR per metric across all tests

This script expects:
  - a "clean" input CSV with at least columns:
      file, condition, enzyme, ligand, zn_present
  - 00_recompute_metrics.py available in PYTHONPATH or same folder,
    providing recompute_metrics_for_row(row_dict, params_dict) -> dict
      containing metric columns:
        contacts_per_atom, polar_contacts_per_atom, clash_rate, polar_close_contacts_per_atom
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from dataclasses import dataclass
from itertools import product
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd

# Local import: must exist (use importlib since module name starts with digit)
import importlib.util
script_dir = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("recompute_metrics", os.path.join(script_dir, "00_recompute_metrics.py"))
recompute_module = importlib.util.module_from_spec(spec)
try:
    spec.loader.exec_module(recompute_module)
    recompute_metrics_for_row = recompute_module.recompute_metrics_for_row
except Exception as e:
    print("ERROR: Could not import recompute_metrics_for_row from 00_recompute_metrics.py", file=sys.stderr)
    raise


# -------------------------
# Stats helpers (correct)
# -------------------------

def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """
    Benjamini-Hochberg FDR.
    Returns q-values aligned to original pvals positions.
    """
    p = np.asarray(pvals, dtype=float)
    q = np.full_like(p, np.nan, dtype=float)

    finite = np.isfinite(p)
    if finite.sum() == 0:
        return q

    p_f = p[finite]
    m = p_f.size
    order = np.argsort(p_f)
    ranked = p_f[order]
    ranks = np.arange(1, m + 1, dtype=float)

    q_sorted = ranked * m / ranks
    # enforce monotonicity
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0.0, 1.0)

    # put back
    q_f = np.empty_like(p_f)
    q_f[order] = q_sorted
    q[finite] = q_f
    return q


def bootstrap_ci_diff(a: np.ndarray,
                      b: np.ndarray,
                      n_boot: int = 10_000,
                      alpha: float = 0.05,
                      seed: int = 0) -> Tuple[float, float]:
    """
    Bootstrap CI for difference in means: mean(a) - mean(b).
    Resample a and b independently with replacement.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    rng = np.random.default_rng(seed)

    na, nb = a.size, b.size
    if na < 2 or nb < 2:
        return (np.nan, np.nan)

    ia = rng.integers(0, na, size=(n_boot, na))
    ib = rng.integers(0, nb, size=(n_boot, nb))

    diffs = a[ia].mean(axis=1) - b[ib].mean(axis=1)
    lo = np.quantile(diffs, alpha / 2.0)
    hi = np.quantile(diffs, 1.0 - alpha / 2.0)
    return (float(lo), float(hi))


def permutation_pvalue_diff(a: np.ndarray,
                            b: np.ndarray,
                            n_perm: int = 10_000,
                            seed: int = 0) -> float:
    """
    Two-sided permutation test for difference in means: mean(a)-mean(b).
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    rng = np.random.default_rng(seed)

    na, nb = a.size, b.size
    if na < 2 or nb < 2:
        return np.nan

    obs = a.mean() - b.mean()
    pooled = np.concatenate([a, b], axis=0)
    n = pooled.size

    # Permute labels by shuffling indices
    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(n)
        pa = pooled[perm[:na]]
        pb = pooled[perm[na:]]
        stat = pa.mean() - pb.mean()
        if abs(stat) >= abs(obs) - 1e-12:
            count += 1
    return (count + 1) / (n_perm + 1)


# -------------------------
# Domain logic
# -------------------------

COGNATE_BY_ENZYME = {
    "ProRS": "PRO",
    "ThrRS": "THR",
}

METRICS = [
    "contacts_per_atom",
    "polar_contacts_per_atom",
    "polar_close_contacts_per_atom",
    "clash_rate",
]


@dataclass(frozen=True)
class Params:
    contact_cutoff: float
    polar_cutoff: float
    polar_close_cutoff: float
    clash_cutoff: float
    zn_engaged_cutoff: float

    def as_dict(self) -> Dict[str, float]:
        return {
            "contact_cutoff": self.contact_cutoff,
            "polar_cutoff": self.polar_cutoff,
            "polar_close_cutoff": self.polar_close_cutoff,
            "clash_cutoff": self.clash_cutoff,
            "zn_engaged_cutoff": self.zn_engaged_cutoff,
        }


def compute_effects_for_df(df_setting: pd.DataFrame,
                           min_n_cognate: int,
                           min_n_noncog: int,
                           n_boot: int,
                           n_perm: int,
                           seed: int) -> pd.DataFrame:
    """
    Given a dataframe with recomputed metric columns for a *single* parameter setting,
    compute condition-level effects for each metric.
    """
    rows = []

    for condition, g in df_setting.groupby("condition", sort=False):
        enzyme = g["enzyme"].iloc[0]
        cog = COGNATE_BY_ENZYME.get(enzyme)
        if cog is None:
            continue

        for metric in METRICS:
            if metric not in g.columns:
                continue

            cog_vals = g.loc[g["ligand"] == cog, metric].dropna().to_numpy()
            non_vals = g.loc[g["ligand"] != cog, metric].dropna().to_numpy()

            n_cog = int(cog_vals.size)
            n_non = int(non_vals.size)

            effect = float(np.mean(cog_vals) - np.mean(non_vals)) if (n_cog > 0 and n_non > 0) else np.nan

            if n_cog >= min_n_cognate and n_non >= min_n_noncog:
                ci_lo, ci_hi = bootstrap_ci_diff(cog_vals, non_vals, n_boot=n_boot, seed=seed)
                p_perm = permutation_pvalue_diff(cog_vals, non_vals, n_perm=n_perm, seed=seed)
            else:
                ci_lo, ci_hi, p_perm = (np.nan, np.nan, np.nan)

            rows.append({
                "metric": metric,
                "condition": condition,
                "enzyme": enzyme,
                "cognate": cog,
                "n_cognate": n_cog,
                "n_noncognate": n_non,
                "effect": effect,
                "ci_low": ci_lo,
                "ci_high": ci_hi,
                "p_perm": p_perm,
            })

    return pd.DataFrame(rows)


def _worker_recompute(args):
    row_dict, params_dict = args
    out = recompute_metrics_for_row(row_dict, params_dict)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_csv", required=True, help="Input clean CSV (e.g., geometry_metrics_clean.csv)")
    ap.add_argument("--out_csv", required=True, help="Output CSV (effects stability table)")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--n_boot", type=int, default=10_000)
    ap.add_argument("--n_perm", type=int, default=10_000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min_n_cognate", type=int, default=3)
    ap.add_argument("--min_n_noncog", type=int, default=3)
    args = ap.parse_args()

    df = pd.read_csv(args.in_csv)

    required = {"file", "condition", "enzyme", "ligand", "zn_present"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Input missing required columns: {sorted(missing)}")

    # Parameter grid (3x3x3x3 + zn cutoff = 3 => 243; if you want 81, keep zn fixed or shrink grid)
    # Your prior run claimed 81 combos; that implies 4 dimensions at 3 levels (3^4=81) and zn cutoff fixed.
    CONTACT_CUTOFFS = [3.5, 4.0, 4.5]
    POLAR_CUTOFFS = [3.5, 4.0, 4.5]
    POLAR_CLOSE_CUTOFFS = [3.2, 3.5, 3.8]
    CLASH_CUTOFFS = [2.0, 2.2, 2.4]
    ZN_ENGAGED_CUTOFF = 3.0  # fixed to preserve 81-combo design

    grid: List[Params] = [
        Params(c, p, pc, cl, ZN_ENGAGED_CUTOFF)
        for c, p, pc, cl in product(CONTACT_CUTOFFS, POLAR_CUTOFFS, POLAR_CLOSE_CUTOFFS, CLASH_CUTOFFS)
    ]
    assert len(grid) == 81, f"Grid size mismatch: expected 81, got {len(grid)}"

    # Prepare lightweight row dicts for recompute
    base_rows = df[["file", "zn_present"]].copy()
    base_rows["zn_present"] = base_rows["zn_present"].astype(bool)
    row_dicts = base_rows.to_dict(orient="records")

    all_results = []

    # Multiprocessing pool
    if args.workers > 1:
        import multiprocessing as mp
        # Use fork on Linux (works with dynamic imports), spawn on others
        ctx = mp.get_context("fork")
        pool = ctx.Pool(processes=args.workers)
        try:
            for i, params in enumerate(grid, start=1):
                params_dict = params.as_dict()
                tasks = [(row_dicts[j], params_dict) for j in range(len(row_dicts))]
                recomputed = pool.map(_worker_recompute, tasks, chunksize=max(1, len(tasks)//(args.workers*8)))
                # Filter out None results (failed parses)
                valid_indices = [j for j, r in enumerate(recomputed) if r is not None]
                recomputed = [recomputed[j] for j in valid_indices]
                df_metrics = pd.DataFrame(recomputed)
                # Keep only corresponding rows from df
                df_subset = df.iloc[valid_indices].reset_index(drop=True)

                # Merge back condition/ligand labels
                df_setting = df_subset[["condition", "enzyme", "ligand"]].reset_index(drop=True).join(df_metrics)

                eff = compute_effects_for_df(
                    df_setting,
                    min_n_cognate=args.min_n_cognate,
                    min_n_noncog=args.min_n_noncog,
                    n_boot=args.n_boot,
                    n_perm=args.n_perm,
                    seed=args.seed,
                )
                # annotate params
                eff["contact_cutoff"] = params.contact_cutoff
                eff["polar_cutoff"] = params.polar_cutoff
                eff["polar_close_cutoff"] = params.polar_close_cutoff
                eff["clash_cutoff"] = params.clash_cutoff
                eff["zn_engaged_cutoff"] = params.zn_engaged_cutoff
                all_results.append(eff)
                print(f"[{i:02d}/81] done", flush=True)
        finally:
            pool.close()
            pool.join()
    else:
        for i, params in enumerate(grid, start=1):
            params_dict = params.as_dict()
            recomputed = [recompute_metrics_for_row(r, params_dict) for r in row_dicts]
            # Filter out None results (failed parses)
            valid_indices = [j for j, r in enumerate(recomputed) if r is not None]
            recomputed = [recomputed[j] for j in valid_indices]
            df_metrics = pd.DataFrame(recomputed)
            df_subset = df.iloc[valid_indices].reset_index(drop=True)
            df_setting = df_subset[["condition", "enzyme", "ligand"]].reset_index(drop=True).join(df_metrics)

            eff = compute_effects_for_df(
                df_setting,
                min_n_cognate=args.min_n_cognate,
                min_n_noncog=args.min_n_noncog,
                n_boot=args.n_boot,
                n_perm=args.n_perm,
                seed=args.seed,
            )
            eff["contact_cutoff"] = params.contact_cutoff
            eff["polar_cutoff"] = params.polar_cutoff
            eff["polar_close_cutoff"] = params.polar_close_cutoff
            eff["clash_cutoff"] = params.clash_cutoff
            eff["zn_engaged_cutoff"] = params.zn_engaged_cutoff
            all_results.append(eff)
            print(f"[{i:02d}/81] done", flush=True)

    out = pd.concat(all_results, ignore_index=True)

    # BH-FDR per metric across all tests (matches your claims)
    out["q_fdr"] = np.nan
    for metric in METRICS:
        mask = (out["metric"] == metric) & np.isfinite(out["p_perm"])
        out.loc[mask, "q_fdr"] = bh_fdr(out.loc[mask, "p_perm"].to_numpy())

    # Sanity: CI should bracket effect when finite
    sane = out[np.isfinite(out["effect"]) & np.isfinite(out["ci_low"]) & np.isfinite(out["ci_high"])]
    if len(sane) > 0:
        bad = (sane["effect"] < sane["ci_low"]) | (sane["effect"] > sane["ci_high"])
        if bad.any():
            n_bad = int(bad.sum())
            print(f"WARNING: {n_bad} rows have CIs that do not bracket effect (should be ~0).", file=sys.stderr)

    out.to_csv(args.out_csv, index=False)
    print(f"Wrote: {args.out_csv}  (rows={len(out)})")


if __name__ == "__main__":
    main()
