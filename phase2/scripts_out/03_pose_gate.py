#!/usr/bin/env python3
"""
03_pose_gate.py

Pose-gating analysis:
  - Compute ligand pose RMSD relative to cognate pose within each condition.
  - Keep only noncognates whose pose is within a threshold (e.g., 2Å or 4Å).
  - Recompute cognate vs gated-noncognate effects:
      effect = mean(cognate) - mean(noncognate_gated)
      bootstrap CI on the difference (correct)
      permutation p-values (correct)
      BH-FDR across the 18 tests (or per metric if you prefer; here: per metric)

This expects your CIFs on disk (file column) and a ligand-residue chain inference
consistent with how you stored them.

Dependencies: numpy, pandas
"""

from __future__ import annotations

import argparse
import sys
from typing import Dict, Tuple, List

import numpy as np
import pandas as pd

# Stats helpers (copied from 02_sensitivity_sweep.py to avoid import issues with dataclasses in Python 3.7)

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


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR. Returns q-values aligned to original pvals positions."""
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
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0.0, 1.0)

    q_f = np.empty_like(p_f)
    q_f[order] = q_sorted
    q[finite] = q_f
    return q


def bootstrap_ci_diff(a: np.ndarray,
                      b: np.ndarray,
                      n_boot: int = 10_000,
                      alpha: float = 0.05,
                      seed: int = 0) -> Tuple[float, float]:
    """Bootstrap CI for difference in means: mean(a) - mean(b)."""
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
    """Two-sided permutation test for difference in means."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    rng = np.random.default_rng(seed)

    na, nb = a.size, b.size
    if na < 2 or nb < 2:
        return np.nan

    obs = a.mean() - b.mean()
    pooled = np.concatenate([a, b], axis=0)
    n = pooled.size

    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(n)
        pa = pooled[perm[:na]]
        pb = pooled[perm[na:]]
        stat = pa.mean() - pb.mean()
        if abs(stat) >= abs(obs) - 1e-12:
            count += 1
    return (count + 1) / (n_perm + 1)


# If you have a robust ligand superposition module, import it.
# Otherwise, this fallback extracts non-H ligand coordinates from CIF atom_site lines.
def extract_ligand_coords_from_cif(cif_file: str, ligand_chain: str = "B") -> np.ndarray:
    """
    CIF atom_site parser for AF3 structures:
    - collects atoms from the ligand chain (typically 'B' for amino acid ligands)
    - excludes hydrogens by element symbol
    - returns (N,3) float array
    """
    coords = []
    with open(cif_file, "r") as f:
        lines = f.readlines()

    # crude atom_site loop parse
    in_loop = False
    headers = []
    for i, line in enumerate(lines):
        if line.startswith("loop_"):
            in_loop = True
            headers = []
            continue
        if in_loop and line.startswith("_atom_site."):
            headers.append(line.strip())
            continue
        if in_loop and headers and not line.startswith("_atom_site.") and line.strip():
            # data lines until blank or new loop
            parts = line.split()
            if len(parts) < len(headers):
                continue
            hmap = {headers[j]: parts[j] for j in range(len(headers))}

            # Check if this atom is in the ligand chain
            chain = hmap.get("_atom_site.label_asym_id", "")
            if chain != ligand_chain:
                continue

            elem = hmap.get("_atom_site.type_symbol", "")
            if elem.upper() == "H":
                continue

            try:
                x = float(hmap["_atom_site.Cartn_x"])
                y = float(hmap["_atom_site.Cartn_y"])
                z = float(hmap["_atom_site.Cartn_z"])
            except Exception:
                continue

            coords.append([x, y, z])

        if in_loop and line.startswith("loop_") and headers:
            # next loop
            in_loop = True

    if len(coords) < 3:
        return np.empty((0, 3), dtype=float)
    return np.asarray(coords, dtype=float)


def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """
    RMSD after optimal rigid alignment (Kabsch).
    Requires same shape (N,3).
    """
    if P.shape != Q.shape or P.shape[0] < 3:
        return np.nan
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    C = Pc.T @ Qc
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    U = V @ D @ Wt
    P_rot = Pc @ U
    diff = P_rot - Qc
    return float(np.sqrt((diff * diff).sum() / P.shape[0]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_csv", required=True, help="geometry_metrics_clean.csv")
    ap.add_argument("--out_csv", required=True, help="effects_after_gate.csv")
    ap.add_argument("--gates", default="2.0,4.0", help="Comma-separated RMSD gates in Å (e.g., 2.0,4.0)")
    ap.add_argument("--n_boot", type=int, default=10_000)
    ap.add_argument("--n_perm", type=int, default=10_000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min_n_cognate", type=int, default=3)
    ap.add_argument("--min_n_noncog", type=int, default=3)
    args = ap.parse_args()

    gates = [float(x) for x in args.gates.split(",") if x.strip()]
    df = pd.read_csv(args.in_csv)

    required = {"file", "condition", "enzyme", "ligand"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Input missing required columns: {sorted(missing)}")

    # Pre-extract ligand coords once per file (cache)
    coord_cache: Dict[str, np.ndarray] = {}
    for f in df["file"].unique():
        try:
            coord_cache[f] = extract_ligand_coords_from_cif(f)
        except Exception:
            coord_cache[f] = np.empty((0, 3), dtype=float)

    rows = []

    for condition, g in df.groupby("condition", sort=False):
        enzyme = g["enzyme"].iloc[0]
        cog = COGNATE_BY_ENZYME.get(enzyme)
        if cog is None:
            continue

        g_cog = g[g["ligand"] == cog]
        if len(g_cog) == 0:
            continue

        # If multiple cognate files, choose first as reference
        cog_file = g_cog["file"].iloc[0]
        cog_coords = coord_cache.get(cog_file, np.empty((0, 3), dtype=float))
        if cog_coords.shape[0] < 3:
            continue

        # RMSD per row vs cognate pose
        rmsd_list = []
        for idx, row in g.iterrows():
            coords = coord_cache.get(row["file"], np.empty((0, 3), dtype=float))
            if coords.shape != cog_coords.shape:
                rmsd = np.nan
            else:
                rmsd = kabsch_rmsd(coords, cog_coords)
            rmsd_list.append(rmsd)

        g = g.copy()
        g["pose_rmsd_to_cognate"] = rmsd_list

        for gate in gates:
            g_gated = g[(g["ligand"] == cog) | ((g["ligand"] != cog) & (g["pose_rmsd_to_cognate"] <= gate))].copy()

            for metric in METRICS:
                if metric not in g_gated.columns:
                    continue

                cog_vals = g_gated.loc[g_gated["ligand"] == cog, metric].dropna().to_numpy()
                non_vals = g_gated.loc[g_gated["ligand"] != cog, metric].dropna().to_numpy()

                n_cog = int(cog_vals.size)
                n_non = int(non_vals.size)

                effect = float(np.mean(cog_vals) - np.mean(non_vals)) if (n_cog > 0 and n_non > 0) else np.nan

                if n_cog >= args.min_n_cognate and n_non >= args.min_n_noncog:
                    ci_lo, ci_hi = bootstrap_ci_diff(cog_vals, non_vals, n_boot=args.n_boot, seed=args.seed)
                    p_perm = permutation_pvalue_diff(cog_vals, non_vals, n_perm=args.n_perm, seed=args.seed)
                else:
                    ci_lo, ci_hi, p_perm = (np.nan, np.nan, np.nan)

                rows.append({
                    "gate_A": gate,
                    "condition": condition,
                    "enzyme": enzyme,
                    "metric": metric,
                    "cognate": cog,
                    "n_cognate": n_cog,
                    "n_noncognate": n_non,
                    "effect": effect,
                    "ci_low": ci_lo,
                    "ci_high": ci_hi,
                    "p_perm": p_perm,
                })

    out = pd.DataFrame(rows)

    # BH-FDR per metric across pose-gated tests
    out["q_fdr"] = np.nan
    for metric in METRICS:
        mask = (out["metric"] == metric) & np.isfinite(out["p_perm"])
        out.loc[mask, "q_fdr"] = bh_fdr(out.loc[mask, "p_perm"].to_numpy())

    # Sanity: CI should bracket effect when finite
    sane = out[np.isfinite(out["effect"]) & np.isfinite(out["ci_low"]) & np.isfinite(out["ci_high"])]
    if len(sane) > 0:
        bad = (sane["effect"] < sane["ci_low"]) | (sane["effect"] > sane["ci_high"])
        if bad.any():
            print(f"WARNING: {int(bad.sum())} rows have CIs that do not bracket effect.", file=sys.stderr)

    out.to_csv(args.out_csv, index=False)
    print(f"Wrote: {args.out_csv} (rows={len(out)})")


if __name__ == "__main__":
    main()
