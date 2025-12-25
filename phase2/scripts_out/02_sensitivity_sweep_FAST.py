#!/usr/bin/env python3
"""
02_sensitivity_sweep_FAST.py

OPTIMIZED sensitivity sweep:
  - Parse CIF files ONCE and cache atom coordinates
  - Pre-compute ALL pairwise distances ONCE
  - For each parameter combo, just apply thresholds to cached distances
  - ~81x faster than re-parsing CIFs for each combo

Statistics (from user's corrected code):
  - effect = mean(cognate) - mean(noncognate)
  - bootstrap CI on the *difference* (independent resampling of both groups)
  - permutation p-values on the same statistic
  - BH-FDR per metric across all tests
"""

from __future__ import annotations

import argparse
import sys
from itertools import product
from typing import Dict, List, Tuple, Optional
from pathlib import Path

import numpy as np
import pandas as pd

# Try to use gemmi for fast CIF parsing
try:
    import gemmi
    HAS_GEMMI = True
except ImportError:
    HAS_GEMMI = False


# -------------------------
# Stats helpers (from user's corrected code)
# -------------------------

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


# -------------------------
# CIF parsing and caching
# -------------------------

def parse_cif_atoms(cif_path: str, protein_chain: str = 'A', ligand_chain: str = 'B',
                    zn_chain: Optional[str] = None) -> Dict:
    """
    Parse CIF file and extract atom coordinates.
    Returns dict with protein_coords, protein_elements, ligand_coords, ligand_elements, zn_coords.
    """
    if HAS_GEMMI:
        return _parse_cif_gemmi(cif_path, protein_chain, ligand_chain, zn_chain)
    else:
        return _parse_cif_manual(cif_path, protein_chain, ligand_chain, zn_chain)


def _parse_cif_gemmi(cif_path: str, protein_chain: str, ligand_chain: str,
                     zn_chain: Optional[str]) -> Dict:
    structure = gemmi.read_structure(cif_path)
    model = structure[0]

    protein_coords = []
    protein_elements = []
    ligand_coords = []
    ligand_elements = []
    zn_coords = []

    for chain in model:
        if chain.name == protein_chain:
            for residue in chain:
                for atom in residue:
                    if atom.element.name != 'H':
                        protein_coords.append([atom.pos.x, atom.pos.y, atom.pos.z])
                        protein_elements.append(atom.element.name)
        elif chain.name == ligand_chain:
            for residue in chain:
                for atom in residue:
                    if atom.element.name != 'H':
                        ligand_coords.append([atom.pos.x, atom.pos.y, atom.pos.z])
                        ligand_elements.append(atom.element.name)
        elif zn_chain and chain.name == zn_chain:
            for residue in chain:
                for atom in residue:
                    if atom.element.name == 'Zn':
                        zn_coords.append([atom.pos.x, atom.pos.y, atom.pos.z])

    return {
        'protein_coords': np.array(protein_coords) if protein_coords else np.empty((0, 3)),
        'protein_elements': protein_elements,
        'ligand_coords': np.array(ligand_coords) if ligand_coords else np.empty((0, 3)),
        'ligand_elements': ligand_elements,
        'zn_coords': np.array(zn_coords) if zn_coords else np.empty((0, 3)),
    }


def _parse_cif_manual(cif_path: str, protein_chain: str, ligand_chain: str,
                      zn_chain: Optional[str]) -> Dict:
    atoms_by_chain = {}

    with open(cif_path) as f:
        lines = f.readlines()

    in_loop = False
    headers = []
    for line in lines:
        if line.startswith("loop_"):
            in_loop = True
            headers = []
            continue
        if in_loop and line.startswith("_atom_site."):
            headers.append(line.strip())
            continue
        if in_loop and headers and not line.startswith("_atom_site.") and line.strip():
            parts = line.split()
            if len(parts) < len(headers):
                continue
            hmap = {headers[j]: parts[j] for j in range(len(headers))}

            chain = hmap.get("_atom_site.label_asym_id", "")
            elem = hmap.get("_atom_site.type_symbol", "")
            if elem.upper() == "H":
                continue

            try:
                x = float(hmap["_atom_site.Cartn_x"])
                y = float(hmap["_atom_site.Cartn_y"])
                z = float(hmap["_atom_site.Cartn_z"])
            except:
                continue

            if chain not in atoms_by_chain:
                atoms_by_chain[chain] = {'coords': [], 'elements': []}
            atoms_by_chain[chain]['coords'].append([x, y, z])
            atoms_by_chain[chain]['elements'].append(elem)

    def get_chain_data(ch):
        if ch in atoms_by_chain:
            return np.array(atoms_by_chain[ch]['coords']), atoms_by_chain[ch]['elements']
        return np.empty((0, 3)), []

    p_coords, p_elem = get_chain_data(protein_chain)
    l_coords, l_elem = get_chain_data(ligand_chain)

    zn_coords = np.empty((0, 3))
    if zn_chain and zn_chain in atoms_by_chain:
        zn_data = atoms_by_chain[zn_chain]
        zn_idx = [i for i, e in enumerate(zn_data['elements']) if e == 'Zn']
        if zn_idx:
            zn_coords = np.array(zn_data['coords'])[zn_idx]

    return {
        'protein_coords': p_coords,
        'protein_elements': p_elem,
        'ligand_coords': l_coords,
        'ligand_elements': l_elem,
        'zn_coords': zn_coords,
    }


def compute_distance_matrix(coords1: np.ndarray, coords2: np.ndarray) -> np.ndarray:
    """Compute all pairwise distances between two sets of coordinates."""
    if coords1.shape[0] == 0 or coords2.shape[0] == 0:
        return np.empty((0, 0))
    # Use broadcasting for speed
    diff = coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :]
    return np.sqrt((diff ** 2).sum(axis=2))


# -------------------------
# Metric computation from cached distances
# -------------------------

POLAR_ATOMS = {'N', 'O', 'S'}

def compute_metrics_from_cache(cache: Dict, contact_cutoff: float, polar_cutoff: float,
                               polar_close_cutoff: float, clash_cutoff: float,
                               zn_engaged_cutoff: float) -> Dict:
    """
    Compute all metrics using cached coordinates and distance matrices.
    """
    dist_matrix = cache['dist_matrix']
    protein_elements = cache['protein_elements']
    ligand_elements = cache['ligand_elements']
    zn_coords = cache['zn_coords']
    ligand_coords = cache['ligand_coords']

    n_ligand = len(ligand_elements)
    if n_ligand == 0:
        return None

    # Contacts: pairs within cutoff
    contacts = (dist_matrix <= contact_cutoff).sum()

    # Clashes: pairs below cutoff
    clashes = (dist_matrix < clash_cutoff).sum()

    # Polar contacts: at least one polar atom, within cutoff
    polar_mask_protein = np.array([e in POLAR_ATOMS for e in protein_elements])
    polar_mask_ligand = np.array([e in POLAR_ATOMS for e in ligand_elements])

    # Either protein OR ligand atom is polar
    either_polar = polar_mask_protein[:, np.newaxis] | polar_mask_ligand[np.newaxis, :]
    polar_contacts = ((dist_matrix <= polar_cutoff) & either_polar).sum()

    # Polar close: both protein AND ligand atoms are polar (N/O/S pairs)
    both_polar = polar_mask_protein[:, np.newaxis] & polar_mask_ligand[np.newaxis, :]
    polar_close = ((dist_matrix < polar_close_cutoff) & both_polar).sum()

    # Normalized metrics
    contacts_per_atom = contacts / n_ligand
    polar_contacts_per_atom = polar_contacts / n_ligand
    clash_rate = clashes / n_ligand
    polar_close_per_atom = polar_close / n_ligand

    # Zn metrics
    zn_min_dist = None
    zn_engaged = False
    if zn_coords.shape[0] > 0 and ligand_coords.shape[0] > 0:
        # Distance to heteroatoms only
        hetero_mask = np.array([e in POLAR_ATOMS for e in ligand_elements])
        if hetero_mask.any():
            hetero_coords = ligand_coords[hetero_mask]
            zn_dists = compute_distance_matrix(zn_coords, hetero_coords)
            zn_min_dist = zn_dists.min()
            zn_engaged = (zn_min_dist <= zn_engaged_cutoff)

    return {
        'contacts_per_atom': contacts_per_atom,
        'polar_contacts_per_atom': polar_contacts_per_atom,
        'clash_rate': clash_rate,
        'polar_close_contacts_per_atom': polar_close_per_atom,
        'zn_min_dist_hetero': zn_min_dist,
        'zn_engaged': zn_engaged,
    }


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


def compute_effects_for_df(df_setting: pd.DataFrame,
                           min_n_cognate: int,
                           min_n_noncog: int,
                           n_boot: int,
                           n_perm: int,
                           seed: int) -> pd.DataFrame:
    """Compute condition-level effects for each metric."""
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_csv", required=True, help="Input clean CSV")
    ap.add_argument("--out_csv", required=True, help="Output CSV")
    ap.add_argument("--n_boot", type=int, default=10_000)
    ap.add_argument("--n_perm", type=int, default=10_000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min_n_cognate", type=int, default=3)
    ap.add_argument("--min_n_noncog", type=int, default=3)
    args = ap.parse_args()

    df = pd.read_csv(args.in_csv)
    print(f"Loaded {len(df)} rows", flush=True)

    required = {"file", "condition", "enzyme", "ligand", "zn_present"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Input missing required columns: {sorted(missing)}")

    # Parameter grid (3^4 = 81 combos)
    CONTACT_CUTOFFS = [3.5, 4.0, 4.5]
    POLAR_CUTOFFS = [3.5, 4.0, 4.5]
    POLAR_CLOSE_CUTOFFS = [3.2, 3.5, 3.8]
    CLASH_CUTOFFS = [2.0, 2.2, 2.4]
    ZN_ENGAGED_CUTOFF = 3.0  # fixed

    grid = list(product(CONTACT_CUTOFFS, POLAR_CUTOFFS, POLAR_CLOSE_CUTOFFS, CLASH_CUTOFFS))
    assert len(grid) == 81

    # PHASE 1: Parse all CIF files ONCE and cache coordinates + distance matrices
    print("Phase 1: Parsing CIF files and caching distances...", flush=True)
    cache_list = []

    for idx, row in df.iterrows():
        cif_path = row['file']
        if not Path(cif_path).exists():
            cache_list.append(None)
            continue

        try:
            # Determine Zn chain
            zn_chain = None
            if row.get('zn_present', False):
                if 'chains' in row:
                    chains_str = str(row['chains'])
                    for ch in ['C', 'D', 'E']:
                        if f'{ch}:1:1' in chains_str:
                            zn_chain = ch
                            break

            atoms = parse_cif_atoms(cif_path, 'A', 'B', zn_chain)

            if atoms['ligand_coords'].shape[0] == 0:
                cache_list.append(None)
                continue

            # Pre-compute distance matrix
            dist_matrix = compute_distance_matrix(atoms['protein_coords'], atoms['ligand_coords'])

            cache_list.append({
                'dist_matrix': dist_matrix,
                'protein_elements': atoms['protein_elements'],
                'ligand_elements': atoms['ligand_elements'],
                'zn_coords': atoms['zn_coords'],
                'ligand_coords': atoms['ligand_coords'],
            })
        except Exception as e:
            print(f"Error parsing {cif_path}: {e}", file=sys.stderr)
            cache_list.append(None)

        if (idx + 1) % 50 == 0:
            print(f"  Parsed {idx+1}/{len(df)}", flush=True)

    valid_count = sum(1 for c in cache_list if c is not None)
    print(f"Phase 1 complete: {valid_count}/{len(df)} valid structures", flush=True)

    # PHASE 2: For each parameter combo, compute metrics from cache
    print("Phase 2: Computing metrics for all parameter combos...", flush=True)
    all_results = []

    for i, (c_cut, p_cut, pc_cut, cl_cut) in enumerate(grid, start=1):
        # Compute metrics for all rows using cached distances
        metrics_list = []
        for j, cache in enumerate(cache_list):
            if cache is None:
                metrics_list.append(None)
            else:
                metrics = compute_metrics_from_cache(
                    cache, c_cut, p_cut, pc_cut, cl_cut, ZN_ENGAGED_CUTOFF
                )
                metrics_list.append(metrics)

        # Filter valid and build dataframe
        valid_indices = [j for j, m in enumerate(metrics_list) if m is not None]
        valid_metrics = [metrics_list[j] for j in valid_indices]

        if not valid_metrics:
            continue

        df_metrics = pd.DataFrame(valid_metrics)
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

        eff["contact_cutoff"] = c_cut
        eff["polar_cutoff"] = p_cut
        eff["polar_close_cutoff"] = pc_cut
        eff["clash_cutoff"] = cl_cut
        eff["zn_engaged_cutoff"] = ZN_ENGAGED_CUTOFF
        all_results.append(eff)
        print(f"[{i:02d}/81] done", flush=True)

    out = pd.concat(all_results, ignore_index=True)

    # BH-FDR per metric
    out["q_fdr"] = np.nan
    for metric in METRICS:
        mask = (out["metric"] == metric) & np.isfinite(out["p_perm"])
        out.loc[mask, "q_fdr"] = bh_fdr(out.loc[mask, "p_perm"].to_numpy())

    # Sanity check: CI should bracket effect
    sane = out[np.isfinite(out["effect"]) & np.isfinite(out["ci_low"]) & np.isfinite(out["ci_high"])]
    if len(sane) > 0:
        bad = (sane["effect"] < sane["ci_low"]) | (sane["effect"] > sane["ci_high"])
        if bad.any():
            n_bad = int(bad.sum())
            print(f"WARNING: {n_bad} rows have CIs that do not bracket effect.", file=sys.stderr)

    out.to_csv(args.out_csv, index=False)
    print(f"Wrote: {args.out_csv}  (rows={len(out)})", flush=True)


if __name__ == "__main__":
    main()
