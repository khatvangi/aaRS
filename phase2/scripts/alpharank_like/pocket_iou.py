#!/usr/bin/env python3
"""Compute pocket IoU metrics for AF3Score outputs."""
import argparse
import json
import logging
import math
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

PROJECT_ROOT = Path("/storage/kiran-stuff/aaRS/phase2")
RESULTS_DIR = PROJECT_ROOT / "results"
METRICS_IN = RESULTS_DIR / "metrics_long.csv"

try:
    from Bio import PDB
    BIO_AVAILABLE = True
except ImportError:
    BIO_AVAILABLE = False


def setup_logger(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


def load_json(path: Path) -> Optional[Dict]:
    try:
        return json.loads(path.read_text())
    except Exception as e:
        logging.warning("Could not read %s: %s", path, e)
        return None


def load_structure(path: Path):
    if not BIO_AVAILABLE:
        raise RuntimeError("Biopython not available")
    if path.suffix.lower() == ".cif":
        parser = PDB.MMCIFParser(QUIET=True)
    else:
        parser = PDB.PDBParser(QUIET=True)
    return parser.get_structure("model", path.as_posix())


EXCLUDE_LIG = {"HOH", "WAT", "NA", "CL", "MG", "K", "ZN"}


def get_ligand_residues(structure, codes: List[str]) -> List:
    ligs = []
    wanted = {c.upper() for c in codes if c}
    model = structure[0]
    for chain in model:
        for res in chain:
            if res.id[0] == " ":
                continue
            resname = res.get_resname().upper()
            if resname in EXCLUDE_LIG:
                continue
            if wanted and resname not in wanted:
                continue
            ligs.append(res)
    if ligs:
        return ligs
    # fallback: largest hetero group
    fallback = []
    max_atoms = 0
    for chain in model:
        for res in chain:
            if res.id[0] == " ":
                continue
            resname = res.get_resname().upper()
            if resname in EXCLUDE_LIG:
                continue
            atoms = [a for a in res.get_atoms()]
            if len(atoms) > max_atoms:
                max_atoms = len(atoms)
                fallback = [res]
    return fallback


def heavy_atom_coords(residue) -> np.ndarray:
    coords = []
    for atom in residue.get_atoms():
        element = atom.element.strip().upper()
        if element == "H":
            continue
        coords.append(atom.get_coord())
    return np.array(coords)


def pocket_residues(structure, ligand_residues: List, cutoff: float = 6.0) -> Set[Tuple[str, int, str, str]]:
    ligand_coords = []
    for res in ligand_residues:
        c = heavy_atom_coords(res)
        if c.size:
            ligand_coords.append(c)
    if not ligand_coords:
        return set()
    ligand_coords = np.concatenate(ligand_coords, axis=0)
    pocket: Set[Tuple[str, int, str, str]] = set()
    cutoff_sq = cutoff * cutoff
    model = structure[0]
    for chain in model:
        for res in chain:
            if res.id[0] != " ":
                continue
            res_atoms = heavy_atom_coords(res)
            if not res_atoms.size:
                continue
            dists = np.sum((res_atoms[:, None, :] - ligand_coords[None, :, :]) ** 2, axis=2)
            if np.any(dists <= cutoff_sq):
                pocket.add((chain.id, res.id[1], res.id[2].strip() if isinstance(res.id[2], str) else "", res.get_resname().upper()))
    return pocket


def compute_iou(a: Set, b: Set) -> float:
    if not a and not b:
        return math.nan
    if not a or not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else math.nan


def infer_era(prefix: str, job_name: str) -> str:
    text = f"{prefix}_{job_name}".lower()
    if "anc" in text or "ancient" in text:
        return "ancient"
    if "modern" in text:
        return "modern"
    return "unknown"


def compute_iou(a: Set, b: Set) -> float:
    if not a and not b:
        return math.nan
    if not a or not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else math.nan


def pockets_for_structure(struct_path: Path, lig1: str, lig2: Optional[str]) -> Tuple[Set, Set]:
    structure = load_structure(struct_path)
    p1 = pocket_residues(structure, get_ligand_residues(structure, [lig1]))
    p2 = set()
    if lig2:
        p2 = pocket_residues(structure, get_ligand_residues(structure, [lig2]))
    return p1, p2


def compute_competition_metrics(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    df = df.copy()
    for col, default in [
        ("pocket_size_1", np.nan),
        ("pocket_size_2", np.nan),
        ("pref_size", np.nan),
        ("iou_1", np.nan),
        ("iou_2", np.nan),
        ("delta_iou", np.nan),
        ("win", np.nan),
        ("pocket_status", "unprocessed"),
        ("pocket_size", np.nan),
        ("pocket_hash", ""),
        ("iou_pref", np.nan),
    ]:
        if col not in df.columns:
            df[col] = default
    cache_comp: Dict[str, Tuple[Set, Set]] = {}
    cache_std: Dict[str, Set] = {}
    exclude = EXCLUDE_LIG | {"CA","MN","FE","CO","NI","CU","SO4","PO4"}
    comp = df[df["run_type"] == "competition"]
    for job, grp in comp.groupby("job_name"):
        lig1 = grp["ligand_1"].iloc[0]
        lig2 = grp["ligand_2"].iloc[0]
        # choose top replicate as preference
        top = grp.sort_values("ranking_score", ascending=False).iloc[0]
        pref_struct = Path(top["structure_file"])
        try:
            pref_p1, pref_p2 = pockets_for_structure(pref_struct, lig1, lig2)
        except Exception as e:
            logging.warning("Pref pocket failed for %s: %s", pref_struct, e)
            df.loc[grp.index, "pocket_status"] = "pref_missing"
            continue
        pref = pref_p1
        pref_size = len(pref)
        if pref_size == 0:
            df.loc[grp.index, "pocket_status"] = "pref_missing"
            continue
        for idx, row in grp.iterrows():
            struct_path = Path(row["structure_file"])
            if struct_path.as_posix() in cache_comp:
                p1, p2 = cache_comp[struct_path.as_posix()]
            else:
                try:
                    p1, p2 = pockets_for_structure(struct_path, lig1, lig2)
                    cache_comp[struct_path.as_posix()] = (p1, p2)
                except Exception as e:
                    logging.warning("Pocket parse failed for %s: %s", struct_path, e)
                    df.at[idx, "pocket_status"] = "parse_error"
                    continue
            df.at[idx, "pocket_size_1"] = len(p1)
            df.at[idx, "pocket_size_2"] = len(p2)
            df.at[idx, "pref_size"] = pref_size
            if len(p1) == 0 or len(p2) == 0:
                df.at[idx, "pocket_status"] = "missing_ligand"
                continue
            iou1 = compute_iou(p1, pref)
            iou2 = compute_iou(p2, pref)
            df.at[idx, "iou_1"] = iou1
            df.at[idx, "iou_2"] = iou2
            delta = iou1 - iou2
            df.at[idx, "delta_iou"] = delta
            df.at[idx, "win"] = 1.0 if delta > 0 else 0.0
            df.at[idx, "pocket_status"] = "ok"
    # STANDARD runs: pocket and IoU to cognate pref
    std = df[df["run_type"] == "standard"]
    # compute pockets per replicate
    for idx, row in std.iterrows():
        lig = row.get("ligand_code")
        struct_path = Path(row["structure_file"])
        if not isinstance(lig, str) or not lig:
            df.at[idx, "pocket_status"] = "missing_ligand"
            continue
        if struct_path.as_posix() in cache_std:
            pocket = cache_std[struct_path.as_posix()]
        else:
            try:
                structure = load_structure(struct_path)
                lig_res = get_ligand_residues(structure, [lig])
                if not lig_res:
                    df.at[idx, "pocket_status"] = "missing_ligand"
                    continue
                pocket = pocket_residues(structure, lig_res)
                cache_std[struct_path.as_posix()] = pocket
            except Exception as e:
                logging.warning("Pocket parse failed for %s: %s", struct_path, e)
                df.at[idx, "pocket_status"] = "parse_error"
                continue
        df.at[idx, "pocket_size"] = len(pocket)
        df.at[idx, "pocket_hash"] = ";".join([f"{c}:{r}:{i}:{n}" for (c,r,i,n) in sorted(list(pocket))])
        df.at[idx, "pocket_status"] = "processed" if len(pocket) else "empty_pocket"
    # build prefs per (era,enzyme)
    cognate_map = {"prors": "PRO", "thrrs": "THR"}
    for (era, enzyme), grp in std.groupby(["era", "enzyme"]):
        cognate = cognate_map.get(str(enzyme).lower(), None)
        if cognate is None:
            continue
        cog_rows = grp[grp["ligand_code"] == cognate]
        if cog_rows.empty:
            continue
        # choose top-ranked cognate replicate
        cog_top = cog_rows.sort_values("ranking_score", ascending=False).iloc[0]
        pref_set = cache_std.get(Path(cog_top["structure_file"]).as_posix(), set())
        pref_size = len(pref_set)
        for idx, row in grp.iterrows():
            pocket = cache_std.get(Path(row["structure_file"]).as_posix(), set())
            if not pocket or not pref_set:
                df.at[idx, "iou_pref"] = math.nan
                df.at[idx, "pocket_status"] = df.at[idx, "pocket_status"] or "missing_pref"
                continue
            df.at[idx, "pref_size"] = pref_size
            df.at[idx, "iou_pref"] = compute_iou(pocket, pref_set)
    return df
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description="Compute pocket IoU metrics.")
    ap.add_argument("-i", "--metrics_in", default=str(METRICS_IN), help="Existing metrics_long.csv")
    ap.add_argument("-o", "--out", default=str(RESULTS_DIR / "metrics_long.csv"), help="Output CSV path")
    ap.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    args = ap.parse_args()
    setup_logger(args.verbose)
    metrics = pd.read_csv(args.metrics_in)
    updated = compute_competition_metrics(metrics)
    args_out = Path(args.out)
    args_out.parent.mkdir(parents=True, exist_ok=True)
    updated.to_csv(args_out, index=False)
    logging.info("Updated metrics with pocket IoU at %s", args_out)


if __name__ == "__main__":
    main()
