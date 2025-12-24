#!/usr/bin/env python3
import json
import csv
from pathlib import Path

# Adjust this to where your AF3 ligand outputs live
BASE = Path("af3_output")

# Map job labels to their summary paths
JOBS = {
    # edit domain
    "deep_editing_pro": BASE / "deep_editing_pro_20251208_175120" / "deep_editing_pro_summary_confidences.json",

    # fill these once LUCA / modern / Thr runs are done:
    # "fulllength_LUCA_Pro": BASE / "fulllength_LUCA_Pro_YYYYMMDD_HHMMSS" / "fulllength_LUCA_Pro_summary_confidences.json",
    # "fulllength_LUCA_Thr": ...
    # "fulllength_modern_Pro": ...
    # "fulllength_modern_Thr": ...
    # "deep_editing_thr": ...
}

OUT_CSV = Path("summary/af3_ligand_ipTM_summary.csv")


def parse_summary(path: Path, job_name: str) -> dict:
    with path.open() as f:
        s = json.load(f)

    chain_iptm = s.get("chain_iptm", [])
    chain_ptm = s.get("chain_ptm", [])
    chain_pair_iptm = s.get("chain_pair_iptm", [])

    n_chain = len(chain_ptm)
    if n_chain == 0:
        raise ValueError(f"{job_name}: no chains found in summary.")

    # Heuristic: ligand is the *last* chain (small molecule, very low ptm)
    ligand_idx = n_chain - 1

    # Protein/RNA chains are all others
    macromol_indices = list(range(0, ligand_idx))

    # Compute max protein–ligand ipTM across macromolecular chains
    protein_ligand_iptm = None
    if chain_pair_iptm and len(chain_pair_iptm) > ligand_idx:
        vals = []
        for i in macromol_indices:
            if i < len(chain_pair_iptm) and ligand_idx < len(chain_pair_iptm[i]):
                vals.append(chain_pair_iptm[i][ligand_idx])
        if vals:
            protein_ligand_iptm = max(vals)

    # We can also compute mean protein cp-iptm if desired
    protein_protein_iptm = None
    if len(macromol_indices) > 1:
        vals_pp = []
        for i in macromol_indices:
            for j in macromol_indices:
                if i == j:
                    continue
                if i < len(chain_pair_iptm) and j < len(chain_pair_iptm[i]):
                    vals_pp.append(chain_pair_iptm[i][j])
        if vals_pp:
            protein_protein_iptm = sum(vals_pp) / len(vals_pp)

    row = {
        "job": job_name,
        "n_chains": n_chain,
        "iptm": s.get("iptm"),
        "ptm": s.get("ptm"),
        "fraction_disordered": s.get("fraction_disordered"),
        "has_clash": s.get("has_clash"),
        "ranking_score": s.get("ranking_score"),
        "ligand_chain_index": ligand_idx,
        "ligand_chain_ptm": chain_ptm[ligand_idx] if ligand_idx < len(chain_ptm) else None,
        "chain_iptm_protein_0": chain_iptm[0] if len(chain_iptm) > 0 else None,
        "chain_iptm_protein_1": chain_iptm[1] if len(chain_iptm) > 1 else None,
        "protein_ligand_iptm_max": protein_ligand_iptm,
        "protein_protein_iptm_mean": protein_protein_iptm,
    }
    return row


def main():
    OUT_CSV.parent.mkdir(exist_ok=True, parents=True)
    rows = []

    for job_name, path in JOBS.items():
        if not path.exists():
            print(f"[WARN] Missing {path} for job {job_name}, skipping.")
            continue
        print(f"[INFO] Parsing {job_name} from {path}")
        row = parse_summary(path, job_name)
        rows.append(row)

    if not rows:
        print("[WARN] No rows parsed; nothing to write.")
        return

    fieldnames = [
        "job", "n_chains",
        "iptm", "ptm",
        "fraction_disordered", "has_clash", "ranking_score",
        "ligand_chain_index", "ligand_chain_ptm",
        "chain_iptm_protein_0", "chain_iptm_protein_1",
        "protein_ligand_iptm_max",
        "protein_protein_iptm_mean",
    ]

    with OUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] Wrote {OUT_CSV} with {len(rows)} rows.")


if __name__ == "__main__":
    main()

