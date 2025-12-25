#!/usr/bin/env python3
"""
build_competition_matrix.py

Parse all COMPETITION AF3 runs and extract winner based on protein-ligand iptm.
The ligand with higher chain_pair_iptm[protein][ligand] is the winner.

Output: competition_tournament.csv with columns:
  - protein_construct, enzyme, era, cognate_ligand, competitor_ligand,
    cognate_iptm, competitor_iptm, winner, margin
"""

import json
import re
from pathlib import Path
import pandas as pd

# Find all competition runs
PHASE2_DIR = Path("/storage/kiran-stuff/aaRS/phase2")

def find_competition_runs():
    """Find all unique COMPETITION directories with summary_confidences.json"""
    runs = {}

    # Search patterns
    for pattern in ["**/COMPETITION*/*_summary_confidences.json",
                    "**/competition*/*_summary_confidences.json"]:
        for f in PHASE2_DIR.glob(pattern):
            # Skip seed-specific files, get the top-level summary
            if "seed-" in str(f.parent):
                continue

            # Get job name from filename
            job_name = f.stem.replace("_summary_confidences", "")

            # Skip duplicates (prefer first found)
            if job_name.upper() in [k.upper() for k in runs.keys()]:
                continue

            runs[job_name] = f

    return runs


def parse_competition_json(conf_path):
    """Parse competition summary_confidences.json and corresponding input JSON"""

    with open(conf_path) as f:
        conf = json.load(f)

    # Find the input JSON
    job_name = conf_path.stem.replace("_summary_confidences", "")
    job_dir = conf_path.parent

    # Try to find input JSON in parent directories
    input_json = None
    for parent in [job_dir.parent, job_dir.parent.parent, job_dir.parent.parent.parent]:
        candidate = parent / f"{job_name}.json"
        if candidate.exists():
            input_json = candidate
            break

    if not input_json:
        # Try _data.json in same directory
        data_json = job_dir / f"{job_name}_data.json"
        if data_json.exists():
            input_json = data_json

    # Parse input to get chain assignments
    chains = {}
    if input_json and input_json.exists():
        with open(input_json) as f:
            inp = json.load(f)

        for seq in inp.get("sequences", []):
            if "protein" in seq:
                chain_id = seq["protein"]["id"]
                if isinstance(chain_id, list):
                    chain_id = chain_id[0]
                chains[chain_id] = "protein"
            elif "ligand" in seq:
                chain_id = seq["ligand"]["id"]
                if isinstance(chain_id, list):
                    chain_id = chain_id[0]
                ccd = seq["ligand"]["ccdCodes"][0]
                chains[chain_id] = ccd

    # Get chain_pair_iptm
    chain_pair_iptm = conf.get("chain_pair_iptm", [])

    return chains, chain_pair_iptm, conf


def extract_matchup_info(job_name):
    """Extract enzyme, era, cognate, competitor from job name"""
    # Patterns:
    # COMPETITION_modern_thrrs_THR_vs_ILE
    # COMPETITION_anc_thrrs_THR_vs_ILE
    # COMPETITION_FULL_anc_prors_PRO_vs_THR

    job_upper = job_name.upper()

    # Extract era
    if "ANC" in job_upper or "ANCESTRAL" in job_upper:
        era = "ancestral"
    elif "MODERN" in job_upper:
        era = "modern"
    else:
        era = "unknown"

    # Extract enzyme
    if "THRRS" in job_upper:
        enzyme = "ThrRS"
    elif "PRORS" in job_upper:
        enzyme = "ProRS"
    else:
        enzyme = "unknown"

    # Extract cognate and competitor from "X_vs_Y" pattern (case insensitive)
    match = re.search(r"([A-Z]{3})_VS_([A-Z]{3})", job_upper)
    if match:
        cognate = match.group(1)
        competitor = match.group(2)
    else:
        cognate = competitor = "unknown"

    return enzyme, era, cognate, competitor


def determine_winner(chains, chain_pair_iptm):
    """
    Determine winner from chain_pair_iptm matrix.
    Returns (cognate_iptm, competitor_iptm, winner)

    Assumes:
    - Chain A = protein (index 0)
    - Chain B = first ligand (index 1)
    - Chain C = second ligand (index 2)
    """

    if not chain_pair_iptm or len(chain_pair_iptm) < 3:
        return None, None, "error"

    # Protein-ligand iptm scores
    # chain_pair_iptm[0] = protein row
    # [0][1] = protein-chainB, [0][2] = protein-chainC

    try:
        lig1_iptm = chain_pair_iptm[0][1]  # First ligand (cognate in name)
        lig2_iptm = chain_pair_iptm[0][2]  # Second ligand (competitor in name)

        if lig1_iptm is None or lig2_iptm is None:
            return lig1_iptm, lig2_iptm, "error"

        if lig1_iptm > lig2_iptm:
            winner = "cognate"
        elif lig2_iptm > lig1_iptm:
            winner = "competitor"
        else:
            winner = "tie"

        return lig1_iptm, lig2_iptm, winner

    except (IndexError, TypeError):
        return None, None, "error"


def main():
    runs = find_competition_runs()
    print(f"Found {len(runs)} unique competition runs")

    results = []

    for job_name, conf_path in sorted(runs.items()):
        print(f"\nProcessing: {job_name}")

        try:
            chains, chain_pair_iptm, conf = parse_competition_json(conf_path)
            enzyme, era, cognate, competitor = extract_matchup_info(job_name)
            cog_iptm, comp_iptm, winner = determine_winner(chains, chain_pair_iptm)

            print(f"  Chains: {chains}")
            print(f"  {cognate} iptm: {cog_iptm}, {competitor} iptm: {comp_iptm} -> {winner}")

            results.append({
                "job_name": job_name,
                "protein_construct": job_name.replace("COMPETITION_", "").replace("FULL_", "").split("_" + cognate)[0],
                "enzyme": enzyme,
                "era": era,
                "cognate_ligand": cognate,
                "competitor_ligand": competitor,
                "cognate_iptm": cog_iptm,
                "competitor_iptm": comp_iptm,
                "margin": (cog_iptm - comp_iptm) if cog_iptm and comp_iptm else None,
                "winner": winner,
                "cognate_wins": 1 if winner == "cognate" else 0,
                "global_iptm": conf.get("iptm"),
                "ptm": conf.get("ptm"),
                "ranking_score": conf.get("ranking_score"),
                "source_file": str(conf_path)
            })

        except Exception as e:
            print(f"  ERROR: {e}")
            continue

    # Create DataFrame
    df = pd.DataFrame(results)

    # Save
    outpath = Path("/storage/kiran-stuff/aaRS/phase3_balanced/competition_tournament.csv")
    df.to_csv(outpath, index=False)
    print(f"\n\nSaved {len(df)} results to {outpath}")

    # Summary
    print("\n" + "="*60)
    print("TOURNAMENT SUMMARY")
    print("="*60)

    if len(df) > 0:
        # Group by enzyme and era
        for (enzyme, era), grp in df.groupby(["enzyme", "era"]):
            n_matchups = len(grp)
            n_wins = grp["cognate_wins"].sum()
            win_rate = n_wins / n_matchups if n_matchups > 0 else 0

            print(f"\n{era.capitalize()} {enzyme}:")
            print(f"  Matchups: {n_matchups}")
            print(f"  Cognate wins: {n_wins} ({win_rate:.0%})")

            for _, row in grp.iterrows():
                marker = "✓" if row["winner"] == "cognate" else "✗"
                print(f"    {marker} {row['cognate_ligand']} vs {row['competitor_ligand']}: "
                      f"{row['cognate_iptm']:.2f} vs {row['competitor_iptm']:.2f} (Δ={row['margin']:.2f})")

    return df


if __name__ == "__main__":
    df = main()
