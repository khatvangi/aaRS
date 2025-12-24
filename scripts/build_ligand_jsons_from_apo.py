#!/usr/bin/env python3
import json
from pathlib import Path

# ------------------------------------------------------------
# CONFIG: tell the script where your *apo* JSONs live and
# what ligand each new input should carry.
#
# FILL THESE PATHS with the actual apo JSONs that worked.
# ------------------------------------------------------------
JOBS = [
    {
        "apo_json": "phase2/apo_af3_inputs/proRS_LUCA.json",
        "out_json": "phase2/af3_ligand_inputs_fixed/proRS_LUCA_Pro.json",
        "name": "proRS_LUCA_Pro",
        "ligand_code": "PRO",
    },
    {
        "apo_json": "phase2/apo_af3_inputs/proRS_LUCA.json",
        "out_json": "phase2/af3_ligand_inputs_fixed/proRS_LUCA_Thr.json",
        "name": "proRS_LUCA_Thr",
        "ligand_code": "THR",
    },
    {
        "apo_json": "phase2/apo_af3_inputs/proRS_modern.json",
        "out_json": "phase2/af3_ligand_inputs_fixed/proRS_modern_Pro.json",
        "name": "proRS_modern_Pro",
        "ligand_code": "PRO",
    },
    {
        "apo_json": "phase2/apo_af3_inputs/proRS_modern.json",
        "out_json": "phase2/af3_ligand_inputs_fixed/proRS_modern_Thr.json",
        "name": "proRS_modern_Thr",
        "ligand_code": "THR",
    },
    {
        "apo_json": "phase2/apo_af3_inputs/proRS_editing.json",
        "out_json": "phase2/af3_ligand_inputs_fixed/proRS_editing_Pro.json",
        "name": "proRS_editing_Pro",
        "ligand_code": "PRO",
    },
    {
        "apo_json": "phase2/apo_af3_inputs/proRS_editing.json",
        "out_json": "phase2/af3_ligand_inputs_fixed/proRS_editing_Thr.json",
        "name": "proRS_editing_Thr",
        "ligand_code": "THR",
    },
]


def add_ligand_block(job_dict: dict, ligand_code: str, new_name: str | None = None) -> dict:
    """
    Take a single AF3 job dict (one 'name', 'dialect', 'sequences' etc.)
    and append the ligand block; optionally update the 'name'.
    """
    if "sequences" not in job_dict or not isinstance(job_dict["sequences"], list):
        raise ValueError("Job dict has no 'sequences' list; unexpected format.")

    ligand_entry = {
        "ligand": {
            "id": "L",
            "ccdCodes": [ligand_code],
        }
    }

    # Avoid duplicating ligand if you re-run the script
    already_has_ligand = any("ligand" in s for s in job_dict["sequences"])
    if not already_has_ligand:
        job_dict["sequences"].append(ligand_entry)

    if new_name is not None:
        job_dict["name"] = new_name

    # Ensure dialect/version are set correctly
    job_dict.setdefault("dialect", "alphafold3")
    job_dict.setdefault("version", 1)

    return job_dict


def process_one_job(apo_path: Path, out_path: Path, ligand_code: str, new_name: str | None):
    print(f"[INFO] Reading apo JSON: {apo_path}")
    with apo_path.open() as f:
        data = json.load(f)

    # Case 1: single job dict
    if isinstance(data, dict) and "sequences" in data:
        job = add_ligand_block(data, ligand_code, new_name=new_name)
        out = job

    # Case 2: top-level list of jobs
    elif isinstance(data, list):
        out_list = []
        for i, job in enumerate(data):
            if not isinstance(job, dict):
                raise ValueError(f"Entry {i} is not a dict in list JSON.")
            name_suffix = f"_{i}" if new_name is None else new_name
            out_list.append(add_ligand_block(job, ligand_code, new_name=name_suffix))
        out = out_list

    # Case 3: dict with 'foldJobs' list (AF server style)
    elif isinstance(data, dict) and "foldJobs" in data:
        jobs = data["foldJobs"]
        if not isinstance(jobs, list):
            raise ValueError("'foldJobs' is not a list.")
        for i, job in enumerate(jobs):
            if not isinstance(job, dict):
                raise ValueError(f"foldJobs[{i}] is not a dict.")
            name_suffix = f"_{i}" if new_name is None else new_name
            jobs[i] = add_ligand_block(job, ligand_code, new_name=name_suffix)
        out = data
    else:
        raise ValueError("Unrecognized apo JSON structure.")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        json.dump(out, f, indent=2)

    print(f"[OK] Wrote ligand JSON: {out_path}")


def main():
    base = Path(".")
    for job in JOBS:
        apo = base / job["apo_json"]
        out = base / job["out_json"]
        ligand_code = job["ligand_code"]
        new_name = job.get("name")
        if not apo.exists():
            print(f"[WARN] Apo JSON does not exist, skipping: {apo}")
            continue
        process_one_job(apo, out, ligand_code=ligand_code, new_name=new_name)


if __name__ == "__main__":
    main()

