Pipelines to run the two mini-analyses locally (tested on the `/storage/kiran-stuff` layout you have on boron).

## 0. One-time setup
- PoincaréMSA repo is already at `/storage/kiran-stuff/PoincareMSA`.
- Clone SWORD2 next to it if it is not present:
  - `cd /storage/kiran-stuff && git clone https://github.com/DSIMB/SWORD2.git`
- Create the two conda envs:
  - `conda env create -f /storage/kiran-stuff/PoincareMSA/env_poincare.yml -n env_poincare`
  - `conda env create -f /storage/kiran-stuff/SWORD2/environment.yml -n sword2`
  - After the SWORD2 env is created, run its post-install once:  
    `cd /storage/kiran-stuff/SWORD2 && conda activate sword2 && bash install.sh`

## 1. PoincaréMSA on ProRS
- Input MSA: `aaRS/phase1/data/interim/ProRS_aligned.fasta` (already aligned; copied as `.mfasta` by the script).
- Output root: `aaRS/poincare_runs/proRS/`.
- Run:
  - `bash aaRS/pipelines/run_poincare_proRS.sh`
- What happens:
  - Activates `env_poincare`.
  - Copies the aligned fasta to `proRS_full.mfasta` under the run folder.
  - Runs `scripts/prepare_data/create_projection.sh` with gap threshold `0.9` to produce PSSMs in `fasta0.9/`.
  - Trains the 2D Poincaré embedding with `knn=5`, `gamma=2`, `batchsize=4`, `epochs=1000`.
  - Outputs CSV + PNG into `aaRS/poincare_runs/proRS/poincare_knn5_gamma2/` (CSV columns: `pm1`, `pm2`, `proteins_id`).

## 2. SWORD2 on AF3 ProRS models (batch)
- Repo path assumed at `/storage/kiran-stuff/SWORD2`.
- Input structures (editable list) live under `aaRS/phase2/outputs/*` (AlphaFold3 `.cif` files).
- Output root: `aaRS/sword2_results/`.
- Run:
  - `bash aaRS/pipelines/run_sword2_batch.sh`
- What happens:
  - Activates `sword2` env.
  - Converts listed `.cif` files to chain-A PDBs under `aaRS/sword2_results/pdb_inputs/`.
  - Runs `SWORD2.py` with `--disable-energies --disable-plots` into `aaRS/sword2_results/runs/` (each chain gets its own folder with `SWORD2_summary.json`).

## Notes / tweaks
- Adjust hyperparameters by exporting env vars before running, e.g. `EPOCHS=600` or `KNN=8` for Poincaré, or edit the `STRUCTURES` array in `run_sword2_batch.sh`.
- If CUDA is visible and you want to force CPU for PoincaréMSA, export `CUDA_VISIBLE_DEVICES=""` before running the script.
- Both scripts are idempotent: rerunning will overwrite outputs in the same run folder.
