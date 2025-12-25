Phase 3 Balanced AF3 Experiment
================================

Design: 6 conditions x 12 ligands x 5 seeds = 360 jobs

CONDITIONS
----------
Modern_ThrRS_Full_Zn
Modern_ThrRS_Full_noZn
Ancestral_ThrRS_Full_Zn
Ancestral_ThrRS_Full_noZn
Modern_ProRS_Full
Ancestral_ProRS_Full

LIGANDS
-------
THR SER PRO ALA GLY VAL ILE LEU CYS MET TYR GLU

SETUP STEPS
-----------

1. Update configs/conditions.yaml with real FASTA paths

2. Generate jobs:
   python3 make_balanced_jobs.py

3. Verify job count (should be 360):
   wc -l manifests/jobs.csv

4. Update run_balanced_af3.sh with your AF3 command

5. Run:
   bash run_balanced_af3.sh

DIRECTORY STRUCTURE
-------------------
configs/
  conditions.yaml    - condition definitions + FASTA paths
  ligands.txt        - 12 ligands

inputs/json/         - AF3 input JSONs (360 files)

manifests/
  jobs.csv           - job manifest

outputs/af3/         - AF3 output directories
outputs/af3/logs/    - per-job logs
