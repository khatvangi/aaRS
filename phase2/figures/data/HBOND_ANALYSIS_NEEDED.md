# H-Bond Analysis Documentation

## Requirements Not Met

BioPython is required for H-bond analysis but is not currently installed.

### To Install BioPython:

```bash
pip install biopython
# or
conda install -c conda-forge biopython
```

## Analysis Plan

### Structures to Analyze:

**For Figure 5D (Zinc Trap):**
1. Modern ThrRS + Zn + THR
   - Expect: 2-3 Zn coordination bonds (hydroxyl + amino)
   - Expect: 5-8 total H-bonds

2. Modern ThrRS + Zn + SER
   - Expect: Similar Zn coordination to THR (the trap!)
   - Should show identical coordination geometry

3. Modern ThrRS + Zn + ILE
   - Expect: 0-1 Zn coordination bonds (hydrophobic)
   - Fewer H-bonds overall

**For Figure 7B (Evolutionary Comparison):**
4. Ancestral ThrRS + THR (no Zn)
5. Ancestral ThrRS + Zn + THR
6. Ancestral ProRS + PRO
7. Ancestral ProRS editing + THR
8. Ancestral ProRS editing + PRO

### Expected Outputs:

1. **hbond_analysis.csv** - Summary table:
   - job_name, total_hbonds, has_zinc, zn_coordination_bonds, avg_distance

2. **hbond_analysis_detailed.json** - Full details:
   - Individual H-bond partners
   - Zn coordination geometry
   - Atom-by-atom breakdown

3. **Figure 7B** - Bar chart comparing:
   - H-bond counts across conditions
   - Zn coordination differences
   - Ancestral vs Modern evolution

### H-Bond Criteria:

- Distance: < 3.5 Å (donor-acceptor)
- Angle: > 120° (if hydrogens present)
- Zn coordination: < 2.8 Å

### Key Questions:

1. Why does SER coordinate Zn like THR?
   → Both have hydroxyl groups in same position

2. Why does ILE fail to coordinate?
   → Hydrophobic side chain cannot coordinate metal

3. Did H-bond network change during evolution?
   → Compare ancestral vs modern counts

## Alternative Analysis (Without BioPython)

If BioPython cannot be installed, manual analysis options:

1. **PyMOL**: Use distance measurements
   ```
   distance hbonds, (chain A), (chain B), 3.5
   ```

2. **ChimeraX**: Use built-in H-bond finding
   ```
   hbonds restrict both
   ```

3. **HBPLUS**: Standalone H-bond calculation tool
   ```
   hbplus structure.pdb
   ```

## Manual Inspection Priority:

Focus on these key structures:
1. modern_thrrs_ecoli_zn_THR vs SER (the trap!)
2. modern_thrrs_ecoli_zn_ILE (rejected)
3. anc_prors_edit_THR vs PRO (editing selectivity)
