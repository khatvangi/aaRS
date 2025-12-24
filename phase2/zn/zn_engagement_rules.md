# Zn Engagement Rules

**Cutoff:** zn_min_dist_hetero ≤ 3.0 Å (Zn to ligand O/N/S atoms)

---

## Classification

| Category | Count | Percentage |
|----------|-------|------------|
| **Engaged** | 23 | 48.9% |
| **Floating** | 24 | 51.1% |

---

## Distance Distribution

- **Mean:** 14.51 ± 13.62 Å
- **Median:** 3.74 Å
- **Range:** [1.70, 32.47] Å

---

## Engagement by Condition

| condition                   |   n_engaged |   n_floating |   pct_engaged |   mean_zn_dist |
|:----------------------------|------------:|-------------:|--------------:|---------------:|
| Modern_ThrRS_Domain_Zn      |           2 |            0 |           100 |        2.11579 |
| Unknown_Other_Domain_Zn     |           2 |            0 |           100 |        2.13761 |
| Modern_ThrRS_Full-length_Zn |          19 |            1 |            95 |        2.36895 |
| Ancestral_ThrRS_Catalytic   |           0 |           20 |             0 |       28.9348  |
| Ancestral_ThrRS_Domain      |           0 |            2 |             0 |       21.8728  |
| Modern_ThrRS_Full-length    |           0 |            1 |             0 |        3.74448 |

---

## Interpretation

### Engaged Zn (≤3.0 Å)
- Zn is within coordination distance of ligand heteroatoms
- Likely represents true Zn-ligand interaction
- **USE for mechanistic interpretation**

### Floating Zn (>3.0 Å)
- Zn is distant from ligand (mean ~10-15 Å)
- AF3 modeling artifact - Zn placed for protein but not ligand
- **EXCLUDE from biological claims**

---

## Recommendation

**Only use structures with `zn_engaged=True` for:**
- Claims about Zn-mediated binding
- Comparisons of Zn effect on selectivity
- Mechanistic interpretation of metal coordination

**Flag `zn_floating` structures as:**
- AF3 structural predictions with incomplete Zn coordination
- Not suitable for Zn-ligand binding analysis
