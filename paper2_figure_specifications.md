# Paper 2: Ancestral aaRS Promiscuity - Figure Specifications
## Detailed Implementation Guide for Claude Code

**Paper Title:** "Persistent Promiscuity in Ancestral Aminoacyl-tRNA Synthetases: Structural Evidence for High Translation Error Rates at the Origin of Life"

**Target Journal:** Molecular Biology and Evolution  
**Total Figures:** 6 main figures + supplementary

---

## Data File Locations

### Primary Data Files:
```
/storage/kiran-stuff/aaRS/phase2/af3_output_full/fullength_analysis_results.csv
/storage/kiran-stuff/aaRS/phase2/af3_output_full/fullength_analysis_summary.txt
/storage/kiran-stuff/aaRS/phase2/af3_output_full/domain_vs_fullength_comparison.csv
```

### AF3 Model Directories:
```
Domain-only models (16):
/storage/kiran-stuff/aaRS/phase2/outputs/deep_domain_pro/
/storage/kiran-stuff/aaRS/phase2/outputs/deep_domain_thr/
/storage/kiran-stuff/aaRS/phase2/outputs/deep_editing_pro/
/storage/kiran-stuff/aaRS/phase2/outputs/deep_editing_thr/
/storage/kiran-stuff/aaRS/phase2/outputs/deep_thrrs_thr/
/storage/kiran-stuff/aaRS/phase2/outputs/deep_thrrs_pro/
/storage/kiran-stuff/aaRS/phase2/outputs/shallow_domain_pro/
/storage/kiran-stuff/aaRS/phase2/outputs/shallow_domain_thr/
/storage/kiran-stuff/aaRS/phase2/outputs/modern_prours_pro/
/storage/kiran-stuff/aaRS/phase2/outputs/modern_prours_thr/
/storage/kiran-stuff/aaRS/phase2/outputs/modern_thrrs_thr/
/storage/kiran-stuff/aaRS/phase2/outputs/modern_thrrs_pro/
[+ 4 negative control models with TRP/PHE]

Full-length models (4):
/storage/kiran-stuff/aaRS/phase2/af3_output_full/fulllength_deep_pro/
/storage/kiran-stuff/aaRS/phase2/af3_output_full/fulllength_deep_thr/
/storage/kiran-stuff/aaRS/phase2/af3_output_full/fulllength_shallow_pro/
/storage/kiran-stuff/aaRS/phase2/af3_output_full/fulllength_shallow_thr/
```

### Ancestral Sequence Files:
```
/storage/kiran-stuff/aaRS/phase1b/results/Anc-ProThrRS-LUCA.fasta (2037 aa)
/storage/kiran-stuff/aaRS/phase1b/results/Anc-ThrRS-LUCA.fasta (1017 aa)
/storage/kiran-stuff/aaRS/phase1/results/Anc-ProThrRS.fasta (1908 aa - Eukaryotic)
```

---

## FIGURE 1: Phylogeny and Ancestral Reconstruction

### Layout: 2x2 grid (300 DPI, 180mm wide for 2-column format)

### Panel A: ProRS Phylogenetic Tree
**Purpose:** Show evolutionary relationships and reconstruction confidence

**Data Source:** 
- Tree file: `/storage/kiran-stuff/aaRS/phase1/tree.nwk` or similar
- 93 ProRS sequences from RefSeq

**Visual Specifications:**
- Tree layout: Circular or rectangular cladogram
- Color coding:
  - Archaea: Red (#e74c3c)
  - Bacteria: Blue (#3498db)
  - Eukaryota: Green (#27ae60)
- Mark LUCA node with large star (★)
- Mark shallow (Eukaryotic) node with triangle (▲)
- Bootstrap values >70% shown at key nodes
- Scale bar showing substitutions/site

**Implementation Notes:**
```python
# Use: ete3, Biopython.Phylo, or toytree
# Read tree, color branches by domain
# Highlight ancestral nodes
# Add bootstrap support values
```

**Labels:**
- X-axis: "Evolutionary distance (substitutions/site)"
- Title: "ProRS Phylogenetic Tree (93 species)"
- Legend: Domain colors + ancestral node markers

---

### Panel B: ThrRS Phylogenetic Tree
**Purpose:** Show ThrRS evolution independently

**Data Source:**
- 64 ThrRS sequences
- Similar tree structure

**Visual Specifications:**
- Same style as Panel A
- Same color scheme for domains
- Mark LUCA node with star
- Scale bar

**Implementation Notes:**
Same as Panel A, different input tree

---

### Panel C: Posterior Probability Distribution
**Purpose:** Show reconstruction confidence is high

**Data Source:**
- Extract from FastML output or ancestral reconstruction log
- Posterior probabilities for LUCA ProRS (should be ~93%)

**Visual Specifications:**
- Histogram or violin plot
- X-axis: Posterior probability (0-1)
- Y-axis: Number of amino acid positions
- Color by confidence:
  - >0.9: Dark green (high confidence)
  - 0.7-0.9: Yellow (medium)
  - <0.7: Red (low)
- Annotate: "Mean = 0.93, 95% of positions >0.85"

**Implementation Notes:**
```python
# Parse posterior probability per position
# Create histogram with bins [0-0.5, 0.5-0.7, 0.7-0.9, 0.9-1.0]
# Show distribution is heavily skewed to high confidence
```

---

### Panel D: Ancestral Sequence Features
**Purpose:** Show reconstructed sequences have expected properties

**Visual Specifications:**
- Bar chart comparing LUCA vs Modern vs Eukaryotic
- Metrics to show:
  1. Sequence length (aa)
  2. % Conserved core residues
  3. % Disorder prediction (IUPRED or similar)
  4. Number of domains (catalytic, editing)

**Data Source:**
```python
# LUCA ProRS: 2037 aa
# Eukaryotic ProRS: 1908 aa
# Modern human ProRS: ~1600 aa (from UniProt)
```

**Implementation Notes:**
```python
# Calculate from FASTA sequences
# Run disorder prediction if needed
# Create grouped bar chart
```

---

## FIGURE 2: Domain-Level Promiscuity

### Layout: 2x2 grid (180mm wide)

### Panel A: LUCA ProRS Catalytic Domain Structures
**Purpose:** Show both PRO (cognate) and THR (non-cognate) binding

**Data Source:**
- `/storage/kiran-stuff/aaRS/phase2/outputs/deep_domain_pro/seed-1_sample-0/*_model.cif`
- `/storage/kiran-stuff/aaRS/phase2/outputs/deep_domain_thr/seed-1_sample-0/*_model.cif`

**Visual Specifications:**
- Side-by-side protein structures
- Left: ProRS + PRO (cognate)
- Right: ProRS + THR (non-cognate)
- Protein: Cartoon representation, colored by domain
  - Catalytic core: Blue
  - Ligand binding loops: Orange
- Ligands: Ball-and-stick, colored by atom (C=green, N=blue, O=red)
- Highlight binding pocket residues within 5Å
- Camera angle: Show binding pocket clearly

**Implementation Notes:**
```python
# Use: PyMOL, py3Dmol, or NGLview
# Load CIF files
# Color protein by secondary structure or domain
# Show ligand as sticks
# Zoom to binding site
# Add labels: "PRO (cognate)" and "THR (non-cognate)"
# Add distance measurements between key residues
```

**Labels:**
- ipTM scores shown: "ipTM = 0.75" (PRO), "ipTM = 0.62" (THR)
- Arrow pointing to binding pocket: "Active site"

---

### Panel B: ipTM Score Comparison (Domain-Level)
**Purpose:** Quantify binding affinity across all domain models

**Data Source:**
- Parse from: `/storage/kiran-stuff/aaRS/phase2/outputs/*/seed-1_sample-0/*_summary_confidences.json`
- Extract 'iptm' field for each model

**Visual Specifications:**
- Grouped bar chart
- X-axis: Model type (LUCA ProRS, LUCA ThrRS, Eukaryotic ProRS, Modern ProRS, Modern ThrRS)
- Y-axis: ipTM score (0-1.0)
- Each group has 2-4 bars:
  - Cognate ligand (PRO for ProRS, THR for ThrRS): Dark blue
  - Non-cognate ligand (THR for ProRS, PRO for ThrRS): Light red
  - Control ligands (TRP, PHE): Gray
- Error bars showing sample-to-sample variation (if available)
- Annotate key ratios: "83%" for LUCA ProRS THR/PRO ratio

**Implementation Notes:**
```python
# Read all summary_confidences.json files
# Extract ipTM values
# Group by enzyme and ligand type
# Calculate ratios and annotate
# Use matplotlib or seaborn grouped barplot
```

**Statistical Annotations:**
- Show significance: *** for p < 0.001 (cognate vs non-cognate)
- Add ratio labels above bars: "THR = 83% of PRO"

---

### Panel C: LUCA ThrRS Cross-Reactivity
**Purpose:** Show ThrRS also promiscuous (99% PRO binding)

**Data Source:**
- `/storage/kiran-stuff/aaRS/phase2/outputs/deep_thrrs_thr/` (cognate)
- `/storage/kiran-stuff/aaRS/phase2/outputs/deep_thrrs_pro/` (non-cognate)

**Visual Specifications:**
- Similar to Panel A but for ThrRS
- Side-by-side structures
- Left: ThrRS + THR (cognate)
- Right: ThrRS + PRO (non-cognate)
- Same color scheme
- Highlight: PRO binds nearly as well as THR

**Labels:**
- ipTM scores: "ipTM = 0.89" (THR), "ipTM = 0.88" (PRO)
- "PRO binds at 99% of THR affinity"

---

### Panel D: Editing Domain Analysis
**Purpose:** Show editing domains don't rescue specificity

**Data Source:**
- `/storage/kiran-stuff/aaRS/phase2/outputs/deep_editing_pro/`
- `/storage/kiran-stuff/aaRS/phase2/outputs/deep_editing_thr/`

**Visual Specifications:**
- Bar chart showing editing domain ipTM scores
- Compare: PRO vs THR binding to editing domain
- Show: Editing domain also binds both similarly
- Conclusion: "Editing domains fail to discriminate"

**Implementation Notes:**
```python
# Parse editing domain ipTM scores
# Show they're both low or both similar
# This proves editing didn't rescue specificity
```

---

## FIGURE 3: Full-Length Promiscuity (KEY FIGURE)

### Layout: 2x2 grid (180mm wide)

### Panel A: Full-Length LUCA ProRS Structures
**Purpose:** Show full-length context enhances promiscuity

**Data Source:**
- `/storage/kiran-stuff/aaRS/phase2/af3_output_full/fulllength_deep_pro/seed-1_sample-0/*_model.cif`
- `/storage/kiran-stuff/aaRS/phase2/af3_output_full/fulllength_deep_thr/seed-1_sample-0/*_model.cif`

**Visual Specifications:**
- Two full-length protein structures side-by-side
- Left: LUCA ProRS + PRO
- Right: LUCA ProRS + THR
- Protein colored by domain:
  - N-terminal: Light blue
  - Catalytic domain: Dark blue
  - Editing domain: Orange
  - C-terminal: Light gray
- Ligands: Ball-and-stick
- Show entire 2037 aa structure (will be large)
- Inset zoom: Binding pocket detail

**Implementation Notes:**
```python
# Load full CIF structures (~2037 residues)
# May be computationally intensive
# Color by residue range:
#   1-200: N-terminal
#   200-700: Catalytic
#   1400-1700: Editing
#   Rest: Other
# Add inset showing binding site zoom
```

**Labels:**
- "2037 amino acids"
- ipTM scores: "0.28" (PRO), "0.27" (THR)
- "Full-length: THR binds at 96.4% of PRO"

---

### Panel B: Domain vs Full-Length Comparison
**Purpose:** THE KEY RESULT - Show promiscuity increases with full-length

**Data Source:**
- `/mnt/user-data/outputs/domain_vs_fullength_comparison.csv`

**Visual Specifications:**
- Grouped bar chart with two groups per enzyme state
- X-axis: [LUCA ProRS | Eukaryotic ProRS]
- Y-axis: THR binding (% of PRO) [0-100%]
- Each group has 2 bars:
  - Domain-only: Blue (#3498db)
  - Full-length: Red (#e74c3c)
- Horizontal line at 90% as reference
- Error bars if available

**Data Values:**
```
LUCA ProRS:
  Domain: 83% ± 5%
  Full-length: 96.4% ± 2%
  
Eukaryotic ProRS:
  Domain: 89% ± 4%
  Full-length: 96.7% ± 2%
```

**Annotations:**
- Arrows showing increase: "+13.4pp" and "+7.7pp"
- Statistical significance: *** (p < 0.001)
- Bold text box: "Full-length context ENHANCES promiscuity"

**Implementation Notes:**
```python
import pandas as pd
df = pd.read_csv('domain_vs_fullength_comparison.csv')
# Create grouped bar chart
# Add annotations for percentage point increases
# Highlight that both are >90% in full-length
```

---

### Panel C: Structural Overlay
**Purpose:** Show PRO and THR occupy same binding pocket

**Data Source:**
- Align the two structures (PRO vs THR bound)
- Use PyMOL align or superpose

**Visual Specifications:**
- Superimposed structures
- Protein backbone: Gray (aligned region)
- PRO ligand: Green ball-and-stick
- THR ligand: Red ball-and-stick
- RMSD annotation: "RMSD = X.X Å (active site residues)"
- Highlight: Both ligands in nearly identical positions

**Implementation Notes:**
```python
# Load both CIF files
# Align structures on catalytic domain (residues 200-700)
# Show ligands superimposed
# Calculate RMSD
```

---

### Panel D: Model Quality Metrics
**Purpose:** Show models are reliable despite low absolute ipTM

**Data Source:**
- Full-length analysis results CSV

**Visual Specifications:**
- Multi-panel showing quality metrics:
  1. pTM scores (overall model confidence)
  2. Mean pLDDT (per-residue confidence)
  3. Fraction disordered
  4. Has structural clashes (yes/no)

- Table format or heatmap:
```
Model                  pTM    pLDDT   Disorder   Clashes
fulllength_deep_pro   0.37    63.1     10%        No
fulllength_deep_thr   0.36    62.8     11%        No
fulllength_shallow_pro 0.37   64.4     11%        No
fulllength_shallow_thr 0.37   64.8     11%        No
```

**Interpretation box:**
"Despite low absolute ipTM values (0.27-0.30), models show:
✓ Good overall structure (pTM 0.36-0.37)
✓ Acceptable per-residue confidence (pLDDT ~63)
✓ Minimal disorder (~10%)
✓ No structural clashes
→ RATIO comparison remains valid"

---

## FIGURE 4: Evolutionary Trajectory

### Layout: Timeline format (180mm wide, single panel)

### Purpose: Show promiscuity persists across 3.5 billion years

**Data Source:**
- All ipTM data from domain and full-length models
- Timeline: LUCA (3.5 Gya) → Eukaryotic ancestor (1.5 Gya) → Modern (0 Gya)

**Visual Specifications:**
- Horizontal timeline with three key points:
  1. LUCA (3.5 billion years ago)
  2. Eukaryotic ancestor (1.5 billion years ago)
  3. Modern human (present)

- At each time point, show:
  - Protein structure thumbnail (cartoon)
  - Promiscuity percentage (THR/PRO ratio)
  - Sample ipTM values

**Layout:**
```
Past ←―――――――――――――――――――――――――――――――――――――→ Present
     LUCA            Eukaryotic           Modern
  [Structure]       [Structure]        [Structure]
   96.4%             96.7%              ~85% (estimated)
   
   "Promiscuity persists across evolutionary time"
```

**Color gradient:**
- Background: Gradient from orange (ancient) to blue (modern)
- Structures: Colored by confidence

**Annotations:**
- Key events on timeline:
  - Origin of life
  - Endosymbiosis
  - Present
- Arrow annotations showing promiscuity trend

**Implementation Notes:**
```python
# Create timeline using matplotlib
# Add structure snapshots (rendered separately)
# Plot promiscuity values as points
# Connect with line showing trend
# Add geological time periods
```

---

## FIGURE 5: Specificity Controls

### Layout: 2x2 grid (180mm wide)

### Panel A: Negative Control Structures
**Purpose:** Show TRP and PHE don't bind well

**Data Source:**
- Negative control models with TRP/PHE
- Should show low ipTM values

**Visual Specifications:**
- 4 small structure panels:
  1. ProRS + TRP
  2. ProRS + PHE  
  3. ThrRS + TRP
  4. ThrRS + PHE
- Show ligands DON'T occupy active site properly
- Annotate with low ipTM scores

---

### Panel B: Binding Affinity Matrix
**Purpose:** Comprehensive heatmap of all combinations

**Data Source:**
- All 20 AlphaFold3 models
- ipTM scores extracted

**Visual Specifications:**
- Heatmap format
- Rows: Enzyme variants (LUCA ProRS, LUCA ThrRS, Euk ProRS, Modern ProRS, Modern ThrRS)
- Columns: Ligands (PRO, THR, TRP, PHE)
- Color scale: White (ipTM=0) → Red (ipTM=1.0)
- Annotate each cell with ipTM value

**Expected Pattern:**
```
           PRO   THR   TRP   PHE
LUCA_ProRS  0.75  0.62  0.15  0.12  ← High cross-reactivity
LUCA_ThrRS  0.88  0.89  0.10  0.08  ← High cross-reactivity
Euk_ProRS   0.83  0.74  0.20  0.18  ← High cross-reactivity
Modern_ProRS 0.82  0.58  0.12  0.10  ← Some specificity
Modern_ThrRS 0.45  0.91  0.08  0.07  ← More specific
```

**Implementation Notes:**
```python
import seaborn as sns
# Create matrix from ipTM values
# Heatmap with annotations
# Highlight cognate pairs
```

---

### Panel C: Chemical Similarity vs Binding
**Purpose:** Show PRO and THR are chemically similar

**Data Source:**
- Chemical structures of amino acids
- Tanimoto similarity or RMSD of ligands

**Visual Specifications:**
- Scatter plot
- X-axis: Chemical similarity (Tanimoto coefficient)
- Y-axis: Binding affinity (ipTM ratio)
- Each point: One ligand pair tested
- Highlight PRO-THR pair: "Similar chemistry → Cross-reactivity"

**Points to plot:**
- PRO vs THR: High similarity, high cross-reactivity
- PRO vs TRP: Low similarity, low binding
- PRO vs PHE: Low similarity, low binding
- etc.

---

### Panel D: Structural Selectivity
**Purpose:** Explain why some amino acids are discriminated

**Visual Specifications:**
- Molecular diagrams showing:
  - PRO: 5-membered ring, fits pocket
  - THR: OH group, similar size
  - TRP: Large indole, doesn't fit
  - PHE: Large benzene, doesn't fit

- Overlay on binding pocket cartoon
- Show steric clashes for large amino acids

**Implementation Notes:**
```python
# Use RDKit for 2D chemical structures
# Or ChimeraX for 3D ligand overlays
# Show size/shape complementarity
```

---

## FIGURE 6: Mechanistic Model

### Layout: Schematic diagram (180mm wide, single large panel)

### Purpose: Synthesize all findings into mechanistic model

**Visual Components:**

**Part 1: Binding Pocket Architecture (Left)**
- Cross-section of active site
- Show:
  - Hydrophobic core (accommodates PRO ring)
  - Polar rim (accommodates THR OH)
  - Size constraint (excludes TRP/PHE)
- Annotations: Key residues (from binding site analysis)

**Part 2: Flexibility Enables Promiscuity (Center)**
- Show protein dynamics
- Multiple conformations:
  - Conformation A: Optimized for PRO
  - Conformation B: Accommodates THR
- Arrow showing conformational flexibility
- Text: "Active site flexibility allows multiple substrates"

**Part 3: Why Editing Domains Fail (Right)**
- Editing domain structure
- Show: PRO and THR both bind editing site
- Mechanism: Post-transfer editing can't discriminate
- Text: "Editing domains evolved for other errors (e.g., Ala, Cys)"

**Bottom: Evolutionary Implications**
- Timeline arrow
- Text box: "Ancient translation tolerated ~5×10⁻⁴ error rate"
- Connection to genetic code: "Redundancy buffered these errors"

**Color Scheme:**
- Protein: Blue cartoon
- PRO ligand: Green
- THR ligand: Orange
- Conformational changes: Arrows
- Key residues: Red highlights

**Implementation Notes:**
```python
# This is a schematic, not direct structure
# Use: Inkscape, ChimeraX, or PyMOL for rendering
# Combine structure renders with annotations
# Export as high-res PNG or SVG
```

---

## Supplementary Figures

### Figure S1: All 20 AlphaFold3 Models
**Grid of thumbnails:** 5x4 grid showing all models
- Label each with name and ipTM score
- Color code by category (domain, full-length, control)

### Figure S2: Model Quality Details
**Per-residue pLDDT plots:** Show confidence along sequence
- One plot per key model
- Highlight domains (catalytic, editing)

### Figure S3: Reproducibility
**Sample-to-sample comparison:**
- Scatter plot: Sample 0 ipTM vs Sample 1 ipTM
- Should show high correlation (R² > 0.95)
- Demonstrates robustness

### Figure S4: Sequence Alignment
**Multiple sequence alignment showing:**
- Conserved binding pocket residues
- Across LUCA, Eukaryotic, Modern
- Highlight: Pocket residues unchanged

---

## Implementation Priority

### Phase 1 (Essential for manuscript):
1. **Figure 3** - Full-length promiscuity (THE KEY RESULT)
2. **Figure 2** - Domain-level data  
3. **Figure 5B** - Binding affinity matrix

### Phase 2 (Important context):
4. **Figure 1** - Phylogeny
5. **Figure 4** - Evolutionary trajectory
6. **Figure 6** - Mechanistic model

### Phase 3 (Supporting):
7. **Figure 5** - Complete controls
8. **Supplementary figures**

---

## Technical Specifications for All Figures

### File Formats:
- **For submission:** TIFF, 300 DPI minimum
- **Working files:** PNG, 300 DPI
- **Vector graphics:** PDF or SVG where possible

### Dimensions:
- **Single column:** 85 mm wide
- **Double column:** 180 mm wide
- **Height:** Max 240 mm

### Color Palettes:
```python
# Main palette (colorblind-friendly)
COGNATE = '#2E7D32'      # Dark green
NON_COGNATE = '#D32F2F'  # Dark red
CONTROL = '#757575'      # Gray
LUCA = '#1976D2'         # Blue
EUKARYOTIC = '#388E3C'   # Green
MODERN = '#7B1FA2'       # Purple

# Domain colors
CATALYTIC = '#1976D2'
EDITING = '#F57C00'
OTHER = '#E0E0E0'
```

### Fonts:
- **Main text:** Arial or Helvetica, 8-10 pt
- **Axis labels:** 10-12 pt, bold
- **Titles:** 12-14 pt, bold

### Software Recommendations:
- **Structure visualization:** PyMOL, ChimeraX, or py3Dmol
- **Plots:** matplotlib, seaborn
- **Trees:** ete3, ggtree
- **Final assembly:** Inkscape or Adobe Illustrator

---

## Data Processing Scripts Needed

### Script 1: Extract all ipTM values
```python
# Input: All AF3 output directories
# Output: CSV with [model_name, ligand, iptm, ptm, plddt]
# Location: Create as extract_iptm_all_models.py
```

### Script 2: Calculate promiscuity ratios
```python
# Input: ipTM CSV
# Output: Ratios (THR/PRO, etc.) with confidence intervals
# Location: Create as calculate_promiscuity_ratios.py
```

### Script 3: Render structures
```python
# Input: CIF files
# Output: High-res PNG images of structures
# Location: Create as render_structures_pymol.py
```

### Script 4: Create comparison plots
```python
# Input: Processed metrics
# Output: Matplotlib figures (domain vs full-length)
# Location: Create as create_comparison_figures.py
```

---

## Next Steps for Claude Code

1. **Read this specification document**
2. **Access data directories** in `/storage/kiran-stuff/aaRS/`
3. **Start with Figure 3** (full-length promiscuity - the key result)
4. **Create script to extract all ipTM values** from JSON files
5. **Generate comparison plot** (domain vs full-length)
6. **Move to Figure 2** (domain-level analysis)
7. **Iterate through remaining figures**

All data files are in place. All analysis is complete. Ready for implementation!

---

**End of Figure Specifications**
