# 🎯 Complete Script Package for Manuscript Figure Generation

## 📦 What You Have

All scripts are ready to run on your computer. Here's what each does:

### **Core Scripts**

1. **`run_figure_pipeline.sh`** ⭐ START HERE
   - Master script that runs everything
   - One command to generate all figures
   - Checks dependencies, extracts data, makes figures

2. **`extract_manuscript_data.py`**
   - Parses your weekend's work
   - Extracts ipTM scores from AF3 outputs
   - Reads domain annotations
   - Creates clean data file for figures

3. **`generate_all_figures.py`**
   - Creates all three main figures
   - Uses your actual data
   - Generates both PDF and PNG versions

4. **`generate_pymol_script.py`**
   - Creates PyMOL rendering script
   - Includes manual rendering guide
   - For Figure 1 Panel D structures

### **Documentation**

5. **`README_figures.md`**
   - Detailed usage instructions
   - Troubleshooting guide
   - Customization tips

6. **`Manuscript_Checklist.md`**
   - Complete publication checklist
   - 3-week timeline
   - All tasks organized

7. **`Introduction_Template.md`**
   - Structured introduction draft
   - Key references included
   - Ready to customize

---

## 🚀 Quick Start (3 Commands)

```bash
# 1. Go to your project directory
cd /storage/kiran-stuff/aaRS

# 2. Copy all scripts there
cp /path/to/downloaded/scripts/* .

# 3. Run the pipeline
bash run_figure_pipeline.sh
```

That's it! All figures will be in `manuscript_figures/` directory.

---

## 📊 What Gets Generated

```
manuscript_figures/
├── extracted_data.json              # Your parsed data
│
├── Figure1_phylogeny_domains.pdf    # Phylogeny + domains
├── Figure1_phylogeny_domains.png    # (preview)
│
├── Figure2_af3_results.pdf          # AF3 binding results  
├── Figure2_af3_results.png          # (preview)
│
├── Figure3_domain_evolution.pdf     # Domain evolution
├── Figure3_domain_evolution.png     # (preview)
│
├── render_structures.pml            # PyMOL script
├── PyMOL_manual_guide.md            # Manual rendering guide
│
└── structures/                      # For rendered PNGs
    ├── deep_domain_pro_render.png
    ├── deep_domain_thr_render.png
    └── ... (6 total)
```

---

## 🎨 Figure Contents

### Figure 1: Foundation (Phylogeny + Domains)
**What it shows:**
- Panel A: Your ProRS/ThrRS phylogenetic trees
- Panel B: Reconstruction quality (93% metric)
- Panel C: Domain architecture evolution timeline
- Panel D: Structure renders (you'll add from PyMOL)

**Data sources:**
- `phase1b/results/*_deep.treefile`
- `domain_analysis_complete/complete_summary.tbl`
- Your ancestral FASTA sequences

---

### Figure 2: Core Results (AF3 Binding)
**What it shows:**
- Panel A: ipTM heatmap (all enzyme-ligand pairs)
- Panel B: Promiscuity trajectory (3.5 Gya → present)
- Panel C: Discrimination comparison (ΔipTM bars)
- Panel D: Validation controls (Phe/Trp negatives)

**Key findings visualized:**
- LUCA ProRS: THR binds at 83% of PRO affinity
- LUCA ThrRS: PRO binds at 99% of THR affinity  
- Shallow fusion: 89% cross-reactivity
- Modern: 98% - promiscuity persisted!

**Data sources:**
- `phase2/outputs/*/summary_confidences.json`
- Automatically extracted from all 16 models

---

### Figure 3: Mechanism (Domain Evolution)
**What it shows:**
- Panel A: Domain presence matrix (which domains in which ancestor)
- Panel B: Editing domain binding (weak for both PRO and THR)
- Panel C: Domain confidence (Pfam E-values)
- Panel D: Evolutionary model schematic

**Key insight:**
- Editing domain present but non-discriminating
- Lost in fusion event
- Supports two-stage discrimination model

**Data sources:**
- Domain annotations from Pfam scan
- Editing domain AF3 models

---

## 🔧 Customization Guide

### Change Colors

Edit `generate_all_figures.py`, line ~20:
```python
COLORS = {
    'ProRS': '#2E86AB',      # Your color here
    'ThrRS': '#F77F00',      # Your color here
    'Fusion': '#A23B72',     # Your color here
}
```

### Adjust Figure Size

Edit individual figure functions:
```python
def generate_figure1():
    fig = plt.figure(figsize=(7, 9))  # width, height in inches
```

### Add Your ipTM Scores Manually

If auto-extraction fails, edit `extract_manuscript_data.py`:
```python
def create_summary_table():
    summary = {
        'catalytic_binding': {
            'LUCA_ProRS': {
                'PRO': 0.75,   # Your actual value
                'THR': 0.62,   # Your actual value
            }
        }
    }
```

### High-Resolution Export

For submission, change DPI in `generate_all_figures.py`:
```python
plt.savefig(output_file, dpi=600)  # Change from 300 to 600
```

---

## 🐛 Common Issues & Solutions

### Issue: "No module named matplotlib"
**Solution:**
```bash
pip install --user matplotlib seaborn pandas numpy biopython ete3
```

### Issue: "File not found: phase1b/results/..."
**Solution:** 
Edit `extract_manuscript_data.py` and update file paths to match your directory structure.

### Issue: "Empty figure panels"
**Solution:** 
Run `extract_manuscript_data.py` first. It creates `extracted_data.json` that the figure script needs.

### Issue: "Font 'Arial' not found"
**Solution:** 
Edit `generate_all_figures.py`, line ~18:
```python
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
```

---

## 📝 Step-by-Step Workflow

### Day 1: Generate Draft Figures
```bash
cd /storage/kiran-stuff/aaRS
bash run_figure_pipeline.sh
# Review PNG previews
```

### Day 2: Render Structures
```bash
# Option 1: Automated
pymol -c manuscript_figures/render_structures.pml

# Option 2: Manual (recommended)
# Open PyMOL GUI and follow manuscript_figures/PyMOL_manual_guide.md
```

### Day 3: Polish and Finalize
```bash
# 1. Update any values in extract_manuscript_data.py
# 2. Adjust colors/layout in generate_all_figures.py
# 3. Re-run:
bash run_figure_pipeline.sh

# 4. Generate high-res versions (600 DPI)
# Edit generate_all_figures.py: dpi=600
python3 generate_all_figures.py
```

---

## ✅ Pre-Submission Checklist

Before sending figures to journal:

- [ ] All panels have real data (no placeholders)
- [ ] Structures rendered and inserted in Figure 1D
- [ ] Colors are publication-quality (not garish)
- [ ] Text is readable (minimum 6pt font)
- [ ] Exported at 600 DPI for print quality
- [ ] File sizes reasonable (<10 MB per figure)
- [ ] Consistent style across all figures
- [ ] Figure legends written (see Introduction_Template.md)
- [ ] Referenced correctly in manuscript text

---

## 🎯 Integration with Writing

### As You Write the Results Section:

```markdown
As shown in Figure 1A, phylogenetic reconstruction of ProRS and 
ThrRS sequences spanning 3.5 billion years revealed... (cite Figure 1A)

Structural modeling demonstrated that ancestral catalytic domains 
bound non-cognate substrates at 83-99% of cognate affinity 
(Figure 2A; Table 1)...

Domain architecture analysis revealed that the editing domain, 
while present in LUCA ProRS (Figure 3A), exhibited weak binding 
affinity to both substrates (Figure 3B)...
```

### Create Supporting Tables:

Use the extracted data to make supplementary tables:
- Table S1: All ipTM scores (from `extracted_data.json`)
- Table S2: Domain annotations (from Pfam output)
- Table S3: Sequence statistics (from FASTA files)

---

## 💡 Pro Tips

1. **Version control**: Keep dated copies as you revise
   ```bash
   cp manuscript_figures manuscript_figures_v1_$(date +%Y%m%d)
   ```

2. **Backup originals**: Keep PDFs as master copies, edit in Illustrator if needed

3. **Consistent style**: Use same colors/fonts across all figures

4. **Accessibility**: Make sure figures work in grayscale (for printing)

5. **Test on different screens**: Check readability on phone/tablet

---

## 🎓 Learning Resources

### Understanding the Code

- **Matplotlib**: https://matplotlib.org/stable/gallery/
- **Seaborn**: https://seaborn.pydata.org/examples/
- **BioPython**: https://biopython.org/wiki/Phylo
- **PyMOL**: https://pymolwiki.org/

### Publication Figure Standards

- **Nature**: https://www.nature.com/nature/for-authors/final-submission
- **Science**: https://www.science.org/content/page/instructions-preparing-initial-manuscript
- **PLOS**: https://journals.plos.org/plosone/s/figures

---

## 🚀 Next Steps After Figures

1. **Write figure legends** (see Introduction_Template.md for examples)
2. **Create supplementary figures** (use same scripts, different data subsets)
3. **Generate supplementary tables** (export from extracted_data.json)
4. **Draft Results section** referencing your figures
5. **Assemble complete manuscript** using provided templates

---

## 📞 Need Help?

If something doesn't work:

1. Check the error message carefully
2. Verify file paths in `extract_manuscript_data.py`
3. Try running scripts individually to isolate the issue
4. Check that all dependencies are installed
5. Share the error with Claude for debugging

**Most common fix:** Update the `BASE_DIR` path in each script to match your actual directory structure.

---

## 🎉 You're Ready!

You have everything you need:
- ✅ Complete scripts
- ✅ Documentation
- ✅ Templates
- ✅ Checklist

**Just run the pipeline and you'll have publication-ready figures!**

```bash
cd /storage/kiran-stuff/aaRS
bash run_figure_pipeline.sh
```

Good luck with your manuscript! 🧬📄🚀
