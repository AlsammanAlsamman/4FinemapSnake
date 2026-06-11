# Quick Reference: Running Colocalization Pipeline

## 🚀 Start Here

### Option 1: Run All Loci (Recommended)
```bash
# Full genome-wide colocalization analysis
cd /s/nath-lab/alsamman/____MyCodes____/4FinemapSnake
bash submit.sh results_Asian_IL12AB/18_coloc/all_coloc.done --cores 16 --jobs 4
```

### Option 2: Test on Single Locus First
```bash
# Quick test before full run
bash submit.sh results_Asian_IL12AB/disc/18_coloc/IL12A/coloc.done --cores 4
```

### Option 3: Dry-Run (Preview What Will Execute)
```bash
# See what would run without executing
bash submit.sh --dry-run results_Asian_IL12AB/disc/18_coloc/IL12A/coloc.done
```

---

## 📊 Monitor Progress

### Check Job Status
```bash
# View SLURM jobs
squeue -u alsammana -j <jobid>

# Check completed files
find results_Asian_IL12AB/18_coloc -name "*.done" | wc -l
```

### View Results Streaming
```bash
# Monitor top hits as they complete
watch -n 5 'tail -20 results_Asian_IL12AB/disc/18_coloc/*/coloc_results.tsv'
```

---

## 📈 Analyze Results

### After Pipeline Completes

1. **Merge Genome-Wide Results**
```bash
bash submit.sh results_Asian_IL12AB/18_coloc/coloc_merge.done --cores 2
```

2. **View Top Colocalization Signals**
```bash
head -50 results_Asian_IL12AB/18_coloc/coloc_genome_wide_results.tsv | cut -f1-9 | column -t
```

3. **Summary Statistics**
```bash
# Open Excel workbook with summary by tissue
open results_Asian_IL12AB/18_coloc/coloc_summary_statistics.xlsx
```

4. **Extract High-Confidence Signals (PP.H4 ≥ 0.75)**
```bash
awk '$6 >= 0.75' results_Asian_IL12AB/18_coloc/coloc_genome_wide_results.tsv | wc -l
```

---

## 🔍 Troubleshooting

| Issue | Solution |
|-------|----------|
| "No rule to produce coloc.done" | Verify Snakefile includes rule 18: `grep include Snakefile` |
| Missing upstream files | Run finemap/SuSiE first: `bash submit.sh results_Asian_IL12AB/disc/10_summary/IL12A/summary.done` |
| Memory errors | Increase `mem_mb=128000` in `rules/18_run_coloc.smk` line 49 |
| Low PP.H4 values | Check GWAS signal strength; verify SNP selection criteria |
| Script errors | Check log: `tail -100 results_Asian_IL12AB/disc/18_coloc/*/coloc_run.log` |

---

## 📚 Output Structure

```
results_Asian_IL12AB/
├── 18_coloc/
│   ├── disc/
│   │   └── IL12A/
│   │       └── 18_coloc/
│   │           ├── coloc_results.tsv           ← All results
│   │           ├── coloc_summary_pp.tsv        ← Significant hits (PP.H4 ≥ 0.75)
│   │           ├── coloc_credible_set95.tsv    ← 95% credible sets
│   │           ├── coloc_diagnostics.json      ← Analysis metadata
│   │           └── coloc.done                  ← Completion marker
│   │
│   ├── all_coloc.done                          ← All loci complete marker
│   ├── coloc_merge.done                        ← Merge complete marker
│   │
│   ├── coloc_genome_wide_results.tsv           ← Merged results
│   ├── coloc_genome_wide_credible_set95.tsv    ← Merged credible sets
│   └── coloc_summary_statistics.xlsx           ← Excel workbook
```

---

## 🧬 Understanding Output Columns

### coloc_results.tsv
- `gene_id`: eGene identifier
- `gene_name`: eGene symbol (e.g., IL12A)
- `dataset`: eQTL source (e.g., Asian_4pop_B_CELL_NAIVE)
- `pp_h0-h4`: Posterior probabilities for each hypothesis
- `nsnps_shared`: Number of variants in common between GWAS & eQTL

### coloc_summary_pp.tsv
- Filtered to PP.H4 ≥ 0.75 (high-confidence colocalization)
- Same columns as results

### coloc_credible_set95.tsv
- `SNP.PP.H4`: Posterior probability specific to this variant
- `cum_pp`: Cumulative posterior (variants sorted until ≥0.95)
- Identifies specific causal variants in colocalized loci

---

## 💡 Interpretation Guide

### Colocalization Signal Strength
- **PP.H4 ≥ 0.75**: Strong evidence for shared causal variant
- **PP.H4 ≥ 0.50**: Moderate evidence (commonly used threshold)
- **PP.H4 0.10-0.50**: Weak evidence (exploratory signals)
- **PP.H4 < 0.10**: No evidence (independent signals likely)

### Key Findings Pattern
```
Gene        | Dataset           | PP.H4  | Interpretation
------------|-------------------|--------|------------------
IL12A       | Asian_4pop_NK     | 0.92   | ✓ Strong signal - NK cell eQTL colocalized with GWAS
IL12A       | Immunex_mDC       | 0.68   | • Moderate signal - Dendritic cell colocalization
IL12A       | eQTLGen           | 0.04   | ✗ No signal - Different tissue doesn't colocalize
IL12B       | Asian_4pop_TH1    | 0.81   | ✓ Strong signal - Helper T cell specific effect
```

---

## 🛠️ Custom Configurations

### Change PP.H4 Threshold
Edit `configs/analysis.yml`, line 157:
```yaml
coloc_params:
  pp_h4_threshold: 0.50      # Change from 0.75
```

### Add/Remove eQTL Datasets
Edit `configs/analysis.yml`, lines 21-130, add/remove datasets:
```yaml
resources:
  eqtl_datasets:
    Asian_4pop_B_CELL_NAIVE: ...
    # Add new dataset here or remove existing ones
```

### Adjust Resource Allocation
Edit `rules/18_run_coloc.smk`, lines 48-51:
```r
resources:
  mem_mb=64000,      # Increase if memory errors
  time="02:00:00",   # Increase if timeout
  cores=4            # Increase for parallelization
```

---

## 📞 Support

For issues, check:
1. **logs/**: Snakemake workflow logs
2. **coloc_run.log**: Per-locus R script logs
3. **coloc_diagnostics.json**: Analysis statistics

Last Updated: May 19, 2026
Status: ✅ Production Ready
