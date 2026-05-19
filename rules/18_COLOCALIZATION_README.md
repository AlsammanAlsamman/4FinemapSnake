# Colocalization Analysis Pipeline (Rule 18)

## Overview

This colocalization pipeline performs **Bayesian statistical colocalization** between GWAS signals and eQTL signals using the `coloc` package (ABF method). The analysis is per-locus, using isolated regions from upstream finemapping/fine-structure analysis (finemap, SuSiE, COJO).

### Key Features

1. **Region-Isolated Analysis**: Uses only SNPs selected by finemap/SuSiE/COJO, ensuring focus on likely causal variants
2. **Multi-Dataset Support**: Analyzes colocalization with 43 eQTL datasets (4pop, eQTLGen, Immunex)
3. **Bayesian Methodology**: Implements ABF (Approximate Bayes Factor) model:
   - H0: No causal variants in either trait
   - H1: Causal variant in GWAS only
   - H2: Causal variant in eQTL only
   - H3: Separate causal variants in both traits
   - H4: **Shared causal variant** (colocalization)
4. **Credible Set Integration**: Maps posterior probabilities to 95% credible sets
5. **Robust Error Handling**: Manages missing SE, AF, MAF data through estimation

---

## Input Requirements

### From Upstream Rules
- **Rule 10** (`export_summary.tsv`): Summary table with finemap/SuSiE/COJO results
  - Contains: SNP, CHR, BP, P, BETA, SE, AF, finemap_pip, finemap_cs, susie_pip, susie_cs, cojo
- **Rule 14** (`selected_snps.tsv`): SNPs selected by any finemapping method
- **Rule 14** (`eqtl_subsets/all_eqtls.tsv`): Combined eQTL data for the locus

### Configuration (analysis.yml)
```yaml
resources:
  eqtl_datasets:
    <DATASET_NAME>:
      path: "resources/eQTLs/<FILE>.tsv"
      chr_col: 1
      pos_col: 2
      rsid_col: 3
      gene_id_col: 4
      gene_name_col: 5
      beta_col: 6
      se_col: null  # Can be null; estimated from p-value
      pval_col: 7
      af_col: null  # Optional
      fdr_col: null # Optional
      type: "4pop"  # or "immunex" or "meta"

coloc_params:
  cojo_keep_label: "iter1"
  pp_h4_threshold: 0.75           # Posterior probability cutoff
  pp_h4_credible_set: 0.95        # Cumulative PP for credible set
  analysis_type: "quant"          # "quant" or "cc"
  output_all_hits: true           # Output non-significant hits too
  include_credible_set: true      # Annotate credible sets
```

---

## Available eQTL Datasets

### 4-Population Immune Cell Types (13 datasets)
Asian GWAS meta-analysis eQTLs from 4-population immune cells:
- `B_CELL_NAIVE`, `CD4_NAIVE`, `CD8_NAIVE`
- `M2`, `MONOCYTES`, `NK`
- `TFH`, `TH1`, `TH17`, `TH2`, `THSTAR`
- `TREG_MEM`, `TREG_NAIVE`

**Characteristics**:
- **Type**: Population-specific (Asian)
- **Format**: CHR, BP, SNP, GeneID, GeneName, Beta, P-value
- **SE**: Estimated from beta and p-value
- **Sample**: ~100-500 individuals per cell type

### eQTLGen Whole-Blood (1 dataset)
European meta-analysis of whole-blood eQTLs:
- `eQTLGen`

**Characteristics**:
- **Type**: Whole-blood meta-analysis
- **Format**: CHR, BP, SNP, GeneID, GeneName, Beta, SE, P-value, FDR
- **SE**: Provided
- **Samples**: ~31,500 individuals

### Immunex Immune Cell Types (28 datasets)
European Immunological Genome (ImmGen) project cell types:
- `immunex_CD16p_Mono`, `immunex_CL_Mono`, `immunex_CM_CD8`, `immunex_DN_B`
- `immunex_EM_CD8`, `immunex_Fr_I_nTreg`, `immunex_Fr_II_eTreg`, `immunex_Fr_III_T`
- `immunex_Int_Mono`, `immunex_LDG`, `immunex_mDC`, `immunex_Mem_CD4`
- `immunex_Mem_CD8`, `immunex_Naive_B`, `immunex_Naive_CD4`, `immunex_Naive_CD8`
- `immunex_NC_Mono`, `immunex_Neu`, `immunex_NK`, `immunex_pDC`
- `immunex_Plasmablast`, `immunex_SM_B`, `immunex_TEMRA_CD8`, `immunex_Tfh`
- `immunex_Th1`, `immunex_Th17`, `immunex_Th2`, `immunex_USM_B`

**Characteristics**:
- **Type**: Fine-grained immune cell types (European)
- **Format**: CHR, BP, SNP, GeneID, GeneName, Beta, P-value
- **SE**: Estimated from beta and p-value
- **Samples**: ~400-2000 individuals per cell type

---

## Output Files

### Per-Locus Outputs (e.g., `results_Asian_IL12AB/disc/18_coloc/IL12A_B/`)

1. **`coloc_results.tsv`** - Full results for all gene-SNP pairs analyzed
   - Columns: gene_id, gene_name, dataset, n_variants, pp_h0, pp_h1, pp_h2, pp_h3, pp_h4, etc.
   - Sorted by PP.H4 (descending)
   - All pairs included if `output_all_hits: true`

2. **`coloc_summary_pp.tsv`** - Significant hits only (PP.H4 >= threshold)
   - Subset of full results meeting posterior probability threshold
   - Primary table for biological interpretation

3. **`coloc_credible_set95.tsv`** - 95% credible set variants
   - Individual variants contributing to credible set
   - Columns: gene_id, gene_name, dataset, SNP, SNP.PP.H4, cum_pp, pp_h4
   - Sorted by cumulative posterior probability

4. **`coloc_diagnostics.json`** - Analysis diagnostics
   - n_genes_analyzed, n_genes_success, n_genes_failed
   - Timestamp, R version, coloc package version

5. **`coloc_run.log`** - Detailed execution log

### Genome-Wide Outputs

6. **`coloc_genome_wide_results.tsv`** - Merged results across all loci
   - Combined from all per-locus `coloc_results.tsv` files

7. **`coloc_genome_wide_credible_set95.tsv`** - Merged credible sets

8. **`coloc_summary_statistics.xlsx`** - Excel workbook with:
   - Genome-Wide Results sheet (top 10,000 hits)
   - Summary Stats sheet (statistics by dataset)

---

## Statistical Model

### Coloc.ABF Method

The colocalization analysis uses the Approximate Bayes Factor (ABF) approach implemented in the `coloc` package.

#### Key Equations

**For each shared variant and each hypothesis:**

$$ABF_h = \frac{p(Z_{GWAS} | H_h)}{p(Z_{GWAS} | H_0)}$$

where $Z$ represents the z-score at each SNP.

**Prior probabilities:**
- Prior causal GWAS: $p_1 = 10^{-4}$ (1 in 10,000 SNPs)
- Prior causal eQTL: $p_2 = 10^{-4}$
- Prior shared causal: $p_{12} = 10^{-5}$ (1 in 100,000)

**Posterior probabilities via Bayes' rule:**

$$PP_h = \frac{ABF_h \times P_h}{\sum_{i=0}^{4} ABF_i \times P_i}$$

### Interpretation

- **PP.H4 > 0.75**: Strong evidence of colocalization
- **PP.H4 > 0.50**: Moderate evidence (often used in sensitivity analyses)
- **PP.H4 < 0.10**: Weak evidence (non-colocalization)
- **PP.H3 >> PP.H4**: Suggests distinct causal variants in each trait

### Credible Sets

The 95% credible set includes variants whose cumulative posterior probability reaches 0.95:

$$CS_{95\%} = \{SNPs: \sum_{i \in CS} PP_i \leq 0.95\}$$

---

## Key Bioinformatics Details

### 1. Standard Error Handling

**When SE is provided**: Use directly
**When SE is missing**: Estimate from beta and p-value using:

$$Z = \left|\frac{\beta}{SE}\right| = \left|\text{qnorm}(P/2)\right|$$
$$SE = \frac{|\beta|}{Z}$$

This approach is mathematically sound for two-tailed tests and used in standard GWAS analyses.

### 2. MHC Region Exclusion

SNPs in the MHC region (chr6:28477897-33448354, GRCh37) are excluded due to:
- Extremely high LD
- Complex haplotype structure
- Multiple independent signals
- Violation of coloc assumptions (single causal variant)

### 3. Minor Allele Frequency Filtering

Currently not applied (set to `null`), but can be configured to filter eQTLs:
- Removes rare allele variants
- Improves signal-to-noise ratio
- Recommend: MAF >= 1% for 500+ sample sizes

### 4. SNP Matching

Only SNPs present in both:
- GWAS summary statistics (from summary.tsv)
- eQTL dataset (from all_eqtls.tsv)

are included in the analysis.

### 5. Duplicate Handling

If a SNP appears multiple times in a dataset:
- First occurrence is kept
- Necessary due to multi-allelic site representation

---

## Running the Pipeline

### Basic Usage

```bash
# Run colocalization for all loci in a target
snakemake \
  --snakefile Snakefile \
  --config analysis_config=configs/analysis.yml \
  --allowed-rules run_coloc,run_coloc_all \
  --cores 8 \
  disc  # target name
```

### Run Specific Locus

```bash
snakemake \
  --snakefile Snakefile \
  --config analysis_config=configs/analysis.yml \
  results_Asian_IL12AB/disc/18_coloc/IL12A_B/coloc.done \
  --cores 8
```

### Generate Genome-Wide Merged Results

```bash
snakemake \
  --snakefile Snakefile \
  --config analysis_config=configs/analysis.yml \
  results_Asian_IL12AB/18_coloc/coloc_merge.done \
  --cores 4
```

---

## Troubleshooting

### Issue: "No variants shared between GWAS and eQTL"

**Cause**: SNP name mismatch (rsID vs coordinate format)
**Solution**: Ensure both datasets use same SNP identifier scheme in preprocessing

### Issue: "coloc.abf failed: Invalid input"

**Cause**: Missing or infinite SE values; extreme p-values
**Solution**: Check SE estimation; verify p-value range (>= 1e-300)

### Issue: Low PP.H4 values across all genes

**Cause**: Weak GWAS signal or weak eQTL signal at locus
**Solution**: Verify SNP selection criteria; check if locus has true finemapping signal

### Issue: Memory errors for large eQTL datasets

**Cause**: Attempting to analyze too many genes simultaneously
**Solution**: Increase `mem_mb` in rule resources; or pre-filter eQTL dataset to relevant genes

---

## Configuration Best Practices

1. **PP.H4 Threshold**: 
   - 0.75: Conservative (high-confidence colocalizations only)
   - 0.50: Standard (often used in literature)
   - 0.10-0.25: Exploratory (report all weak signals)

2. **Credible Set Coverage**: 
   - 0.95: Standard (95% credible set)
   - 0.99: More restrictive
   - 0.50: Very restrictive, only top variants

3. **Sample Size in GWAS**: 
   - Always populate `sample_size_gwas` from meta-analysis N
   - Critical for ABF calculations
   - If missing, coloc uses default calculations

4. **Multiple Testing**: 
   - With 43 datasets × ~10,000-30,000 genes per locus
   - Apply multiple testing correction if using liberal PP.H4 threshold
   - Consider FDR or Bonferroni

---

## Output Interpretation Example

```
gene_id    gene_name  dataset       pp_h0   pp_h1   pp_h2   pp_h3   pp_h4
ENSG123    IL12A      4pop_CD4      0.001   0.005   0.002   0.012   0.980  ← Strong colocalization
ENSG123    IL12A      eQTLGen       0.100   0.450   0.100   0.300   0.050  ← Weak signal
ENSG124    IL6        4pop_NK       0.500   0.250   0.200   0.040   0.010  ← No colocalization
```

### Interpretation:
- **ENSG123 with 4pop_CD4**: Strong evidence of shared causal variant (98% posterior)
- **ENSG123 with eQTLGen**: Ambiguous; possibly different causal variants
- **ENSG124**: No colocalization; distinct signals or no signal

---

## References

1. **Coloc Package**: Giambartolomei et al. (2014). Nature Genetics 46: 881-885
   - "Bayesian test for colocalisation between pairs of genetic associations"
   - https://github.com/chr1swallace/coloc

2. **Bayesian ABF Method**: Wakefield J (2009). Genet Epidemiol 33:79-86
   - "Bayes factors for genome-wide association studies"

3. **Credible Sets**: Maller et al. (2012). Genome Biol 13:R39
   - "Bayesian refinement of association signals for 14 loci in 3 common diseases"

4. **eQTL Colocalization**: Ongen et al. (2017). PLOS Computational Biology 13(2): e1005357
   - "Fast and efficient QTL mapper for thousands of molecular phenotypes"
