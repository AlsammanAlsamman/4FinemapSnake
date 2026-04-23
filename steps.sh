#!/usr/bin/env bash
# Pipeline execution guide — run steps one by one using submit.sh.
# Each command targets a specific rule via --until.
# Dry-run first to validate graph, then execute.

# ─── DRY-RUN (validate DAG before any real execution) ────────────────────────
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --dry-run --snakefile Snakefile"

# ─── STEP 1: Extract locus GWAS and standardize columns ──────────────────────
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile results/asian/01_extract/ARID3A/extracted.done results/disc/01_extract/ARID3A/extracted.done"

# ─── STEP 2: Harmonize GWAS to reference panel PLINK alleles ─────────────────
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile results/asian/02_harmonize/ARID3A/harmonized.done results/disc/02_harmonize/ARID3A/harmonized.done"

# ─── STEP 3: Filter low-MAF SNPs and enforce reference SNP order ─────────────
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile results/asian/03_match/ARID3A/matched.done results/disc/03_match/ARID3A/matched.done"

# ─── STEP 4: LD-score vs -log10(P) diagnostic plots ─────────────────────────
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile results/asian/04_ldscore_plot/ARID3A/ldscore_plot.done results/disc/04_ldscore_plot/ARID3A/ldscore_plot.done"

# ─── STEP 5: Build LD matrix ─────────────────────────────────────────────────
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile results/asian/05_ld_matrix/ARID3A/ld_matrix.done results/disc/05_ld_matrix/ARID3A/ld_matrix.done"

# ─── STEP 6: QC / symmetrize LD matrix ───────────────────────────────────────
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile results/asian/06_ld_qc/ARID3A/ld_qc.done results/disc/06_ld_qc/ARID3A/ld_qc.done"

# ─── STEP 7: FINEMAP ─────────────────────────────────────────────────────────
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile results/asian/07_finemap/ARID3A/finemap.done results/disc/07_finemap/ARID3A/finemap.done"

# ─── STEP 8: SuSiE ───────────────────────────────────────────────────────────
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile results/asian/08_susier/ARID3A/susier.done results/disc/08_susier/ARID3A/susier.done"

# ─── STEP 9: COJO iterative conditional analysis ─────────────────────────────
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile results/asian/09_cojo/ARID3A/cojo.done results/disc/09_cojo/ARID3A/cojo.done"

# ─── STEP 9b: COJO iterative conditional analysis (GCTA-COJO) ───────────────
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile results/asian/09_cojo_gcta/ARID3A/cojo.done results/disc/09_cojo_gcta/ARID3A/cojo.done"

# ─── STEP 10: Manhattan plots ────────────────────────────────────────────────
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile results/asian/10_manhattan/ARID3A/manhattan_all.png results/disc/10_manhattan/ARID3A/manhattan_all.png"

# ─── STEP 11: Export per-locus summary tables (TSV + Excel) ──────────────────
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile results/summary/ARID3A_summary.tsv"

# ─── FULL PIPELINE (all steps in one shot) ───────────────────────────────────
# bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile"
