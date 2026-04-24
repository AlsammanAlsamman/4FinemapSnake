#!/usr/bin/env bash
# Pipeline execution guide — run one step at a time with simple rule-based commands.
# Each command runs ALL configured datasets/loci from configs/analysis.yml.

# Dry-run (validate DAG before execution)
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile --dry-run"

# Step 1: Extract locus GWAS and reference panel slice
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile --allowed-rules extract_locus,extract_locus_all extract_locus_all"

# Step 2: Harmonize GWAS to reference panel alleles
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile --allowed-rules harmonize_with_refpanel,harmonize_with_refpanel_all harmonize_with_refpanel_all"

# Step 3: Filter and match SNPs to reference order
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile --allowed-rules filter_and_match_snps,filter_and_match_snps_all filter_and_match_snps_all"

# Step 4: Build LD matrix
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile --allowed-rules build_ld_matrix,build_ld_matrix_all build_ld_matrix_all"

# Step 5: LD-score vs signal diagnostics
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile --allowed-rules ldscore_diagnostics,ldscore_diagnostics_all ldscore_diagnostics_all"

# Step 6: QC / symmetrize LD matrix
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile --allowed-rules qc_fix_ld_matrix,qc_fix_ld_matrix_all qc_fix_ld_matrix_all"

# Step 7: FINEMAP
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile --allowed-rules run_finemap,run_finemap_all run_finemap_all"

# Step 8: SuSiE
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile --allowed-rules run_susier,run_susier_all run_susier_all"

# Step 9: COJO iterative conditional analysis (GCTA-COJO)
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile --allowed-rules run_cojo_iterative_gcta,run_cojo_iterative_gcta_all run_cojo_iterative_gcta_all"

# Step 10: Manhattan plots
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile --allowed-rules plot_manhattan,plot_manhattan_all plot_manhattan_all"

# Step 11: Export per-locus summary tables (TSV + Excel)
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile --allowed-rules export_summary,export_summary_all export_summary_all"

# Step 12: Export COJO iteration SNP table (TSV + Excel)
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile --allowed-rules export_cojo_iteration_table,export_cojo_iteration_table_all export_cojo_iteration_table_all"

# Step 13: GWAS locus-zoom + LD R² triangle plot (PDF + PNG)
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile --allowed-rules plot_locus_ld,plot_locus_ld_all plot_locus_ld_all"

# Step 14: Export FINEMAP all-SNP table (TSV + Excel)
bash --login -c "cd $(pwd) && ml slurm && bash submit.sh --snakefile Snakefile --allowed-rules export_finemap_table,export_finemap_table_all export_finemap_table_all"