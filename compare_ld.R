library(data.table)

target <- "disc"
locus <- "ARID3A"

matched_path <- sprintf("results/%s/03_match/%s/matched.tsv", target, locus)
ld_a_path <- sprintf("results/%s/05_ld_matrix/%s/ld_matrix.tsv", target, locus)
ld_b_path <- sprintf("results/%s/06_ld_qc/%s/ld_matrix_qc.tsv", target, locus)

# Read matched SNPs
matched <- fread(matched_path)
# In matched.tsv, look for P columns. Let's assume 'P' or similar. 
# In some pipelines it might be 'p' or 'P_value'. Check column names first.
# For now, try 'P'.
sig_snps <- matched[P < 5e-4, SNP]

compute_ld_stats <- function(ld_path, sig_snps, label) {
  if (!file.exists(ld_path)) {
    cat(sprintf("%s: File not found\n", label))
    return(NULL)
  }
  
  # Read LD matrix
  ld_mat <- fread(ld_path)
  ld_snps <- colnames(ld_mat)[-1]
  if (is.null(ld_snps)) ld_snps <- colnames(ld_mat)
  
  overlap_sig <- intersect(sig_snps, ld_snps)
  
  mat_num <- as.matrix(ld_mat[, -1, with=FALSE])
  rownames(mat_num) <- ld_mat[[1]]
  
  valid_sig <- intersect(overlap_sig, rownames(mat_num))
  
  if (length(valid_sig) > 0) {
    ld_scores <- rowSums(mat_num[valid_sig, , drop=FALSE]^2)
    ld_score_norm <- ld_scores / ncol(mat_num)
    
    cat(sprintf("Source: %s\n", label))
    cat(sprintf("  overlap_sig: %d\n", length(valid_sig)))
    cat(sprintf("  unique ld_score_norm: %d\n", length(unique(ld_score_norm))))
    cat(sprintf("  ld_score_norm min: %f, median: %f, max: %f\n", 
                min(ld_score_norm), median(ld_score_norm), max(ld_score_norm)))
  } else {
    cat(sprintf("Source: %s - No significant SNPs found in LD matrix\n", label))
  }
}

compute_ld_stats(ld_a_path, sig_snps, "05_ld_matrix")
compute_ld_stats(ld_b_path, sig_snps, "06_ld_qc")
