#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx == length(args)) {
    stop(sprintf("Missing required argument: %s", flag))
  }
  args[idx + 1]
}

input_tsv <- get_arg("--input-tsv")
ref_bim <- get_arg("--ref-bim")
out_matrix <- get_arg("--out-matrix")
diag_json <- get_arg("--diag-json")
done_file <- get_arg("--done-file")

mkdir_for <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

mkdir_for(out_matrix)
mkdir_for(diag_json)
mkdir_for(done_file)

# Read matched GWAS table.
df <- read.table(input_tsv, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (!("SNP" %in% colnames(df))) {
  stop("Input file must include column SNP")
}

# Read matched PLINK subset BIM (CHR SNP CM BP A1 A2).
bim <- read.table(ref_bim, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
if (ncol(bim) < 6) {
  stop("Reference BIM must have at least 6 columns")
}
colnames(bim)[1:6] <- c("CHR", "SNP", "CM", "BP", "A1", "A2")

ref_snps <- as.character(bim$SNP)
gwas_snps <- as.character(df$SNP)

if (length(ref_snps) == 0) {
  stop("Reference BIM has zero SNPs")
}

set_equal <- setequal(ref_snps, gwas_snps)
if (!set_equal) {
  only_ref <- head(setdiff(ref_snps, gwas_snps), 5)
  only_gwas <- head(setdiff(gwas_snps, ref_snps), 5)
  stop(sprintf(
    "SNP mismatch between matched GWAS TSV and refpanel_matched BIM. Examples only_in_bim=%s, only_in_gwas=%s",
    paste(only_ref, collapse = ","),
    paste(only_gwas, collapse = ",")
  ))
}

# Reorder GWAS rows to PLINK subset order for deterministic output.
row_idx <- match(ref_snps, gwas_snps)
if (any(is.na(row_idx))) {
  stop("Failed to align GWAS rows to BIM SNP order")
}
df <- df[row_idx, , drop = FALSE]

# Build LD proxy from distance if explicit genotypes are unavailable.
if ("BP" %in% colnames(df)) {
  bp <- suppressWarnings(as.numeric(df$BP))
} else {
  bp <- seq_len(nrow(df)) - 1
}

bad_bp <- !is.finite(bp)
if (any(bad_bp)) {
  bp[bad_bp] <- (seq_len(nrow(df)) - 1)[bad_bp]
}

if (length(bp) <= 1) {
  mat <- matrix(1.0, nrow = length(bp), ncol = length(bp))
  span <- 0.0
  tau <- 0.0
} else {
  span <- max(bp) - min(bp)
  tau <- max(span / 20.0, 5000.0)
  dist_mat <- abs(outer(bp, bp, `-`))
  mat <- exp(-dist_mat / tau)
  diag(mat) <- 1.0
}

snps <- as.character(df$SNP)
out <- data.frame(SNP = snps, mat, check.names = FALSE)
colnames(out)[-1] <- snps
write.table(out, file = out_matrix, sep = "\t", quote = FALSE, row.names = FALSE)

diag <- list(
  n_snps = as.integer(length(snps)),
  n_ref_bim_snps = as.integer(length(ref_snps)),
  snp_set_equal_to_ref_bim = isTRUE(set_equal),
  matrix_shape = c(as.integer(nrow(mat)), as.integer(ncol(mat))),
  ld_model = "exp_decay_by_bp",
  bp_span = as.numeric(span),
  tau = as.numeric(tau),
  max_abs_asymmetry = if (length(mat) == 0) 0.0 else max(abs(mat - t(mat)))
)

json_lines <- c(
  "{",
  sprintf("  \"n_snps\": %d,", diag$n_snps),
  sprintf("  \"n_ref_bim_snps\": %d,", diag$n_ref_bim_snps),
  sprintf("  \"snp_set_equal_to_ref_bim\": %s,", tolower(as.character(diag$snp_set_equal_to_ref_bim))),
  sprintf("  \"matrix_shape\": [%d, %d],", diag$matrix_shape[1], diag$matrix_shape[2]),
  sprintf("  \"ld_model\": \"%s\",", diag$ld_model),
  sprintf("  \"bp_span\": %.12g,", diag$bp_span),
  sprintf("  \"tau\": %.12g,", diag$tau),
  sprintf("  \"max_abs_asymmetry\": %.12g", diag$max_abs_asymmetry),
  "}"
)
writeLines(json_lines, diag_json)
writeLines("ok", done_file)
