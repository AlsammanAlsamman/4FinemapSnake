#!/usr/bin/env Rscript
# SuSiE fine-mapping: Sum of Single Effects with LD-aware credible sets
suppressPackageStartupMessages({
  library(susieR)
  library(data.table)
  library(jsonlite)
})

# Parse command line arguments
parse_args <- function(x) {
  out <- list()
  i <- 1
  while (i <= length(x)) {
    key <- gsub("^--", "", x[i])
    out[[key]] <- x[i + 1]
    i <- i + 2
  }
  out
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))

input_tsv <- opt[["input-tsv"]]
ld_matrix_file <- opt[["ld-matrix"]]
out_tsv <- opt[["out-tsv"]]
diag_json <- opt[["diag-json"]]
done_file <- opt[["done-file"]]
L <- as.integer(opt[["L"]])
coverage <- as.numeric(opt[["coverage"]])

# Create output directories
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(diag_json), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(done_file), recursive = TRUE, showWarnings = FALSE)

# Read matched GWAS data
df <- fread(input_tsv, sep = "\t", check.names = FALSE, data.table = FALSE)

if (!"SNP" %in% colnames(df)) {
  stop("Matched table must include SNP column")
}

# Compute z-scores from BETA and SE
if (all(c("BETA", "SE") %in% colnames(df))) {
  beta <- as.numeric(df$BETA)
  se <- as.numeric(df$SE)
  z_scores <- beta / se
} else {
  stop("Matched table must include BETA and SE columns for z-score computation")
}

# Read LD matrix
ld_matrix <- as.matrix(fread(ld_matrix_file, sep = "\t", header = TRUE, data.table = FALSE))
rownames(ld_matrix) <- colnames(ld_matrix)

# Ensure SNP order matches LD matrix
if (length(df$SNP) != nrow(ld_matrix)) {
  stop(sprintf("SNP count mismatch: matched=%d, LD matrix=%d", length(df$SNP), nrow(ld_matrix)))
}

# Fit SuSiE model
# susie() expects z-scores and LD matrix (correlation matrix)
# It will compute L causal effect estimates assuming exactly L causal variants
susie_fit <- susie(z = z_scores, LD = ld_matrix, L = L, verbose = FALSE, standardize = FALSE)

# Extract posterior inclusion probabilities (PIPs)
pips <- susie_fit$pip

# Create output table with top variants by PIP
result_df <- data.frame(
  SNP = df$SNP,
  P = df$P,
  BETA = df$BETA,
  SE = df$SE,
  PIP = pips,
  stringsAsFactors = FALSE
)

# Sort by PIP descending
result_df <- result_df[order(-result_df$PIP), ]

# Select top L variants (or fewer if not enough with non-zero PIP)
n_select <- min(L, sum(result_df$PIP > 0))
result_df_selected <- result_df[seq_len(n_select), ]

# Write credible set
write.table(result_df_selected, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Extract credible sets from susie_fit
cs_list <- susie_fit$sets$cs
cs_sizes <- if (is.null(cs_list)) 0 else sapply(cs_list, length)

# Calculate 95% credible set size
cs_95_n <- 0
cs_95_snps <- character(0)
if (!is.null(cs_list) && length(cs_list) > 0) {
  cumsum_pip <- 0
  for (i in 1:length(cs_list)) {
    cs_snps <- cs_list[[i]]
    cumsum_pip <- cumsum_pip + sum(pips[cs_snps])
    cs_95_snps <- c(cs_95_snps, df$SNP[cs_snps])
    if (cumsum_pip >= coverage) {
      cs_95_n <- length(cs_95_snps)
      break
    }
  }
}

# Get lead SNP (highest PIP)
lead_idx <- which.max(pips)
lead_snp <- df$SNP[lead_idx]
lead_pip <- pips[lead_idx]

# Create diagnostics
diag_list <- list(
  method = "susie_ld_aware",
  L = L,
  coverage = coverage,
  n_selected = nrow(result_df_selected),
  n_credible_set_95 = cs_95_n,
  lead_snp = lead_snp,
  lead_pip = round(lead_pip, 4),
  n_variants_total = nrow(df),
  n_causal_assumed = L
)

# Write diagnostics JSON
writeLines(toJSON(diag_list, pretty = TRUE), con = diag_json)

# Write done marker
writeLines("ok", con = done_file)

cat(sprintf("SuSiE complete: %d/%d variants selected, lead SNP: %s (PIP=%.4f)\n",
            nrow(result_df_selected), nrow(df), lead_snp, lead_pip))

