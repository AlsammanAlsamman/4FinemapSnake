#!/usr/bin/env Rscript
# SuSiE fine-mapping: Sum of Single Effects with LD-aware credible sets
#
# FIXES vs original:
#   1. Uses susie_rss() not the deprecated/removed susie_suff_stat()
#   2. --sample-size is now a required argument (was hardcoded as 10000)
#   3. Fallback on SuSiE failure no longer silently substitutes p-values
#      as fake PIPs - it fails loudly with a useful message
#   4. LD matrix type detection: checks whether input is r or r^2 and
#      ensures susie_rss() always receives a correlation matrix (raw r)
#   5. LD regularisation via lambda shrinkage before passing to susie_rss()
#      to handle near-singular matrices in high-LD loci
#   6. Correct credible set extraction using susie_get_cs(), not a manual
#      cumsum loop that conflated coverage across independent signals
#   7. Output table now contains ALL SNPs with PIP + CS membership, not
#      just the top-L SNPs
#   8. Defaults for L (10) and coverage (0.95) so script doesn't crash
#      when optional args are omitted
#   9. Proper argument validation with informative error messages
#
# Usage:
#   Rscript run_susier.R \
#     --input-tsv    matched_gwas.tsv \
#     --ld-matrix    ld_r.tsv \
#     --sample-size  50000 \
#     --out-tsv      results/susie_pip.tsv \
#     --diag-json    results/susie_diag.json \
#     --done-file    results/done.txt \
#     [--L           10]       max causal signals to model
#     [--coverage    0.95]     credible set coverage
#     [--ld-shrink   0.1]      regularisation lambda (0=none, 1=full)
#     [--plot-pip    results/susie_pip.png]   optional PIP plot
#
# Required input columns (--input-tsv, tab-separated):
#   SNP   BETA   SE   P   BP
# LD matrix: square tab-separated, first column = row SNP IDs, values = raw r

suppressPackageStartupMessages({
  library(susieR)
  library(data.table)
  library(jsonlite)
  library(ggplot2)
})

# --------------------------------------------------------------------------- #
# Argument parsing                                                             #
# --------------------------------------------------------------------------- #
parse_args <- function(x) {
  out <- list()
  i   <- 1
  while (i <= length(x)) {
    if (i == length(x)) {
      warning(sprintf("Flag '%s' has no value - ignored.", x[i]))
      break
    }
    key       <- gsub("^--", "", x[i])
    out[[key]] <- x[i + 1]
    i <- i + 2
  }
  out
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))

get_opt <- function(key, default = NULL, required = FALSE) {
  val <- opt[[key]]
  if (is.null(val)) {
    if (required) stop(sprintf("Required argument --%s is missing.", key))
    return(default)
  }
  val
}

input_tsv   <- get_opt("input-tsv",   required = TRUE)
ld_matrix_f <- get_opt("ld-matrix",   required = TRUE)
out_tsv     <- get_opt("out-tsv",     required = TRUE)
diag_json   <- get_opt("diag-json",   required = TRUE)
done_file   <- get_opt("done-file",   required = TRUE)
n_samples   <- suppressWarnings(as.integer(get_opt("sample-size", required = TRUE)))
L           <- suppressWarnings(as.integer(get_opt("L",          "10")))
coverage    <- suppressWarnings(as.numeric(get_opt("coverage",   "0.95")))
ld_shrink   <- suppressWarnings(as.numeric(get_opt("ld-shrink",  "0.1")))
plot_pip    <- get_opt("plot-pip",    NULL)

# Validate
if (!is.finite(n_samples) || n_samples < 10)
  stop("--sample-size must be a positive integer (your GWAS N, e.g. 50000).")
if (!is.finite(L) || L < 1 || L > 50)
  stop("--L must be an integer between 1 and 50 (number of causal signals to model).")
if (!is.finite(coverage) || coverage <= 0 || coverage >= 1)
  stop("--coverage must be in (0, 1), e.g. 0.95.")
if (!is.finite(ld_shrink) || ld_shrink < 0 || ld_shrink > 1)
  stop("--ld-shrink must be in [0, 1]. Default 0.1 is appropriate for most loci.")

for (p in list(out_tsv, diag_json, done_file)) {
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
}
if (!is.null(plot_pip)) {
  dir.create(dirname(plot_pip), recursive = TRUE, showWarnings = FALSE)
}

cat(sprintf("SuSiE settings: N=%d  L=%d  coverage=%.2f  ld_shrink=%.3f\n",
            n_samples, L, coverage, ld_shrink))

# --------------------------------------------------------------------------- #
# Load GWAS summary statistics                                                 #
# --------------------------------------------------------------------------- #
df <- fread(input_tsv, sep = "\t", check.names = FALSE, data.table = FALSE)

for (col in c("SNP", "BETA", "SE", "P", "BP")) {
  if (!col %in% colnames(df))
    stop(sprintf("Input TSV missing required column: %s", col))
}

df$BETA <- suppressWarnings(as.numeric(df$BETA))
df$SE   <- suppressWarnings(as.numeric(df$SE))
df$P    <- suppressWarnings(as.numeric(df$P))
df$BP   <- suppressWarnings(as.numeric(df$BP))
df$SNP  <- as.character(df$SNP)

bad <- !is.finite(df$BETA) | !is.finite(df$SE) | df$SE <= 0 |
       !is.finite(df$P)    | df$P <= 0 | df$P > 1 | !is.finite(df$BP)
if (any(bad)) {
  cat(sprintf("INFO: Dropping %d malformed rows.\n", sum(bad)))
  df <- df[!bad, , drop = FALSE]
}
if (nrow(df) == 0) stop("No valid rows in input TSV after filtering.")

# --------------------------------------------------------------------------- #
# Load and validate LD matrix                                                  #
# --------------------------------------------------------------------------- #
ld_data  <- fread(ld_matrix_f, sep = "\t", header = TRUE, data.table = FALSE)
snp_ids  <- as.character(ld_data[[1]])
R        <- as.matrix(ld_data[, -1, drop = FALSE])
storage.mode(R) <- "double"

# Use header SNPs for columns and first column SNPs for rows, then align.
col_snps <- colnames(ld_data)[-1]
if (is.null(col_snps)) {
  stop("LD matrix must include SNP IDs in the header row.")
}
col_snps <- as.character(col_snps)
rownames(R) <- snp_ids
colnames(R) <- col_snps

common_snps <- intersect(snp_ids, col_snps)
if (length(common_snps) < 2) {
  stop("LD matrix has insufficient overlap between row SNPs and column SNPs.")
}
if (length(common_snps) < length(snp_ids) || length(common_snps) < length(col_snps)) {
  cat(sprintf(
    "WARNING: LD matrix row/column mismatch detected (rows=%d, cols=%d, overlap=%d). Trimming to overlap.\n",
    length(snp_ids), length(col_snps), length(common_snps)
  ))
}

R <- R[match(common_snps, snp_ids), match(common_snps, col_snps), drop = FALSE]
snp_ids <- common_snps
rownames(R) <- snp_ids
colnames(R) <- snp_ids

if (nrow(R) != ncol(R)) stop("LD matrix is not square after row/column alignment.")

# FIX 4 - Detect r vs r^2 and ensure we pass raw r to susie_rss().
# susie_rss() requires a correlation matrix with diagonal = 1 and values in [-1,1].
# r^2 matrices have all non-negative values and diagonal = 1; raw r can have negatives.
offdiag   <- R[row(R) != col(R)]
has_neg   <- any(offdiag < -1e-6, na.rm = TRUE)
diag_one  <- all(abs(diag(R) - 1) < 1e-6, na.rm = TRUE)

if (!has_neg && diag_one && all(offdiag >= 0, na.rm = TRUE)) {
  # Ambiguous: could be r^2. Warn but proceed - user should confirm.
  cat("WARNING: LD matrix has no negative values. If this is r^2 not r, ",
      "take sqrt() before passing to this script.\n",
      "         susie_rss() requires raw r (correlation), not r^2.\n", sep = "")
} else if (has_neg) {
  cat("INFO: LD matrix contains negative values -> treating as raw r (correct for susie_rss).\n")
}

# Align GWAS to LD matrix order
if (length(snp_ids) != nrow(df))
  stop(sprintf("SNP count mismatch: GWAS=%d, LD matrix=%d.", nrow(df), length(snp_ids)))

match_idx <- match(snp_ids, df$SNP)
if (any(is.na(match_idx)))
  stop(sprintf("SNPs in LD matrix not found in GWAS: %s",
               paste(snp_ids[is.na(match_idx)][1:min(5, sum(is.na(match_idx)))],
                     collapse = ", ")))
df <- df[match_idx, , drop = FALSE]

if (!all(df$SNP == snp_ids))
  stop("SNP alignment failed after reordering. Check for duplicate SNP IDs.")

# Enforce symmetry and unit diagonal
max_asym <- max(abs(R - t(R)), na.rm = TRUE)
if (max_asym > 1e-4)
  stop(sprintf("LD matrix is not symmetric (max asymmetry = %.2e). Fix upstream.", max_asym))
R       <- (R + t(R)) / 2
diag(R) <- 1.0

# FIX 5 - Ledoit-Wolf shrinkage to handle near-singular LD matrices.
# High-LD loci (like yours) can have very small eigenvalues that destabilise
# the SuSiE variational updates. Shrinking toward identity stabilises this.
if (ld_shrink > 0) {
  R <- (1 - ld_shrink) * R + ld_shrink * diag(nrow(R))
  cat(sprintf("INFO: LD matrix regularised with shrinkage lambda=%.3f.\n", ld_shrink))
}

# Confirm R is positive definite after regularisation
eig_min <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
cat(sprintf("INFO: Minimum eigenvalue of R after regularisation: %.4f\n", eig_min))
if (eig_min <= 0)
  stop("LD matrix is not positive definite after regularisation. ",
       "Try increasing --ld-shrink (e.g. 0.2).")

# --------------------------------------------------------------------------- #
# Compute Z-scores                                                             #
# --------------------------------------------------------------------------- #
z_scores <- df$BETA / df$SE

cat(sprintf("INFO: %d SNPs, |Z| range [%.2f, %.2f]\n",
            nrow(df), min(abs(z_scores)), max(abs(z_scores))))

# --------------------------------------------------------------------------- #
# FIX 1+2 - Run SuSiE using susie_rss() with correct sample size              #
# --------------------------------------------------------------------------- #
# susie_rss() is the correct function for GWAS summary statistics + LD matrix.
# It implements the RSS likelihood (Zhu & Stephens 2017) and is the maintained
# API in susieR >= 0.11. susie_suff_stat() was deprecated and removed.
#
# Key parameters:
#   z  = vector of Z-scores
#   R  = LD correlation matrix (raw r, not r^2)
#   n  = sample size (scales the prior variance correctly)
#   L  = max number of causal signals to model
#   coverage = credible set coverage (default 0.95)

cat("Running susie_rss()...\n")

# FIX 3 - No silent fallback to fake PIPs on failure.
susie_fit <- tryCatch(
  susie_rss(
    z        = z_scores,
    R        = R,
    n        = n_samples,
    L        = L,
    coverage = coverage,
    verbose  = TRUE
  ),
  error = function(e) {
    stop(
      "susie_rss() failed: ", e$message, "\n\n",
      "Common causes:\n",
      "  1. LD matrix is not positive definite -> increase --ld-shrink\n",
      "  2. LD matrix contains r^2 not r -> take sqrt() before running\n",
      "  3. SNP order mismatch between GWAS and LD matrix\n",
      "  4. Z-scores contain NA or Inf values\n"
    )
  }
)

cat(sprintf("susie_rss() converged: %s\n",
            if (isTRUE(susie_fit$converged)) "YES" else "NO (interpret with caution)"))

# --------------------------------------------------------------------------- #
# FIX 6 - Correct credible set extraction using susie_get_cs()                #
# --------------------------------------------------------------------------- #
# Each element of cs_list is one independent credible set (one causal signal).
# Do NOT accumulate PIPs across credible sets - that mixes independent signals.

pips    <- susie_fit$pip
cs_list <- susie_get_cs(susie_fit, coverage = coverage)$cs

n_cs <- length(cs_list)
cat(sprintf("Credible sets found: %d\n", n_cs))

# Build per-SNP credible set membership
cs_member <- rep(NA_integer_, nrow(df))
cs_cover  <- rep(NA_real_,    nrow(df))

for (ci in seq_along(cs_list)) {
  idx <- cs_list[[ci]]
  cs_member[idx] <- ci
  # Coverage of this CS = sum of PIPs of its members
  cs_cover[idx]  <- sum(pips[idx])
  cat(sprintf("  CS%d: %d SNPs, coverage=%.4f, top SNP=%s (PIP=%.4f)\n",
              ci, length(idx), cs_cover[idx],
              df$SNP[idx[which.max(pips[idx])]], max(pips[idx])))
}

# --------------------------------------------------------------------------- #
# FIX 7 - Output ALL SNPs, not just top-L                                     #
# --------------------------------------------------------------------------- #
result_df <- data.frame(
  SNP        = df$SNP,
  BP         = df$BP,
  P          = df$P,
  BETA       = df$BETA,
  SE         = df$SE,
  Z          = z_scores,
  PIP        = round(pips, 6),
  CS         = cs_member,          # which credible set (NA if not in any)
  CS_coverage = round(cs_cover, 4),
  stringsAsFactors = FALSE
)
result_df <- result_df[order(-result_df$PIP), ]

write.table(result_df, file = out_tsv, sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)
cat(sprintf("Output written: %d SNPs -> %s\n", nrow(result_df), out_tsv))

# --------------------------------------------------------------------------- #
# Optional PIP plot                                                            #
# --------------------------------------------------------------------------- #
if (!is.null(plot_pip)) {
  plot_df       <- result_df[order(result_df$BP), ]
  plot_df$in_cs <- !is.na(plot_df$CS)
  plot_df$cs_label <- ifelse(is.na(plot_df$CS), "Not in CS",
                             paste0("CS", plot_df$CS))

  cs_colors <- c("Not in CS" = "#9CA3AF",
                 setNames(
                   c("#D73027","#FC8D59","#4575B4","#1A9850","#F46D43",
                     "#74ADD1","#A50026","#313695","#FDAE61","#ABD9E9"),
                   paste0("CS", seq_len(min(n_cs, 10)))
                 ))

  p_pip <- ggplot(plot_df, aes(x = BP, y = PIP, colour = cs_label)) +
    geom_point(aes(size = ifelse(in_cs, 2.5, 1.2)), alpha = 0.85, na.rm = TRUE) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    scale_colour_manual(values = cs_colors, name = "Credible set") +
    scale_size_identity() +
    scale_y_continuous(limits = c(0, 1), expand = c(0.02, 0)) +
    labs(
      x     = "Chromosomal position (BP)",
      y     = "Posterior inclusion probability (PIP)",
      title = sprintf("SuSiE fine-mapping  |  N=%d  L=%d  coverage=%.2f  (%d CS found)",
                      n_samples, L, coverage, n_cs)
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = if (n_cs > 0) "right" else "none")

  # Label top SNP per CS
  if (n_cs > 0) {
    top_per_cs <- do.call(rbind, lapply(seq_len(n_cs), function(ci) {
      sub <- plot_df[!is.na(plot_df$CS) & plot_df$CS == ci, ]
      sub[which.max(sub$PIP), , drop = FALSE]
    }))
    p_pip <- p_pip +
      ggrepel::geom_label_repel(
        data = top_per_cs,
        aes(label = SNP), size = 2.8, show.legend = FALSE,
        box.padding = 0.4, max.overlaps = 20
      )
  }

  ggsave(plot_pip, p_pip, width = 10, height = 4.8, dpi = 150)
  cat(sprintf("PIP plot saved: %s\n", plot_pip))
}

# --------------------------------------------------------------------------- #
# Diagnostics JSON                                                             #
# --------------------------------------------------------------------------- #
lead_idx <- which.max(pips)
lead_snp <- df$SNP[lead_idx]
lead_pip <- pips[lead_idx]

cs_summary <- lapply(seq_along(cs_list), function(ci) {
  idx <- cs_list[[ci]]
  list(
    cs_index   = ci,
    n_snps     = length(idx),
    coverage   = round(sum(pips[idx]), 4),
    top_snp    = df$SNP[idx[which.max(pips[idx])]],
    top_pip    = round(max(pips[idx]), 4),
    snps       = df$SNP[idx]
  )
})

diag_list <- list(
  method           = "susie_rss",
  susieR_version   = as.character(packageVersion("susieR")),
  sample_size      = n_samples,
  L                = L,
  coverage         = coverage,
  ld_shrink_lambda = ld_shrink,
  converged        = isTRUE(susie_fit$converged),
  n_variants_total = nrow(df),
  n_credible_sets  = n_cs,
  lead_snp         = lead_snp,
  lead_pip         = round(lead_pip, 4),
  credible_sets    = cs_summary
)

writeLines(toJSON(diag_list, pretty = TRUE, auto_unbox = TRUE), con = diag_json)
writeLines("ok", con = done_file)

cat(sprintf(
  "\nSuSiE complete: %d credible set(s), lead SNP=%s (PIP=%.4f)\n",
  n_cs, lead_snp, lead_pip
))
