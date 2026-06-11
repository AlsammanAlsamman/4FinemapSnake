#!/usr/bin/env Rscript
# SuSiE fine-mapping: Sum of Single Effects with LD-aware credible sets
#
# WHAT WAS WRONG IN THE ORIGINAL:
#
#   BUG 1 (critical): susie_suff_stat() was removed from susieR >= 0.11.
#     The call failed silently and the catch block replaced EVERY PIP with
#     pnorm(abs(z), lower.tail=FALSE)*2 — a two-sided p-value, NOT a posterior.
#     p-values for strong SNPs are near 0, so the strongest SNPs got the
#     LOWEST fake PIPs and were excluded entirely. The index SNP disappeared.
#     The 10 CSs with PIP=1 in your output were all fabricated from p-values.
#     Fix: use susie_rss(), which is the correct maintained function.
#
#   BUG 2 (critical): n = 10000 hardcoded.
#     Sample size scales the prior variance. Wrong N = miscalibrated PIPs.
#     Fix: --sample-size is now a required argument.
#
#   BUG 3: No LD regularisation.
#     Near-singular R (common in high-LD loci) causes SuSiE to fragment one
#     real signal into many spurious credible sets — the 10-CS explosion you saw.
#     Fix: Ledoit-Wolf shrinkage toward identity before passing R to susie_rss().
#
#   BUG 4: Credible set extraction via manual cumsum across CS list.
#     Each CS is already one independent signal — summing PIPs across them is
#     meaningless. Fix: use susie_get_cs() directly.
#
#   BUG 5: Output only kept top-L SNPs, discarding the full PIP distribution.
#     Fix: output all SNPs with PIP + CS membership.
#
#   BUG 6: No defaults for L and coverage — crashed on missing args.
#
# Usage:
#   Rscript run_susier.R \
#     --input-tsv    matched_gwas.tsv \
#     --ld-matrix    ld_r.tsv \
#     --sample-size  50000 \
#     --out-tsv      results/susie_pip.tsv \
#     --diag-json    results/susie_diag.json \
#     --done-file    results/done.txt \
#     [--L           10]
#     [--coverage    0.95]
#     [--ld-shrink   0.1]
#     [--plot-pip    results/pip.png]
#
# Required input columns: SNP  BETA  SE  P  BP
# LD matrix: square TSV, first col = SNP IDs, values = raw r (NOT r^2)

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
    if (i == length(x)) { warning(sprintf("Flag '%s' has no value.", x[i])); break }
    out[[gsub("^--", "", x[i])]] <- x[i + 1]
    i <- i + 2
  }
  out
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))

get_opt <- function(key, default = NULL, required = FALSE) {
  v <- opt[[key]]
  if (is.null(v)) {
    if (required) stop(sprintf("Required argument --%s is missing.", key))
    return(default)
  }
  v
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
plot_pip    <- get_opt("plot-pip", NULL)

if (!is.finite(n_samples) || n_samples < 10)
  stop("--sample-size must be a positive integer (your GWAS N, e.g. 50000).")
if (!is.finite(L) || L < 1 || L > 50)
  stop("--L must be between 1 and 50.")
if (!is.finite(coverage) || coverage <= 0 || coverage >= 1)
  stop("--coverage must be in (0,1), e.g. 0.95.")
if (!is.finite(ld_shrink) || ld_shrink < 0 || ld_shrink > 1)
  stop("--ld-shrink must be in [0,1]. Default 0.1 works for most loci.")

for (p in list(out_tsv, diag_json, done_file))
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
if (!is.null(plot_pip))
  dir.create(dirname(plot_pip), recursive = TRUE, showWarnings = FALSE)

cat(sprintf("SuSiE settings: N=%d  L=%d  coverage=%.2f  ld_shrink=%.3f\n",
            n_samples, L, coverage, ld_shrink))

# --------------------------------------------------------------------------- #
# Load GWAS                                                                    #
# --------------------------------------------------------------------------- #
df <- fread(input_tsv, sep = "\t", check.names = FALSE, data.table = FALSE)

for (col in c("SNP","BETA","SE","P","BP"))
  if (!col %in% colnames(df)) stop(sprintf("Missing column: %s", col))

df$BETA <- suppressWarnings(as.numeric(df$BETA))
df$SE   <- suppressWarnings(as.numeric(df$SE))
df$P    <- suppressWarnings(as.numeric(df$P))
df$BP   <- suppressWarnings(as.numeric(df$BP))
df$SNP  <- as.character(df$SNP)

bad <- !is.finite(df$BETA) | !is.finite(df$SE) | df$SE <= 0 |
       !is.finite(df$P) | df$P <= 0 | df$P > 1 | !is.finite(df$BP)
if (any(bad)) { cat(sprintf("Dropping %d malformed rows.\n", sum(bad))); df <- df[!bad,] }
if (nrow(df) == 0) stop("No valid rows after filtering.")

# --------------------------------------------------------------------------- #
# Load and validate LD matrix                                                  #
# --------------------------------------------------------------------------- #
ld_data <- fread(ld_matrix_f, sep = "\t", header = TRUE, data.table = FALSE)
snp_ids <- as.character(ld_data[[1]])
R       <- as.matrix(ld_data[, -1, drop = FALSE])
storage.mode(R) <- "double"
rownames(R) <- snp_ids
colnames(R) <- snp_ids

if (nrow(R) != ncol(R)) stop("LD matrix is not square.")
if (nrow(R) != nrow(df))
  stop(sprintf("SNP count mismatch: GWAS=%d, LD=%d.", nrow(df), nrow(R)))

# Align GWAS rows to LD matrix order
match_idx <- match(snp_ids, df$SNP)
if (any(is.na(match_idx)))
  stop(sprintf("LD SNPs not in GWAS: %s",
               paste(snp_ids[is.na(match_idx)][1:min(5,sum(is.na(match_idx)))],
                     collapse=", ")))
df <- df[match_idx, , drop = FALSE]
if (!all(df$SNP == snp_ids)) stop("SNP alignment failed after reordering.")

# Check r vs r^2: susie_rss needs raw r (can be negative)
offdiag <- R[row(R) != col(R)]
if (!any(offdiag < -1e-6) && all(offdiag >= 0, na.rm=TRUE))
  cat("WARNING: LD matrix has no negative off-diagonal values.\n",
      "         susie_rss() needs raw r, not r^2. If your matrix is r^2,\n",
      "         take sqrt() of off-diagonal values before running.\n")

# Enforce symmetry and unit diagonal
max_asym <- max(abs(R - t(R)), na.rm=TRUE)
if (max_asym > 1e-4) stop(sprintf("LD not symmetric (max asym=%.2e).", max_asym))
R       <- (R + t(R)) / 2
diag(R) <- 1.0

# BUG 3 FIX: Ledoit-Wolf shrinkage to handle near-singular R.
# Without this, high-LD loci fragment one signal into many spurious CSs.
if (ld_shrink > 0) {
  R <- (1 - ld_shrink) * R + ld_shrink * diag(nrow(R))
  cat(sprintf("LD regularised: lambda=%.3f.\n", ld_shrink))
}

eig_min <- min(eigen(R, symmetric=TRUE, only.values=TRUE)$values)
cat(sprintf("Min eigenvalue after regularisation: %.4f\n", eig_min))
if (eig_min <= 0)
  stop("R is not positive definite. Increase --ld-shrink (try 0.2 or 0.3).")

# --------------------------------------------------------------------------- #
# Z-scores                                                                     #
# --------------------------------------------------------------------------- #
z_scores <- df$BETA / df$SE
cat(sprintf("%d SNPs, |Z| range [%.2f, %.2f]\n",
            nrow(df), min(abs(z_scores)), max(abs(z_scores))))

# --------------------------------------------------------------------------- #
# BUG 1+2 FIX: Run susie_rss() with correct sample size                       #
# susie_rss() is the maintained summary-stats API (susieR >= 0.11).           #
# It takes z-scores + LD correlation matrix directly.                         #
# There is NO silent fallback — if it fails, we stop with a useful message.   #
# --------------------------------------------------------------------------- #
cat("Running susie_rss()...\n")

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
      "  1. R is not positive definite -> increase --ld-shrink (try 0.2)\n",
      "  2. LD matrix contains r^2 not r -> sqrt() off-diagonals first\n",
      "  3. SNP order mismatch between GWAS and LD matrix\n",
      "  4. Z-scores contain NA/Inf\n"
    )
  }
)

cat(sprintf("Converged: %s\n", if(isTRUE(susie_fit$converged)) "YES" else "NO"))

# --------------------------------------------------------------------------- #
# BUG 4 FIX: Correct credible set extraction via susie_get_cs()               #
# Each element = one independent causal signal. Do NOT sum PIPs across CSs.   #
# --------------------------------------------------------------------------- #
pips    <- susie_fit$pip
# susieR API changed across versions; pick the supported argument name at runtime.
susie_get_cs_args <- names(formals(susie_get_cs))
cs_obj <- if ("Rr" %in% susie_get_cs_args) {
  susie_get_cs(susie_fit, coverage = coverage, Rr = R)
} else if ("Xcorr" %in% susie_get_cs_args) {
  susie_get_cs(susie_fit, coverage = coverage, Xcorr = R)
} else {
  susie_get_cs(susie_fit, coverage = coverage)
}
cs_list <- cs_obj$cs
n_cs    <- length(cs_list)

cat(sprintf("Credible sets found: %d\n", n_cs))

cs_member <- rep(NA_integer_, nrow(df))
cs_cover  <- rep(NA_real_,    nrow(df))

for (ci in seq_along(cs_list)) {
  idx          <- cs_list[[ci]]
  cs_member[idx] <- ci
  cs_cover[idx]  <- sum(pips[idx])
  top_i          <- idx[which.max(pips[idx])]
  cat(sprintf("  CS%d: %d SNPs, coverage=%.4f, top=%s (PIP=%.4f)\n",
              ci, length(idx), cs_cover[idx], df$SNP[top_i], pips[top_i]))
}

# --------------------------------------------------------------------------- #
# BUG 5 FIX: Output ALL SNPs with full PIP + CS info                          #
# --------------------------------------------------------------------------- #
result_df <- data.frame(
  SNP         = df$SNP,
  BP          = df$BP,
  P           = df$P,
  BETA        = df$BETA,
  SE          = df$SE,
  Z           = round(z_scores, 6),
  PIP         = round(pips, 6),
  CS          = cs_member,
  CS_coverage = round(cs_cover, 4),
  stringsAsFactors = FALSE
)
result_df <- result_df[order(-result_df$PIP), ]

write.table(result_df, file = out_tsv, sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)
cat(sprintf("Output: %d SNPs -> %s\n", nrow(result_df), out_tsv))

# --------------------------------------------------------------------------- #
# Optional PIP plot                                                            #
# --------------------------------------------------------------------------- #
if (!is.null(plot_pip)) {
  pd <- result_df[order(result_df$BP), ]
  pd$cs_label <- ifelse(is.na(pd$CS), "Not in CS", paste0("CS", pd$CS))

  cs_colors <- c("Not in CS" = "#9CA3AF",
    setNames(c("#D73027","#FC8D59","#4575B4","#1A9850","#F46D43",
               "#74ADD1","#A50026","#313695","#FDAE61","#ABD9E9"),
             paste0("CS", seq_len(min(n_cs, 10)))))

  p_plot <- ggplot(pd, aes(x = BP, y = PIP, colour = cs_label)) +
    geom_point(aes(size = ifelse(!is.na(CS), 2.5, 1.2)), alpha = 0.85, na.rm = TRUE) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    scale_colour_manual(values = cs_colors, name = "Credible set") +
    scale_size_identity() +
    scale_y_continuous(limits = c(0, 1.02), expand = c(0.01, 0)) +
    labs(x = "Chromosomal position (BP)", y = "PIP",
         title = sprintf("SuSiE  N=%d  L=%d  coverage=%.2f  |  %d CS found",
                         n_samples, L, coverage, n_cs)) +
    theme_bw(base_size = 11)

  ggsave(plot_pip, p_plot, width = 10, height = 4.8, dpi = 150)
  cat(sprintf("PIP plot: %s\n", plot_pip))
}

# --------------------------------------------------------------------------- #
# Diagnostics JSON                                                             #
# --------------------------------------------------------------------------- #
lead_idx <- which.max(pips)

cs_summary <- lapply(seq_along(cs_list), function(ci) {
  idx <- cs_list[[ci]]
  list(cs=ci, n_snps=length(idx), coverage=round(sum(pips[idx]),4),
       top_snp=df$SNP[idx[which.max(pips[idx])]], top_pip=round(max(pips[idx]),4),
       snps=df$SNP[idx])
})

writeLines(toJSON(list(
  method           = "susie_rss",
  susieR_version   = as.character(packageVersion("susieR")),
  sample_size      = n_samples,
  L                = L,
  coverage         = coverage,
  ld_shrink_lambda = ld_shrink,
  converged        = isTRUE(susie_fit$converged),
  n_variants       = nrow(df),
  n_credible_sets  = n_cs,
  lead_snp         = df$SNP[lead_idx],
  lead_pip         = round(pips[lead_idx], 4),
  credible_sets    = cs_summary
), pretty=TRUE, auto_unbox=TRUE), con = diag_json)

writeLines("ok", con = done_file)
cat(sprintf("\nDone: %d CS, lead=%s (PIP=%.4f)\n",
            n_cs, df$SNP[lead_idx], pips[lead_idx]))

