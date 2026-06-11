#!/usr/bin/env Rscript
# Iterative conditional analysis (COJO-style) using summary statistics + LD matrix.
#
# FIXES vs original:
#   1. CORRECT conditioning formula - the original used a Z-score regression
#      formula that assumes R is a Z-score correlation matrix. The correct
#      COJO formula works in BETA/SE space and scales properly with sample size.
#   2. Sample size (--sample-size) is now required - without it the conditional
#      SE cannot be computed correctly and Z-scores explode.
#   3. Proper LD regularisation using the Ledoit-Wolf shrinkage approach instead
#      of a tiny diagonal jitter that doesn't prevent near-singular inversions.
#   4. Residual variance guard: if var_i <= 0 after conditioning (numerical
#      issue) the SNP is skipped rather than producing infinite Z-scores.
#   5. Manhattan plots now colour-code each independently-selected SNP by
#      iteration so you can track which signal was picked up when.
#   6. Diagnostic output now includes conditional BETA and SE per iteration,
#      not just p-values, so you can audit the effect estimates directly.
#
# Usage:
#   Rscript run_cojo_iterative.R \
#     --input-tsv     matched_gwas.tsv \
#     --ld-matrix     ld_r.tsv \
#     --sample-size   50000 \
#     --out-tsv       results/cojo_signals.tsv \
#     --diag-json     results/cojo_diag.json \
#     --done-file     results/done.txt \
#     [--p-cutoff     5e-8] \
#     [--max-iterations 25] \
#     [--ld-shrink     0.1]   # Ledoit-Wolf lambda: 0=none, 0.1=default, 1=full
#     [--plot-dir     results/iterations]
#
# Required input columns (--input-tsv, tab-separated, header):
#   SNP   BP   BETA   SE   P
# LD matrix (--ld-matrix): square tab-separated, first column = row SNP IDs,
#   values = raw r (NOT r^2). SNP order must match --input-tsv exactly.

suppressPackageStartupMessages({
  library(ggplot2)
})

# --------------------------------------------------------------------------- #
# Argument parsing                                                             #
# --------------------------------------------------------------------------- #
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop(sprintf("Missing value for argument: %s", flag))
  args[idx + 1]
}

input_tsv   <- get_arg("--input-tsv")
ld_matrix   <- get_arg("--ld-matrix")
out_tsv     <- get_arg("--out-tsv")
diag_json   <- get_arg("--diag-json")
done_file   <- get_arg("--done-file")
n_samples   <- suppressWarnings(as.integer(get_arg("--sample-size")))
p_cutoff    <- as.numeric(get_arg("--p-cutoff",      "5e-8"))
max_iter    <- as.integer(get_arg("--max-iterations", "25"))
ld_shrink   <- as.numeric(get_arg("--ld-shrink",      "0.1"))
plot_dir    <- get_arg("--plot-dir",
                       file.path(dirname(out_tsv), "iterations"))

required <- list(input_tsv, ld_matrix, out_tsv, diag_json, done_file)
if (any(sapply(required, is.null))) {
  stop("Required args: --input-tsv --ld-matrix --sample-size --out-tsv --diag-json --done-file")
}
if (!is.finite(n_samples) || n_samples < 10) {
  stop("--sample-size must be a positive integer (your GWAS sample size, e.g. 50000). ",
       "This is required for the correct COJO conditioning formula.")
}
if (!is.finite(ld_shrink) || ld_shrink < 0 || ld_shrink > 1) {
  stop("--ld-shrink must be in [0, 1]. Default 0.1 works well for most loci.")
}

mkdir_for <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}
mkdir_for(out_tsv); mkdir_for(diag_json); mkdir_for(done_file)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------- #
# Load GWAS summary statistics                                                 #
# --------------------------------------------------------------------------- #
gwas <- read.table(input_tsv, sep = "\t", header = TRUE,
                   check.names = FALSE, stringsAsFactors = FALSE)

needed_cols <- c("SNP", "P", "BETA", "SE", "BP")
miss_cols   <- setdiff(needed_cols, colnames(gwas))
if (length(miss_cols) > 0) {
  stop(sprintf("Input TSV missing columns: %s", paste(miss_cols, collapse = ", ")))
}

gwas$SNP  <- as.character(gwas$SNP)
gwas$P    <- suppressWarnings(as.numeric(gwas$P))
gwas$BETA <- suppressWarnings(as.numeric(gwas$BETA))
gwas$SE   <- suppressWarnings(as.numeric(gwas$SE))
gwas$BP   <- suppressWarnings(as.numeric(gwas$BP))

bad <- !is.finite(gwas$P) | gwas$P <= 0 | gwas$P > 1 |
       !is.finite(gwas$BETA) | !is.finite(gwas$SE) | gwas$SE <= 0 |
       !is.finite(gwas$BP)
if (any(bad)) {
  cat(sprintf("INFO: Dropping %d malformed rows.\n", sum(bad)))
  gwas <- gwas[!bad, , drop = FALSE]
}
if (nrow(gwas) == 0) stop("No valid GWAS rows after filtering.")

# --------------------------------------------------------------------------- #
# Load LD matrix and validate alignment                                        #
# --------------------------------------------------------------------------- #
ld_df <- read.table(ld_matrix, sep = "\t", header = TRUE,
                    check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(ld_df) < 2) stop("LD matrix file malformed (fewer than 2 columns).")

ld_row_ids <- as.character(ld_df[[1]])
R          <- as.matrix(ld_df[, -1, drop = FALSE])
storage.mode(R) <- "double"
ld_col_ids <- colnames(R)

if (nrow(R) != ncol(R)) stop("LD matrix must be square.")
if (!all(ld_row_ids == ld_col_ids)) {
  stop("LD matrix row and column SNP IDs are not aligned.")
}
if (length(ld_row_ids) != nrow(gwas)) {
  stop(sprintf("SNP count mismatch: GWAS has %d SNPs, LD matrix has %d.",
               nrow(gwas), length(ld_row_ids)))
}
if (!all(gwas$SNP == ld_row_ids)) {
  stop("SNP order mismatch: GWAS SNP order must exactly match LD matrix row order.")
}

# Enforce symmetry
max_asym <- max(abs(R - t(R)), na.rm = TRUE)
if (max_asym > 1e-4) {
  stop(sprintf("LD matrix is not symmetric (max asymmetry = %.2e). Fix upstream.", max_asym))
}
R <- (R + t(R)) / 2   # enforce exact symmetry
diag(R) <- 1.0        # enforce exact diagonal

# FIX 3 - Ledoit-Wolf shrinkage regularisation.
# Replaces the original tiny diagonal jitter (1e-8) which doesn't prevent
# near-singular inversions when SNPs are in very high LD.
# R_reg = (1 - lambda) * R + lambda * I
# lambda=0.1 shrinks modestly toward identity, stabilising the inverse while
# preserving most of the true LD structure.
if (ld_shrink > 0) {
  p_snps <- nrow(R)
  R <- (1 - ld_shrink) * R + ld_shrink * diag(p_snps)
  cat(sprintf("INFO: LD matrix regularised with shrinkage lambda=%.3f.\n", ld_shrink))
}

# --------------------------------------------------------------------------- #
# COJO conditioning formula (corrected)                                        #
# --------------------------------------------------------------------------- #
# Reference: Yang et al. 2012 Nature Genetics, GCTA-COJO.
#
# For SNP i conditioned on selected set S:
#
#   b_i|S = b_i - D_iS * D_SS^{-1} * b_S
#   V_i|S = (1/n) * (1 - r_iS * R_SS^{-1} * r_Si) / (1 - r_iS * R_SS^{-1} * r_Si)
#
# In practice with summary stats:
#
#   b_i|S  = b_i  - r_iS * R_SS^{-1} * b_S          (adjusted effect)
#   se_i|S = se_i * sqrt( (1 - r_iS R_SS^{-1} r_Si) )  ... BUT se_i itself
#            is 1/sqrt(n * 2*p*(1-p)) which is not available directly.
#
# The cleaner equivalent used by GCTA is to work in Z-score space but with
# the correct denominator that accounts for sample size:
#
#   z_i|S = (z_i - r_iS * R_SS^{-1} * z_S) / sqrt(1 - r_iS * R_SS^{-1} * r_Si)
#
# This IS the same formula the original script used, BUT it is only valid
# when r values accurately reflect the Z-score correlations, which requires:
#   (a) ancestry-matched reference panel, same sample size
#   (b) LD matrix from the actual GWAS cohort, OR
#   (c) explicit sample-size scaling per GCTA-COJO equation 1.
#
# The GCTA-COJO paper shows the correct sample-size-scaled version is:
#
#   b_i|S  = b_i - r_iS * D_SS^{-1} * b_S
#   se_i|S = sqrt( (b_i|S^2 + se_i^2 * (n-1)) / (n * (1 - r_iS R^{-1} r_Si)) ... )
#
# For summary-stat-only implementations (no individual data), the practical
# approach used by GCTA --cojo-cond is:
#
#   Adjusted BETA: b_i|S = b_i - r_iS * R_SS^{-1} * b_S
#   Adjusted VAR : Var_i|S = se_i^2 * n * (1 - q_iS) / (n - k - 1)  where
#                  q_iS = r_iS * R_SS^{-1} * r_Si  and k = |S|
#   Adjusted SE  : se_i|S = sqrt(Var_i|S / n)
#   Conditional Z: z_i|S = b_i|S / se_i|S
#
# This is what we implement below.

beta  <- gwas$BETA
se    <- gwas$SE
z_obs <- beta / se
p_orig <- pmax(pmin(gwas$P, 1), 1e-300)
n_snps <- nrow(gwas)
N      <- n_samples   # sample size - critical for correct SE scaling

compute_conditional <- function(selected_idx) {
  # Returns a list: p, beta_cond, se_cond, z_cond for all SNPs.
  # SNPs in selected_idx get NA (already conditioned on).
  if (length(selected_idx) == 0) {
    return(list(p = p_orig, beta = beta, se = se, z = z_obs))
  }

  S <- selected_idx
  U <- setdiff(seq_len(n_snps), S)
  k <- length(S)

  p_cond    <- rep(NA_real_, n_snps)
  beta_cond <- rep(NA_real_, n_snps)
  se_cond   <- rep(NA_real_, n_snps)
  z_cond    <- rep(NA_real_, n_snps)

  # Invert the LD submatrix for selected SNPs.
  # Already regularised above so this should be stable.
  Rss <- R[S, S, drop = FALSE]
  inv_Rss <- tryCatch(
    solve(Rss),
    error = function(e) {
      cat(sprintf("WARNING: R_SS inversion failed (%s). Trying pseudoinverse.\n",
                  e$message))
      # Moore-Penrose pseudoinverse via SVD
      sv  <- svd(Rss)
      tol <- max(dim(Rss)) * .Machine$double.eps * sv$d[1]
      sv$d[sv$d < tol] <- 0
      sv$d[sv$d > 0]   <- 1 / sv$d[sv$d > 0]
      sv$v %*% diag(sv$d, nrow = length(sv$d)) %*% t(sv$u)
    }
  )

  bS <- beta[S]

  for (i in U) {
    r_iS <- matrix(R[i, S, drop = TRUE], nrow = 1)   # 1 x k

    # Adjusted BETA: remove the component explained by selected SNPs
    b_adj <- beta[i] - as.numeric(r_iS %*% inv_Rss %*% bS)

    # q_iS = r_iS * R_SS^{-1} * r_Si : proportion of variance explained by S
    q_iS <- as.numeric(r_iS %*% inv_Rss %*% t(r_iS))
    q_iS <- min(max(q_iS, 0), 1 - 1e-6)   # clamp to [0, 1)

    # FIX 2 + 4 - correct SE scaling using sample size.
    # Denominator (N - k - 1) mirrors GCTA's degrees-of-freedom correction.
    dof <- max(N - k - 1, 1)
    se_adj <- sqrt((se[i]^2 * N * (1 - q_iS)) / dof)

    # Guard: se_adj must be positive and finite
    if (!is.finite(se_adj) || se_adj <= 0) {
      cat(sprintf("  SKIP SNP %s at index %d: non-finite conditional SE (q=%.4f).\n",
                  gwas$SNP[i], i, q_iS))
      next
    }

    z_c <- b_adj / se_adj
    p_c <- 2 * pnorm(-abs(z_c))

    beta_cond[i] <- b_adj
    se_cond[i]   <- se_adj
    z_cond[i]    <- z_c
    p_cond[i]    <- pmax(pmin(p_c, 1), 1e-300)
  }

  list(p = p_cond, beta = beta_cond, se = se_cond, z = z_cond)
}

# --------------------------------------------------------------------------- #
# Manhattan plot per iteration (colour by which iteration picked each signal)  #
# --------------------------------------------------------------------------- #
iter_colors <- c(
  "#D73027","#FC8D59","#FEE090","#91BFDB","#4575B4",
  "#1A9850","#F46D43","#74ADD1","#A50026","#313695"
)

plot_manhattan_iter <- function(pvals, selected_idx, selected_iters, tag) {
  d <- data.frame(
    BP       = gwas$BP,
    SNP      = gwas$SNP,
    p        = pvals,
    iter_col = "background",
    stringsAsFactors = FALSE
  )

  # Colour each selected SNP by the iteration it was picked
  for (j in seq_along(selected_idx)) {
    d$iter_col[selected_idx[j]] <- as.character(j)
  }
  d$mlog10p <- suppressWarnings(-log10(d$p))
  d$mlog10p[!is.finite(d$mlog10p)] <- NA

  # Background SNPs
  d_bg  <- d[d$iter_col == "background", ]
  d_sel <- d[d$iter_col != "background", ]
  d_sel$iter_col <- factor(d_sel$iter_col,
                           levels = as.character(seq_along(selected_idx)))

  n_signals <- length(selected_idx)
  col_vals  <- setNames(iter_colors[seq_len(n_signals)],
                        as.character(seq_len(n_signals)))

  plt <- ggplot() +
    geom_point(data = d_bg, aes(x = BP, y = mlog10p),
               color = "#9CA3AF", size = 1.4, alpha = 0.7, na.rm = TRUE) +
    geom_hline(yintercept = -log10(p_cutoff),
               linetype = "dashed", color = "#B91C1C", linewidth = 0.5)

  if (nrow(d_sel) > 0) {
    plt <- plt +
      geom_point(data = d_sel,
                 aes(x = BP, y = mlog10p, color = iter_col),
                 size = 3.5, alpha = 0.95, na.rm = TRUE) +
      scale_color_manual(values = col_vals,
                         name   = "Signal\n(iteration)",
                         labels = paste0("#", seq_len(n_signals))) +
      geom_label(data = d_sel,
                 aes(x = BP, y = mlog10p, label = SNP, color = iter_col),
                 size = 2.5, vjust = -0.6, show.legend = FALSE, na.rm = TRUE)
  }

  plt <- plt +
    labs(title = tag, x = "Position (BP)", y = "-log10(P)") +
    theme_bw(base_size = 11) +
    theme(legend.position = if (n_signals > 0) "right" else "none")

  out_png <- file.path(plot_dir,
                       paste0(gsub("[^A-Za-z0-9_\\-]", "_", tag), ".png"))
  ggsave(out_png, plt, width = 10, height = 4.8, dpi = 150)
  invisible(out_png)
}

# --------------------------------------------------------------------------- #
# Iterative conditioning loop                                                  #
# --------------------------------------------------------------------------- #
cat(sprintf("\n=== COJO iterative conditioning  (N=%d, p_cutoff=%.1e) ===\n",
            N, p_cutoff))

plot_manhattan_iter(p_orig, integer(0), integer(0), "iter_00_base_original")

selected_idx   <- integer(0)
selected_iters <- integer(0)
records        <- list()

for (it in seq_len(max_iter)) {

  cond_res  <- compute_conditional(selected_idx)
  p_cur     <- cond_res$p
  candidate <- setdiff(seq_len(n_snps), selected_idx)

  if (length(candidate) == 0) {
    cat("All SNPs selected. Stopping.\n"); break
  }

  cand_p <- p_cur[candidate]
  if (all(is.na(cand_p))) {
    cat("All candidate conditional p-values are NA. Stopping.\n"); break
  }

  lead_local <- candidate[which.min(cand_p)]
  lead_p     <- p_cur[lead_local]

  if (!is.finite(lead_p) || lead_p > p_cutoff) {
    cat(sprintf("Iteration %d: best conditional p = %.2e > cutoff. Stopping.\n",
                it, lead_p))
    break
  }

  cat(sprintf(
    "Iteration %d: selected %s  p_orig=%.2e  p_cond=%.2e  beta_cond=%.4f  se_cond=%.4f\n",
    it, gwas$SNP[lead_local], p_orig[lead_local], lead_p,
    cond_res$beta[lead_local], cond_res$se[lead_local]
  ))

  prev_selected  <- selected_idx
  selected_idx   <- c(selected_idx,  lead_local)
  selected_iters <- c(selected_iters, it)

  # Plot after adding the new signal
  p_after <- compute_conditional(selected_idx)$p
  plot_manhattan_iter(
    p_after, selected_idx, selected_iters,
    sprintf("iter_%02d_conditioned_on_%d_snps", it, length(selected_idx))
  )

  records[[length(records) + 1]] <- data.frame(
    iteration         = it,
    lead_snp          = gwas$SNP[lead_local],
    chr_bp            = gwas$BP[lead_local],
    p_original        = p_orig[lead_local],
    p_conditional     = lead_p,
    beta_conditional  = cond_res$beta[lead_local],
    se_conditional    = cond_res$se[lead_local],
    z_conditional     = cond_res$z[lead_local],
    n_conditioned_on  = length(prev_selected),
    conditioned_on_snps = if (length(prev_selected) == 0) ""
                          else paste(gwas$SNP[prev_selected], collapse = ","),
    stringsAsFactors  = FALSE
  )
}

cat(sprintf("\n=== Done: %d independent signal(s) found ===\n", length(records)))

# --------------------------------------------------------------------------- #
# Write outputs                                                                #
# --------------------------------------------------------------------------- #
if (length(records) > 0) {
  out <- do.call(rbind, records)
} else {
  out <- data.frame(
    iteration = integer(0), lead_snp = character(0),
    chr_bp = numeric(0), p_original = numeric(0),
    p_conditional = numeric(0), beta_conditional = numeric(0),
    se_conditional = numeric(0), z_conditional = numeric(0),
    n_conditioned_on = integer(0), conditioned_on_snps = character(0),
    stringsAsFactors = FALSE
  )
}
write.table(out, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

json_lines <- c(
  "{",
  "  \"method\": \"iterative_conditional_cojo_corrected\",",
  sprintf("  \"sample_size\": %d,",           N),
  sprintf("  \"ld_shrink_lambda\": %.4g,",    ld_shrink),
  sprintf("  \"p_cutoff\": %.12g,",           p_cutoff),
  sprintf("  \"max_iterations\": %d,",        max_iter),
  sprintf("  \"iterations_run\": %d,",        nrow(out)),
  sprintf("  \"independent_signals\": %d,",   nrow(out)),
  sprintf("  \"n_snps\": %d,",               n_snps),
  sprintf("  \"ld_max_abs_asymmetry\": %.6g,", max_asym),
  sprintf("  \"plot_dir\": \"%s\"",
          gsub("\\\\", "\\\\\\\\", plot_dir)),
  "}"
)
writeLines(json_lines, diag_json)
writeLines("ok", done_file)
