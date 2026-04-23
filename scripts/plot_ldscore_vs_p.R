#!/usr/bin/env Rscript
# Five-panel fine-mapping diagnostic plot
#   Panel 1 : |Z-score| track by position
#   Panel 2 : Normalized LD score track by position
#   Panel 3 : chi-sq vs LD score (LDSC regression)
#   Panel 4 : Causal score scatter  (|Z| vs LD score, bubble = -log10p)
#   Panel 5 : R^2 / D-prime quadrant map
#
# Required input columns (--input-tsv, tab-separated, header):
#   SNP   BETA   SE   P   CHR   BP
# Optional columns (used when present):
#   R2_index   DP_index   MAF
# LD matrix (--ld-matrix): square tab-separated, row names = SNP ids, values = r
#   (raw r, NOT r^2 - the script squares them internally)
#
# Usage:
#   Rscript plot_ldscore_vs_p.R \
#     --input-tsv    results.tsv \
#     --ld-matrix    ld.tsv \
#     --plot-png     out/plots.png \
#     --diag-tsv     out/diag.tsv \
#     --done-file    out/done.txt \
#     [--require-correlation true] \
#     [--min-correlation 0.05] \
#     [--index-snp rs12345]

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

PVALUE_CUTOFF <- 5e-4

# ?? Argument parsing ??????????????????????????????????????????????????????????
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  args[idx + 1]
}

input_tsv        <- get_arg("--input-tsv")
ld_matrix_file   <- get_arg("--ld-matrix")
plot_png         <- get_arg("--plot-png")
diag_tsv         <- get_arg("--diag-tsv")
done_file        <- get_arg("--done-file")
require_corr_str <- get_arg("--require-correlation", "false")
min_corr         <- as.numeric(get_arg("--min-correlation", "0.05"))
index_snp_arg    <- get_arg("--index-snp", NULL)   # optional: force a specific index SNP

if (any(sapply(list(input_tsv, ld_matrix_file, plot_png, diag_tsv, done_file), is.null))) {
  stop("Required: --input-tsv --ld-matrix --plot-png --diag-tsv --done-file")
}

require_corr <- tolower(trimws(require_corr_str)) %in% c("1", "true", "yes", "y", "on")

mkdir_for <- function(path) dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
mkdir_for(plot_png); mkdir_for(diag_tsv); mkdir_for(done_file)

# ?? Load GWAS data ????????????????????????????????????????????????????????????
gwas <- read.table(input_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

for (col in c("SNP", "BETA", "SE", "P")) {
  if (!col %in% colnames(gwas)) stop(sprintf("Column '%s' not found in input TSV", col))
}

gwas$BETA <- suppressWarnings(as.numeric(gwas$BETA))
gwas$SE   <- suppressWarnings(as.numeric(gwas$SE))
gwas$P    <- suppressWarnings(as.numeric(gwas$P))

# Remove malformed rows
gwas <- gwas[
  !is.na(gwas$BETA) & !is.na(gwas$SE) & gwas$SE > 0 &
  !is.na(gwas$P)    & gwas$P > 0 & gwas$P <= 1,
  , drop = FALSE
]

# Focus diagnostics on stronger association signals
gwas_sig <- gwas[gwas$P < PVALUE_CUTOFF, , drop = FALSE]

# ?? Empty-result fallback ?????????????????????????????????????????????????????
write_empty <- function(msg) {
  p_empty <- ggplot() +
    annotate("text", x = 0, y = 0, label = msg, size = 5) +
    theme_void()
  ggsave(plot_png, plot = p_empty, width = 14, height = 10, dpi = 150)
  write.table(
    data.frame(
      correlation_ldnorm_znormsq = NA_real_, n_snps = 0L,
      ld_norm_denominator = NA_real_,        z_mean = NA_real_,
      z_sd = NA_real_,                       pvalue_cutoff = PVALUE_CUTOFF
    ),
    file = diag_tsv, sep = "\t", row.names = FALSE, quote = FALSE
  )
  writeLines("ok", done_file)
  cat(msg, "\n")
  quit(save = "no", status = 0)
}

if (nrow(gwas_sig) == 0) {
  write_empty(sprintf("No SNPs pass P < %.1e", PVALUE_CUTOFF))
}

# ?? Z-score + normalisation ???????????????????????????????????????????????????
gwas_sig$zscore <- gwas_sig$BETA / gwas_sig$SE

z_mean <- mean(gwas_sig$zscore, na.rm = TRUE)
z_sd   <- sd(gwas_sig$zscore,   na.rm = TRUE)

gwas_sig$zscore_norm <- if (!is.finite(z_sd) || z_sd == 0) {
  0
} else {
  (gwas_sig$zscore - z_mean) / z_sd
}

# chi-sq = Z^2  (NOT abs(Z^2) - redundant, and confusing for signed diagnostics)
gwas_sig$chisq     <- gwas_sig$zscore^2
gwas_sig$chisq_norm <- gwas_sig$zscore_norm^2
gwas_sig$absz       <- abs(gwas_sig$zscore)
gwas_sig$logp       <- -log10(gwas_sig$P)

# ?? Load LD matrix ????????????????????????????????????????????????????????????
ld_raw <- tryCatch(
  read.table(ld_matrix_file, header = TRUE, sep = "\t",
             stringsAsFactors = FALSE, check.names = FALSE, row.names = 1),
  error = function(e) stop(sprintf("Cannot read LD matrix: %s", e$message))
)

ld_mat <- as.matrix(ld_raw)

# Detect whether matrix contains raw r or already r^2.
#
# BUG FIX 1 - the old test (all values >= 0 AND diagonal ~ 1) is ambiguous:
# a raw-r matrix from a single positively-correlated locus satisfies it too,
# so the squaring step was silently skipped even when it should have run.
#
# Correct test: r^2 is symmetric AND all off-diagonal values are in [0,1].
# Raw r has negative values whenever two SNPs are in repulsion LD.
# If *any* off-diagonal value is negative, the input is definitely raw r.
# If all values are non-negative we cannot be certain, so we check whether
# the off-diagonal values cluster below 1 the way r^2 does vs the way r does
# (raw r off-diagonals can exceed 0.9 across large stretches; r^2 rarely does
# for a typical locus). As a conservative default we require the user to
# confirm via a flag, but we also emit a clear warning.
diag_vals    <- diag(ld_mat)
offdiag_vals <- ld_mat[row(ld_mat) != col(ld_mat)]
offdiag_vals <- offdiag_vals[!is.na(offdiag_vals)]

has_negatives   <- any(offdiag_vals < -1e-6)
diag_all_one    <- all(abs(diag_vals[!is.na(diag_vals)] - 1) < 1e-6)
offdiag_max     <- max(abs(offdiag_vals), na.rm = TRUE)

# If negative values exist -> must be raw r. Square it.
# If all non-negative AND diagonal == 1 AND max off-diagonal <= 1 ->
#   ambiguous, but we lean toward r^2 only if values are strictly <= 1
#   and the matrix is symmetric (transpose equals itself within tolerance).
is_symmetric <- isTRUE(all.equal(ld_mat, t(ld_mat), tolerance = 1e-6,
                                  check.attributes = FALSE))
already_r2 <- !has_negatives && diag_all_one &&
              offdiag_max <= 1 + 1e-6 && is_symmetric

if (has_negatives) {
  cat("INFO: LD matrix contains negative values -> treating as raw r, squaring.\n")
} else if (already_r2) {
  cat("INFO: LD matrix looks like r^2 (all non-negative, diagonal=1, symmetric) -> no squaring.\n")
} else {
  cat("INFO: LD matrix ambiguous, assuming raw r -> squaring. Pass pre-squared r^2 matrix if wrong.\n")
  already_r2 <- FALSE
}

if (!already_r2) {
  ld_mat <- ld_mat^2
}

# Raw LD scores (sum of r^2 per SNP, excluding self).
# Subtract actual diagonal (not hard-coded 1) to handle imperfect diagonals.
ld_score <- rowSums(ld_mat, na.rm = TRUE) - diag(ld_mat)
ld_score[ld_score < 0] <- 0   # numerical safety

# BUG FIX 2 - the old denominator was ncol(ld_mat)-1, i.e. the full matrix
# size. When the LD matrix covers many more SNPs than the GWAS results that
# survive the p-value filter and the merge, every ld_score_norm value is
# crushed near zero because the denominator is enormous. This made panel 2
# look almost flat (only 2 visible points) and compressed panels 3 & 4 onto
# the far left of the x-axis - even though the raw LD scores were fine.
#
# Fix: normalise AFTER the merge, using only the SNPs that are actually
# present in the merged data frame. We store raw ld_score here and compute
# ld_score_norm later once we know the post-merge SNP set.
ld_df <- data.frame(
  SNP      = rownames(ld_mat),
  ld_score = ld_score,
  stringsAsFactors = FALSE
)

# ?? Merge GWAS + LD scores ????????????????????????????????????????????????????
# Keep ALL columns from gwas_sig (including R2_index, DP_index, MAF if present).
df <- merge(gwas_sig, ld_df, by = "SNP")
n  <- nrow(df)

if (n == 0) stop("No SNPs remain after merging GWAS with LD matrix - check that SNP IDs match.")

# BUG FIX 2 (continued): normalise ld_score using only the n SNPs that are
# actually in the merged data. This keeps the denominator meaningful - it
# reflects "out of the SNPs we are plotting, how many does each SNP tag?"
# rather than being diluted by thousands of unrelated SNPs in the full matrix.
ld_denominator   <- max(n - 1, 1)
df$ld_score_norm <- df$ld_score / ld_denominator
cat(sprintf("INFO: ld_score normalised by n_merged - 1 = %d (was full matrix ncol-1 = %d)\n",
            ld_denominator, max(ncol(ld_mat) - 1, 1)))

# ?? Optional columns: R2_index, DP_index, MAF ?????????????????????????????????
has_r2  <- "R2_index" %in% colnames(df)
has_dp  <- "DP_index" %in% colnames(df)
has_maf <- "MAF"      %in% colnames(df)

# ?? Identify index SNP ????????????????????????????????????????????????????????
if (!is.null(index_snp_arg) && index_snp_arg %in% df$SNP) {
  index_snp <- index_snp_arg
} else {
  index_snp <- df$SNP[which.min(df$P)]
}

df$is_index <- df$SNP == index_snp

# ?? Causal score = |Z| * (1 - R2) * (1 / LD_score_norm) ?????????????????????
# If R2_index is absent, approximate as 1 - (ld_score_norm / max(ld_score_norm))
# so tag SNPs in dense LD blocks are still penalised.
if (has_r2) {
  df$r2_use <- pmin(pmax(as.numeric(df$R2_index), 0), 1)
} else {
  max_ld <- max(df$ld_score_norm, na.rm = TRUE)
  df$r2_use <- if (max_ld > 0) df$ld_score_norm / max_ld else 0
}

safe_ld <- pmax(df$ld_score_norm, 0.01)   # avoid /0
df$causal_score <- df$absz * (1 - df$r2_use) * (1 / safe_ld)

# ?? Correlation (chi-sq ~ LD score) for LDSC diagnostics ?????????????????????
corr <- if (n > 1) cor(df$ld_score_norm, df$chisq_norm, use = "complete.obs") else 0.0

# ?? Shared aesthetics ?????????????????????????????????????????????????????????
pal_logp <- scale_color_gradient(low = "#74add1", high = "#d73027",
                                  name = expression(-log[10](P)))
theme_base <- theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.key.size = unit(0.4, "cm"))


# ?? Panel 1: |Z-score| by position ???????????????????????????????????????????
if ("BP" %in% colnames(df)) {
  df$plot_x <- as.numeric(df$BP)
  x_lab1 <- "Chromosomal position (BP)"
} else {
  df$plot_x <- seq_len(n)
  x_lab1 <- "SNP index"
}

index_row <- df[df$is_index, ]
p1 <- ggplot(df, aes(x = plot_x, y = absz)) +
  geom_point(aes(color = logp, size = logp), alpha = 0.8) +
  geom_point(data = index_row, aes(x = plot_x, y = absz),
             shape = 23, fill = "#d73027", color = "white", size = 4) +
  geom_hline(yintercept = qnorm(5e-8 / 2, lower.tail = FALSE),
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  pal_logp +
  scale_size_continuous(range = c(1.2, 4), name = expression(-log[10](P))) +
  labs(x = x_lab1, y = "|Z-score|", title = "Panel 1 - |Z-score| track") +
  theme_base

# ?? Panel 2: LD score by position ????????????????????????????????????????????
# Use a point/line track instead of bars: on BP-scale x values, bar widths in
# data units become visually negligible (e.g., width 0.8 bp on Mb axis).
p2 <- ggplot(df, aes(x = plot_x, y = ld_score_norm)) +
  geom_line(color = "#7f7f7f", linewidth = 0.35, alpha = 0.6) +
  geom_point(aes(color = ld_score_norm), size = 2.2, alpha = 0.9) +
  geom_point(data = index_row, aes(x = plot_x, y = ld_score_norm),
       shape = 23, fill = "#d73027", color = "white", size = 4) +
  scale_color_gradient(low = "#AFA9EC", high = "#534AB7",
           name = expression(LD~score~(l^2))) +
  labs(x = x_lab1, y = expression(Normalized~LD~score~(l^2)),
    title = "Panel 2 - Normalized LD score track") +
  theme_base

# ?? Panel 3: chi-sq vs LD score (LDSC regression) ????????????????????????????
p3 <- ggplot(df, aes(x = ld_score_norm, y = chisq_norm)) +
  geom_point(aes(color = logp, size = logp), alpha = 0.75) +
  geom_point(data = index_row, aes(x = ld_score_norm, y = chisq_norm),
             shape = 23, fill = "#d73027", color = "white", size = 4) +
  pal_logp +
  scale_size_continuous(range = c(1.2, 4), name = expression(-log[10](P))) +
  labs(
    x     = expression(Normalized~LD~score~(l^2)),
    y     = expression(Normalized~chi^2~(Z[norm]^2)),
    title = sprintf("Panel 3 - LDSC regression  (r=%.3f, n=%d)", corr, n)
  ) +
  theme_base

if (n > 1) {
  p3 <- p3 + geom_smooth(method = "lm", formula = y ~ x,
                          se = FALSE, linewidth = 0.9,
                          color = "#1D9E75", linetype = "dashed")
}

# ?? Panel 4: causal score scatter - |Z| vs LD score, bubble = -log10p ????????
p4 <- ggplot(df, aes(x = ld_score_norm, y = absz)) +
  geom_point(aes(color = r2_use, size = logp), alpha = 0.8) +
  geom_point(data = index_row, aes(x = ld_score_norm, y = absz),
             shape = 23, fill = "#d73027", color = "white", size = 5) +
  scale_color_gradient2(
    low = "#1D9E75", mid = "#EF9F27", high = "#d73027", midpoint = 0.5,
    name = if (has_r2) expression(R^2~(index)) else expression(LD~proxy~R^2)
  ) +
  scale_size_continuous(range = c(1.5, 5), name = expression(-log[10](P))) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
           label = "Causal zone", size = 3, color = "grey40", fontface = "italic") +
  labs(
    x     = expression(Normalized~LD~score~(l^2)),
    y     = "|Z-score|",
    title = "Panel 4 - Causal score scatter  (top-left = causal candidate)"
  ) +
  theme_base

# ?? Panel 5: R2 vs D-prime quadrant map (only if both columns present) ????????
if (has_r2 && has_dp) {
  df$r2_col <- as.numeric(df$R2_index)
  df$dp_col <- as.numeric(df$DP_index)

  p5 <- ggplot(df, aes(x = r2_col, y = dp_col)) +
    # quadrant shading
    annotate("rect", xmin = 0,   xmax = 0.5, ymin = 0.5, ymax = 1.05,
             fill = "#EF9F27", alpha = 0.08) +   # high D', low R2 = rare
    annotate("rect", xmin = 0.5, xmax = 1.05, ymin = 0.5, ymax = 1.05,
             fill = "#378ADD", alpha = 0.08) +   # high D', high R2 = tag
    annotate("rect", xmin = 0,   xmax = 0.5, ymin = -0.05, ymax = 0.5,
             fill = "#1D9E75", alpha = 0.08) +   # low D', low R2 = independent
    # quadrant dividers
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    # quadrant labels
    annotate("text", x = 0.25, y = 0.97, label = "Rare / haplotype",
             size = 2.8, color = "#BA7517", fontface = "italic") +
    annotate("text", x = 0.75, y = 0.97, label = "Tag SNP",
             size = 2.8, color = "#185FA5", fontface = "italic") +
    annotate("text", x = 0.25, y = 0.03, label = "Independent signal",
             size = 2.8, color = "#0F6E56", fontface = "italic") +
    annotate("text", x = 0.75, y = 0.03, label = "(rare in practice)",
             size = 2.8, color = "grey50", fontface = "italic") +
    geom_point(aes(color = logp, size = absz), alpha = 0.8) +
    geom_point(data = index_row, aes(x = r2_col, y = dp_col),
               shape = 23, fill = "#d73027", color = "white", size = 5) +
    pal_logp +
    scale_size_continuous(range = c(1.5, 5), name = "|Z-score|") +
    coord_cartesian(xlim = c(-0.02, 1.05), ylim = c(-0.02, 1.05)) +
    labs(
      x     = expression(R^2~with~index~SNP),
      y     = expression(D*minute~with~index~SNP),
      title = "Panel 5 - R2 / Dprime quadrant map"
    ) +
    theme_base
} else {
  # Fallback: show causal score ranked bar if R2/DP not available
  df_ranked <- df[order(df$causal_score, decreasing = TRUE), ]
  df_ranked$rank <- seq_len(nrow(df_ranked))
  df_top <- head(df_ranked, 30)

  p5 <- ggplot(df_top, aes(x = reorder(SNP, causal_score), y = causal_score)) +
    geom_col(aes(fill = logp), alpha = 0.85) +
    scale_fill_gradient(low = "#74add1", high = "#d73027",
                        name = expression(-log[10](P))) +
    coord_flip() +
    labs(
      x     = "SNP",
      y     = expression(Causal~score~("|Z|" %.% (1-R^2) %.% LD^{-1})),
      title = "Panel 5 - Causal score ranking  (top 30 SNPs)\n[Add R2_index + DP_index columns for quadrant map]"
    ) +
    theme_base
}

# ?? Assemble & save ???????????????????????????????????????????????????????????
combined <- gridExtra::arrangeGrob(p1, p2, p3, p4, p5, ncol = 2)
ggsave(plot_png, plot = combined, width = 14, height = 18, dpi = 150)

# ?? Diagnostics TSV ???????????????????????????????????????????????????????????
write.table(
  data.frame(
    index_snp                  = index_snp,
    correlation_ldnorm_znormsq = round(corr, 6),
    n_snps                     = n,
    ld_norm_denominator        = ld_denominator,   # = n_merged - 1 (post-fix)
    z_mean                     = round(z_mean, 6),
    z_sd                       = round(z_sd, 6),
    pvalue_cutoff              = PVALUE_CUTOFF
  ),
  file = diag_tsv, sep = "\t", row.names = FALSE, quote = FALSE
)

# ?? Correlation guard ?????????????????????????????????????????????????????????
if (require_corr && n >= 2 && corr < min_corr) {
  stop(sprintf(
    "chi-sq vs LD-score correlation too low: %.4f < %.4f (threshold)", corr, min_corr
  ))
}

writeLines("ok", done_file)
cat(sprintf(
  "Done. n_snps=%d  index_snp=%s  cor(LD_norm, chi2_norm)=%.4f\n",
  n, index_snp, corr
))
