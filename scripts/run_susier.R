#!/usr/bin/env Rscript
# SuSiE fine-mapping with bigsnpr LD computation from PLINK files
#
# KEY CHANGES vs previous version:
#   - Accepts PLINK binary (.bed/.bim/.fam) via --plink-prefix instead of
#     a pre-computed LD matrix TSV.
#   - Uses bigsnpr::snp_cor() to compute RAW r (not r^2) directly from
#     genotypes. This was the root cause of the 10-CS explosion.
#   - Window (--window-kb) is applied to BOTH the GWAS data AND the PLINK
#     bim before LD computation, so the LD matrix covers exactly the SNPs
#     that enter SuSiE.
#   - SNP matching: GWAS <-> PLINK bim by rsID, with optional allele
#     harmonisation (A1/A2 flip detection) if A1/A2 columns are present.
#   - MAF filter applied inside the window before computing LD.
#   - LD correction pipeline:
#       1. Clip values to [-1, 1]
#       2. Enforce symmetry
#       3. Nearest PD projection (Matrix::nearPD) if any negative eigenvalues
#       4. Ledoit-Wolf shrinkage toward identity
#   - All diagnostic logs, JSON and plots from the previous version kept.
#
# Usage:
#   Rscript run_susier.R \
#     --input-tsv      matched_gwas.tsv \
#     --plink-prefix   /data/plink/cohort_chr22 \
#     --sample-size    159219 \
#     --out-tsv        results/susie_pip.tsv \
#     --diag-json      results/susie_diag.json \
#     --diag-log       results/susie_diag.log \
#     --plot-overview  results/susie_overview.png \
#     --plot-pip       results/susie_pip.png \
#     --plot-ld        results/susie_ld.png \
#     --plot-zscore    results/susie_zscore.png \
#     --done-file      results/done.txt \
#     [--window-kb     250]
#     [--L             10]
#     [--coverage      0.95]
#     [--ld-shrink     0.1]
#     [--min-maf       0.01]
##
# Required input-tsv columns: SNP  BETA  SE  P  BP
# Optional columns: A1  A2  (enables allele harmonisation)

suppressPackageStartupMessages({
  library(bigsnpr)
  library(susieR)
  library(data.table)
  library(jsonlite)
  library(ggplot2)
  library(gridExtra)
  library(Matrix)
})

# --------------------------------------------------------------------------- #
# Argument parsing                                                             #
# --------------------------------------------------------------------------- #
parse_args <- function(x) {
  out <- list(); i <- 1
  while (i <= length(x)) {
    if (i == length(x)) { warning(sprintf("Flag '%s' has no value.", x[i])); break }
    out[[gsub("^--", "", x[i])]] <- x[i + 1]; i <- i + 2
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

input_tsv     <- get_opt("input-tsv",     required = TRUE)
plink_prefix  <- get_opt("plink-prefix",  required = TRUE)
out_tsv       <- get_opt("out-tsv",       required = TRUE)
diag_json     <- get_opt("diag-json",     required = TRUE)
diag_log      <- get_opt("diag-log",      required = TRUE)
done_file     <- get_opt("done-file",     required = TRUE)
plot_overview <- get_opt("plot-overview", NULL)
plot_pip      <- get_opt("plot-pip",      NULL)
plot_ld       <- get_opt("plot-ld",       NULL)
plot_zscore   <- get_opt("plot-zscore",   NULL)
n_samples     <- suppressWarnings(as.integer(get_opt("sample-size", required = TRUE)))
L             <- suppressWarnings(as.integer(get_opt("L",           "10")))
coverage      <- suppressWarnings(as.numeric(get_opt("coverage",    "0.95")))
ld_shrink     <- suppressWarnings(as.numeric(get_opt("ld-shrink",   "0.1")))
window_kb     <- suppressWarnings(as.numeric(get_opt("window-kb",   "250")))
min_maf       <- suppressWarnings(as.numeric(get_opt("min-maf",     "0.01")))

if (!is.finite(n_samples) || n_samples < 10)
  stop("--sample-size must be a positive integer.")
if (!is.finite(L) || L < 1 || L > 50)
  stop("--L must be 1-50.")
if (!is.finite(coverage) || coverage <= 0 || coverage >= 1)
  stop("--coverage must be in (0, 1).")
if (!is.finite(ld_shrink) || ld_shrink < 0 || ld_shrink > 1)
  stop("--ld-shrink must be in [0, 1].")
if (!is.finite(window_kb) || window_kb <= 0)
  stop("--window-kb must be positive.")

for (ext in c(".bed", ".bim", ".fam")) {
  f <- paste0(plink_prefix, ext)
  if (!file.exists(f)) stop(sprintf("PLINK file not found: %s", f))
}

for (p in list(out_tsv, diag_json, diag_log, done_file))
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
for (p in list(plot_overview, plot_pip, plot_ld, plot_zscore))
  if (!is.null(p)) dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------- #
# Logging helpers                                                              #
# --------------------------------------------------------------------------- #
log_lines <- character(0)
LOG     <- function(...) { msg <- paste0(...); cat(msg, "\n"); log_lines <<- c(log_lines, msg) }
LOGSEP  <- function() LOG(paste(rep("=", 70), collapse = ""))
LOGWARN <- function(...) LOG("  [WARNING] ", ...)
LOGOK   <- function(...) LOG("  [OK]      ", ...)
LOGFAIL <- function(...) LOG("  [FAIL]    ", ...)
LOGINFO <- function(...) LOG("  [INFO]    ", ...)

LOGSEP()
LOG("SuSiE + bigsnpr DIAGNOSTIC RUN")
LOG("Timestamp    : ", format(Sys.time()))
LOG("Input TSV    : ", input_tsv)
LOG("PLINK prefix : ", plink_prefix)
LOG(sprintf("Settings     : N=%d  L=%d  cov=%.2f  shrink=%.3f  window_kb=%.0f  min_maf=%.3f",
            n_samples, L, coverage, ld_shrink, window_kb, min_maf))
LOGSEP()

# --------------------------------------------------------------------------- #
# SECTION 1: Load GWAS summary statistics                                      #
# --------------------------------------------------------------------------- #
LOG("SECTION 1: GWAS summary statistics")

df_all <- fread(input_tsv, sep = "\t", check.names = FALSE, data.table = FALSE)

for (col in c("SNP", "BETA", "SE", "P", "BP"))
  if (!col %in% colnames(df_all))
    stop(sprintf("Required column missing from input TSV: %s", col))

df_all$SNP  <- as.character(df_all$SNP)
df_all$BETA <- suppressWarnings(as.numeric(df_all$BETA))
df_all$SE   <- suppressWarnings(as.numeric(df_all$SE))
df_all$P    <- suppressWarnings(as.numeric(df_all$P))
df_all$BP   <- suppressWarnings(as.numeric(df_all$BP))

has_alleles <- all(c("A1", "A2") %in% colnames(df_all))
if (has_alleles) {
  df_all$A1 <- toupper(as.character(df_all$A1))
  df_all$A2 <- toupper(as.character(df_all$A2))
  LOGINFO("A1/A2 columns found -- allele harmonisation will be performed.")
} else {
  LOGWARN("No A1/A2 columns. Allele flip check will be skipped.")
}

bad <- !is.finite(df_all$BETA) | !is.finite(df_all$SE) | df_all$SE <= 0 |
       !is.finite(df_all$P) | df_all$P <= 0 | df_all$P > 1 | !is.finite(df_all$BP)
if (any(bad)) {
  LOGWARN(sprintf("Dropping %d malformed rows.", sum(bad)))
  df_all <- df_all[!bad, , drop = FALSE]
}

LOGINFO(sprintf("Total SNPs loaded    : %d", nrow(df_all)))
LOGINFO(sprintf("BP range             : %d -- %d  (%.2f Mb)",
                min(df_all$BP), max(df_all$BP),
                (max(df_all$BP) - min(df_all$BP)) / 1e6))

# --------------------------------------------------------------------------- #
# SECTION 2: Window definition                                                 #
# --------------------------------------------------------------------------- #
LOGSEP()
LOG("SECTION 2: Index SNP and window")

# Index SNP = min-P SNP in the input TSV.
# The GWAS input already contains only the locus region, so the strongest
# association in the table is by definition the index SNP for this locus.
index_row <- df_all[which.min(df_all$P), , drop = FALSE]

index_snp <- index_row$SNP[1]
index_bp  <- index_row$BP[1]
index_p   <- index_row$P[1]
bp_lo     <- index_bp - window_kb * 1000
bp_hi     <- index_bp + window_kb * 1000

LOGINFO(sprintf("Index SNP  : %s  P=%.2e  BP=%d", index_snp, index_p, index_bp))
LOGINFO(sprintf("Window BP  : %d -- %d  (+/- %.0f kb)", bp_lo, bp_hi, window_kb))

df <- df_all[df_all$BP >= bp_lo & df_all$BP <= bp_hi, , drop = FALSE]
df <- df[order(df$BP), ]

LOGINFO(sprintf("GWAS SNPs in window  : %d", nrow(df)))
if (nrow(df) == 0) stop("No GWAS SNPs in window. Adjust --window-kb.")
if (nrow(df) > 2000) LOGWARN(sprintf("Window has %d SNPs -- consider narrowing.", nrow(df)))

# --------------------------------------------------------------------------- #
# SECTION 3: PLINK bim loading and intersection                                #
# --------------------------------------------------------------------------- #
LOGSEP()
LOG("SECTION 3: PLINK bim and SNP intersection")

bim_file <- paste0(plink_prefix, ".bim")
bim <- fread(bim_file, header = FALSE, data.table = FALSE,
             col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
bim$SNP <- as.character(bim$SNP)
bim$A1  <- toupper(bim$A1)
bim$A2  <- toupper(bim$A2)

LOGINFO(sprintf("PLINK bim total SNPs : %d", nrow(bim)))

# Filter bim to window
bim_win <- bim[bim$BP >= bp_lo & bim$BP <= bp_hi, , drop = FALSE]
LOGINFO(sprintf("PLINK SNPs in window : %d", nrow(bim_win)))

# --------------------------------------------------------------------------- #
# SECTION 4: SNP matching and allele harmonisation                             #
# --------------------------------------------------------------------------- #
LOGSEP()
LOG("SECTION 4: SNP matching and allele harmonisation")

common_snps  <- intersect(df$SNP, bim_win$SNP)
n_gwas_only  <- sum(!df$SNP %in% bim_win$SNP)
n_plink_only <- sum(!bim_win$SNP %in% df$SNP)

LOGINFO(sprintf("GWAS window SNPs     : %d", nrow(df)))
LOGINFO(sprintf("PLINK window SNPs    : %d", nrow(bim_win)))
LOGINFO(sprintf("Matched by SNP ID    : %d", length(common_snps)))
LOGINFO(sprintf("GWAS only (no PLINK) : %d", n_gwas_only))
LOGINFO(sprintf("PLINK only (no GWAS) : %d", n_plink_only))

if (length(common_snps) < 10)
  stop(sprintf("Only %d SNPs matched. Check SNP ID format in GWAS vs PLINK bim.", length(common_snps)))
if (n_gwas_only > nrow(df) * 0.25) {
  LOGWARN(sprintf(">25%% GWAS SNPs not in PLINK (%d SNPs).", n_gwas_only))
} else {
  LOGOK(sprintf("%.1f%% GWAS SNPs matched in PLINK.", 100 * length(common_snps) / nrow(df)))
}

# Subset and sort both by BP order
df      <- df[df$SNP %in% common_snps, , drop = FALSE]
df      <- df[order(df$BP), ]
bim_win <- bim_win[bim_win$SNP %in% common_snps, , drop = FALSE]
bim_win <- bim_win[order(bim_win$BP), ]

# Align df rows to bim_win order (crucial for LD matrix ordering)
df <- df[match(bim_win$SNP, df$SNP), , drop = FALSE]
stopifnot(all(df$SNP == bim_win$SNP))
LOGOK(sprintf("SNP order verified: %d SNPs aligned in BP order.", nrow(df)))

# Allele harmonisation
n_flip  <- 0L
n_ambig <- 0L
flip_idx <- rep(FALSE, nrow(df))

if (has_alleles) {
  LOG("  Harmonising alleles (GWAS A1/A2 vs PLINK bim A1/A2) ...")
  complement <- c(A = "T", T = "A", C = "G", G = "C")

  for (k in seq_len(nrow(df))) {
    g_a1 <- df$A1[k]; g_a2 <- df$A2[k]
    p_a1 <- bim_win$A1[k]; p_a2 <- bim_win$A2[k]
    g_a1c <- complement[g_a1]; g_a2c <- complement[g_a2]

    # Perfect match
    if (isTRUE(g_a1 == p_a1) && isTRUE(g_a2 == p_a2)) next

    # Allele swap: effect coded on opposite allele
    if (isTRUE(g_a1 == p_a2) && isTRUE(g_a2 == p_a1)) {
      flip_idx[k] <- TRUE; n_flip <- n_flip + 1L; next
    }

    # Strand complement
    if (!is.na(g_a1c) && isTRUE(g_a1c == p_a1) && !is.na(g_a2c) && isTRUE(g_a2c == p_a2)) next

    # Strand complement + swap
    if (!is.na(g_a1c) && isTRUE(g_a1c == p_a2) && !is.na(g_a2c) && isTRUE(g_a2c == p_a1)) {
      flip_idx[k] <- TRUE; n_flip <- n_flip + 1L; next
    }

    # Ambiguous palindromic SNP (A/T or C/G)
    if ((g_a1 %in% c("A","T") && g_a2 %in% c("A","T")) ||
        (g_a1 %in% c("C","G") && g_a2 %in% c("C","G"))) {
      n_ambig <- n_ambig + 1L; next
    }

    # True mismatch
    LOGWARN(sprintf("Allele mismatch at %s: GWAS %s/%s vs PLINK %s/%s -- SNP kept but verify.",
                    df$SNP[k], g_a1, g_a2, p_a1, p_a2))
  }

  # Apply flips: negate BETA so effect is on the same allele as PLINK A1
  if (any(flip_idx)) df$BETA[flip_idx] <- -df$BETA[flip_idx]
  LOGINFO(sprintf("Allele flips applied : %d", n_flip))
  LOGINFO(sprintf("Ambiguous SNPs       : %d  (A/T or C/G -- not flipped)", n_ambig))
  if (n_ambig > nrow(df) * 0.1)
    LOGWARN("Many ambiguous SNPs. Consider removing A/T and C/G SNPs.")
} else {
  n_ambig <- 0L
  LOGWARN("Allele harmonisation skipped. Z-score sign consistency may be affected.")
}

n_snps <- nrow(df)
LOGINFO(sprintf("SNPs entering LD computation: %d", n_snps))

# --------------------------------------------------------------------------- #
# SECTION 5: bigsnpr -- compute raw r LD matrix                               #
# --------------------------------------------------------------------------- #
LOGSEP()
LOG("SECTION 5: LD computation with bigsnpr")

# bigsnpr stores a .bk/.rds backing file; use a temp directory
tmp_dir     <- tempdir()
backing_file <- file.path(tmp_dir, paste0("bigsnp_", basename(plink_prefix)))
rds_path    <- paste0(backing_file, ".rds")

if (!file.exists(rds_path)) {
  LOG("  Reading PLINK bed file into bigsnpr backing store ...")
  rds_created <- snp_readBed(paste0(plink_prefix, ".bed"), backingfile = backing_file)
  snp_obj <- snp_attach(rds_created)
} else {
  LOG("  Attaching existing bigsnpr backing store ...")
  snp_obj <- snp_attach(rds_path)
}

G        <- snp_obj$genotypes
bim_full <- snp_obj$map
fam      <- snp_obj$fam
n_plink_samples <- nrow(fam)

LOGINFO(sprintf("PLINK samples        : %d", n_plink_samples))
LOGINFO(sprintf("PLINK SNPs total     : %d", nrow(bim_full)))

# Map windowed SNPs to indices in the full bigsnpr object
# bigsnpr map uses 'marker.ID' as the SNP ID column
bim_full$SNP_id <- as.character(bim_full$marker.ID)
plink_idx <- match(df$SNP, bim_full$SNP_id)

if (any(is.na(plink_idx)))
  stop(sprintf("%d windowed SNPs not found in full PLINK bim after reading. Check IDs: %s",
               sum(is.na(plink_idx)),
               paste(df$SNP[is.na(plink_idx)][1:min(5, sum(is.na(plink_idx)))], collapse = ", ")))

LOGINFO(sprintf("PLINK column indices : %d to %d", min(plink_idx), max(plink_idx)))

# MAF filter within the window
LOG("  Computing MAF for windowed SNPs ...")
maf_vals <- snp_MAF(G, ind.col = plink_idx)
low_maf  <- maf_vals < min_maf

LOGINFO(sprintf("MAF >= %.3f          : %d SNPs  (removed: %d)", min_maf, sum(!low_maf), sum(low_maf)))
if (any(low_maf)) {
  LOGWARN(sprintf("Removing %d low-MAF SNPs from LD computation.", sum(low_maf)))
  plink_idx <- plink_idx[!low_maf]
  df        <- df[!low_maf, , drop = FALSE]
  flip_idx  <- flip_idx[!low_maf]
  n_snps    <- nrow(df)
  LOGINFO(sprintf("SNPs after MAF filter: %d", n_snps))
}

# Compute raw r correlation matrix
# snp_cor() computes Pearson correlation between columns (SNPs) of G
# This gives raw r in [-1, 1], which is exactly what susie_rss() requires.
LOG(sprintf("  Computing %d x %d raw r matrix (snp_cor) ...", n_snps, n_snps))
t_ld <- proc.time()

R_bigsnp <- snp_cor(
  Gna     = G,
  ind.col = plink_idx,
  ncores  = max(1L, parallel::detectCores() - 1L)
)

t_ld_elapsed <- (proc.time() - t_ld)[["elapsed"]]
LOG(sprintf("  LD computation done: %.1f seconds", t_ld_elapsed))

R_raw <- as.matrix(R_bigsnp)
rownames(R_raw) <- df$SNP
colnames(R_raw) <- df$SNP

# --------------------------------------------------------------------------- #
# SECTION 6: LD matrix diagnostics and correction                              #
# --------------------------------------------------------------------------- #
LOGSEP()
LOG("SECTION 6: LD matrix validation and correction")

# 6a: Check raw r properties
offdiag    <- R_raw[row(R_raw) != col(R_raw)]
n_negative <- sum(offdiag < -1e-6, na.rm = TRUE)
n_above1   <- sum(abs(offdiag) > 1 + 1e-6, na.rm = TRUE)
oq         <- quantile(offdiag, c(0.01, 0.25, 0.5, 0.75, 0.99), na.rm = TRUE)

LOGINFO(sprintf("Off-diagonal negatives : %d  (should be >0 for raw r)", n_negative))
LOGINFO(sprintf("Off-diagonal |val| > 1 : %d  (should be 0)", n_above1))
LOGINFO(sprintf("Off-diag quantiles [1,25,50,75,99]%%: %.3f %.3f %.3f %.3f %.3f",
                oq[1], oq[2], oq[3], oq[4], oq[5]))

if (n_negative > 0) {
  LOGOK("Negative values present -- bigsnpr computed raw r correctly.")
} else {
  LOGWARN("No negative off-diagonal values. Unexpected for raw r from real genotypes.")
  LOGWARN("Check that G contains dosage genotypes (0/1/2), not something else.")
}

# 6b: Clip values to [-1, 1] (numerical safety)
if (n_above1 > 0) {
  R_raw <- pmax(pmin(R_raw, 1.0), -1.0)
  LOGWARN(sprintf("%d values clipped to [-1, 1].", n_above1))
}

# 6c: Enforce exact symmetry
max_asym <- max(abs(R_raw - t(R_raw)), na.rm = TRUE)
LOGINFO(sprintf("Max asymmetry          : %.2e", max_asym))
R <- (R_raw + t(R_raw)) / 2
diag(R) <- 1.0
LOGOK("Symmetry enforced: R = (R + R') / 2, diagonal set to 1.")

# 6d: Eigenvalue spectrum
LOG("  Computing eigenvalue spectrum ...")
eig_vals <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
eig_min  <- min(eig_vals)
eig_max  <- max(eig_vals)
cond_num <- eig_max / max(abs(eig_min), 1e-12)
eff_rank <- sum(eig_vals > 0.01 * eig_max)
pct_top5 <- 100 * sum(head(sort(eig_vals, decreasing = TRUE), 5)) / sum(eig_vals)

LOGINFO(sprintf("Min eigenvalue (pre-correction) : %.6f", eig_min))
LOGINFO(sprintf("Max eigenvalue                  : %.6f", eig_max))
LOGINFO(sprintf("Condition number                : %.1f", cond_num))
LOGINFO(sprintf("Effective rank                  : %d / %d", eff_rank, n_snps))
LOGINFO(sprintf("Top-5 eigenvectors variance     : %.1f%%", pct_top5))

if (eig_min < 0) {
  LOGWARN(sprintf("Negative eigenvalues (min=%.4f). Applying Matrix::nearPD projection.", eig_min))
}
if (cond_num > 1e5) {
  LOGWARN(sprintf("Condition number %.1e is very high. LD shrinkage is important.", cond_num))
}

# 6e: Nearest positive-definite projection if R has negative eigenvalues
nearPD_applied <- FALSE
if (eig_min < 0) {
  R_pd <- as.matrix(nearPD(R, corr = TRUE, keepDiag = TRUE, do2eigen = TRUE,
                            ensureSymmetry = TRUE)$mat)
  diag(R_pd) <- 1.0
  eig_min_pd <- min(eigen(R_pd, symmetric = TRUE, only.values = TRUE)$values)
  LOGINFO(sprintf("Min eigenvalue after nearPD     : %.6f", eig_min_pd))
  R <- R_pd
  nearPD_applied <- TRUE
  LOGOK("nearPD projection applied successfully.")
} else {
  LOGOK(sprintf("R is positive semi-definite (min eig = %.4f). nearPD not needed.", eig_min))
}

# 6f: Ledoit-Wolf shrinkage toward identity
# R_shrunk = (1 - lambda) * R + lambda * I
# Guarantees strict positive definiteness, stabilises SuSiE variational updates
if (ld_shrink > 0) {
  R <- (1 - ld_shrink) * R + ld_shrink * diag(n_snps)
  LOG(sprintf("  Ledoit-Wolf shrinkage applied: lambda = %.3f", ld_shrink))
}

eig_min_post <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
LOGINFO(sprintf("Min eigenvalue (post-correction): %.6f", eig_min_post))

if (eig_min_post <= 0) {
  stop(sprintf(
    "R is not positive definite after correction (min eig = %.4e). ",
    "Increase --ld-shrink (try 0.2 or 0.3)."
  ))
}
LOGOK(sprintf("R is positive definite after all corrections (min eig = %.4f).", eig_min_post))

# --------------------------------------------------------------------------- #
# SECTION 7: Z-score diagnostics                                               #
# --------------------------------------------------------------------------- #
LOGSEP()
LOG("SECTION 7: Z-score distribution check")

z_scores  <- df$BETA / df$SE
top10_idx <- order(abs(z_scores), decreasing = TRUE)[1:min(10, n_snps)]
top_signs <- sign(z_scores[top10_idx])
n_pos_top <- sum(top_signs > 0)
n_neg_top <- sum(top_signs < 0)

LOGINFO(sprintf("N SNPs               : %d", n_snps))
LOGINFO(sprintf("|Z| range            : %.2f -- %.2f", min(abs(z_scores)), max(abs(z_scores))))
LOGINFO(sprintf("Z mean / SD          : %.4f / %.4f", mean(z_scores), sd(z_scores)))
LOGINFO(sprintf("Top-10 Z signs       : %d positive, %d negative", n_pos_top, n_neg_top))

if (n_pos_top > 0 && n_neg_top > 0) {
  LOGWARN("Mixed Z-score signs in top 10 SNPs. Causes:")
  LOGWARN("  (a) Two signals with opposite effects  (b) Allele mismatch  (c) Window too wide")
} else {
  LOGOK("Top Z-scores have consistent sign.")
}

se_med    <- median(df$SE)
n_implied <- 1 / (se_med^2 * 0.42)
LOGINFO(sprintf("Median SE            : %.4f", se_med))
LOGINFO(sprintf("Implied N (MAF~0.3)  : %.0f  (provided: %d)", n_implied, n_samples))
if (abs(log10(n_implied) - log10(n_samples)) > 0.5) {
  LOGWARN(sprintf("Implied N (~%.0f) differs from provided N (%d) by >3x.", n_implied, n_samples))
  LOGWARN("Verify --sample-size matches the GWAS N for this locus.")
} else {
  LOGOK("Sample size consistent with SE magnitudes.")
}

# --------------------------------------------------------------------------- #
# SECTION 8: Run SuSiE                                                         #
# --------------------------------------------------------------------------- #
LOGSEP()
LOG("SECTION 8: susie_rss()")
LOG(sprintf("  N=%d  L=%d  coverage=%.2f  n_snps=%d", n_samples, L, coverage, n_snps))

t0 <- proc.time()
susie_fit <- tryCatch(
  susie_rss(
    z        = z_scores,
    R        = R,
    n        = n_samples,
    L        = L,
    coverage = coverage,
    verbose  = FALSE
  ),
  error = function(e) { LOGFAIL(paste("susie_rss() error:", e$message)); stop(e) }
)
t_susie <- (proc.time() - t0)[["elapsed"]]

LOG(sprintf("  Converged : %s   Elapsed : %.1f s",
            if (isTRUE(susie_fit$converged)) "YES" else "NO", t_susie))
if (!isTRUE(susie_fit$converged))
  LOGWARN("SuSiE did not converge. Try increasing max_iter or adjusting ld_shrink.")

pips <- susie_fit$pip

# Credible sets
cs_args <- names(formals(susie_get_cs))
cs_obj  <- if ("Rr"    %in% cs_args) susie_get_cs(susie_fit, coverage = coverage, Rr    = R) else
            if ("Xcorr" %in% cs_args) susie_get_cs(susie_fit, coverage = coverage, Xcorr = R) else
            susie_get_cs(susie_fit, coverage = coverage)

cs_list   <- cs_obj$cs
n_cs      <- length(cs_list)
cs_member <- rep(NA_integer_, n_snps)
cs_cover  <- rep(NA_real_,    n_snps)

# --------------------------------------------------------------------------- #
# SECTION 9: Credible set summary                                              #
# --------------------------------------------------------------------------- #
LOGSEP()
LOG(sprintf("SECTION 9: Credible sets  (%d found)", n_cs))

if (n_cs > 5)
  LOGWARN(sprintf("%d CSs found. If window is correct and LD is raw r, run COJO first.", n_cs))

for (ci in seq_along(cs_list)) {
  idx <- cs_list[[ci]]
  cs_member[idx] <- ci
  cs_cover[idx]  <- sum(pips[idx])
  ti <- idx[which.max(pips[idx])]
  flag <- if (abs(z_scores[ti]) < 4) " [WEAK Z]" else
          if (length(idx) > 20) " [LARGE CS]" else ""
  LOG(sprintf("  CS%02d: %3d SNPs  cov=%.4f  top=%s  PIP=%.4f  BP=%d  Z=%.2f  P=%.2e%s",
              ci, length(idx), cs_cover[idx], df$SNP[ti],
              pips[ti], df$BP[ti], z_scores[ti], df$P[ti], flag))
}

idx_in_cs <- cs_member[which(df$SNP == index_snp)]
idx_in_cs <- if (length(idx_in_cs) == 0) NA_integer_ else idx_in_cs[1]

if (is.na(idx_in_cs)) {
  LOGWARN(sprintf("Index SNP %s NOT in any credible set.", index_snp))
  idx_pip <- pips[df$SNP == index_snp]
  if (length(idx_pip) > 0) LOGINFO(sprintf("  Index SNP PIP = %.4f", idx_pip[1]))
} else {
  LOGOK(sprintf("Index SNP %s is in CS%d.", index_snp, idx_in_cs))
}

# --------------------------------------------------------------------------- #
# SECTION 10: Output table                                                     #
# --------------------------------------------------------------------------- #
result_df <- data.frame(
  SNP            = df$SNP,
  BP             = df$BP,
  P              = df$P,
  BETA           = df$BETA,
  SE             = df$SE,
  Z              = round(z_scores, 6),
  PIP            = round(pips, 6),
  CS             = cs_member,
  CS_coverage    = round(cs_cover, 4),
  allele_flipped = flip_idx,
  stringsAsFactors = FALSE
)
result_df <- result_df[order(-result_df$PIP), ]
write.table(result_df, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
LOGINFO(sprintf("Results written: %d SNPs -> %s", nrow(result_df), out_tsv))

# --------------------------------------------------------------------------- #
# SECTION 11: Plots                                                            #
# --------------------------------------------------------------------------- #
LOGSEP()
LOG("SECTION 11: Plots")

cs_pal <- c("Not in CS" = "#9CA3AF",
  setNames(c("#D73027","#FC8D59","#4575B4","#1A9850","#F46D43",
             "#74ADD1","#A50026","#313695","#FDAE61","#ABD9E9"),
           paste0("CS", seq_len(min(n_cs, 10)))))

pd <- result_df[order(result_df$BP), ]
pd$cs_label <- ifelse(is.na(pd$CS), "Not in CS", paste0("CS", pd$CS))
pd$logp     <- -log10(pd$P)
pd$zsign    <- ifelse(pd$Z > 0, "Positive Z", "Negative Z")

window_actual_kb <- (max(df$BP) - min(df$BP)) / 1000

th <- theme_bw(base_size = 10) +
  theme(legend.key.size  = unit(0.35, "cm"),
        plot.title       = element_text(size = 9, face = "bold"),
        plot.subtitle    = element_text(size = 8, colour = "grey40"))

pA <- ggplot(pd, aes(x = BP / 1e6, y = PIP, colour = cs_label)) +
  geom_point(aes(size = ifelse(!is.na(CS), 2.5, 1.2)), alpha = 0.85, na.rm = TRUE) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  scale_colour_manual(values = cs_pal, name = "CS") + scale_size_identity() +
  scale_y_continuous(limits = c(0, 1.02), expand = c(0.01, 0)) +
  labs(x = "Position (Mb)", y = "PIP",
       title = sprintf("PIP track  |  %d CS  |  N=%d  L=%d  shrink=%.2f",
                       n_cs, n_samples, L, ld_shrink),
       subtitle = sprintf("Window: %.0f kb  |  %d SNPs  |  index: %s  |  LD: bigsnpr raw r",
                          window_actual_kb, n_snps, index_snp)) + th

pB <- ggplot(pd, aes(x = BP / 1e6, y = logp, colour = PIP)) +
  geom_point(size = 1.5, alpha = 0.85, na.rm = TRUE) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed",
             colour = "red", linewidth = 0.4) +
  scale_colour_gradient2(low = "#74add1", mid = "#fee090", high = "#d73027",
                         midpoint = 0.5, limits = c(0, 1), name = "PIP") +
  labs(x = "Position (Mb)", y = expression(-log[10](P)),
       title = sprintf("Manhattan  |  %d SNPs in +/-%.0f kb window", n_snps, window_kb),
       subtitle = "Red dashed = 5e-8. Points coloured by SuSiE PIP.") + th

pC <- ggplot(pd, aes(x = BP / 1e6, y = Z, colour = zsign, size = abs(Z))) +
  geom_point(alpha = 0.8, na.rm = TRUE) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  scale_colour_manual(values = c("Positive Z" = "#4575B4", "Negative Z" = "#D73027"),
                      name = "Sign") +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  labs(x = "Position (Mb)", y = "Z-score",
       title = "Z-score track",
       subtitle = "Mixed signs = possible 2nd signal or allele mismatch") + th

eig_df <- data.frame(rank  = seq_along(eig_vals),
                     value = sort(eig_vals, decreasing = TRUE))
pD <- ggplot(eig_df[eig_df$rank <= min(60, nrow(eig_df)), ], aes(x = rank, y = value)) +
  geom_line(colour = "#4575B4", linewidth = 0.7) +
  geom_point(size = 1.2, colour = "#4575B4") +
  geom_hline(yintercept = 0.01 * max(eig_df$value), linetype = "dashed",
             colour = "red", linewidth = 0.4) +
  labs(x = "Eigenvalue rank", y = "Eigenvalue (pre-shrink)",
       title = sprintf("LD eigenvalue spectrum  |  eff.rank=%d/%d  cond=%.0f",
                       eff_rank, n_snps, cond_num),
       subtitle = "Red = 1% threshold. Computed from bigsnpr raw r.") + th

if (!is.null(plot_overview)) {
  g <- arrangeGrob(pA, pB, pC, pD, ncol = 2,
    top = sprintf("SuSiE + bigsnpr  |  %s  |  %d SNPs  |  +/-%.0f kb  |  LD: raw r from PLINK",
                  basename(input_tsv), n_snps, window_kb))
  ggsave(plot_overview, g, width = 14, height = 10, dpi = 150)
  LOG(sprintf("  Overview : %s", plot_overview))
}
if (!is.null(plot_pip)) {
  ggsave(plot_pip, pA, width = 10, height = 4.5, dpi = 150)
  LOG(sprintf("  PIP      : %s", plot_pip))
}
if (!is.null(plot_zscore)) {
  ggsave(plot_zscore, pC, width = 10, height = 4.5, dpi = 150)
  LOG(sprintf("  Z-score  : %s", plot_zscore))
}
if (!is.null(plot_ld)) {
  idx_sub <- if (n_snps > 150) sort(sample(seq_len(n_snps), 150)) else seq_len(n_snps)
  R_sub   <- R[idx_sub, idx_sub]
  melt_df <- data.frame(
    i   = rep(seq_along(idx_sub), each  = length(idx_sub)),
    j   = rep(seq_along(idx_sub), times = length(idx_sub)),
    val = as.vector(R_sub)
  )
  pLD <- ggplot(melt_df, aes(x = j, y = i, fill = val)) +
    geom_raster() +
    scale_fill_gradient2(low = "#313695", mid = "white", high = "#D73027",
                         midpoint = 0, limits = c(-1, 1), name = "r") +
    scale_y_reverse() +
    labs(title = sprintf("LD heatmap (raw r from bigsnpr)  |  %s%d SNPs",
                         if (n_snps > 150) "subsampled 150 of " else "", n_snps),
         subtitle = sprintf("Blue = negative r, Red = positive r  |  shrinkage lambda=%.2f applied",
                            ld_shrink),
         x = "SNP index (BP order)", y = "SNP index (BP order)") +
    th + theme(axis.text = element_blank(), axis.ticks = element_blank())
  ggsave(plot_ld, pLD, width = 7, height = 6, dpi = 150)
  LOG(sprintf("  LD heatmap: %s", plot_ld))
}

# --------------------------------------------------------------------------- #
# SECTION 12: Problem summary and recommendations                              #
# --------------------------------------------------------------------------- #
LOGSEP()
LOG("SUMMARY OF PROBLEMS DETECTED:")

problems <- character(0)
if (!isTRUE(susie_fit$converged))
  problems <- c(problems, "SuSiE did not converge.")
if (n_cs > 5)
  problems <- c(problems, sprintf("%d CSs found. Run COJO first if multiple independent signals expected.", n_cs))
if (is.na(idx_in_cs))
  problems <- c(problems, sprintf("Index SNP %s NOT in any CS.", index_snp))
if (n_negative == 0)
  problems <- c(problems, "No negative LD values (unexpected for raw r from real genotypes).")
if (eig_min < 0)
  problems <- c(problems, sprintf("R had negative eigenvalues (%.4f). nearPD applied.", eig_min))
if (n_pos_top > 0 && n_neg_top > 0)
  problems <- c(problems, "Mixed Z-score signs in top 10 SNPs.")
if (abs(log10(n_implied) - log10(n_samples)) > 0.5)
  problems <- c(problems, sprintf("N mismatch: implied ~%.0f, provided %d.", n_implied, n_samples))
if (n_gwas_only > nrow(df_all[df_all$BP >= bp_lo & df_all$BP <= bp_hi, ]) * 0.25)
  problems <- c(problems, sprintf(">25%% GWAS SNPs missing from PLINK file.", n_gwas_only))
if (n_ambig > n_snps * 0.1)
  problems <- c(problems, sprintf(">10%% ambiguous palindromic SNPs (%d).", n_ambig))

if (length(problems) == 0) {
  LOGOK("No obvious problems detected.")
} else {
  for (pr in problems) LOG(sprintf("  !! %s", pr))
}

LOGSEP()
LOG("RECOMMENDATIONS:")
LOG(sprintf("  1. LD computed by bigsnpr from PLINK genotypes (raw r, not r^2)."))
LOG(sprintf("     PLINK prefix used: %s", plink_prefix))
LOG(sprintf("     Verify this reference matches your GWAS ancestry (Asian = use Asian panel)."))
LOG(sprintf("  2. Window: +/-%.0f kb around %s (BP=%d) = %d SNPs.",
            window_kb, index_snp, index_bp, n_snps))
if (n_cs > 3)
  LOG("  3. Multiple CSs: confirm with GCTA-COJO, then fine-map each signal window separately.")
LOG(sprintf("  4. Sample size: provided N=%d, implied from SE ~%.0f.", n_samples, n_implied))
LOG(sprintf("  5. Allele flips applied: %d. Ambiguous SNPs: %d.", n_flip, n_ambig))
LOG(sprintf("  6. LD corrections applied: clip->symmetry%s->shrinkage(%.2f).",
            if (nearPD_applied) "->nearPD" else "", ld_shrink))
LOGSEP()

# --------------------------------------------------------------------------- #
# Write JSON diagnostics                                                       #
# --------------------------------------------------------------------------- #
cs_summary <- lapply(seq_along(cs_list), function(ci) {
  idx <- cs_list[[ci]]; ti <- idx[which.max(pips[idx])]
  list(cs = ci, n_snps = length(idx), coverage = round(sum(pips[idx]), 4),
       top_snp = df$SNP[ti], top_pip = round(pips[ti], 4),
       top_bp = df$BP[ti], top_z = round(z_scores[ti], 4),
       top_p = df$P[ti], snps = df$SNP[idx])
})

writeLines(toJSON(list(
  method               = "susie_rss_bigsnpr",
  susieR_version       = as.character(packageVersion("susieR")),
  bigsnpr_version      = as.character(packageVersion("bigsnpr")),
  timestamp            = format(Sys.time()),
  plink_prefix         = plink_prefix,
  n_samples            = n_samples, L = L, coverage = coverage,
  ld_shrink            = ld_shrink, window_kb = window_kb, min_maf = min_maf,
  n_snps_gwas_total    = nrow(df_all),
  n_snps_window        = n_snps,
  n_snps_maf_removed   = sum(low_maf),
  n_allele_flips       = n_flip,
  n_ambiguous_snps     = n_ambig,
  index_snp            = index_snp, index_bp = index_bp, index_p = index_p,
  index_snp_in_cs      = !is.na(idx_in_cs),
  converged            = isTRUE(susie_fit$converged),
  n_credible_sets      = n_cs,
  n_snps_in_any_cs     = sum(!is.na(cs_member)),
  lead_pip_snp         = df$SNP[which.max(pips)],
  lead_pip             = round(max(pips), 4),
  eig_min_pre          = round(eig_min, 6),
  eig_min_post         = round(eig_min_post, 6),
  nearPD_applied       = nearPD_applied,
  condition_number     = round(cond_num, 1),
  effective_rank       = eff_rank,
  ld_has_negatives     = n_negative > 0,
  n_gwas_only          = n_gwas_only,
  z_mixed_signs        = n_pos_top > 0 && n_neg_top > 0,
  implied_n            = round(n_implied, 0),
  ld_computation_sec   = round(t_ld_elapsed, 1),
  problems_detected    = problems,
  credible_sets        = cs_summary
), pretty = TRUE, auto_unbox = TRUE), con = diag_json)

writeLines(log_lines, con = diag_log)
writeLines("ok", con = done_file)
cat(sprintf("\nDone. %d CS found. Lead PIP SNP: %s (PIP=%.4f). Log: %s\n",
            n_cs, df$SNP[which.max(pips)], max(pips), diag_log))
