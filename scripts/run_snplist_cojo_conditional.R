#!/usr/bin/env Rscript
# SNP-list conditional analysis using GCTA-COJO --cojo-cond
# For each SNP in the user-provided list that falls within the locus, condition
# the locus summary stats on that SNP and measure impact.

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(openxlsx)
  library(dplyr)
})

option_list <- list(
  make_option("--input-tsv",    type = "character", help = "Matched GWAS TSV (SNP CHR BP A1 A2 P BETA SE REF_MAF)"),
  make_option("--bfile-prefix", type = "character", help = "PLINK bfile prefix for reference panel"),
  make_option("--snp-list",     type = "character", help = "One SNP ID per line"),
  make_option("--loci-file",    type = "character", help = "Loci file (loci chr start end)"),
  make_option("--locus",        type = "character", help = "Locus name (matches loci file column 1)"),
  make_option("--out-dir",      type = "character", help = "Output directory"),
  make_option("--sample-size",  type = "integer",   help = "GWAS sample size"),
  make_option("--gcta-bin",     type = "character", default = "gcta", help = "Path/name of GCTA binary"),
  make_option("--done-file",    type = "character", help = "Marker file written on success")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ── Validate required args ────────────────────────────────────────────────────
required_args <- c("input-tsv", "bfile-prefix", "snp-list", "loci-file",
                   "locus", "out-dir", "sample-size", "done-file")
for (a in required_args) {
  if (is.null(opt[[a]])) stop(paste("Missing required argument: --", a, sep = ""))
}

dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)
plots_dir <- file.path(opt$`out-dir`, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) message("[snplist_cojo] ", ...)

# ── Read matched GWAS ─────────────────────────────────────────────────────────
gwas <- read.table(opt$`input-tsv`, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
gwas$P    <- as.numeric(gwas$P)
gwas$BETA <- as.numeric(gwas$BETA)
gwas$SE   <- as.numeric(gwas$SE)
gwas$BP   <- as.integer(gwas$BP)
gwas$CHR  <- as.integer(gwas$CHR)
gwas <- gwas[!is.na(gwas$P) & !is.na(gwas$BETA) & !is.na(gwas$SE) & gwas$SE > 0, ]

# ── Get locus boundaries ──────────────────────────────────────────────────────
loci_df <- read.table(opt$`loci-file`, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(loci_df)[1] <- "loci"
locus_row <- loci_df[loci_df$loci == opt$locus, ]
if (nrow(locus_row) == 0) stop(paste("Locus not found in loci file:", opt$locus))

locus_chr   <- as.integer(locus_row$chr[1])
locus_start <- as.integer(locus_row$start[1])
locus_end   <- as.integer(locus_row$end[1])
log_msg("Locus: ", opt$locus, "  chr", locus_chr, ":", locus_start, "-", locus_end)

# ── Restrict GWAS to locus ────────────────────────────────────────────────────
gwas_locus <- gwas[gwas$CHR == locus_chr &
                   gwas$BP  >= locus_start &
                   gwas$BP  <= locus_end, ]
gwas_locus <- gwas_locus[order(gwas_locus$BP), ]
log_msg("GWAS SNPs in locus: ", nrow(gwas_locus))

# ── Load SNP list and filter to locus SNPs ────────────────────────────────────
snp_list_all <- tryCatch(
  read.table(opt$`snp-list`, header = FALSE, stringsAsFactors = FALSE)[, 1],
  error = function(e) character(0)
)
cond_snps <- intersect(snp_list_all, gwas_locus$SNP)
log_msg("Conditioning SNPs found in locus: ", length(cond_snps), " / ", length(snp_list_all), " total")

if (length(cond_snps) == 0) {
  log_msg("No SNPs from snp_list within locus. Writing empty outputs and quitting.")
  writeLines(c(
    paste0("No SNPs from ", opt$`snp-list`, " were found within locus ", opt$locus),
    paste0("Locus bounds: chr", locus_chr, ":", locus_start, "-", locus_end),
    paste0("SNPs checked: ", length(snp_list_all))
  ), file.path(opt$`out-dir`, "no_snps_found.txt"))
  writeLines("ok", opt$`done-file`)
  quit(status = 0)
}

# ── Build GCTA .ma file ───────────────────────────────────────────────────────
freq_col <- if ("REF_MAF" %in% colnames(gwas)) gwas$REF_MAF else rep(0.5, nrow(gwas))
ma_df <- data.frame(
  SNP  = gwas$SNP,
  A1   = gwas$A1,
  A2   = gwas$A2,
  freq = as.numeric(freq_col),
  b    = gwas$BETA,
  se   = gwas$SE,
  p    = gwas$P,
  N    = as.integer(opt$`sample-size`),
  stringsAsFactors = FALSE
)
ma_df <- ma_df[!is.na(ma_df$p) & !is.na(ma_df$b) & !is.na(ma_df$se) & ma_df$se > 0, ]
ma_file <- file.path(opt$`out-dir`, "gwas.ma")
write.table(ma_df, ma_file, quote = FALSE, row.names = FALSE, sep = " ")
log_msg("Wrote .ma file: ", nrow(ma_df), " SNPs")

# ── AUC helper (area above genome-wide threshold) ─────────────────────────────
GW_THRESHOLD <- -log10(5e-8)   # ≈ 7.30

calc_auc <- function(pvals) {
  logp <- -log10(pmax(as.numeric(pvals), 1e-300))
  sum(pmax(0, logp - GW_THRESHOLD), na.rm = TRUE)
}

auc_before <- calc_auc(gwas_locus$P)

# ── Manhattan plot helper ─────────────────────────────────────────────────────
plot_manhattan <- function(df, title, out_file, highlight_snps = NULL) {
  df <- df[order(df$BP), ]
  df$logP <- -log10(pmax(as.numeric(df$P), 1e-300))
  df$category <- "ns"
  df$category[df$logP > GW_THRESHOLD] <- "sig"
  if (!is.null(highlight_snps)) {
    df$category[df$SNP %in% highlight_snps] <- "highlight"
  }

  top_idx <- which.max(df$logP)
  top_row <- df[top_idx, , drop = FALSE]
  x_span <- max(df$BP, na.rm = TRUE) - min(df$BP, na.rm = TRUE)
  if (!is.finite(x_span) || x_span <= 0) x_span <- 1
  y_span <- max(df$logP, na.rm = TRUE) - min(df$logP, na.rm = TRUE)
  if (!is.finite(y_span) || y_span <= 0) y_span <- 1

  label_df <- data.frame(
    SNP = as.character(top_row$SNP[1]),
    x_point_mb = as.numeric(top_row$BP[1]) / 1e6,
    y_point = as.numeric(top_row$logP[1]),
    x_label_mb = (as.numeric(top_row$BP[1]) + 0.03 * x_span) / 1e6,
    y_label = as.numeric(top_row$logP[1]) + 0.08 * y_span,
    stringsAsFactors = FALSE
  )

  p <- ggplot(df, aes(x = BP / 1e6, y = logP, color = category)) +
    geom_point(size = 1.8, alpha = 0.85) +
    geom_hline(yintercept = GW_THRESHOLD, linetype = "dashed", color = "red", linewidth = 0.6) +
    scale_color_manual(
      values  = c(ns = "grey60", sig = "steelblue", highlight = "darkorange"),
      labels  = c(ns = "Not significant", sig = "Significant", highlight = "SNP list"),
      name    = "",
      drop    = FALSE
    ) +
    scale_x_continuous(labels = function(x) sprintf("%.2f", x)) +
    labs(
      title = title,
      x     = paste0("Position (Mb, chr", locus_chr, ")"),
      y     = expression(-log[10](italic(P)))
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")

  p <- p +
    geom_segment(
      data = label_df,
      aes(x = x_label_mb, y = y_label, xend = x_point_mb, yend = y_point),
      inherit.aes = FALSE,
      linewidth = 0.4,
      color = "black"
    ) +
    geom_label(
      data = label_df,
      aes(x = x_label_mb, y = y_label, label = SNP),
      inherit.aes = FALSE,
      size = 3,
      label.size = 0.25,
      fill = "white",
      color = "black"
    )

  ggsave(out_file, plot = p, width = 10, height = 4, dpi = 150)
}

# Unconditional locus plot
plot_manhattan(
  gwas_locus,
  title       = paste0(opt$locus, "  —  unconditional"),
  out_file    = file.path(plots_dir, "locus_unconditional.png"),
  highlight_snps = cond_snps
)
log_msg("Unconditional plot done")

# ── Run GCTA cojo-cond per SNP ────────────────────────────────────────────────
run_cojo_cond <- function(snp_id, out_prefix) {
  cond_file <- paste0(out_prefix, ".condsnp")
  writeLines(snp_id, cond_file)
  cmd <- sprintf(
    "%s --bfile %s --cojo-file %s --cojo-cond %s --out %s --thread-num 1 2>&1",
    opt$`gcta-bin`,
    opt$`bfile-prefix`,
    ma_file,
    cond_file,
    out_prefix
  )
  ret <- system(cmd)
  cojo_out <- paste0(out_prefix, ".cma.cojo")
  if (!file.exists(cojo_out)) {
    log_msg("WARNING: no .cma.cojo output for SNP ", snp_id, " (exit code ", ret, ")")
    return(NULL)
  }
  res <- tryCatch(
    read.table(cojo_out, header = TRUE, stringsAsFactors = FALSE),
    error = function(e) { log_msg("Parse error for ", snp_id, ": ", e$message); NULL }
  )
  return(res)
}

cond_results <- list()

for (snp_id in cond_snps) {
  log_msg("Conditioning on: ", snp_id)
  snp_safe   <- gsub("[^A-Za-z0-9_.-]", "_", snp_id)
  out_prefix <- file.path(opt$`out-dir`, paste0("gcta_cond_", snp_safe))
  res <- run_cojo_cond(snp_id, out_prefix)
  if (is.null(res)) next
  # GCTA cma.cojo columns: Chr SNP bp refA freq b se p n freq_geno bJ bJ_se pJ LD_r
  if (!"pC" %in% colnames(res)) {
    log_msg("No pC column for SNP: ", snp_id, " (cols: ", paste(colnames(res), collapse=","), ")")
    next
  }
  res$pC    <- as.numeric(res$pC)
  res$bC    <- as.numeric(res$bC)
  res$bC_se <- as.numeric(res$bC_se)

  # Merge conditional p/beta back to locus GWAS
  cond_df <- merge(gwas_locus, res[, c("SNP", "pC", "bC", "bC_se"), drop = FALSE],
                   by = "SNP", all.x = TRUE)
  cond_df <- cond_df[order(cond_df$BP), ]
  # SNPs where pC is NA are collinear with the conditioning SNP (r²>0.9 cutoff)
  # GCTA cannot estimate a conditional effect — treat them as non-significant (P=1)
  cond_df$pC[is.na(cond_df$pC)] <- 1.0
  cond_df$bC[is.na(cond_df$bC)] <- 0.0
  cond_df$bC_se[is.na(cond_df$bC_se)] <- NA_real_
  cond_results[[snp_id]] <- cond_df

  # Conditional Manhattan (use pJ where available, fallback to original P)
  cond_plot_df         <- cond_df
  cond_plot_df$P       <- cond_df$pC   # NA already replaced with 1 above
  plot_manhattan(
    cond_plot_df,
    title       = paste0(opt$locus, "  —  conditioned on ", snp_id),
    out_file    = file.path(plots_dir, paste0("locus_conditional_", snp_safe, ".png")),
    highlight_snps = snp_id
  )
}

log_msg("Conditional runs complete (", length(cond_results), " succeeded)")

# ── Build SNP table ───────────────────────────────────────────────────────────
snp_table <- gwas_locus[, c("SNP", "CHR", "BP", "A1", "A2", "P", "BETA", "SE")]
colnames(snp_table)[colnames(snp_table) == "P"]    <- "P_unconditional"
colnames(snp_table)[colnames(snp_table) == "BETA"] <- "BETA_unconditional"
colnames(snp_table)[colnames(snp_table) == "SE"]   <- "SE_unconditional"

for (snp_id in names(cond_results)) {
  snp_safe <- gsub("[^A-Za-z0-9_.-]", "_", snp_id)
  cr <- cond_results[[snp_id]]
  merge_cols <- cr[, c("SNP", "pC", "bC", "bC_se"), drop = FALSE]
  colnames(merge_cols) <- c("SNP",
                             paste0("P_cond_",    snp_safe),
                             paste0("BETA_cond_", snp_safe),
                             paste0("SE_cond_",   snp_safe))
  snp_table <- merge(snp_table, merge_cols, by = "SNP", all.x = TRUE)
}
snp_table <- snp_table[order(snp_table$P_unconditional), ]

# ── Impact summary table ──────────────────────────────────────────────────────
impact_rows <- lapply(names(cond_results), function(snp_id) {
  cr       <- cond_results[[snp_id]]
  p_after  <- ifelse(is.na(cr$pC), cr$P, cr$pC)
  auc_after <- calc_auc(p_after)
  pct_dec   <- if (auc_before > 0) (auc_before - auc_after) / auc_before * 100 else 0
  data.frame(
    conditioning_SNP  = snp_id,
    AUC_before        = round(auc_before, 4),
    AUC_after         = round(auc_after, 4),
    pct_AUC_decrease  = round(pct_dec, 2),
    n_sig_before      = sum(cr$P    < 5e-8, na.rm = TRUE),
    n_sig_after       = sum(p_after < 5e-8, na.rm = TRUE),
    stringsAsFactors  = FALSE
  )
})

if (length(impact_rows) > 0) {
  impact_df <- do.call(rbind, impact_rows)
  impact_df <- impact_df[order(-impact_df$pct_AUC_decrease), ]
} else {
  impact_df <- data.frame(
    conditioning_SNP = character(), AUC_before = numeric(), AUC_after = numeric(),
    pct_AUC_decrease = numeric(), n_sig_before = integer(), n_sig_after = integer()
  )
}

# ── Impact bar plot ───────────────────────────────────────────────────────────
if (nrow(impact_df) > 0) {
  impact_df$conditioning_SNP <- factor(
    impact_df$conditioning_SNP,
    levels = impact_df$conditioning_SNP[order(impact_df$pct_AUC_decrease)]
  )
  p_impact <- ggplot(impact_df, aes(x = conditioning_SNP,
                                     y = pct_AUC_decrease,
                                     fill = pct_AUC_decrease)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_gradient(low = "steelblue", high = "firebrick", name = "% decrease") +
    labs(
      title = paste0(opt$locus, "  —  conditioning impact (% AUC decrease)"),
      x     = "Conditioning SNP",
      y     = "% decrease in area above genome-wide threshold"
    ) +
    theme_bw(base_size = 12)
  plot_h <- max(4, nrow(impact_df) * 0.45 + 2)
  ggsave(file.path(plots_dir, "impact_plot.png"),
         plot = p_impact, width = 8, height = plot_h, dpi = 150)
  log_msg("Impact plot saved")
}

# ── Write Excel workbook ──────────────────────────────────────────────────────
header_style <- createStyle(
  fontColour = "#FFFFFF", fgFill = "#2C5F8A",
  halign = "CENTER", textDecoration = "Bold", border = "Bottom"
)

wb <- createWorkbook()

## Sheet 1: p-values only (placed before SNP_Table)
pval_cols <- grep("^P_", colnames(snp_table), value = TRUE)
pval_table <- if (length(pval_cols) > 0) {
  snp_table[, c("SNP", pval_cols), drop = FALSE]
} else {
  snp_table[, "SNP", drop = FALSE]
}

addWorksheet(wb, "PValues_Only")
writeData(wb, "PValues_Only", pval_table)
addStyle(wb, "PValues_Only", header_style,
         rows = 1, cols = seq_len(ncol(pval_table)), gridExpand = TRUE)
setColWidths(wb, "PValues_Only", cols = seq_len(ncol(pval_table)), widths = "auto")

if (nrow(pval_table) > 0 && ncol(pval_table) > 1) {
  pval_num_style <- createStyle(numFmt = "0.00E+00")
  sig_p_style <- createStyle(fgFill = "#FDE2E1")

  # Apply scientific notation with 2 decimals to all p-value cells.
  addStyle(
    wb,
    "PValues_Only",
    pval_num_style,
    rows = 2:(nrow(pval_table) + 1),
    cols = 2:ncol(pval_table),
    gridExpand = TRUE,
    stack = TRUE
  )

  # Highlight significant p-values (p < 5e-8) in light red.
  for (col_idx in 2:ncol(pval_table)) {
    vals <- suppressWarnings(as.numeric(pval_table[[col_idx]]))
    hit_rows <- which(!is.na(vals) & vals < 5e-8)
    if (length(hit_rows) > 0) {
      addStyle(
        wb,
        "PValues_Only",
        sig_p_style,
        rows = hit_rows + 1,
        cols = col_idx,
        gridExpand = FALSE,
        stack = TRUE
      )
    }
  }
}

addWorksheet(wb, "SNP_Table")
writeData(wb, "SNP_Table", snp_table)
addStyle(wb, "SNP_Table", header_style,
         rows = 1, cols = seq_len(ncol(snp_table)), gridExpand = TRUE)
setColWidths(wb, "SNP_Table", cols = seq_len(ncol(snp_table)), widths = "auto")

addWorksheet(wb, "Impact_Summary")
writeData(wb, "Impact_Summary", impact_df)
addStyle(wb, "Impact_Summary", header_style,
         rows = 1, cols = seq_len(ncol(impact_df)), gridExpand = TRUE)
setColWidths(wb, "Impact_Summary", cols = seq_len(ncol(impact_df)), widths = "auto")

out_xlsx <- file.path(opt$`out-dir`, "snplist_conditional_results.xlsx")
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
log_msg("Excel written: ", out_xlsx)

writeLines("ok", opt$`done-file`)
log_msg("Finished successfully — locus: ", opt$locus)
