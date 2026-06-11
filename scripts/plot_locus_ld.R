#!/usr/bin/env Rscript
# plot_locus_ld.R
#
# Produces a two-panel PDF / PNG for a GWAS locus:
#   (top)    Locus-zoom scatter  –log10(P) vs position, points coloured by R²
#            with the lead variant.
#   (bottom) LD R² triangle heatmap (Haploview-style rotated 45°) for the
#            same SNP set, x-axis aligned with the GWAS panel.
#
# Only SNPs within --window-bp of the lead variant are shown.
# LD matrix values are interpreted as Pearson R (not R²); they are squared
# internally.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--matched-tsv", type = "character",
              help = "Path to matched.tsv (GWAS, tab-delimited)"),
  make_option("--ld-matrix",   type = "character",
              help = "Path to LD matrix (R values, first col = row SNP IDs)"),
  make_option("--out-pdf",     type = "character", default = NULL,
              help = "Output PDF path [optional]"),
  make_option("--out-png",     type = "character", default = NULL,
              help = "Output PNG path [optional]"),
  make_option("--window-bp",   type = "integer",   default = 100000L,
              help = "Half-window in bp around lead SNP [default: %default]"),
  make_option("--done-file",   type = "character", default = NULL,
              help = "Sentinel done file to write on success")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt[["matched-tsv"]]) || is.null(opt[["ld-matrix"]])) {
  stop("--matched-tsv and --ld-matrix are required.", call. = FALSE)
}
if (is.null(opt[["out-pdf"]]) && is.null(opt[["out-png"]])) {
  stop("At least one of --out-pdf or --out-png is required.", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Load GWAS
# ---------------------------------------------------------------------------
message("Loading matched GWAS ...")
gwas <- fread(opt[["matched-tsv"]], sep = "\t", header = TRUE)

# Normalise column names (case-insensitive)
cn        <- tolower(names(gwas))
names(gwas) <- cn
# Expected header: SNP CHR BP A1 A2 P BETA SE REF_MAF
# Rename if needed
col_map <- c(snp = "snp", chr = "chr", bp = "bp",
             a1 = "a1",   a2 = "a2",   p  = "p",
             beta = "beta", se = "se",  ref_maf = "ref_maf")
for (new_nm in names(col_map)) {
  old_nm <- col_map[[new_nm]]
  if (!new_nm %in% names(gwas) && old_nm %in% names(gwas))
    setnames(gwas, old_nm, new_nm)
}

required_cols <- c("snp", "chr", "bp", "p")
missing_cols  <- setdiff(required_cols, names(gwas))
if (length(missing_cols))
  stop("matched.tsv missing columns: ", paste(missing_cols, collapse = ", "))

gwas$bp <- as.integer(gwas$bp)
gwas$p  <- as.numeric(gwas$p)

# ---------------------------------------------------------------------------
# Load LD matrix  (R values)
# ---------------------------------------------------------------------------
message("Loading LD matrix ...")
ld_raw  <- fread(opt[["ld-matrix"]], sep = "\t", header = TRUE)
ld_snps <- ld_raw[[1]]          # row IDs (first column)
ld_cols <- names(ld_raw)[-1]    # column IDs

ld_mat           <- as.matrix(ld_raw[, -1, with = FALSE])
rownames(ld_mat) <- ld_snps
colnames(ld_mat) <- ld_cols

r2_mat <- ld_mat ^ 2           # convert R → R²

# ---------------------------------------------------------------------------
# Identify lead SNP and window
# ---------------------------------------------------------------------------
lead_idx <- which.min(gwas$p)
lead_snp <- gwas$snp[lead_idx]
lead_bp  <- gwas$bp[lead_idx]
chr_val  <- gwas$chr[lead_idx]
window   <- opt[["window-bp"]]

message(sprintf("Lead SNP : %s  |  chr%s : %d  |  P = %g",
                lead_snp, chr_val, lead_bp, gwas$p[lead_idx]))
message(sprintf("Window   : ± %d bp  (%.3f Mb)", window, window / 1e6))

# ---------------------------------------------------------------------------
# Filter GWAS to window, intersect with LD matrix
# ---------------------------------------------------------------------------
gwas_win <- gwas[abs(gwas$bp - lead_bp) <= window, ]
keep     <- intersect(gwas_win$snp, rownames(r2_mat))

if (length(keep) < 2L)
  stop(sprintf(
    "Only %d SNP(s) overlap between GWAS window and LD matrix. Increase --window-bp.",
    length(keep)))

gwas_win <- gwas_win[gwas_win$snp %in% keep, ]
gwas_win <- gwas_win[order(gwas_win$bp), ]
gwas_win$idx <- seq_len(nrow(gwas_win))   # 1-based index sorted by bp
n_snps <- nrow(gwas_win)

message(sprintf("SNPs in window ∩ LD matrix: %d", n_snps))

# Subset R² matrix to window SNPs, in bp order
snp_ord <- gwas_win$snp
r2_sub  <- r2_mat[snp_ord, snp_ord, drop = FALSE]

# ---------------------------------------------------------------------------
# R² of each SNP with lead variant  (for GWAS point colouring)
# ---------------------------------------------------------------------------
if (lead_snp %in% rownames(r2_sub)) {
  r2_lead <- as.numeric(r2_sub[lead_snp, snp_ord])
} else {
  warning("Lead SNP not found in LD matrix; colouring by distance proxy.")
  r2_lead <- rep(NA_real_, n_snps)
}
gwas_win$r2_lead <- r2_lead
gwas_win$r2_lead[gwas_win$snp == lead_snp] <- 1.0   # force self = 1

# ---------------------------------------------------------------------------
# Shared Haploview-style colour scale
# ---------------------------------------------------------------------------
ld_pal    <- c("#FFFFFF", "#4169E1", "#00BFFF", "#00CC44",
               "#FFDD00", "#FF6600", "#CC0000")
ld_values <- rescale(c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0))

gwas_colour_scale <- scale_color_gradientn(
  colours  = ld_pal,
  values   = ld_values,
  limits   = c(0, 1),
  na.value = "grey60",
  name     = bquote(italic(r)^2 ~ "with lead SNP"),
  guide    = guide_colorbar(barwidth = 0.7, barheight = 5)
)

ld_fill_scale <- scale_fill_gradientn(
  colours  = ld_pal,
  values   = ld_values,
  limits   = c(0, 1),
  name     = bquote(italic(r)^2),
  guide    = guide_colorbar(barwidth = 0.7, barheight = 5)
)

# ---------------------------------------------------------------------------
# Shared x-axis ticks  (index space → Mb labels)
# ---------------------------------------------------------------------------
tick_bp  <- pretty(gwas_win$bp, n = 8)
# Map each tick bp to nearest SNP index
tick_idx <- sapply(tick_bp, function(b)
  gwas_win$idx[which.min(abs(gwas_win$bp - b))])
tick_lbl <- sprintf("%.3f", tick_bp / 1e6)

x_scale_top <- scale_x_continuous(
  limits = c(0.5, n_snps + 0.5),
  breaks = tick_idx,
  labels = tick_lbl,
  expand = c(0, 0)
)
x_scale_bot <- scale_x_continuous(
  limits = c(0.5, n_snps + 0.5),
  breaks = tick_idx,
  labels = tick_lbl,
  expand = c(0, 0),
  name   = sprintf("Position on chr%s (Mb)", chr_val)
)

# ---------------------------------------------------------------------------
# GWAS panel  (top)
# ---------------------------------------------------------------------------
gwas_win$logp   <- -log10(gwas_win$p)
is_lead_row     <- gwas_win$snp == lead_snp

theme_gwas <- theme_classic(base_size = 11) +
  theme(
    axis.title.x    = element_blank(),
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    plot.margin     = margin(6, 6, 0, 6, unit = "pt"),
    legend.position = "right"
  )

p_gwas <- ggplot(gwas_win, aes(x = idx, y = logp, colour = r2_lead)) +
  # Significance lines
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed",
             colour = "red",  linewidth = 0.4) +
  geom_hline(yintercept = -log10(1e-5), linetype = "dotted",
             colour = "blue", linewidth = 0.3) +
  # All non-lead SNPs
  geom_point(data  = gwas_win[!is_lead_row, ],
             size  = 1.8, alpha = 0.85) +
  # Lead SNP (diamond, purple)
  geom_point(data  = gwas_win[ is_lead_row, ],
             size  = 4, shape = 18, colour = "purple") +
  # Lead SNP label
  annotate("text",
           x     = gwas_win$idx[is_lead_row],
           y     = gwas_win$logp[is_lead_row],
           label = lead_snp,
           vjust = -0.9, hjust = 0.5,
           size  = 2.8, colour = "purple") +
  gwas_colour_scale +
  x_scale_top +
  scale_y_continuous(
    name   = expression(-log[10](italic(P))),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  theme_gwas

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Helper: build Haploview-style rotated-45 triangle polygon data frame
#
# Correct geometry:
#   xc = (i + j) / 2   -- midpoint, aligns with SNP index in Manhattan above
#   yc = (j - i) / 2   -- POSITIVE, diagonal yc=0, apex yc=(n-1)/2
#
# scale_y_reverse(limits=c(max_y, -0.5)):
#   renders max_y (apex) at BOTTOM and -0.5 (above diagonal) at TOP
#   => flat edge at top, apex pointing down = correct Haploview orientation
# ---------------------------------------------------------------------------
make_ld_triangle <- function(r2_matrix) {
  n_s      <- nrow(r2_matrix)
  pair_idx <- which(row(r2_matrix) <= col(r2_matrix), arr.ind = TRUE)
  idx_i    <- pair_idx[, 1]
  idx_j    <- pair_idx[, 2]
  r2_v     <- r2_matrix[pair_idx]
  n_p      <- length(r2_v)
  xc <- (idx_i + idx_j) / 2
  yc <- (idx_j - idx_i) / 2
  data.frame(
    group = rep(seq_len(n_p), each = 4L),
    r2    = rep(r2_v, each = 4L),
    x     = as.vector(rbind(xc - 0.5, xc,        xc + 0.5, xc       )),
    y     = as.vector(rbind(yc,        yc - 0.5,  yc,        yc + 0.5))
  )
}

theme_ld <- theme_classic(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    plot.margin     = margin(0, 6, 6, 6, unit = "pt"),
    legend.position = "right",
    axis.title.y    = element_blank(),
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.line.y     = element_blank()
  )

# ---------------------------------------------------------------------------
# PLOT 1 -- Full pairwise LD triangle (Haploview style, all SNP pairs)
# ---------------------------------------------------------------------------
message(sprintf("Building full LD triangle (%d x %d) ...", n_snps, n_snps))

poly_full  <- make_ld_triangle(r2_sub)
max_y_full <- (n_snps - 1) / 2 + 0.5

p_ld_full <- ggplot(poly_full, aes(x = x, y = y, group = group, fill = r2)) +
  geom_polygon(colour = NA) +
  ld_fill_scale +
  x_scale_bot +
  scale_y_reverse(limits = c(max_y_full, -0.5), expand = c(0, 0)) +
  labs(subtitle = sprintf(
    "Full pairwise LD heatmap  |  %d SNPs  |  chr%s: %s-%s Mb",
    n_snps, chr_val,
    formatC(min(gwas_win$bp) / 1e6, digits = 3, format = "f"),
    formatC(max(gwas_win$bp) / 1e6, digits = 3, format = "f")
  )) +
  theme_ld +
  theme(plot.subtitle = element_text(size = 8, colour = "grey30"))

# ---------------------------------------------------------------------------
# PLOT 2 -- Index SNP r2 bar chart
#
# One bar per SNP showing r2 with the lead/index SNP only.
# Bars coloured by the same Haploview palette.
# x-axis aligned with the Manhattan panel above.
# The lead SNP bar is highlighted with a purple fill and black border.
# This shows LD decay around the index SNP directly and unambiguously.
# ---------------------------------------------------------------------------
message(sprintf("Building index-SNP r2 bar chart (index = %s) ...", lead_snp))

bar_df <- data.frame(
  idx     = gwas_win$idx,
  r2_lead = gwas_win$r2_lead,
  is_lead = gwas_win$snp == lead_snp
)

lead_bp_val <- gwas_win$bp[gwas_win$snp == lead_snp]
lead_p_val  <- gwas_win$p[gwas_win$snp  == lead_snp]

p_ld_index <- ggplot(bar_df, aes(x = idx, y = r2_lead, fill = r2_lead)) +
  geom_col(width = 0.8, alpha = 0.9) +
  geom_col(data = bar_df[bar_df$is_lead, ],
           aes(x = idx, y = r2_lead),
           fill = "purple", colour = "black", linewidth = 0.6, width = 0.8) +
  ld_fill_scale +
  x_scale_bot +
  scale_y_continuous(
    limits = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    expand = c(0, 0), name = bquote(italic(r)^2 ~ "vs index SNP")
  ) +
  labs(subtitle = sprintf(
    "r2 vs index SNP: %s  |  chr%s:%d  |  P = %s",
    lead_snp, chr_val, lead_bp_val,
    formatC(lead_p_val, format = "e", digits = 2)
  )) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    plot.margin     = margin(4, 6, 6, 6, unit = "pt"),
    legend.position = "none",
    plot.subtitle   = element_text(size = 8, colour = "grey30")
  )

# ---------------------------------------------------------------------------
# Shared title / subtitle
# ---------------------------------------------------------------------------
locus_title    <- sprintf("Locus zoom  |  Lead SNP: %s  |  Window: +/-%g kb",
                          lead_snp, window / 1e3)
locus_subtitle <- sprintf("chr%s: %s - %s Mb",
                           chr_val,
                           formatC(min(gwas_win$bp) / 1e6, digits = 3, format = "f"),
                           formatC(max(gwas_win$bp) / 1e6, digits = 3, format = "f"))
title_theme <- theme(
  plot.title    = element_text(size = 12, face = "bold"),
  plot.subtitle = element_text(size = 10, colour = "grey20")
)

# ---------------------------------------------------------------------------
# Assemble Plot 1: Manhattan + full LD triangle
# ---------------------------------------------------------------------------
combined_full <- (p_gwas / p_ld_full) +
  plot_layout(heights = c(2, 1.5), guides = "keep") +
  plot_annotation(
    title    = paste0(locus_title, " | Plot 1: Full pairwise LD triangle"),
    subtitle = locus_subtitle,
    theme    = title_theme
  )

# ---------------------------------------------------------------------------
# Assemble Plot 2: Manhattan + index SNP r2 bar chart
# ---------------------------------------------------------------------------
combined_index <- (p_gwas / p_ld_index) +
  plot_layout(heights = c(2, 0.8), guides = "keep") +
  plot_annotation(
    title    = paste0(locus_title, " | Plot 2: r2 vs index SNP only"),
    subtitle = locus_subtitle,
    theme    = title_theme
  )

# ---------------------------------------------------------------------------
# Save -- Plot 1 uses original paths; Plot 2 inserts "_index" before extension
# ---------------------------------------------------------------------------
save_plot <- function(p, out_path, w = 11, h = 8) {
  if (is.null(out_path) || nchar(out_path) == 0) return(invisible(NULL))
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_path, p, width = w, height = h, units = "in", dpi = 150)
  message("Saved: ", out_path)
}

add_suffix <- function(path, suffix) {
  if (is.null(path) || nchar(path) == 0) return(NULL)
  sub("(\\.[^.]+)$", paste0(suffix, "\\1"), path)
}

save_plot(combined_full,  opt[["out-pdf"]])
save_plot(combined_full,  opt[["out-png"]])
save_plot(combined_index, add_suffix(opt[["out-pdf"]], "_index"))
save_plot(combined_index, add_suffix(opt[["out-png"]], "_index"))

# Done sentinel
done_file <- opt[["done-file"]]
if (!is.null(done_file) && nchar(done_file) > 0) {
  dir.create(dirname(done_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(format(Sys.time()), done_file)
}
message("Done. Two plots saved:")
message("  Plot 1 (full LD triangle) : ", opt[["out-png"]])
message("  Plot 2 (index SNP r2 bars): ", add_suffix(opt[["out-png"]], "_index"))
