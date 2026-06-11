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
# LD triangle panel  (bottom) -- Haploview style
#
# The triangle hangs BELOW the Manhattan with the flat edge (diagonal,
# self-correlations = 1) at the TOP and the apex pointing DOWN.
# Each SNP column in the triangle aligns exactly with its dot above.
#
# Geometry for pair (i, j), bp-sorted indices, i <= j:
#
#   xc =  (i + j) / 2     exact midpoint -> aligns with SNP index axis
#   yc = -(j - i) / 2     negative so triangle grows downward
#                          self-pair i=j: yc=0 (top/diagonal)
#                          most distant: yc=-(N-1)/2 (apex)
#
# FIX 1: yc sign flipped (was positive -> pointed up, now negative -> down)
# FIX 2: xc formula fixed (was (i+j-1)/2 -> shifted left 0.5; now (i+j)/2)
# FIX 3: scale_y_reverse() so y=0 is at top, apex at bottom of panel
# FIX 4: y limits include full extent of all diamonds without clipping
# ---------------------------------------------------------------------------
message(sprintf("Building LD triangle polygons for %d x %d matrix ...", n_snps, n_snps))

pair_idx <- which(row(r2_sub) <= col(r2_sub), arr.ind = TRUE)  # upper tri + diag
idx_i    <- pair_idx[, 1]
idx_j    <- pair_idx[, 2]
r2_vals  <- r2_sub[pair_idx]
n_pairs  <- length(r2_vals)

# Correct Haploview geometry:
#
# xc = (i + j) / 2      midpoint x-coordinate, aligns exactly with SNP index axis
# yc = (j - i) / 2      POSITIVE, grows upward in data space
#                        diagonal (i==j): yc=0
#                        most distant pair: yc=(n-1)/2  (apex)
#
# scale_y_reverse() then flips the panel:
#   yc=0      (diagonal) renders at the TOP    <- flat edge touching Manhattan
#   yc=max    (apex)     renders at the BOTTOM <- triangle points down
#
# BUGS FIXED:
#   - yc was negative (-( j-i)/2) which double-flipped with scale_y_reverse
#     making the triangle point UP again
#   - scale_y_reverse(limits=c(0.5, min_y)) was backwards: in ggplot2
#     scale_y_reverse limits=c(a,b) renders a at BOTTOM, b at TOP,
#     so c(0.5, negative) put the apex at top and diagonal at bottom
#   - Correct limits: c(max_y + 0.5, -0.5) puts max_y at bottom, 0 at top

xc <- (idx_i + idx_j) / 2          # exact midpoint -> aligns with SNP indices
yc <- (idx_j - idx_i) / 2          # positive, apex is largest value
max_y <- max(yc) + 0.5             # top of apex diamond

# Diamond corners: Left, Top(toward diagonal), Right, Bottom(toward apex)
poly_df <- data.frame(
  group = rep(seq_len(n_pairs), each = 4L),
  r2    = rep(r2_vals,          each = 4L),
  x     = as.vector(rbind(xc - 0.5,  xc,         xc + 0.5,  xc        )),
  y     = as.vector(rbind(yc,         yc - 0.5,   yc,         yc + 0.5 ))
)

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

p_ld <- ggplot(poly_df, aes(x = x, y = y, group = group, fill = r2)) +
  geom_polygon(colour = NA) +
  ld_fill_scale +
  x_scale_bot +
  # scale_y_reverse: limits=c(a,b) -> a rendered at BOTTOM, b at TOP
  # We want: max_y (apex) at BOTTOM, -0.5 (just above diagonal) at TOP
  scale_y_reverse(
    limits = c(max_y, -0.5),
    expand = c(0, 0)
  ) +
  theme_ld

# ---------------------------------------------------------------------------
# Combine and save
# ---------------------------------------------------------------------------
combined <- (p_gwas / p_ld) +
  plot_layout(heights = c(2, 1.5), guides = "keep") +
  plot_annotation(
    title    = sprintf("Locus: %s  |  Lead SNP: %s  |  Window: ±%g kb",
                       lead_snp, lead_snp, window / 1e3),
    subtitle = sprintf("chr%s: %s - %s Mb",
                       chr_val,
                       formatC(min(gwas_win$bp) / 1e6, digits = 3, format = "f"),
                       formatC(max(gwas_win$bp) / 1e6, digits = 3, format = "f")),
    theme    = theme(plot.title    = element_text(size = 12, face = "bold"),
                     plot.subtitle = element_text(size = 10))
  )

for (out_path in c(opt[["out-pdf"]], opt[["out-png"]])) {
  if (!is.null(out_path) && nchar(out_path) > 0) {
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    ggsave(out_path, combined, width = 11, height = 8, units = "in", dpi = 150)
    message("Saved: ", out_path)
  }
}

# Done sentinel
done_file <- opt[["done-file"]]
if (!is.null(done_file) && nchar(done_file) > 0) {
  dir.create(dirname(done_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(format(Sys.time()), done_file)
}
message("Done.")
