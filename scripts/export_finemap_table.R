#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--matched-tsv", type = "character"),
  make_option("--afreq", type = "character"),
  make_option("--finemap-tsv", type = "character"),
  make_option("--out-tsv", type = "character"),
  make_option("--out-xlsx", type = "character"),
  make_option("--diag-json", type = "character"),
  make_option("--done-file", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

required_args <- c("matched-tsv", "afreq", "finemap-tsv", "out-tsv", "out-xlsx", "diag-json", "done-file")
for (a in required_args) {
  if (is.null(opt[[a]]) || opt[[a]] == "") {
    stop(sprintf("Missing required argument: --%s", a), call. = FALSE)
  }
}

mkdir_for <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

fmt_sci <- function(x) {
  x <- as.numeric(x)
  out <- rep("", length(x))
  ok <- is.finite(x)
  out[ok] <- sprintf("%.3e", x[ok])
  out
}

fmt_fixed <- function(x, d) {
  x <- as.numeric(x)
  out <- rep("", length(x))
  ok <- is.finite(x)
  out[ok] <- sprintf(paste0("%.", d, "f"), x[ok])
  out
}

or_ci <- function(beta, se) {
  b <- suppressWarnings(as.numeric(beta))
  s <- suppressWarnings(as.numeric(se))
  or_v <- exp(b)
  lo <- exp(b - 1.96 * s)
  hi <- exp(b + 1.96 * s)
  list(or = or_v, lo = lo, hi = hi)
}

mkdir_for(opt[["out-tsv"]])
mkdir_for(opt[["out-xlsx"]])
mkdir_for(opt[["diag-json"]])
mkdir_for(opt[["done-file"]])

matched <- fread(opt[["matched-tsv"]], sep = "\t", data.table = FALSE)
af <- fread(opt[["afreq"]], sep = "\t", data.table = FALSE)
finemap <- fread(opt[["finemap-tsv"]], sep = "\t", data.table = FALSE)

need_matched <- c("SNP", "CHR", "BP", "A1", "A2", "P", "BETA", "SE")
miss_matched <- setdiff(need_matched, colnames(matched))
if (length(miss_matched) > 0) {
  stop(sprintf("matched.tsv missing columns: %s", paste(miss_matched, collapse = ", ")), call. = FALSE)
}

need_af <- c("ID", "REF", "ALT", "ALT_FREQS")
miss_af <- setdiff(need_af, colnames(af))
if (length(miss_af) > 0) {
  stop(sprintf("afreq missing columns: %s", paste(miss_af, collapse = ", ")), call. = FALSE)
}

matched$SNP <- as.character(matched$SNP)
matched$A1 <- toupper(as.character(matched$A1))
matched$A2 <- toupper(as.character(matched$A2))

af$ID <- as.character(af$ID)
af$REF <- toupper(as.character(af$REF))
af$ALT <- toupper(as.character(af$ALT))
af$ALT_FREQS <- suppressWarnings(as.numeric(af$ALT_FREQS))

# Find FINEMAP PIP column case-insensitively.
fm_names_l <- tolower(colnames(finemap))
pip_idx <- which(fm_names_l == "pip")
if (length(pip_idx) == 0) {
  finemap$pip <- NA_real_
} else {
  finemap$pip <- suppressWarnings(as.numeric(finemap[[pip_idx[1]]]))
}

if (!("SNP" %in% colnames(finemap))) {
  stop("finemap TSV must include SNP column", call. = FALSE)
}
finemap$SNP <- as.character(finemap$SNP)

# Keep the maximum PIP per SNP if duplicates exist.
finemap_agg <- aggregate(pip ~ SNP, data = finemap[, c("SNP", "pip")], FUN = function(x) max(x, na.rm = TRUE))
finemap_agg$pip[!is.finite(finemap_agg$pip)] <- NA_real_

x <- merge(matched[, need_matched], af[, need_af], by.x = "SNP", by.y = "ID", all.x = TRUE)
x <- merge(x, finemap_agg, by = "SNP", all.x = TRUE)

x$af_ea_from_afreq <- ifelse(x$A1 == x$ALT, x$ALT_FREQS,
                             ifelse(x$A1 == x$REF, 1 - x$ALT_FREQS, NA_real_))

x$af <- x$af_ea_from_afreq
x$ea <- x$A1
x$oa <- ifelse(x$A1 == x$REF, x$ALT,
               ifelse(x$A1 == x$ALT, x$REF, x$A2))
x$allele_align_ok <- x$A1 == x$REF | x$A1 == x$ALT

or_stats <- or_ci(x$BETA, x$SE)
x$or <- or_stats$or
x$or_l95 <- or_stats$lo
x$or_u95 <- or_stats$hi

out_num <- data.frame(
  af = round(as.numeric(x$af), 3),
  snp = x$SNP,
  chr = x$CHR,
  bp = x$BP,
  ea = x$ea,
  oa = x$oa,
  p = as.numeric(x$P),
  or = round(as.numeric(x$or), 2),
  or_l95 = round(as.numeric(x$or_l95), 2),
  or_u95 = round(as.numeric(x$or_u95), 2),
  pip = as.numeric(x$pip),
  allele_align_ok = as.logical(x$allele_align_ok),
  af_ea_from_afreq = as.numeric(x$af_ea_from_afreq)
)

ord <- order(out_num$p, out_num$bp, na.last = TRUE)
out_num <- out_num[ord, , drop = FALSE]

# TSV: human-formatted strings.
out_tsv_df <- out_num
out_tsv_df$af <- fmt_fixed(out_tsv_df$af, 3)
out_tsv_df$p <- fmt_sci(out_tsv_df$p)
out_tsv_df$or <- fmt_fixed(out_tsv_df$or, 2)
out_tsv_df$or_l95 <- fmt_fixed(out_tsv_df$or_l95, 2)
out_tsv_df$or_u95 <- fmt_fixed(out_tsv_df$or_u95, 2)
out_tsv_df$pip <- fmt_fixed(out_tsv_df$pip, 4)
out_tsv_df$af_ea_from_afreq <- fmt_fixed(out_tsv_df$af_ea_from_afreq, 3)

fwrite(out_tsv_df, file = opt[["out-tsv"]], sep = "\t", quote = FALSE, na = "")

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("R package 'openxlsx' is required to write Excel output. Please install or load it.", call. = FALSE)
}
openxlsx::write.xlsx(out_num, file = opt[["out-xlsx"]], asTable = TRUE)

n_total <- nrow(out_num)
n_with_pip <- sum(is.finite(out_num$pip))
n_allele_not_aligned <- sum(!out_num$allele_align_ok & !is.na(out_num$allele_align_ok))

json_lines <- c(
  "{",
  sprintf('  "matched_tsv": "%s",', opt[["matched-tsv"]]),
  sprintf('  "afreq": "%s",', opt[["afreq"]]),
  sprintf('  "finemap_tsv": "%s",', opt[["finemap-tsv"]]),
  sprintf('  "n_rows": %d,', n_total),
  sprintf('  "n_rows_with_pip": %d,', n_with_pip),
  sprintf('  "n_allele_alignment_false": %d', n_allele_not_aligned),
  "}"
)
writeLines(json_lines, con = opt[["diag-json"]])

writeLines(format(Sys.time()), con = opt[["done-file"]])
cat("Done\n")
