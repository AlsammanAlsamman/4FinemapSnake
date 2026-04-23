#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop(sprintf("Missing value for argument: %s", flag))
  args[idx + 1]
}

input_tsv <- get_arg("--input-tsv")
bfile_prefix <- get_arg("--bfile-prefix")
out_tsv <- get_arg("--out-tsv")
diag_json <- get_arg("--diag-json")
done_file <- get_arg("--done-file")
p_cutoff <- as.numeric(get_arg("--p-cutoff", "5e-8"))
max_iter <- as.integer(get_arg("--max-iterations", "25"))
sample_size <- as.numeric(get_arg("--sample-size"))
plot_dir <- get_arg("--plot-dir", file.path(dirname(out_tsv), "iterations_gcta"))
gcta_bin <- get_arg("--gcta-bin", "gcta")

required <- list(input_tsv, bfile_prefix, out_tsv, diag_json, done_file)
if (any(sapply(required, is.null))) {
  stop("Required args: --input-tsv --bfile-prefix --out-tsv --diag-json --done-file")
}
if (!is.finite(sample_size) || sample_size <= 0) {
  stop("--sample-size must be a positive number")
}

mkdir_for <- function(path) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
mkdir_for(out_tsv)
mkdir_for(diag_json)
mkdir_for(done_file)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

tmp_dir <- file.path(dirname(out_tsv), "gcta_tmp")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

# Load matched GWAS that already matches refpanel subset from Step 3.
gwas <- read.table(input_tsv, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
need <- c("SNP", "A1", "A2", "P", "BETA", "SE", "BP")
miss <- setdiff(need, colnames(gwas))
if (length(miss) > 0) {
  stop(sprintf("Matched GWAS missing columns: %s", paste(miss, collapse = ",")))
}

gwas$SNP <- as.character(gwas$SNP)
gwas$A1 <- as.character(gwas$A1)
gwas$A2 <- as.character(gwas$A2)
gwas$P <- suppressWarnings(as.numeric(gwas$P))
gwas$BETA <- suppressWarnings(as.numeric(gwas$BETA))
gwas$SE <- suppressWarnings(as.numeric(gwas$SE))
gwas$BP <- suppressWarnings(as.numeric(gwas$BP))

bad <- !is.finite(gwas$P) | gwas$P <= 0 | gwas$P > 1 |
  !is.finite(gwas$BETA) | !is.finite(gwas$SE) | gwas$SE <= 0 | !is.finite(gwas$BP)
if (any(bad)) {
  gwas <- gwas[!bad, , drop = FALSE]
}
if (nrow(gwas) == 0) stop("No valid GWAS rows after filtering")

freq_col <- if ("REF_MAF" %in% colnames(gwas)) "REF_MAF" else if ("MAF" %in% colnames(gwas)) "MAF" else NULL
if (is.null(freq_col)) {
  gwas$freq_use <- 0.2
} else {
  gwas$freq_use <- suppressWarnings(as.numeric(gwas[[freq_col]]))
  gwas$freq_use[!is.finite(gwas$freq_use) | gwas$freq_use <= 0 | gwas$freq_use >= 1] <- 0.2
}

# Build .ma file required by GCTA-COJO.
ma <- data.frame(
  SNP = gwas$SNP,
  A1 = gwas$A1,
  A2 = gwas$A2,
  freq = gwas$freq_use,
  b = gwas$BETA,
  se = gwas$SE,
  p = gwas$P,
  N = rep(sample_size, nrow(gwas)),
  stringsAsFactors = FALSE
)
ma_file <- file.path(tmp_dir, "cojo_input.ma")
write.table(ma, ma_file, sep = "\t", row.names = FALSE, quote = FALSE)

run_gcta_cond <- function(selected_snps, prefix) {
  cond_file <- file.path(tmp_dir, paste0(prefix, ".cond.snps"))
  writeLines(selected_snps, cond_file)
  out_prefix <- file.path(tmp_dir, prefix)

  cmd_args <- c(
    "--bfile", bfile_prefix,
    "--cojo-file", ma_file,
    "--cojo-cond", cond_file,
    "--out", out_prefix
  )

  status <- system2(gcta_bin, args = cmd_args, stdout = TRUE, stderr = TRUE)
  status_code <- attr(status, "status")
  if (!is.null(status_code) && status_code != 0) {
    stop(sprintf("GCTA failed for %s: %s", prefix, paste(status, collapse = "\n")))
  }

  cma_file <- paste0(out_prefix, ".cma.cojo")
  if (!file.exists(cma_file)) {
    stop(sprintf("Missing GCTA output: %s", cma_file))
  }

  cma <- read.table(cma_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (!("SNP" %in% colnames(cma))) {
    stop(sprintf("GCTA cma missing SNP column: %s", cma_file))
  }
  cma
}

plot_manhattan <- function(pvals, selected, tag) {
  d <- data.frame(BP = gwas$BP, SNP = gwas$SNP, p = pvals, selected = FALSE, stringsAsFactors = FALSE)
  d$selected[d$SNP %in% selected] <- TRUE
  d$mlog10p <- -log10(pmax(pmin(d$p, 1), 1e-300))

  p <- ggplot(d, aes(x = BP, y = mlog10p)) +
    geom_point(color = "#9CA3AF", size = 1.5, alpha = 0.75) +
    geom_point(data = d[d$selected, , drop = FALSE], color = "#D73027", size = 2.3, alpha = 0.95) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "#B91C1C", linewidth = 0.5) +
    labs(title = tag, x = "Position (BP)", y = "-log10(P)") +
    theme_bw(base_size = 11)

  png <- file.path(plot_dir, paste0(gsub("[^A-Za-z0-9_\\-]", "_", tag), ".png"))
  ggsave(png, p, width = 10, height = 4.8, dpi = 150)
}

p_current <- gwas$P
selected <- character(0)
records <- list()
stop_reason <- "p_threshold"

# Base plot at iteration 0 using original GWAS p-values.
plot_manhattan(p_current, selected, "iter_00_base_original")

for (it in seq_len(max_iter)) {
  remaining <- !(gwas$SNP %in% selected)
  if (!any(remaining)) {
    stop_reason <- "no_candidates"
    break
  }

  idx_rem <- which(remaining)
  lead_idx <- idx_rem[which.min(p_current[idx_rem])]
  lead_p <- p_current[lead_idx]

  if (!is.finite(lead_p) || lead_p > p_cutoff) {
    stop_reason <- "p_threshold"
    break
  }

  lead_snp <- gwas$SNP[lead_idx]
  prev_selected <- selected
  selected <- c(selected, lead_snp)

  cma <- run_gcta_cond(selected, sprintf("iter_%02d", it))

  p_col <- if ("pC" %in% colnames(cma)) "pC" else if ("p" %in% colnames(cma)) "p" else NULL
  if (is.null(p_col)) {
    stop("GCTA cma output has no pC/p column")
  }
  beta_col <- if ("bC" %in% colnames(cma)) "bC" else if ("b" %in% colnames(cma)) "b" else NA
  se_col <- if ("bC_se" %in% colnames(cma)) "bC_se" else if ("se" %in% colnames(cma)) "se" else NA

  p_next <- setNames(rep(NA_real_, nrow(gwas)), gwas$SNP)
  p_next[as.character(cma$SNP)] <- suppressWarnings(as.numeric(cma[[p_col]]))
  p_next <- pmax(pmin(p_next, 1), 1e-300)

  # Plot after conditioning on selected SNPs.
  plot_manhattan(p_next, selected, sprintf("iter_%02d_conditioned_on_%d_snps", it, length(selected)))

  cma_row <- cma[cma$SNP == lead_snp, , drop = FALSE]
  beta_cond <- NA_real_
  se_cond <- NA_real_
  z_cond <- NA_real_
  if (nrow(cma_row) > 0 && !is.na(beta_col) && !is.na(se_col)) {
    beta_cond <- suppressWarnings(as.numeric(cma_row[[beta_col]][1]))
    se_cond <- suppressWarnings(as.numeric(cma_row[[se_col]][1]))
    if (is.finite(beta_cond) && is.finite(se_cond) && se_cond > 0) {
      z_cond <- beta_cond / se_cond
    }
  }

  records[[length(records) + 1]] <- data.frame(
    iteration = it,
    lead_snp = lead_snp,
    chr_bp = gwas$BP[lead_idx],
    p_original = gwas$P[lead_idx],
    p_conditional = lead_p,
    beta_conditional = beta_cond,
    se_conditional = se_cond,
    z_conditional = z_cond,
    n_conditioned_on = length(prev_selected),
    conditioned_on_snps = paste(prev_selected, collapse = ","),
    stringsAsFactors = FALSE
  )

  p_current <- p_next

  if (it == max_iter) {
    stop_reason <- "max_iterations"
  }
}

if (length(records) > 0) {
  out <- do.call(rbind, records)
} else {
  out <- data.frame(
    iteration = integer(0),
    lead_snp = character(0),
    chr_bp = numeric(0),
    p_original = numeric(0),
    p_conditional = numeric(0),
    beta_conditional = numeric(0),
    se_conditional = numeric(0),
    z_conditional = numeric(0),
    n_conditioned_on = integer(0),
    conditioned_on_snps = character(0),
    stringsAsFactors = FALSE
  )
}

write.table(out, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

json_lines <- c(
  "{",
  "  \"method\": \"iterative_conditional_gcta_cojo\",",
  sprintf("  \"gcta_bin\": \"%s\",", gcta_bin),
  sprintf("  \"sample_size\": %.12g,", sample_size),
  sprintf("  \"p_cutoff\": %.12g,", p_cutoff),
  sprintf("  \"max_iterations\": %d,", max_iter),
  sprintf("  \"iterations_run\": %d,", nrow(out)),
  sprintf("  \"stop_reason\": \"%s\",", stop_reason),
  sprintf("  \"n_snps\": %d,", nrow(gwas)),
  sprintf("  \"plot_dir\": \"%s\"", gsub("\\\\", "\\\\\\\\", plot_dir)),
  "}"
)
writeLines(json_lines, diag_json)
writeLines("ok", done_file)
