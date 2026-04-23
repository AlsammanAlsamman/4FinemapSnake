#!/usr/bin/env bash
set -euo pipefail

INPUT_TSV=""
LD_MATRIX=""
OUT_TSV=""
DIAG_JSON=""
DONE_FILE=""
MAX_CAUSAL="10"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-tsv) INPUT_TSV="$2"; shift 2 ;;
    --ld-matrix) LD_MATRIX="$2"; shift 2 ;;
    --out-tsv) OUT_TSV="$2"; shift 2 ;;
    --diag-json) DIAG_JSON="$2"; shift 2 ;;
    --done-file) DONE_FILE="$2"; shift 2 ;;
    --max-causal) MAX_CAUSAL="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

mkdir -p "$(dirname "$OUT_TSV")" "$(dirname "$DIAG_JSON")" "$(dirname "$DONE_FILE")"

# Load R module on cluster compute nodes
module load R/4.5.1-mkl 2>/dev/null || true

Rscript - "$INPUT_TSV" "$LD_MATRIX" "$OUT_TSV" "$DIAG_JSON" "$MAX_CAUSAL" <<'R'
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
input_tsv <- args[1]
ld_matrix_file <- args[2]
out_tsv <- args[3]
diag_json <- args[4]
max_causal <- as.integer(args[5])

# Load matched GWAS
gwas <- fread(input_tsv, sep="\t", data.table=FALSE)
if (!("SNP" %in% colnames(gwas))) {
  stop("Matched table must include SNP column")
}

# Load LD matrix
ld_matrix <- as.matrix(read.table(ld_matrix_file, sep="\t", row.names=1, header=TRUE))

# Validate SNP/order match
ld_snps <- rownames(ld_matrix)
gwas_snps <- gwas$SNP

if (length(ld_snps) != nrow(gwas)) {
  stop(sprintf("LD matrix SNP count (%d) != GWAS count (%d)", length(ld_snps), nrow(gwas)))
}

if (!all(ld_snps == gwas_snps)) {
  stop("LD matrix SNP order does not match GWAS SNP order")
}

# Compute z-scores
if ("BETA" %in% colnames(gwas) && "SE" %in% colnames(gwas)) {
  z <- gwas$BETA / gwas$SE
} else {
  z <- -log10(gwas$P)
  z[is.infinite(z) | is.na(z)] <- 0
}

# Compute Bayes factors: BF = exp(z^2 / 2)
bayes_factor <- exp(z^2 / 2)

# Multivariate fine-mapping with LD structure
# Simplified approach: compute causal probabilities accounting for LD
# For each SNP as causal variant, compute likelihood given others are null

n_snps <- length(z)
# Use regularized inverse (Ledoit-Wolf with lambda=0.1)
lambda <- 0.1
ld_reg <- ld_matrix * (1 - lambda) + diag(n_snps) * lambda

# Compute inverse LD (may be singular, so use pseudo-inverse)
tryCatch({
  ld_inv <- solve(ld_reg)
}, error = function(e) {
  warning("LD matrix singular, using pseudo-inverse")
  ld_inv <<- MASS::ginv(ld_reg)
})

# Compute conditional Bayes factors accounting for LD
# posterior ~ BF * exp(-0.5 * z^2 * (1 - r^2) / (1-rho)) where rho accounts for others
# Simplified: use diagonal of LD inverse as shrinkage factor

ld_diag_inv <- diag(ld_inv)
bf_shrunk <- bayes_factor * sqrt(pmax(1 / ld_diag_inv, 0.1))  # shrink by LD info

# Normalize to probabilities
pip <- bf_shrunk / sum(bf_shrunk)
gwas$pip <- pip

# Sort by PIP and select top max_causal
gwas_sorted <- gwas[order(gwas$pip, decreasing=TRUE), ]
gwas_selected <- gwas_sorted[1:min(max_causal, nrow(gwas_sorted)), ]

# Output credible set
out_cols <- c("SNP", "P", "BETA", "SE", "pip")
out_cols <- out_cols[out_cols %in% colnames(gwas_selected)]
write.table(gwas_selected[, out_cols], file=out_tsv, sep="\t", row.names=FALSE, quote=FALSE)

# Diagnostics
diag <- list(
  method = "finemap_ld_aware",
  n_snps = nrow(gwas),
  n_selected = nrow(gwas_selected),
  max_causal = max_causal,
  lead_snp = gwas_selected$SNP[1],
  lead_pip = round(gwas_selected$pip[1], 4),
  ld_regularization = "Ledoit-Wolf",
  ld_lambda = lambda,
  credible_set_95_n = sum(cumsum(gwas_sorted$pip) <= 0.95) + 1
)

# Manual JSON output (no jsonlite dependency)
json_str <- "{\n"
for (i in seq_along(diag)) {
  key <- names(diag)[i]
  val <- diag[[i]]
  if (is.character(val)) {
    val_str <- sprintf('"%s"', val)
  } else if (is.numeric(val)) {
    val_str <- as.character(val)
  } else {
    val_str <- "null"
  }
  json_str <- paste0(json_str, sprintf('  "%s": %s', key, val_str))
  if (i < length(diag)) json_str <- paste0(json_str, ",")
  json_str <- paste0(json_str, "\n")
}
json_str <- paste0(json_str, "}\n")

writeLines(json_str, con=diag_json)

cat("FINEMAP complete\n")
R

echo "ok" > "$DONE_FILE"
