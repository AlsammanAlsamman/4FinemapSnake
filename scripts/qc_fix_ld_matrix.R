#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx == length(args)) {
    stop(sprintf("Missing required argument: %s", flag))
  }
  args[idx + 1]
}

input_matrix <- get_arg("--input-matrix")
out_matrix <- get_arg("--out-matrix")
diag_json <- get_arg("--diag-json")
done_file <- get_arg("--done-file")
tolerance <- as.numeric(get_arg("--tolerance"))

mkdir_for <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

mkdir_for(out_matrix)
mkdir_for(diag_json)
mkdir_for(done_file)

mat_df <- read.table(input_matrix, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(mat_df) < 2) {
  stop("Input matrix file is malformed")
}

row_ids <- mat_df[[1]]
mat <- as.matrix(mat_df[, -1, drop = FALSE])
storage.mode(mat) <- "double"
col_ids <- colnames(mat)

if (nrow(mat) != ncol(mat)) {
  stop("LD matrix must be square")
}
if (length(row_ids) != nrow(mat) || !all(as.character(row_ids) == as.character(col_ids))) {
  stop("Row/column SNP IDs in matrix are not aligned")
}

before <- if (length(mat) == 0) 0.0 else max(abs(mat - t(mat)))
fixed <- (mat + t(mat)) / 2.0
if (nrow(fixed) > 0) {
  diag(fixed) <- 1.0
}
after <- if (length(fixed) == 0) 0.0 else max(abs(fixed - t(fixed)))

out <- data.frame(SNP = as.character(row_ids), fixed, check.names = FALSE)
colnames(out)[-1] <- as.character(row_ids)
write.table(out, file = out_matrix, sep = "\t", quote = FALSE, row.names = FALSE)

diag <- list(
  max_abs_asymmetry_before = as.numeric(before),
  max_abs_asymmetry_after = as.numeric(after),
  tolerance = as.numeric(tolerance),
  corrected = isTRUE(before > 0)
)

json_lines <- c(
  "{",
  sprintf("  \"max_abs_asymmetry_before\": %.12g,", diag$max_abs_asymmetry_before),
  sprintf("  \"max_abs_asymmetry_after\": %.12g,", diag$max_abs_asymmetry_after),
  sprintf("  \"tolerance\": %.12g,", diag$tolerance),
  sprintf("  \"corrected\": %s", tolower(as.character(diag$corrected))),
  "}"
)
writeLines(json_lines, diag_json)
writeLines("ok", done_file)
