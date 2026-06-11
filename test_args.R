#!/usr/bin/env Rscript

# Test argument parsing
args <- commandArgs(trailingOnly=TRUE)

cat("Number of arguments:", length(args), "\n")
cat("Arguments received:\n")
for (i in seq_along(args)) {
    cat(sprintf("  [%d] '%s'\n", i, args[i]))
}

# Parse like in run_coloc.R
test_args <- list(summary_tsv = NULL, out_tsv = NULL)

for (i in seq(1, length(args), by=2)) {
    if (i+1 <= length(args)) {
        key <- sub("^--", "", args[i])
        value <- args[i+1]
        cat(sprintf("Parsing: key='%s', value='%s'\n", key, value))
        
        if (key %in% names(test_args)) {
            test_args[[key]] <- value
            cat(sprintf("  -> Assigned %s = %s\n", key, value))
        } else {
            cat(sprintf("  -> key '%s' not in names\n", key))
        }
    }
}

cat("\nFinal parsed args:\n")
print(test_args)
