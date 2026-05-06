#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--loci-file", type = "character"),
  make_option("--cluster-dir", type = "character"),
  make_option("--target", type = "character"),
  make_option("--out-loci", type = "character"),
  make_option("--done-file", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("loci-file", "cluster-dir", "target", "out-loci", "done-file")
for (k in required) {
  if (is.null(opt[[k]]) || !nzchar(opt[[k]])) {
    stop(paste("Missing required argument --", k, sep = ""))
  }
}

safe_mkdir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

safe_mkdir(opt$`out-loci`)
safe_mkdir(opt$`done-file`)

# Read original loci file
loci <- fread(opt$`loci-file`, sep = "\t", header = TRUE)
setnames(loci, tolower(names(loci)))

if (!all(c("loci", "chr", "start", "end") %in% names(loci))) {
  stop("loci file must contain loci, chr, start, end columns")
}

message(sprintf("Read %d loci from %s", nrow(loci), opt$`loci-file`))

# Find all cluster summary files for this target
cluster_dir <- opt$`cluster-dir`
summary_files <- list.files(
  cluster_dir,
  pattern = "cluster_summary\\.tsv$",
  full.names = TRUE,
  recursive = TRUE
)

message(sprintf("Found %d cluster summary files", length(summary_files)))

# Read all cluster summaries
cluster_summaries <- list()
for (f in summary_files) {
  tryCatch({
    dt <- fread(f, sep = "\t", header = TRUE)
    if (nrow(dt) > 0 && "locus" %in% names(dt)) {
      cluster_summaries[[length(cluster_summaries) + 1]] <- dt
    }
  }, error = function(e) {
    warning(sprintf("Failed to read %s: %s", f, e$message))
  })
}

if (length(cluster_summaries) == 0) {
  message("No cluster summaries found; output loci unchanged")
  fwrite(loci, opt$`out-loci`, sep = "\t")
  writeLines("ok", opt$`done-file`)
  quit(status = 0)
}

all_clusters <- rbindlist(cluster_summaries, use.names = TRUE, fill = TRUE)
if (nrow(all_clusters) == 0) {
  message("No cluster data; output loci unchanged")
  fwrite(loci, opt$`out-loci`, sep = "\t")
  writeLines("ok", opt$`done-file`)
  quit(status = 0)
}

setnames(all_clusters, tolower(names(all_clusters)))

message(sprintf("Read %d cluster rows total", nrow(all_clusters)))

# Process each locus
split_loci <- list()

for (i in seq_len(nrow(loci))) {
  locus_name <- as.character(loci[[i, "loci"]])
  chr <- loci[[i, "chr"]]
  locus_start <- as.integer(loci[[i, "start"]])
  locus_end <- as.integer(loci[[i, "end"]])
  locus_width <- locus_end - locus_start

  # Get clusters for this locus
  locus_clusters <- all_clusters[locus == locus_name, ]

  if (nrow(locus_clusters) == 0) {
    # No clusters: keep unchanged
    split_loci[[length(split_loci) + 1]] <- data.table(
      loci = locus_name,
      chr = chr,
      start = locus_start,
      end = locus_end
    )
  } else if (nrow(locus_clusters) == 1) {
    # Single cluster
    n_sig_snps <- locus_clusters[[1, "n_sig_snps"]]
    lead_snp_p <- locus_clusters[[1, "lead_p"]]

    if (locus_width <= 250000) {
      # Keep unchanged if <= 250kb
      split_loci[[length(split_loci) + 1]] <- data.table(
        loci = locus_name,
        chr = chr,
        start = locus_start,
        end = locus_end
      )
    } else if (is.na(n_sig_snps) || n_sig_snps == 0) {
      # No sig SNPs: keep unchanged
      split_loci[[length(split_loci) + 1]] <- data.table(
        loci = locus_name,
        chr = chr,
        start = locus_start,
        end = locus_end
      )
    } else {
      # Shrink to 250kb window around lead SNP
      # The lead_snp position is not in the cluster summary, so we use cluster center as proxy
      cluster_mid <- (locus_clusters[[1, "cluster_start"]] + locus_clusters[[1, "cluster_end"]]) / 2

      # Center 250kb window on cluster midpoint
      window_half <- 250000 / 2
      win_start <- as.integer(cluster_mid - window_half)
      win_end <- as.integer(cluster_mid + window_half)

      # Clamp to locus boundaries
      win_start <- max(win_start, locus_start)
      win_end <- min(win_end, locus_end)

      split_loci[[length(split_loci) + 1]] <- data.table(
        loci = locus_name,
        chr = chr,
        start = win_start,
        end = win_end
      )
    }
  } else {
    # Multiple clusters: split by midpoints
    # Sort clusters by position
    setorder(locus_clusters, cluster_start)

    n_clusters <- nrow(locus_clusters)
    sub_loci_list <- list()

    # Calculate midpoints between consecutive clusters
    midpoints <- numeric(n_clusters - 1)
    for (j in seq_len(n_clusters - 1)) {
      c_end_curr <- locus_clusters[[j, "cluster_end"]]
      c_start_next <- locus_clusters[[j + 1, "cluster_start"]]
      midpoints[j] <- (c_end_curr + c_start_next) / 2
    }

    # Define boundaries for each sub-locus
    for (j in seq_len(n_clusters)) {
      if (j == 1) {
        sub_start <- locus_start
        sub_end <- as.integer(midpoints[1])
      } else if (j == n_clusters) {
        sub_start <- as.integer(midpoints[n_clusters - 1])
        sub_end <- locus_end
      } else {
        sub_start <- as.integer(midpoints[j - 1])
        sub_end <- as.integer(midpoints[j])
      }

      sub_name <- sprintf("%s_c%d", locus_name, j)
      sub_loci_list[[length(sub_loci_list) + 1]] <- data.table(
        loci = sub_name,
        chr = chr,
        start = sub_start,
        end = sub_end
      )
    }

    split_loci <- c(split_loci, sub_loci_list)
  }
}

# Combine all split loci into one table
result <- rbindlist(split_loci, use.names = TRUE)

message(sprintf("Output %d loci (from original %d)", nrow(result), nrow(loci)))

fwrite(result, opt$`out-loci`, sep = "\t")
writeLines("ok", opt$`done-file`)
