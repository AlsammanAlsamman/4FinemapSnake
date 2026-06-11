#!/usr/bin/env Rscript

#============================================================
# Bayesian Colocalization Analysis using coloc.abf
#============================================================
# This script performs genome-wide colocalization between GWAS
# signals and eQTL signals at isolated loci using Bayesian
# ABF methodology (coloc package).
#
# Key Features:
# - Per-locus analysis using pre-selected causal variants
# - Multi-dataset eQTL support (tissue-specific or whole-blood)
# - Credible set annotation (95% cumulative posterior)
# - Robust error handling and diagnostics
#============================================================

suppressMessages({
    library(data.table)
    library(stringr)
    library(dplyr)
    library(jsonlite)
    library(coloc)
})

options(scipen=10)
options(datatable.fread.datatable=FALSE)
options(warn=1)

#============================================================
# CONSTANTS AND CONFIGURATION
#============================================================

# Thresholds for colocalization signal strength
PRIOR_CAUSAL_GWAS <- 1e-4       # Prior probability of GWAS causal
PRIOR_CAUSAL_EQTL <- 1e-4       # Prior probability of eQTL causal
PRIOR_SHARED_CAUSAL <- 1e-5     # Prior probability of shared causal variant

# Minimum p-value for numerical stability
MIN_PVAL <- 1e-300
MAX_PVAL <- 1

# SE calculation fallbacks
DUMMY_SE <- 0.00001             # Fallback SE when p-value is 0 or 1

# MHC region exclusion (GRCh37/GRCh38)
MHC_CHR <- 6
MHC_START <- 28477897
MHC_END <- 33448354


#============================================================
# UTILITY FUNCTIONS
#============================================================

#' Normalize chromosome notation
#' Converts chr1, 1, CHR1, etc. to consistent "1" format
normalize_chr <- function(chr_str) {
    chr_str <- as.character(chr_str)
    chr_str <- gsub("^[Cc][Hh][Rr]?", "", chr_str)
    chr_str <- trimws(chr_str)
    return(chr_str)
}

#' Check if variant is in MHC region
is_in_mhc <- function(chr, pos) {
    chr_norm <- normalize_chr(chr)
    is_mhc <- (chr_norm == "6") & (pos >= MHC_START) & (pos <= MHC_END)
    return(is_mhc)
}

#' Estimate standard error from beta and p-value
#' Using: Z = beta / SE, where Z ~ N(0,1) under null
#' Therefore: SE = |beta| / |Z|, where |Z| = abs(qnorm(p/2))
estimate_se_from_pval <- function(beta, pval) {
    # Handle edge cases
    pval <- pmax(pval, MIN_PVAL)
    pval <- pmin(pval, MAX_PVAL)
    
    # Calculate z-score magnitude
    z_score <- abs(qnorm(pval / 2))
    
    # Handle zero z-scores
    z_score <- pmax(z_score, 1e-10)
    
    # Calculate SE
    se <- abs(beta) / z_score
    
    # Replace any remaining zeros or NaNs with fallback
    se[is.na(se) | se == 0 | is.infinite(se)] <- DUMMY_SE
    
    return(se)
}

#' Format coloc output with proper naming and calculations
format_coloc_results <- function(gwas_data, eqtl_data, coloc_abf_results, 
                                  gene_id, gene_name, dataset_name) {
    if (is.null(coloc_abf_results) || is.list(coloc_abf_results$summary) == FALSE) {
        return(NULL)
    }
    
    # Extract posterior probabilities
    pp_h0 <- as.numeric(coloc_abf_results$summary["PP.H0.abf"])
    pp_h1 <- as.numeric(coloc_abf_results$summary["PP.H1.abf"])
    pp_h2 <- as.numeric(coloc_abf_results$summary["PP.H2.abf"])
    pp_h3 <- as.numeric(coloc_abf_results$summary["PP.H3.abf"])
    pp_h4 <- as.numeric(coloc_abf_results$summary["PP.H4.abf"])
    
    # Replace NaN with 0
    pp_h0[is.na(pp_h0)] <- 0
    pp_h1[is.na(pp_h1)] <- 0
    pp_h2[is.na(pp_h2)] <- 0
    pp_h3[is.na(pp_h3)] <- 0
    pp_h4[is.na(pp_h4)] <- 0
    
    # Normalize to ensure sum = 1
    pp_sum <- pp_h0 + pp_h1 + pp_h2 + pp_h3 + pp_h4
    if (pp_sum > 0) {
        pp_h0 <- pp_h0 / pp_sum
        pp_h1 <- pp_h1 / pp_sum
        pp_h2 <- pp_h2 / pp_sum
        pp_h3 <- pp_h3 / pp_sum
        pp_h4 <- pp_h4 / pp_sum
    }
    
    # Prepare output data frame
    result <- data.frame(
        gene_id = gene_id,
        gene_name = gene_name,
        dataset = dataset_name,
        n_variants = nrow(coloc_abf_results$results),
        pp_h0 = pp_h0,
        pp_h1 = pp_h1,
        pp_h2 = pp_h2,
        pp_h3 = pp_h3,
        pp_h4 = pp_h4,
        nsnps_gwas = nrow(gwas_data),
        nsnps_eqtl = nrow(eqtl_data),
        nsnps_shared = length(which(!is.na(coloc_abf_results$results$SNP))),
        stringsAsFactors = FALSE
    )
    
    return(result)
}

#' Identify 95% credible set from posterior probabilities
identify_credible_set <- function(coloc_results, credible_threshold = 0.95) {
    if (is.null(coloc_results$results) || nrow(coloc_results$results) == 0) {
        return(NULL)
    }
    
    # Get posterior probabilities for each variant
    variants_df <- as.data.frame(coloc_results$results)
    
    # Sort by posterior probability (descending)
    variants_df <- variants_df[order(variants_df$SNP.PP.H4, decreasing = TRUE), ]
    
    # Calculate cumulative probability
    variants_df$cum_pp <- cumsum(as.numeric(variants_df$SNP.PP.H4))
    
    # Identify credible set (cumulative PP up to threshold)
    credible_set <- variants_df[which(variants_df$cum_pp <= credible_threshold), ]
    
    return(credible_set)
}

#' Prepare GWAS data for coloc.abf analysis
prepare_gwas_data <- function(summary_df, selected_snps_df) {
    # Filter summary to selected SNPs only
    gwas_dat <- summary_df[summary_df$SNP %in% selected_snps_df$SNP, ]
    
    if (nrow(gwas_dat) == 0) {
        warning("No selected SNPs found in summary data")
        return(NULL)
    }
    
    # Ensure required columns
    required_cols <- c("SNP", "CHR", "BP", "P", "BETA")
    missing_cols <- setdiff(required_cols, colnames(gwas_dat))
    
    if (length(missing_cols) > 0) {
        warning(sprintf("GWAS data missing columns: %s", paste(missing_cols, collapse=", ")))
        return(NULL)
    }
    
    # Handle SE column
    if ("SE" %in% colnames(gwas_dat)) {
        gwas_dat$SE <- as.numeric(gwas_dat$SE)
        # Replace zeros with estimated SE
        zero_se_idx <- which(gwas_dat$SE == 0 | is.na(gwas_dat$SE))
        if (length(zero_se_idx) > 0) {
            gwas_dat$SE[zero_se_idx] <- estimate_se_from_pval(
                gwas_dat$BETA[zero_se_idx],
                gwas_dat$P[zero_se_idx]
            )
        }
    } else {
        # Estimate SE from beta and p-value
        gwas_dat$SE <- estimate_se_from_pval(gwas_dat$BETA, gwas_dat$P)
    }
    
    # Ensure numeric columns
    gwas_dat$P <- as.numeric(gwas_dat$P)
    gwas_dat$BETA <- as.numeric(gwas_dat$BETA)
    gwas_dat$SE <- as.numeric(gwas_dat$SE)
    
    # Remove duplicate SNPs (keep first)
    gwas_dat <- gwas_dat[!duplicated(gwas_dat$SNP), ]
    
    # Filter out MHC region
    gwas_dat$CHR_norm <- normalize_chr(gwas_dat$CHR)
    gwas_dat$BP <- as.numeric(gwas_dat$BP)
    mhc_idx <- which(is_in_mhc(gwas_dat$CHR_norm, gwas_dat$BP))
    if (length(mhc_idx) > 0) {
        warning(sprintf("Excluding %d SNPs from MHC region", length(mhc_idx)))
        gwas_dat <- gwas_dat[-mhc_idx, ]
    }
    
    return(gwas_dat[, c("SNP", "P", "BETA", "SE")])
}

#' Prepare eQTL data for a specific gene in coloc.abf analysis
prepare_eqtl_data_for_gene <- function(eqtl_df, gene_id, gene_col="eGeneID") {
    # Filter to specific gene
    eqtl_gene <- eqtl_df[eqtl_df[[gene_col]] == gene_id, ]
    
    if (nrow(eqtl_gene) == 0) {
        return(NULL)
    }
    
    # Handle both column naming schemes
    snp_col <- if ("SNP" %in% colnames(eqtl_gene)) "SNP" else if ("rsID" %in% colnames(eqtl_gene)) "rsID" else NULL
    beta_col <- if ("BETA" %in% colnames(eqtl_gene)) "BETA" else if ("beta_eQTL" %in% colnames(eqtl_gene)) "beta_eQTL" else NULL
    pval_col <- if ("P" %in% colnames(eqtl_gene)) "P" else if ("pval_eQTL" %in% colnames(eqtl_gene)) "pval_eQTL" else NULL
    se_col <- if ("SE" %in% colnames(eqtl_gene)) "SE" else if ("SE_eQTL" %in% colnames(eqtl_gene)) "SE_eQTL" else NULL
    
    if (is.null(snp_col) || is.null(beta_col) || is.null(pval_col)) {
        return(NULL)
    }
    
    # Rename columns to standard names for processing
    eqtl_gene$SNP <- eqtl_gene[[snp_col]]
    eqtl_gene$P <- as.numeric(eqtl_gene[[pval_col]])
    eqtl_gene$BETA <- as.numeric(eqtl_gene[[beta_col]])
    
    # Handle SE column
    if (!is.null(se_col) && se_col %in% colnames(eqtl_gene)) {
        eqtl_gene$SE <- as.numeric(eqtl_gene[[se_col]])
        # Replace zeros with estimated SE
        zero_se_idx <- which(eqtl_gene$SE == 0 | is.na(eqtl_gene$SE))
        if (length(zero_se_idx) > 0) {
            eqtl_gene$SE[zero_se_idx] <- estimate_se_from_pval(
                eqtl_gene$BETA[zero_se_idx],
                eqtl_gene$P[zero_se_idx]
            )
        }
    } else {
        # Estimate SE from beta and p-value
        eqtl_gene$SE <- estimate_se_from_pval(eqtl_gene$BETA, eqtl_gene$P)
    }
    
    # Ensure numeric columns
    eqtl_gene$P <- as.numeric(eqtl_gene$P)
    eqtl_gene$BETA <- as.numeric(eqtl_gene$BETA)
    eqtl_gene$SE <- as.numeric(eqtl_gene$SE)
    
    # Remove duplicate SNPs
    eqtl_gene <- eqtl_gene[!duplicated(eqtl_gene$SNP), ]
    
    return(eqtl_gene[, c("SNP", "P", "BETA", "SE")])
}

#' Run coloc.abf analysis between GWAS and eQTL signals
run_coloc_abf_locus <- function(gwas_data, eqtl_data, 
                                 sample_size_gwas = NULL, 
                                 sample_size_eqtl = NULL,
                                 analysis_type = "quant") {
    if (is.null(gwas_data) || is.null(eqtl_data)) {
        return(NULL)
    }
    
    if (nrow(gwas_data) < 2 || nrow(eqtl_data) < 2) {
        return(NULL)
    }
    
    # Prepare datasets for coloc
    gwas_list <- list(
        beta = gwas_data$BETA,
        se = gwas_data$SE,
        pval = gwas_data$P,
        N = sample_size_gwas,
        type = analysis_type,
        snp = gwas_data$SNP
    )
    
    eqtl_list <- list(
        beta = eqtl_data$BETA,
        se = eqtl_data$SE,
        pval = eqtl_data$P,
        N = sample_size_eqtl,
        type = analysis_type,
        snp = eqtl_data$SNP
    )
    
    # Run coloc.abf with error handling
    coloc_result <- tryCatch({
        coloc::coloc.abf(gwas_list, eqtl_list,
                         p12 = 1e-5,
                         p1 = PRIOR_CAUSAL_GWAS,
                         p2 = PRIOR_CAUSAL_EQTL,
                         p12 = PRIOR_SHARED_CAUSAL)
    }, error = function(e) {
        warning(sprintf("coloc.abf failed: %s", e$message))
        return(NULL)
    })
    
    return(coloc_result)
}


#============================================================
# MAIN ANALYSIS PIPELINE
#============================================================

main <- function() {
    # Parse command-line arguments
    args <- parse_arguments()
    
    # Read input files
    cat("[INFO] Reading input files...\n")
    
    summary_df <- data.table::fread(args$summary_tsv, sep="\t", header=TRUE)
    summary_df <- as.data.frame(summary_df)
    
    selected_snps_df <- data.table::fread(args$selected_snps_tsv, sep="\t", header=TRUE)
    selected_snps_df <- as.data.frame(selected_snps_df)
    
    eqtl_combined_df <- data.table::fread(args$eqtl_combined_tsv, sep="\t", header=TRUE)
    eqtl_combined_df <- as.data.frame(eqtl_combined_df)
    
    cat(sprintf("[INFO] Loaded summary: %d SNPs\n", nrow(summary_df)))
    cat(sprintf("[INFO] Loaded selected SNPs: %d variants\n", nrow(selected_snps_df)))
    cat(sprintf("[INFO] Loaded eQTL combined: %d records\n", nrow(eqtl_combined_df)))
    
    # Prepare GWAS data
    cat("[INFO] Preparing GWAS data...\n")
    gwas_prepared <- prepare_gwas_data(summary_df, selected_snps_df)
    
    if (is.null(gwas_prepared) || nrow(gwas_prepared) < 2) {
        stop("Failed to prepare GWAS data or insufficient SNPs")
    }
    
    cat(sprintf("[INFO] GWAS data prepared: %d SNPs\n", nrow(gwas_prepared)))
    
    # Get unique genes in eQTL data
    # Handle both column naming schemes
    gene_col <- if ("gene_id" %in% colnames(eqtl_combined_df)) "gene_id" else if ("eGeneID" %in% colnames(eqtl_combined_df)) "eGeneID" else NULL
    
    if (is.null(gene_col)) {
        stop("Cannot find gene ID column in eQTL data. Expected 'gene_id' or 'eGeneID'")
    }
    
    gene_list <- unique(eqtl_combined_df[[gene_col]])
    cat(sprintf("[INFO] Found %d unique genes in eQTL data\n", length(gene_list)))
    
    # Initialize results storage
    all_results <- list()
    all_credible_sets <- list()
    diagnostics <- list()
    
    # Main colocalization loop
    cat("[INFO] Starting colocalization analysis...\n")
    n_genes <- length(gene_list)
    n_successes <- 0
    n_failures <- 0
    
    # Determine gene and dataset columns
    gene_col <- if ("gene_id" %in% colnames(eqtl_combined_df)) "gene_id" else "eGeneID"
    gene_name_col <- if ("gene_name" %in% colnames(eqtl_combined_df)) "gene_name" else if ("eGeneName" %in% colnames(eqtl_combined_df)) "eGeneName" else NULL
    dataset_col <- if ("dataset" %in% colnames(eqtl_combined_df)) "dataset" else if ("source_dataset" %in% colnames(eqtl_combined_df)) "source_dataset" else NULL
    
    for (i in seq_along(gene_list)) {
        gene_id <- gene_list[i]
        
        if (i %% 1000 == 0) {
            cat(sprintf("[INFO] Processing gene %d/%d: %s\n", i, n_genes, gene_id))
        }
        
        # Prepare eQTL data for this gene
        eqtl_prepared <- prepare_eqtl_data_for_gene(eqtl_combined_df, gene_id, gene_col)
        
        if (is.null(eqtl_prepared) || nrow(eqtl_prepared) < 2) {
            next
        }
        
        # Get gene name and dataset
        gene_info <- eqtl_combined_df[eqtl_combined_df[[gene_col]] == gene_id, ][1, ]
        gene_name <- if (!is.null(gene_name_col) && gene_name_col %in% colnames(gene_info)) gene_info[[gene_name_col]] else gene_id
        dataset_name <- if (!is.null(dataset_col) && dataset_col %in% colnames(gene_info)) gene_info[[dataset_col]] else "unknown"
        
        # Run coloc.abf
        coloc_result <- run_coloc_abf_locus(
            gwas_prepared, 
            eqtl_prepared,
            sample_size_gwas = as.numeric(args$sample_size_gwas),
            sample_size_eqtl = NULL,
            analysis_type = args$analysis_type
        )
        
        if (is.null(coloc_result)) {
            n_failures <- n_failures + 1
            next
        }
        
        # Format results
        formatted_result <- format_coloc_results(
            gwas_prepared, eqtl_prepared, coloc_result,
            gene_id, gene_name, dataset_name
        )
        
        if (!is.null(formatted_result)) {
            all_results[[length(all_results) + 1]] <- formatted_result
            n_successes <- n_successes + 1
            
            # If PP.H4 meets threshold or output_all_hits enabled, store credible set
            if (formatted_result$pp_h4 >= as.numeric(args$pp_h4_threshold) || 
                args$output_all_hits == "true") {
                
                credible_set <- identify_credible_set(coloc_result, 
                                                      as.numeric(args$pp_h4_credible))
                if (!is.null(credible_set)) {
                    credible_set$gene_id <- gene_id
                    credible_set$gene_name <- gene_name
                    credible_set$dataset <- dataset_name
                    credible_set$pp_h4 <- formatted_result$pp_h4
                    
                    all_credible_sets[[length(all_credible_sets) + 1]] <- credible_set
                }
            }
        }
    }
    
    cat(sprintf("[INFO] Analysis complete. Successes: %d, Failures: %d\n", 
                n_successes, n_failures))
    
    # Combine and write results
    cat("[INFO] Writing results...\n")
    
    if (length(all_results) > 0) {
        results_df <- do.call(rbind, all_results)
        rownames(results_df) <- NULL
        
        # Sort by PP.H4 descending
        results_df <- results_df[order(results_df$pp_h4, decreasing=TRUE), ]
        
        # Write full results
        data.table::fwrite(results_df, args$out_tsv, sep="\t", quote=FALSE)
        cat(sprintf("[INFO] Wrote %d results to %s\n", nrow(results_df), args$out_tsv))
        
        # Write summary (top hits)
        summary_df_final <- results_df[results_df$pp_h4 >= as.numeric(args$pp_h4_threshold), ]
        data.table::fwrite(summary_df_final, args$out_summary_tsv, sep="\t", quote=FALSE)
        cat(sprintf("[INFO] Wrote %d summary results to %s\n", 
                    nrow(summary_df_final), args$out_summary_tsv))
    }
    
    # Write credible sets
    if (length(all_credible_sets) > 0) {
        credible_df <- do.call(rbind, all_credible_sets)
        rownames(credible_df) <- NULL
        
        # Select key columns for output
        key_cols <- c("gene_id", "gene_name", "dataset", "SNP", "SNP.PP.H4", "cum_pp", "pp_h4")
        out_cols <- intersect(key_cols, colnames(credible_df))
        credible_df <- credible_df[, out_cols]
        
        data.table::fwrite(credible_df, args$out_credible_tsv, sep="\t", quote=FALSE)
        cat(sprintf("[INFO] Wrote %d credible set variants to %s\n", 
                    nrow(credible_df), args$out_credible_tsv))
    }
    
    # Write diagnostics
    diagnostics$n_genes_analyzed <- n_genes
    diagnostics$n_genes_success <- n_successes
    diagnostics$n_genes_failed <- n_failures
    diagnostics$analysis_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    diagnostics$r_version <- paste(R.version$major, R.version$minor, sep=".")
    tryCatch({
        diagnostics$coloc_version <- as.character(packageVersion("coloc"))
    }, error = function(e) {
        diagnostics$coloc_version <- "unknown"
    })
    
    jsonlite::write_json(diagnostics, args$out_diagnostics_json, pretty=TRUE)
    cat(sprintf("[INFO] Wrote diagnostics to %s\n", args$out_diagnostics_json))
    
    # Write log file marker
    cat(sprintf("[INFO] Analysis completed successfully at %s\n", 
                format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    cat("Done\n", file=args$out_log_file, append=TRUE)
    
    # Create done marker
    writeLines("COLOC analysis completed successfully", args$done_file)
    cat("[INFO] Analysis marked as complete\n")
}

#' Parse command-line arguments
parse_arguments <- function() {
    parser_args <- commandArgs(trailingOnly=TRUE)
    
    args <- list(
        summary_tsv = NULL,
        selected_snps_tsv = NULL,
        eqtl_combined_tsv = NULL,
        target = NULL,
        locus = NULL,
        sample_size_gwas = 100000,
        pp_h4_threshold = 0.75,
        pp_h4_credible = 0.95,
        analysis_type = "quant",
        output_all_hits = "true",
        out_tsv = NULL,
        out_summary_tsv = NULL,
        out_credible_tsv = NULL,
        out_diagnostics_json = NULL,
        out_log_file = NULL,
        done_file = NULL
    )
    
    for (i in seq(1, length(parser_args), by=2)) {
        if (i+1 <= length(parser_args)) {
            key <- sub("^--", "", parser_args[i])
            key <- gsub("-", "_", key)  # Convert dashes to underscores
            value <- parser_args[i+1]
            
            if (key %in% names(args)) {
                args[[key]] <- value
            }
        }
    }
    
    # Validate required arguments
    required <- c("summary_tsv", "selected_snps_tsv", "eqtl_combined_tsv", 
                  "out_tsv", "sample_size_gwas")
    missing <- required[sapply(args[required], is.null)]
    
    if (length(missing) > 0) {
        stop(sprintf("Missing required arguments: %s", paste(missing, collapse=", ")))
    }
    
    return(args)
}

# Run main analysis
main()
