#!/usr/bin/env Rscript

targets <- c("asian", "disc")
locus <- "ARID3A"

for (t in targets) {
	matched_file <- sprintf("results/%s/03_match/%s/matched.tsv", t, locus)
	ld_file <- sprintf("results/%s/06_ld_qc/%s/ld_matrix_qc.tsv", t, locus)

	if (!file.exists(matched_file) || !file.exists(ld_file)) {
		cat(sprintf("target=%s missing_input=1\n", t))
		next
	}

	m <- read.table(matched_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
	if (!("SNP" %in% colnames(m))) {
		cat(sprintf("target=%s missing_SNP_column=1\n", t))
		next
	}

	m_sig <- m[m$P < 5e-4 & is.finite(m$P), , drop = FALSE]

	ld_raw <- read.table(ld_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
	ld_mat <- as.matrix(ld_raw)

	overlap_all <- intersect(m$SNP, rownames(ld_mat))
	overlap_sig <- intersect(m_sig$SNP, rownames(ld_mat))

	cat(sprintf(
		"target=%s n_matched=%d n_sig=%d n_ld_rows=%d overlap_all=%d overlap_sig=%d\n",
		t, nrow(m), nrow(m_sig), nrow(ld_mat), length(overlap_all), length(overlap_sig)
	))

	if (length(overlap_sig) > 1) {
		ld_sub <- ld_mat[overlap_sig, overlap_sig, drop = FALSE]
		ld_score <- rowSums(ld_sub^2) - 1
		ld_score_norm <- ld_score / (length(overlap_sig) - 1)
		cat(sprintf(
			"target=%s unique_ld_score_norm=%d min=%.8f median=%.8f max=%.8f\n",
			t,
			length(unique(round(ld_score_norm, 10))),
			min(ld_score_norm, na.rm = TRUE),
			median(ld_score_norm, na.rm = TRUE),
			max(ld_score_norm, na.rm = TRUE)
		))
	} else {
		cat(sprintf("target=%s insufficient_overlap_for_norm=1\n", t))
	}
}
