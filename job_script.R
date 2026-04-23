library(data.table)
dt <- fread("inputs/Disc.tsv")
# Ensure chrom is treated correctly (assuming it might have 'chr' prefix or just number)
# Taking chr19:238434-1785772
# We'll filter based on chrom == 19 (or "19" or "chr19") and pos between 238434 and 1785772
# First, let's detect the format of chrom column
sample_chrom <- as.character(dt$chrom[1])
chrom_val <- if (grepl("chr", sample_chrom)) "chr19" else "19"

matches <- dt[chrom == chrom_val & pos >= 238434 & pos <= 1785772]
cat("Count:", nrow(matches), "\n")
if (nrow(matches) > 0) {
    print(head(matches[, .(chrom, pos, snpid, markername, ea, nea, p)]))
}
