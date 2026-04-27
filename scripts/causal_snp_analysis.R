# =============================================================================
# CAUSAL SNP IDENTIFICATION ANALYSIS
# Inputs : snp_master_table.tsv | ld_matrix_r2.tsv | ld_matrix_r.tsv | ld_matrix_Dprime.tsv
# Outputs: all_plots.pdf | causal_snp_summary.xlsx
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) {
    if (!is.null(default)) {
      return(default)
    }
    stop(sprintf("Missing required argument: %s", flag))
  }
  if (idx == length(args)) {
    stop(sprintf("Missing value for argument: %s", flag))
  }
  args[idx + 1]
}

mkdir_for <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

snp_master_tsv <- get_arg("--snp-master")
ld_r2_tsv <- get_arg("--ld-r2")
ld_r_tsv <- get_arg("--ld-r")
ld_dprime_tsv <- get_arg("--ld-dprime")
out_pdf <- get_arg("--out-pdf")
out_xlsx <- get_arg("--out-xlsx")
done_file <- get_arg("--done-file", default = "")

mkdir_for(out_pdf)
mkdir_for(out_xlsx)
if (nzchar(done_file)) {
  mkdir_for(done_file)
}

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
required <- c("tidyverse","ggrepel","pheatmap","cowplot","openxlsx",
              "RColorBrewer","viridis","scales","corrplot","ggridges",
              "patchwork","ggpubr","reshape2","glue")
to_install <- required[!required %in% installed.packages()[,"Package"]]
if(length(to_install)) install.packages(to_install, repos="https://cloud.r-project.org")
invisible(lapply(required, library, character.only=TRUE))

# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────
message("Loading data...")
snp <- read_tsv(snp_master_tsv, na = "NA", show_col_types = FALSE) %>%
  mutate(
    SNP_rsID = as.character(SNP_rsID),
    CS_membership = as.character(CS_membership),
    CS_size = as.numeric(CS_size)
  )

r2_mat     <- read_tsv(ld_r2_tsv)     %>% column_to_rownames("SNP_rsID") %>% as.matrix()
r_mat      <- read_tsv(ld_r_tsv)      %>% column_to_rownames("SNP_rsID") %>% as.matrix()
dprime_mat <- read_tsv(ld_dprime_tsv) %>% column_to_rownames("SNP_rsID") %>% as.matrix()

# ── 2. DERIVED COLUMNS ────────────────────────────────────────────────────────
eqtl_datasets <- snp %>% select(starts_with("pval_eQTL@")) %>% names() %>%
  str_remove("pval_eQTL@")

snp <- snp %>%
  mutate(
    in_CS            = !is.na(CS_membership),
    log10p_norm      = GWAS_log10p / max(GWAS_log10p, na.rm=TRUE),
    ld_score_norm    = (LD_score - min(LD_score,na.rm=TRUE)) /
                       (max(LD_score,na.rm=TRUE) - min(LD_score,na.rm=TRUE)),
    pip_tier         = case_when(
      Composite_PIP > 0.5  ~ "High (>0.5)",
      Composite_PIP > 0.1  ~ "Moderate (0.1–0.5)",
      TRUE                 ~ "Low (<0.1)"
    ) %>% factor(levels=c("High (>0.5)","Moderate (0.1–0.5)","Low (<0.1)")),
    n_eqtl_sig       = rowSums(across(starts_with("pval_eQTL@"),
                                      ~(!is.na(.x) & .x < 1e-5))),
    min_eqtl_pval    = do.call(pmin, c(across(starts_with("pval_eQTL@")),
                                        list(na.rm=TRUE))),
    log10_min_eqtl   = -log10(min_eqtl_pval + 1e-320),
    eqtl_sig_tier    = case_when(
      min_eqtl_pval < 1e-10 ~ "Strong (<1e-10)",
      min_eqtl_pval < 1e-5  ~ "Moderate (<1e-5)",
      !is.na(min_eqtl_pval) ~ "Weak / none",
      TRUE                  ~ "No data"
    ) %>% factor(levels=c("Strong (<1e-10)","Moderate (<1e-5)","Weak / none","No data")),
    label_top        = if_else(Composite_PIP > 0.1 | in_CS, SNP_rsID, NA_character_)
  )

cs_snps   <- snp %>% filter(in_CS)
top_snps  <- snp %>% filter(Composite_PIP > 0.05)

# ── 3. COLOUR PALETTES ────────────────────────────────────────────────────────
pip_cols  <- c("High (>0.5)"="#1D9E75","Moderate (0.1–0.5)"="#EF9F27","Low (<0.1)"="#B4B2A9")
eqtl_cols <- c("Strong (<1e-10)"="#534AB7","Moderate (<1e-5)"="#378ADD",
                "Weak / none"="#B4B2A9","No data"="#E8E6DF")
cs_col    <- "#D85A30"
hi_col    <- "#1D9E75"

theme_set(theme_bw(base_size=11) +
          theme(plot.title=element_text(face="bold", size=12),
                plot.subtitle=element_text(size=9, colour="grey40"),
                strip.background=element_rect(fill="grey95"),
                panel.grid.minor=element_blank()))

pdf(out_pdf, width=14, height=10, onefile=TRUE)

# =============================================================================
# SECTION A  – FINE-MAPPING OVERVIEW
# =============================================================================

# A1: FINEMAP PIP vs SuSiE PIP scatter ----------------------------------------
message("Plot A1: FINEMAP vs SuSiE PIP scatter")
p_a1 <- ggplot(snp, aes(FINEMAP_PIP, SuSiE_PIP)) +
  annotate("rect", xmin=0.5, xmax=1, ymin=0.5, ymax=1,
           fill="#1D9E75", alpha=0.07) +
  geom_point(aes(colour=pip_tier, size=GWAS_log10p), alpha=0.75) +
  geom_abline(linetype="dashed", colour="grey50") +
  geom_vline(xintercept=0.5, linetype="dotted", colour="grey60") +
  geom_hline(yintercept=0.5, linetype="dotted", colour="grey60") +
  geom_text_repel(data=filter(snp, !is.na(label_top)),
                  aes(label=label_top), size=2.8, max.overlaps=20,
                  segment.colour="grey60", segment.size=0.3) +
  scale_colour_manual(values=pip_cols, name="PIP tier") +
  scale_size_continuous(name="–log₁₀(p)", range=c(1,5)) +
  annotate("text", x=0.75, y=0.52, label="Consensus causal zone",
           colour="#0F6E56", size=3, fontface="italic") +
  labs(title="A1 · FINEMAP PIP vs SuSiE PIP",
       subtitle="Top-right = both tools agree on causality | diagonal = perfect agreement | size = GWAS significance",
       x="FINEMAP PIP", y="SuSiE PIP") +
  coord_fixed(xlim=c(0,1), ylim=c(0,1))
print(p_a1)

# A2: Composite PIP ranked lollipop -------------------------------------------
message("Plot A2: Composite PIP lollipop")
top30 <- snp %>% arrange(desc(Composite_PIP)) %>% slice_head(n=30) %>%
  mutate(SNP_rsID = fct_reorder(SNP_rsID, Composite_PIP))
p_a2 <- ggplot(top30, aes(Composite_PIP, SNP_rsID)) +
  geom_segment(aes(x=0, xend=Composite_PIP, y=SNP_rsID, yend=SNP_rsID,
                   colour=eqtl_sig_tier), linewidth=0.8) +
  geom_point(aes(colour=eqtl_sig_tier, shape=in_CS), size=3.5) +
  geom_vline(xintercept=0.5, linetype="dashed", colour="grey40") +
  scale_colour_manual(values=eqtl_cols, name="eQTL evidence") +
  scale_shape_manual(values=c(`TRUE`=19,`FALSE`=1), name="In credible set") +
  labs(title="A2 · Top 30 SNPs by composite PIP",
       subtitle="Color = best eQTL signal across all tissues | filled = credible set member | dashed = PIP 0.5 threshold",
       x="Composite PIP [(FINEMAP + SuSiE) / 2]", y=NULL) +
  theme(axis.text.y=element_text(size=8))
print(p_a2)

# A3: PIP distribution histogram ----------------------------------------------
message("Plot A3: PIP distribution")
snp_long_pip <- snp %>% select(SNP_rsID, FINEMAP_PIP, SuSiE_PIP) %>%
  pivot_longer(-SNP_rsID, names_to="Tool", values_to="PIP") %>%
  mutate(Tool=str_remove(Tool,"_PIP"))
p_a3 <- ggplot(snp_long_pip, aes(PIP, fill=Tool)) +
  geom_histogram(position="identity", alpha=0.55, bins=50, colour="white") +
  geom_vline(xintercept=0.5, linetype="dashed") +
  scale_fill_manual(values=c(FINEMAP="#534AB7", SuSiE="#D85A30")) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  labs(title="A3 · PIP distribution — FINEMAP vs SuSiE",
       subtitle="Most SNPs have low PIP; right tail shows fine-mapped candidates | dashed = 0.5 threshold",
       x="PIP", y="Number of SNPs")
print(p_a3)

# A4: Credible set PIP stacked bar --------------------------------------------
message("Plot A4: Credible set PIPs")
p_a4 <- ggplot(cs_snps, aes(reorder(SNP_rsID, -Composite_PIP), Composite_PIP)) +
  geom_col(aes(fill=Composite_PIP), width=0.7) +
  geom_errorbar(aes(ymin=SuSiE_PIP, ymax=FINEMAP_PIP), width=0.25, colour="grey30") +
  scale_fill_gradient(low="#9FE1CB", high="#0F6E56", name="Composite PIP") +
  labs(title="A4 · Credible set SNPs — composite PIP with FINEMAP/SuSiE range",
       subtitle="Bar = composite PIP | whiskers = range between FINEMAP PIP (top) and SuSiE PIP (bottom)",
       x="SNP", y="PIP") +
  theme(axis.text.x=element_text(angle=35, hjust=1, size=8))
print(p_a4)

# =============================================================================
# SECTION B  – GWAS ARCHITECTURE
# =============================================================================

# B1: Regional association (LocusZoom-style) ----------------------------------
message("Plot B1: Regional association")
r2_lead <- if(nrow(cs_snps)>0){
  lead_id <- cs_snps$SNP_rsID[which.max(cs_snps$Composite_PIP)]
  if(lead_id %in% colnames(r2_mat)){
    tibble(SNP_rsID=rownames(r2_mat), r2_lead=r2_mat[,lead_id])
  } else NULL
} else NULL

snp_reg <- snp
if(!is.null(r2_lead)) snp_reg <- left_join(snp, r2_lead, by="SNP_rsID")

r2_breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
r2_labels <- c("<0.2","0.2–0.4","0.4–0.6","0.6–0.8","≥0.8")
r2_palette <- c("#B4B2A9","#85B7EB","#378ADD","#EF9F27","#E24B4A")

if("r2_lead" %in% names(snp_reg)){
  snp_reg <- snp_reg %>%
    mutate(r2_band=cut(r2_lead, breaks=r2_breaks, labels=r2_labels, include.lowest=TRUE))
  p_b1 <- ggplot(snp_reg, aes(POS/1e6, GWAS_log10p)) +
    geom_point(aes(colour=r2_band,
                   shape=ifelse(in_CS,"Credible set","Other"),
                   size=ifelse(in_CS, 3, 1.8)), alpha=0.85) +
    geom_text_repel(data=filter(snp_reg, !is.na(label_top)),
                    aes(label=label_top), size=2.6, max.overlaps=15,
                    segment.colour="grey50", segment.size=0.3) +
    scale_colour_manual(values=setNames(r2_palette, r2_labels),
                        name=glue("r² with lead SNP\n{if(!is.null(r2_lead)) head(cs_snps$SNP_rsID[which.max(cs_snps$Composite_PIP)],1) else ''}"),
                        na.value="grey80") +
    scale_shape_manual(values=c("Credible set"=18,"Other"=16), name=NULL) +
    scale_size_identity() +
    labs(title="B1 · Regional association plot (LocusZoom-style)",
         subtitle="Color = LD r² with lead SNP | diamond = credible set member",
         x="Genomic position (Mb)", y="–log₁₀(p-value)")
} else {
  p_b1 <- ggplot(snp, aes(POS/1e6, GWAS_log10p, colour=pip_tier)) +
    geom_point(alpha=0.75) +
    scale_colour_manual(values=pip_cols, name="PIP tier") +
    labs(title="B1 · Regional association plot", x="Genomic position (Mb)", y="–log₁₀(p-value)")
}
print(p_b1)

# B2: PIP vs –log10(p) bubble (bubble = LD score) ----------------------------
message("Plot B2: PIP vs -log10p bubble")
p_b2 <- ggplot(snp, aes(Composite_PIP, GWAS_log10p)) +
  geom_point(aes(size=LD_score, colour=eqtl_sig_tier), alpha=0.7) +
  geom_vline(xintercept=0.5, linetype="dashed", colour="grey40") +
  geom_text_repel(data=filter(snp, !is.na(label_top)),
                  aes(label=label_top), size=2.6, max.overlaps=15) +
  scale_colour_manual(values=eqtl_cols, name="eQTL evidence") +
  scale_size_continuous(name="LD score", range=c(1,7)) +
  labs(title="B2 · Composite PIP vs GWAS –log₁₀(p) · bubble = LD score",
       subtitle="High PIP + moderate p = fine-mapped causal | large bubble = high LD background",
       x="Composite PIP", y="–log₁₀(GWAS p-value)")
print(p_b2)

# B3: MAF vs PIP coloured by OR -----------------------------------------------
message("Plot B3: MAF vs PIP")
p_b3 <- ggplot(snp, aes(AF_EA, Composite_PIP, colour=GWAS_OR)) +
  geom_point(aes(size=GWAS_log10p, shape=in_CS), alpha=0.8) +
  geom_vline(xintercept=0.01, linetype="dashed", colour="grey50", linewidth=0.6) +
  scale_colour_gradient2(low="#185FA5", mid="#EF9F27", high="#D85A30",
                         midpoint=1, name="GWAS OR") +
  scale_size_continuous(range=c(1,5), name="–log₁₀(p)") +
  scale_shape_manual(values=c(`TRUE`=18, `FALSE`=16), name="In credible set") +
  scale_x_log10() +
  annotation_logticks(sides="b") +
  labs(title="B3 · Effect allele frequency vs composite PIP — colored by GWAS OR",
       subtitle="Log x-axis | dashed = MAF 0.01 | rare high-effect SNPs (top-left) are particularly noteworthy",
       x="Effect allele frequency (log scale)", y="Composite PIP")
print(p_b3)

# B4: LD score vs PIP (slope diagnostic) -------------------------------------
message("Plot B4: LD score vs PIP")
lm_all <- lm(Composite_PIP ~ LD_score, data=snp)
lm_cs  <- if(nrow(cs_snps)>1) lm(Composite_PIP ~ LD_score, data=cs_snps) else NULL
p_b4 <- ggplot(snp, aes(LD_score, Composite_PIP)) +
  geom_point(aes(colour=pip_tier), alpha=0.5, size=1.5) +
  geom_smooth(method="lm", se=TRUE, colour="grey40", linetype="dashed",
              linewidth=0.8, fill="grey80") +
  {if(!is.null(lm_cs))
    geom_smooth(data=cs_snps, method="lm", se=TRUE, colour=cs_col,
                linewidth=1, fill="#FAECE7")} +
  geom_text_repel(data=filter(snp, !is.na(label_top)),
                  aes(label=label_top), size=2.6, max.overlaps=12) +
  scale_colour_manual(values=pip_cols, name="PIP tier") +
  annotate("text", x=max(snp$LD_score,na.rm=TRUE)*0.7,
           y=max(snp$Composite_PIP,na.rm=TRUE)*0.9,
           label=glue("All SNPs slope = {round(coef(lm_all)[2],5)}"),
           colour="grey40", size=3) +
  labs(title="B4 · LD score vs composite PIP — confounding diagnostic",
       subtitle="Flat/negative slope in credible set (orange) = fine-mapping corrected LD inflation | positive = residual confounding",
       x="LD score (Σr²)", y="Composite PIP")
print(p_b4)

# =============================================================================
# SECTION C  – eQTL INTEGRATION
# =============================================================================

# Build long eQTL table -------------------------------------------------------
eqtl_long <- snp %>%
  select(SNP_rsID, Composite_PIP, in_CS, pip_tier,
         all_of(paste0("pval_eQTL@", eqtl_datasets)),
         all_of(paste0("beta_eQTL@", eqtl_datasets))) %>%
  pivot_longer(
    cols = -c(SNP_rsID, Composite_PIP, in_CS, pip_tier),
    names_to = c(".value","tissue"),
    names_pattern = "(.+)@(.+)"
  ) %>%
  rename(pval=pval_eQTL, beta=beta_eQTL) %>%
  filter(!is.na(pval)) %>%
  mutate(log10p = -log10(pval),
         sig    = pval < 1e-5)

# C1: eQTL significance heatmap (tissue × top SNPs) --------------------------
message("Plot C1: eQTL heatmap")
top_n_snps <- snp %>% arrange(desc(Composite_PIP)) %>%
  slice_head(n=min(40, nrow(snp))) %>% pull(SNP_rsID)

snp_heat <- snp %>%
  filter(SNP_rsID %in% top_n_snps) %>%
  filter(!is.na(SNP_rsID), SNP_rsID != "") %>%
  arrange(desc(Composite_PIP)) %>%
  distinct(SNP_rsID, .keep_all = TRUE)

if (nrow(snp_heat) == 0) {
  stop("No valid SNP_rsID values available for C1 heatmap")
}

heat_dat <- snp_heat %>%
  select(SNP_rsID, all_of(paste0("log10p_eQTL@", eqtl_datasets))) %>%
  column_to_rownames("SNP_rsID") %>%
  rename_with(~str_remove(.,"log10p_eQTL@"))

heat_dat[is.na(heat_dat)] <- 0
heat_mat <- as.matrix(heat_dat)
row_ord  <- snp %>% filter(SNP_rsID %in% rownames(heat_mat)) %>%
  arrange(desc(Composite_PIP)) %>% pull(SNP_rsID)
heat_mat <- heat_mat[intersect(row_ord, rownames(heat_mat)),, drop=FALSE]

pip_annot <- snp %>% filter(SNP_rsID %in% rownames(heat_mat)) %>%
  select(SNP_rsID, Composite_PIP, in_CS) %>%
  mutate(CS = ifelse(in_CS,"Yes","No")) %>%
  column_to_rownames("SNP_rsID") %>%
  select(Composite_PIP, CS)

ann_colors <- list(CS=c(Yes=cs_col, No="grey85"),
                   Composite_PIP=colorRampPalette(c("white","#1D9E75"))(100))

pheatmap(heat_mat,
         color=colorRampPalette(c("white","#AFA9EC","#3C3489"))(100),
         cluster_rows=FALSE, cluster_cols=TRUE,
         annotation_row=pip_annot, annotation_colors=ann_colors,
         fontsize=7, fontsize_row=6, fontsize_col=7,
         angle_col=45, na_col="grey95",
         main="C1 · eQTL –log₁₀(p) heatmap — top 40 SNPs × all tissues\n(rows sorted by composite PIP | columns clustered by eQTL pattern)")

# C2: Number of significant eQTL tissues vs PIP -------------------------------
message("Plot C2: N eQTL tissues vs PIP")
p_c2 <- ggplot(snp, aes(n_eqtl_sig, Composite_PIP)) +
  geom_jitter(aes(colour=pip_tier, shape=in_CS), width=0.3, alpha=0.7, size=2) +
  geom_smooth(method="lm", se=TRUE, colour="grey30", fill="grey80") +
  stat_cor(method="pearson", label.x=max(snp$n_eqtl_sig,na.rm=TRUE)*0.6,
           label.y=max(snp$Composite_PIP,na.rm=TRUE)*0.9) +
  scale_colour_manual(values=pip_cols, name="PIP tier") +
  scale_shape_manual(values=c(`TRUE`=18,`FALSE`=16), name="In credible set") +
  labs(title="C2 · Number of tissues with significant eQTL (p<1e-5) vs composite PIP",
       subtitle="SNPs with broad tissue eQTL signal tend to have higher PIP | Pearson r shown",
       x="Number of tissues with eQTL p < 1e-5", y="Composite PIP")
print(p_c2)

# C3: Best eQTL –log10p vs composite PIP --------------------------------------
message("Plot C3: Min eQTL p vs PIP")
p_c3 <- ggplot(filter(snp, !is.na(log10_min_eqtl)),
               aes(log10_min_eqtl, Composite_PIP)) +
  geom_point(aes(colour=eqtl_sig_tier, size=GWAS_log10p, shape=in_CS), alpha=0.8) +
  geom_smooth(method="lm", se=TRUE, colour="grey30") +
  stat_cor(method="pearson") +
  geom_vline(xintercept=5, linetype="dashed", colour="grey50") +
  geom_vline(xintercept=10, linetype="dotted", colour="grey50") +
  geom_text_repel(data=filter(snp, !is.na(label_top) & !is.na(log10_min_eqtl)),
                  aes(label=label_top), size=2.6, max.overlaps=12) +
  scale_colour_manual(values=eqtl_cols, name="eQTL tier") +
  scale_shape_manual(values=c(`TRUE`=18,`FALSE`=16), name="In credible set") +
  scale_size_continuous(range=c(1,5), name="–log₁₀(GWAS p)") +
  labs(title="C3 · Strongest eQTL signal (any tissue) vs composite PIP",
       subtitle="Dashed = p<1e-5 | dotted = p<1e-10 | SNPs in top-right corner are both causal and functionally supported",
       x="–log₁₀(best eQTL p-value, any tissue)", y="Composite PIP")
print(p_c3)

# C4: eQTL beta vs GWAS beta per tissue (faceted, top SNPs only) -------------
message("Plot C4: eQTL beta vs GWAS beta per tissue")
gwas_beta_join <- snp %>% select(SNP_rsID, GWAS_beta, in_CS)
eqtl_beta_long <- snp %>%
  select(SNP_rsID, all_of(paste0("beta_eQTL@", eqtl_datasets))) %>%
  pivot_longer(-SNP_rsID, names_to="tissue", values_to="eqtl_beta",
               names_prefix="beta_eQTL@") %>%
  filter(!is.na(eqtl_beta)) %>%
  left_join(gwas_beta_join, by="SNP_rsID")

# Pick tissues with most data
tissue_n <- eqtl_beta_long %>% count(tissue, sort=TRUE) %>% slice_head(n=9)
eqtl_beta_top <- eqtl_beta_long %>%
  filter(tissue %in% tissue_n$tissue)

p_c4 <- ggplot(eqtl_beta_top, aes(eqtl_beta, GWAS_beta)) +
  geom_point(aes(colour=in_CS), alpha=0.5, size=1.5) +
  geom_smooth(method="lm", se=TRUE, colour="#D85A30", fill="#FAECE7", linewidth=0.8) +
  stat_cor(method="pearson", size=2.8) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey60") +
  geom_vline(xintercept=0, linetype="dashed", colour="grey60") +
  facet_wrap(~tissue, scales="free_x", ncol=3) +
  scale_colour_manual(values=c(`TRUE`=cs_col,`FALSE`="grey60"), name="In CS") +
  labs(title="C4 · eQTL β vs GWAS β per tissue (top 9 tissues by data coverage)",
       subtitle="Positive correlation → risk allele increases expression | Pearson r per tissue",
       x="eQTL effect size (β)", y="GWAS β") +
  theme(strip.text=element_text(size=7))
print(p_c4)

# C5: Per-tissue eQTL significance violin for CS vs non-CS -------------------
message("Plot C5: CS vs non-CS eQTL per tissue")
eqtl_cs_comp <- eqtl_long %>%
  filter(tissue %in% tissue_n$tissue) %>%
  mutate(group=ifelse(in_CS,"Credible set","Other SNPs"))

p_c5 <- ggplot(eqtl_cs_comp, aes(group, log10p, fill=group)) +
  geom_violin(alpha=0.6, colour=NA) +
  geom_boxplot(width=0.12, outlier.size=0.5, fill="white") +
  facet_wrap(~tissue, ncol=3, scales="free_y") +
  stat_compare_means(method="wilcox.test", size=2.8,
                     label.y.npc=0.95) +
  scale_fill_manual(values=c("Credible set"=cs_col,"Other SNPs"="grey70"), guide="none") +
  labs(title="C5 · eQTL –log₁₀(p) distribution — credible set vs other SNPs per tissue",
       subtitle="Wilcoxon p shown | credible set SNPs should show stronger eQTL signal if functionally relevant",
       x=NULL, y="–log₁₀(eQTL p-value)") +
  theme(strip.text=element_text(size=7), axis.text.x=element_text(size=8))
print(p_c5)

# =============================================================================
# SECTION D  – LD MATRICES
# =============================================================================

# D1: r² heatmap --------------------------------------------------------------
message("Plot D1: r2 heatmap")
pip_row_annot <- data.frame(
  Composite_PIP = cs_snps$Composite_PIP[match(rownames(r2_mat), cs_snps$SNP_rsID)],
  row.names = rownames(r2_mat)
)
pheatmap(r2_mat,
         color=colorRampPalette(c("white","#FAEEDA","#EF9F27","#633806"))(100),
         display_numbers=TRUE, number_format="%.2f", fontsize_number=7,
         annotation_row=pip_row_annot,
         annotation_colors=list(Composite_PIP=colorRampPalette(c("white","#1D9E75"))(100)),
         cluster_rows=TRUE, cluster_cols=TRUE,
         fontsize=8, angle_col=45,
         main="D1 · Pairwise LD r² — credible set SNPs\n(low off-diagonal r² = independent causal signals)")

# D2: r heatmap (signed) ------------------------------------------------------
message("Plot D2: r heatmap")
pheatmap(r_mat,
         color=colorRampPalette(c("#185FA5","white","#D85A30"))(100),
         display_numbers=TRUE, number_format="%.2f", fontsize_number=7,
         cluster_rows=TRUE, cluster_cols=TRUE,
         breaks=seq(-1, 1, length.out=101),
         fontsize=8, angle_col=45,
         main="D2 · Pairwise signed r — credible set SNPs\n(sign indicates allelic phase between SNP pairs)")

# D3: D' heatmap --------------------------------------------------------------
message("Plot D3: D-prime heatmap")
pheatmap(dprime_mat,
         color=colorRampPalette(c("#185FA5","white","#D85A30"))(100),
         display_numbers=TRUE, number_format="%.2f", fontsize_number=7,
         cluster_rows=TRUE, cluster_cols=TRUE,
         breaks=seq(-1, 1, length.out=101),
         fontsize=8, angle_col=45,
         main="D3 · D' (normalized LD) — credible set SNPs\n(|D'|=1 = complete LD | values <1 indicate historical recombination)")

# D4: LD score distribution by PIP tier ---------------------------------------
message("Plot D4: LD score by PIP tier")
p_d4 <- ggplot(snp, aes(LD_score, fill=pip_tier)) +
  geom_density(alpha=0.55, colour=NA) +
  geom_rug(data=cs_snps, aes(LD_score), colour=cs_col, alpha=0.6, linewidth=0.4) +
  scale_fill_manual(values=pip_cols, name="PIP tier") +
  labs(title="D4 · LD score distribution by PIP tier",
       subtitle="High-PIP SNPs at lower LD scores = fine-mapping deprioritised high-LD tags | rug = credible set members",
       x="LD score (Σr²)", y="Density")
print(p_d4)

# D5: LD score distribution eQTL sig vs not -----------------------------------
message("Plot D5: LD score eQTL groups")
snp_eqtl_grp <- snp %>%
  mutate(eqtl_group=case_when(
    min_eqtl_pval < 1e-10 ~ "Strong eQTL",
    min_eqtl_pval < 1e-5  ~ "Moderate eQTL",
    !is.na(min_eqtl_pval) ~ "No sig eQTL",
    TRUE                  ~ "No data"
  ) %>% factor(levels=c("Strong eQTL","Moderate eQTL","No sig eQTL","No data")))

p_d5 <- ggplot(filter(snp_eqtl_grp, eqtl_group!="No data"),
               aes(LD_score, fill=eqtl_group)) +
  geom_density(alpha=0.55, colour=NA) +
  scale_fill_manual(values=c("Strong eQTL"="#534AB7","Moderate eQTL"="#378ADD",
                              "No sig eQTL"="#B4B2A9"), name="eQTL group") +
  labs(title="D5 · LD score distribution — eQTL significance groups",
       subtitle="Lower LD score in eQTL-significant SNPs = independent functional signal, not just LD hitchhiking",
       x="LD score", y="Density")
print(p_d5)

# =============================================================================
# SECTION E  – STATISTICAL MODELS
# =============================================================================

# E1: LASSO / linear model — feature importance -------------------------------
message("Plot E1: Feature importance")
model_df <- snp %>%
  select(Composite_PIP, GWAS_log10p, GWAS_OR, GWAS_SE,
         LD_score, AF_EA, log10_min_eqtl, n_eqtl_sig) %>%
  drop_na()

# Standardise predictors
model_df_sc <- model_df %>%
  mutate(across(-Composite_PIP, scale))

lm_full <- lm(Composite_PIP ~ ., data=model_df_sc)
lm_summ <- summary(lm_full)
coef_df <- as.data.frame(lm_summ$coefficients) %>%
  rownames_to_column("Feature") %>%
  filter(Feature != "(Intercept)") %>%
  rename(beta=Estimate, se=`Std. Error`, p=`Pr(>|t|)`) %>%
  mutate(sig = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                         p < 0.05 ~ "*", TRUE ~ ""),
         ci_lo=beta - 1.96*se, ci_hi=beta + 1.96*se,
         direction=ifelse(beta>0,"Positive","Negative"))

p_e1 <- ggplot(coef_df, aes(reorder(Feature, beta), beta, fill=direction)) +
  geom_col(alpha=0.85, width=0.6) +
  geom_errorbar(aes(ymin=ci_lo, ymax=ci_hi), width=0.2, colour="grey30") +
  geom_text(aes(label=sig, y=ci_hi+0.001), size=4, colour="grey20") +
  geom_hline(yintercept=0, colour="grey30") +
  coord_flip() +
  scale_fill_manual(values=c(Positive="#1D9E75",Negative="#D85A30"), name="Direction") +
  labs(title="E1 · Multiple linear regression — predictors of composite PIP",
       subtitle=glue("Standardised β ± 95% CI | R² = {round(lm_summ$r.squared,3)} | * p<0.05 ** p<0.01 *** p<0.001"),
       x=NULL, y="Standardised β coefficient")
print(p_e1)

# E2: Partial regression plots ------------------------------------------------
message("Plot E2: Partial regression plots")
predictors <- c("GWAS_log10p","LD_score","log10_min_eqtl","n_eqtl_sig","AF_EA","GWAS_OR")
pred_plots <- map(predictors, function(prd){
  dat_sub <- snp %>% select(all_of(c("Composite_PIP","in_CS","pip_tier", prd))) %>% drop_na()
  ggplot(dat_sub, aes(.data[[prd]], Composite_PIP)) +
    geom_point(aes(colour=pip_tier, shape=in_CS), alpha=0.6, size=1.5) +
    geom_smooth(method="lm", se=TRUE, colour="grey30", fill="grey80", linewidth=0.7) +
    stat_cor(method="pearson", size=2.8) +
    scale_colour_manual(values=pip_cols, guide="none") +
    scale_shape_manual(values=c(`TRUE`=18,`FALSE`=16), guide="none") +
    labs(x=prd, y="Composite PIP", title=prd) +
    theme(plot.title=element_text(size=9, face="bold"))
})
print(wrap_plots(pred_plots, ncol=3) +
      plot_annotation(title="E2 · Partial regression — each predictor vs composite PIP",
                      subtitle="Pearson r shown | orange diamond = credible set member"))

# E3: Correlation matrix of all numeric features ------------------------------
message("Plot E3: Feature correlation matrix")
feat_cor_df <- snp %>%
  select(Composite_PIP, FINEMAP_PIP, SuSiE_PIP, GWAS_log10p, GWAS_OR,
         GWAS_SE, LD_score, AF_EA, log10_min_eqtl, n_eqtl_sig) %>%
  drop_na()

cor_mat <- cor(feat_cor_df, method="spearman")
p_e3 <- corrplot(cor_mat, method="color", type="upper", order="hclust",
                 col=colorRampPalette(c("#185FA5","white","#D85A30"))(200),
                 addCoef.col="black", number.cex=0.7, tl.cex=0.75,
                 tl.col="black", tl.srt=45, cl.cex=0.7,
                 title="E3 · Spearman correlation matrix — all numeric features\n",
                 mar=c(0,0,3,0))

# E4: PIP × eQTL joint scoring plot ------------------------------------------
message("Plot E4: Joint score")
snp_score <- snp %>%
  filter(!is.na(log10_min_eqtl)) %>%
  mutate(
    pip_sc    = (Composite_PIP - min(Composite_PIP,na.rm=TRUE)) /
                (max(Composite_PIP,na.rm=TRUE) - min(Composite_PIP,na.rm=TRUE)),
    eqtl_sc   = (log10_min_eqtl - min(log10_min_eqtl,na.rm=TRUE)) /
                (max(log10_min_eqtl,na.rm=TRUE) - min(log10_min_eqtl,na.rm=TRUE)),
    joint_score = (pip_sc + eqtl_sc) / 2
  )
p_e4 <- ggplot(snp_score, aes(pip_sc, eqtl_sc)) +
  geom_point(aes(colour=joint_score, size=joint_score, shape=in_CS), alpha=0.85) +
  geom_text_repel(data=filter(snp_score, !is.na(joint_score) & joint_score > quantile(joint_score, 0.92, na.rm=TRUE)),
                  aes(label=SNP_rsID), size=2.6, max.overlaps=15) +
  scale_colour_viridis_c(option="D", name="Joint\nscore") +
  scale_size_continuous(range=c(1,6), name="Joint\nscore") +
  scale_shape_manual(values=c(`TRUE`=18,`FALSE`=16), name="In CS") +
  geom_hline(yintercept=0.5, linetype="dashed", colour="grey50") +
  geom_vline(xintercept=0.5, linetype="dashed", colour="grey50") +
  labs(title="E4 · Joint causal score — normalised PIP × normalised eQTL signal",
       subtitle="Top-right quadrant = high fine-mapping confidence AND strong functional support",
       x="Normalised composite PIP", y="Normalised best eQTL –log₁₀(p)")
print(p_e4)

# E5: OR winner's curse: LD score vs OR ---------------------------------------
message("Plot E5: LD score vs OR")
p_e5 <- ggplot(snp, aes(LD_score, GWAS_OR)) +
  geom_point(aes(colour=pip_tier, shape=in_CS), alpha=0.7, size=2) +
  geom_smooth(method="lm", se=TRUE, colour="grey30", fill="grey80") +
  stat_cor(method="pearson") +
  geom_hline(yintercept=1, linetype="dashed", colour="grey50") +
  geom_text_repel(data=filter(snp, !is.na(label_top)),
                  aes(label=label_top), size=2.6, max.overlaps=12) +
  scale_colour_manual(values=pip_cols, name="PIP tier") +
  scale_shape_manual(values=c(`TRUE`=18,`FALSE`=16), name="In CS") +
  labs(title="E5 · LD score vs GWAS OR — winner's curse diagnostic",
       subtitle="Positive slope = OR inflated by LD background | use slope to apply shrinkage correction",
       x="LD score (Σr²)", y="GWAS OR")
print(p_e5)

dev.off()
message("✔  all_plots.pdf written")

# =============================================================================
# SECTION F  – EXCEL SUMMARY TABLE
# =============================================================================
message("Building Excel workbook...")

wb <- createWorkbook()
modifyBaseFont(wb, fontName="Arial", fontSize=10)

# Styles
h_style  <- createStyle(fontName="Arial", fontSize=10, fontColour="white",
                         fgFill="#2C4770", halign="center", valign="center",
                         textDecoration="bold", wrapText=TRUE, border="Bottom",
                         borderColour="white")
body_ev  <- createStyle(fontName="Arial", fontSize=9, fgFill="#F7F8FC")
body_od  <- createStyle(fontName="Arial", fontSize=9, fgFill="white")
hi_style <- createStyle(fontName="Arial", fontSize=9, fontColour="#0F6E56",
                         textDecoration="bold")
cs_style <- createStyle(fontName="Arial", fontSize=9, fgFill="#E1F5EE",
                         textDecoration="bold")
num4     <- createStyle(numFmt="0.0000")
num6     <- createStyle(numFmt="0.000000")
sci_fmt  <- createStyle(numFmt="0.00E+00")

add_sheet_data <- function(wb, name, df, freeze_col=1){
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle=h_style)
  n_rows <- nrow(df)
  n_cols <- ncol(df)
  # Batch alternating row styles (2 calls instead of n calls)
  even_rows <- seq(2, n_rows + 1, by = 2)
  odd_rows  <- seq(3, n_rows + 1, by = 2)
  if (length(even_rows) > 0)
    addStyle(wb, name, body_ev, rows = even_rows, cols = seq_len(n_cols), gridExpand = TRUE, stack = TRUE)
  if (length(odd_rows) > 0)
    addStyle(wb, name, body_od, rows = odd_rows,  cols = seq_len(n_cols), gridExpand = TRUE, stack = TRUE)
  freezePane(wb, name, firstActiveCol = freeze_col + 1)
  setColWidths(wb, name, cols = seq_len(n_cols), widths = 14)
  invisible(NULL)
}

# ── Sheet 1: Master table (top 100 by composite PIP) ──────────────────────────
master_out <- snp %>%
  arrange(desc(Composite_PIP)) %>%
  slice_head(n=100) %>%
  select(SNP_rsID, CHR, POS, EA, OA,
         Composite_PIP, FINEMAP_PIP, SuSiE_PIP,
         CS_membership, CS_size,
         GWAS_pval, GWAS_OR, GWAS_beta, GWAS_log10p, GWAS_chi2, GWAS_SE,
         GWAS_CI95_lower, GWAS_CI95_upper,
         LD_score, AF_EA, AF_OA,
         n_eqtl_sig, log10_min_eqtl)
add_sheet_data(wb, "Top100_Master", master_out, freeze_col=5)

# Highlight CS rows
cs_rows <- which(master_out$CS_membership == "1.0") + 1
if(length(cs_rows)>0)
  addStyle(wb, "Top100_Master", cs_style, rows=cs_rows,
           cols=seq_len(ncol(master_out)), gridExpand=TRUE, stack=TRUE)

# ── Sheet 2: Credible set detail ──────────────────────────────────────────────
cs_detail <- snp %>%
  filter(in_CS) %>%
  arrange(desc(Composite_PIP)) %>%
  select(SNP_rsID, CHR, POS, EA, OA,
         FINEMAP_PIP, SuSiE_PIP, Composite_PIP,
         CS_membership, CS_size,
         GWAS_pval, GWAS_OR, GWAS_beta, GWAS_log10p, GWAS_chi2, GWAS_SE,
         GWAS_CI95_lower, GWAS_CI95_upper,
         LD_score, AF_EA, AF_OA,
         n_eqtl_sig, log10_min_eqtl,
         all_of(paste0("pval_eQTL@", eqtl_datasets)),
         all_of(paste0("beta_eQTL@", eqtl_datasets)),
         all_of(paste0("OR_eQTL@",   eqtl_datasets)))
add_sheet_data(wb, "CredibleSet_Detail", cs_detail, freeze_col=5)

# ── Sheet 3: eQTL summary across tissues (credible set SNPs) ──────────────────
eqtl_summary <- eqtl_long %>%
  filter(in_CS) %>%
  group_by(SNP_rsID, tissue) %>%
  summarise(eQTL_pval=min(pval,na.rm=TRUE),
            eQTL_log10p=max(log10p,na.rm=TRUE),
            eQTL_beta=mean(beta,na.rm=TRUE),
            .groups="drop") %>%
  left_join(select(snp, SNP_rsID, Composite_PIP, GWAS_OR, GWAS_beta), by="SNP_rsID") %>%
  mutate(sig_flag=ifelse(eQTL_pval < 1e-5,"*",""),
         direction_concordant=ifelse(sign(eQTL_beta)==sign(GWAS_beta),"Yes","No")) %>%
  arrange(SNP_rsID, eQTL_pval)
add_sheet_data(wb, "eQTL_Tissue_Summary", eqtl_summary, freeze_col=2)

# ── Sheet 4: Feature correlation matrix ───────────────────────────────────────
cor_out <- as.data.frame(round(cor_mat, 4)) %>%
  rownames_to_column("Feature")
add_sheet_data(wb, "Feature_Correlation", cor_out, freeze_col=1)
# Add conditional formatting (3-color scale)
conditionalFormatting(wb, "Feature_Correlation",
  cols=2:(ncol(cor_out)), rows=2:(nrow(cor_out)+1),
  style=c("#185FA5","white","#D85A30"),
  type="colourScale")

# ── Sheet 5: Regression model summary ─────────────────────────────────────────
reg_out <- coef_df %>%
  select(Feature, beta, se, ci_lo, ci_hi, p, sig) %>%
  rename(Std_beta=beta, Std_SE=se, CI95_lower=ci_lo, CI95_upper=ci_hi,
         p_value=p, Significance=sig) %>%
  mutate(across(c(Std_beta,Std_SE,CI95_lower,CI95_upper), ~round(.x,5)),
         p_value=signif(p_value,4)) %>%
  arrange(p_value)

# Add model fit row
model_fit <- tibble(Feature="MODEL FIT",
                    Std_beta=round(lm_summ$r.squared,4),
                    Std_SE=round(lm_summ$adj.r.squared,4),
                    CI95_lower=NA, CI95_upper=NA,
                    p_value=signif(pf(lm_summ$fstatistic[1],
                                      lm_summ$fstatistic[2],
                                      lm_summ$fstatistic[3], lower.tail=FALSE),4),
                    Significance="R² | Adj-R² | F-test p")
reg_out <- bind_rows(reg_out, model_fit)
add_sheet_data(wb, "Regression_Model", reg_out, freeze_col=1)

# ── Sheet 6: LD matrix r² ─────────────────────────────────────────────────────
r2_out <- as.data.frame(round(r2_mat,4)) %>% rownames_to_column("SNP_rsID")
add_sheet_data(wb, "LD_matrix_r2", r2_out, freeze_col=1)
conditionalFormatting(wb, "LD_matrix_r2",
  cols=2:(ncol(r2_out)), rows=2:(nrow(r2_out)+1),
  style=c("white","#EF9F27","#633806"), type="colourScale")

# ── Sheet 7: LD matrix r ──────────────────────────────────────────────────────
r_out <- as.data.frame(round(r_mat,4)) %>% rownames_to_column("SNP_rsID")
add_sheet_data(wb, "LD_matrix_r", r_out, freeze_col=1)
conditionalFormatting(wb, "LD_matrix_r",
  cols=2:(ncol(r_out)), rows=2:(nrow(r_out)+1),
  style=c("#185FA5","white","#D85A30"), type="colourScale")

# ── Sheet 8: LD matrix D' ─────────────────────────────────────────────────────
dp_out <- as.data.frame(round(dprime_mat,4)) %>% rownames_to_column("SNP_rsID")
add_sheet_data(wb, "LD_matrix_Dprime", dp_out, freeze_col=1)
conditionalFormatting(wb, "LD_matrix_Dprime",
  cols=2:(ncol(dp_out)), rows=2:(nrow(dp_out)+1),
  style=c("#185FA5","white","#D85A30"), type="colourScale")

# ── Sheet 9: Joint score ranking ──────────────────────────────────────────────
joint_out <- snp_score %>%
  arrange(desc(joint_score)) %>%
  select(SNP_rsID, CHR, POS, EA, OA,
         FINEMAP_PIP, SuSiE_PIP, Composite_PIP,
         pip_sc, eqtl_sc, joint_score,
         n_eqtl_sig, log10_min_eqtl,
         GWAS_pval, GWAS_OR, GWAS_log10p,
         LD_score, in_CS, CS_membership) %>%
  rename(PIP_score_norm=pip_sc, eQTL_score_norm=eqtl_sc, Joint_causal_score=joint_score)
add_sheet_data(wb, "Joint_Causal_Score", joint_out, freeze_col=5)

# ── Sheet 10: README ──────────────────────────────────────────────────────────
addWorksheet(wb, "README")
readme_txt <- c(
  "CAUSAL SNP IDENTIFICATION — ANALYSIS SUMMARY",
  "",
  "SHEETS:",
  "1. Top100_Master        — Top 100 SNPs by composite PIP with all key metrics",
  "2. CredibleSet_Detail   — All credible set SNPs with full eQTL data per tissue",
  "3. eQTL_Tissue_Summary  — eQTL p-value, beta, direction per tissue for CS SNPs",
  "4. Feature_Correlation  — Spearman correlation matrix (colour-coded) among all numeric features",
  "5. Regression_Model     — Standardised linear regression: Composite_PIP ~ all predictors",
  "6. LD_matrix_r2         — Pairwise r² among credible set SNPs (colour: white→dark amber)",
  "7. LD_matrix_r          — Signed r (blue=negative, red=positive phase)",
  "8. LD_matrix_Dprime     — D' normalised LD (blue=neg, red=pos; <1 = historical recombination)",
  "9. Joint_Causal_Score   — Combined ranking: normalised PIP + normalised eQTL signal",
  "10. README              — This sheet",
  "",
  "KEY METRICS:",
  "Composite_PIP    = (FINEMAP_PIP + SuSiE_PIP) / 2",
  "Joint_score      = (norm_PIP + norm_eQTL_log10p) / 2",
  "n_eqtl_sig       = number of tissues with eQTL p < 1e-5",
  "log10_min_eqtl   = -log10(best eQTL p-value across all tissues)",
  "",
  "DECISION CRITERIA FOR CAUSALITY:",
  "Tier 1 (high confidence) : Composite_PIP > 0.5  AND  n_eqtl_sig >= 3  AND  in credible set",
  "Tier 2 (moderate)        : Composite_PIP > 0.1  AND  n_eqtl_sig >= 1",
  "Tier 3 (GWAS support)    : In credible set but no eQTL signal (possible non-coding or rare variant)",
  "",
  "PLOTS (all_plots.pdf):",
  "A1  FINEMAP vs SuSiE PIP scatter",
  "A2  Composite PIP lollipop (top 30)",
  "A3  PIP distribution histogram",
  "A4  Credible set PIPs with tool range",
  "B1  Regional association (LocusZoom-style)",
  "B2  PIP vs -log10p bubble (LD score = size)",
  "B3  MAF vs PIP colored by OR",
  "B4  LD score vs PIP (slope diagnostic)",
  "C1  eQTL heatmap: tissue x top-40 SNPs",
  "C2  N significant eQTL tissues vs PIP",
  "C3  Best eQTL -log10p vs PIP",
  "C4  eQTL beta vs GWAS beta per tissue",
  "C5  CS vs non-CS eQTL per tissue (Wilcoxon)",
  "D1  r² heatmap (credible set)",
  "D2  Signed r heatmap (credible set)",
  "D3  D-prime heatmap (credible set)",
  "D4  LD score density by PIP tier",
  "D5  LD score density by eQTL group",
  "E1  Feature importance (linear regression)",
  "E2  Partial regression plots",
  "E3  Spearman feature correlation matrix",
  "E4  Joint causal score scatter",
  "E5  LD score vs OR (winner's curse)"
)
writeData(wb, "README", data.frame(Notes=readme_txt), colNames=FALSE)
setColWidths(wb, "README", cols=1, widths=80)
addStyle(wb, "README", createStyle(fontName="Arial", fontSize=10,
         textDecoration="bold", fontColour="#2C4770"),
         rows=1, cols=1)

saveWorkbook(wb, out_xlsx, overwrite=TRUE)
message(glue("✔  {out_xlsx} written"))
message("\n=== DONE ===")
message(glue("Outputs: {out_pdf} | {out_xlsx}"))

if (nzchar(done_file)) {
  writeLines("ok", done_file)
}
