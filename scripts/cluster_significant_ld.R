#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(igraph)
  library(ggplot2)
})

option_list <- list(
  make_option("--matched-tsv", type = "character"),
  make_option("--ld-matrix", type = "character"),
  make_option("--loci-file", type = "character"),
  make_option("--locus", type = "character"),
  make_option("--p-threshold", type = "double", default = 5e-8),
  make_option("--r2-threshold", type = "double", default = 0.6),
  make_option("--distance-kb", type = "double", default = 50),
  make_option("--out-summary", type = "character"),
  make_option("--out-html", type = "character"),
  make_option("--out-diag", type = "character"),
  make_option("--out-dir", type = "character"),
  make_option("--done-file", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("matched-tsv", "ld-matrix", "loci-file", "locus", "out-summary", "out-html", "out-diag", "out-dir", "done-file")
for (k in required) {
  if (is.null(opt[[k]]) || !nzchar(opt[[k]])) {
    stop(paste("Missing required argument --", k, sep = ""))
  }
}

safe_mkdir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

safe_mkdir(opt$`out-summary`)
safe_mkdir(opt$`out-html`)
safe_mkdir(opt$`out-diag`)
safe_mkdir(opt$`done-file`)
write_interactive_plot <- function(plot_df, locus_name, p_threshold, out_html, meta = list()) {
  plot_df <- as.data.frame(plot_df)
  plot_df$bp_mb <- as.numeric(plot_df$bp) / 1e6
  plot_df$logp <- -log10(pmax(pmin(as.numeric(plot_df$p), 1), 1e-300))
  plot_df$cluster_label <- as.character(plot_df$cluster_label)

  uniq <- unique(plot_df$cluster_label)
  uniq <- uniq[order(uniq)]

  palette <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
    "#e377c2", "#17becf", "#bcbd22", "#7f7f7f", "#4daf4a", "#984ea3"
  )

  traces <- list()
  traces[[length(traces) + 1]] <- list(
    x = plot_df$bp_mb,
    y = plot_df$logp,
    mode = "markers",
    type = "scattergl",
    name = "All SNPs (background)",
    marker = list(color = "rgba(130,130,130,0.35)", size = 5),
    hovertemplate = paste0(
      "SNP: %{customdata[0]}<br>",
      "CHR: %{customdata[1]}<br>",
      "BP: %{customdata[2]}<br>",
      "P: %{customdata[3]:.3e}<br>",
      "Cluster: %{customdata[4]}<extra></extra>"
    ),
    customdata = cbind(plot_df$snp, plot_df$chr, plot_df$bp, plot_df$p, plot_df$cluster_label),
    showlegend = TRUE
  )

  cid <- 0
  for (lab in uniq) {
    if (lab == "background") next
    cid <- cid + 1
    d <- plot_df[plot_df$cluster_label == lab, , drop = FALSE]
    if (nrow(d) == 0) next
    color <- palette[((cid - 1) %% length(palette)) + 1]
    traces[[length(traces) + 1]] <- list(
      x = d$bp_mb,
      y = d$logp,
      mode = "markers",
      type = "scattergl",
      name = sprintf("%s (n=%d)", lab, nrow(d)),
      marker = list(color = color, size = 8, line = list(color = "#111111", width = 0.4)),
      hovertemplate = paste0(
        "SNP: %{customdata[0]}<br>",
        "CHR: %{customdata[1]}<br>",
        "BP: %{customdata[2]}<br>",
        "P: %{customdata[3]:.3e}<br>",
        "Cluster: %{customdata[4]}<extra></extra>"
      ),
      customdata = cbind(d$snp, d$chr, d$bp, d$p, d$cluster_label),
      showlegend = TRUE
    )
  }

  threshold_line <- list(
    type = "line",
    xref = "paper", x0 = 0, x1 = 1,
    yref = "y", y0 = -log10(p_threshold), y1 = -log10(p_threshold),
    line = list(color = "#b91c1c", width = 1.2, dash = "dash")
  )

  layout <- list(
    title = list(text = sprintf("%s | Interactive Manhattan with LD Clusters", locus_name)),
    xaxis = list(title = "Position (Mb)"),
    yaxis = list(title = "-log10(P)"),
    hovermode = "closest",
    legend = list(orientation = "v", x = 1.02, y = 1),
    margin = list(l = 70, r = 260, t = 60, b = 60),
    shapes = list(threshold_line),
    plot_bgcolor = "#ffffff",
    paper_bgcolor = "#ffffff"
  )

  payload <- list(data = traces, layout = layout)
  payload_json <- jsonlite::toJSON(payload, auto_unbox = TRUE, digits = 12)
  locus_name_json <- jsonlite::toJSON(locus_name, auto_unbox = TRUE)

  # Build info card HTML from meta list
  fmt_bp <- function(x) if (is.na(x) || !is.finite(x)) "N/A" else format(as.integer(x), big.mark = ",")
  fmt_mb <- function(x) if (is.na(x) || !is.finite(x)) "N/A" else sprintf("%.3f Mb", x / 1e6)
  fmt_kb <- function(x) if (is.na(x) || !is.finite(x)) "N/A" else sprintf("%.1f kb", x / 1e3)
  fmt_sci <- function(x) if (is.na(x) || !is.finite(x)) "N/A" else format(x, scientific = TRUE, digits = 3)
  fmt_n  <- function(x) if (is.na(x) || !is.finite(x)) "N/A" else format(as.integer(x), big.mark = ",")

  locus_chr   <- if (!is.null(meta$chr))          meta$chr          else "N/A"
  locus_start <- if (!is.null(meta$locus_start))  meta$locus_start  else NA
  locus_end   <- if (!is.null(meta$locus_end))    meta$locus_end    else NA
  locus_width <- if (is.finite(locus_start) && is.finite(locus_end)) locus_end - locus_start else NA
  n_total     <- if (!is.null(meta$n_total))      meta$n_total      else nrow(plot_df)
  n_sig       <- if (!is.null(meta$n_sig))        meta$n_sig        else sum(plot_df$p <= p_threshold, na.rm = TRUE)
  n_clusters  <- if (!is.null(meta$n_clusters))   meta$n_clusters   else length(unique(plot_df$cluster_label[plot_df$cluster_label != "background"]))
  r2_thr      <- if (!is.null(meta$r2_threshold)) meta$r2_threshold else "N/A"
  dist_kb_val <- if (!is.null(meta$distance_kb))  meta$distance_kb  else "N/A"
  lead_snp    <- if (!is.null(meta$lead_snp))     meta$lead_snp     else "N/A"
  lead_p      <- if (!is.null(meta$lead_p))       meta$lead_p       else NA

  cards <- list(
    list(label = "Locus",              value = locus_name),
    list(label = "Chromosome",         value = paste0("chr", locus_chr)),
    list(label = "Start",              value = fmt_bp(locus_start)),
    list(label = "End",                value = fmt_bp(locus_end)),
    list(label = "Width",              value = fmt_kb(locus_width)),
    list(label = "Total SNPs",         value = fmt_n(n_total)),
    list(label = paste0("SNPs < ", format(p_threshold, scientific = TRUE)),
                                       value = fmt_n(n_sig)),
    list(label = "LD clusters",        value = fmt_n(n_clusters)),
    list(label = "R\u00b2 threshold",  value = as.character(r2_thr)),
    list(label = "Merge distance",     value = paste0(dist_kb_val, " kb")),
    list(label = "Lead SNP",           value = as.character(lead_snp)),
    list(label = "Lead P",             value = fmt_sci(lead_p))
  )

  card_html <- paste(sapply(cards, function(cc) {
    sprintf(
      "<div class='card'><span class='clabel'>%s</span><span class='cval'>%s</span></div>",
      cc$label, cc$value
    )
  }), collapse = "\n")

  html <- c(
    "<!doctype html>",
    "<html>",
    "<head>",
    "  <meta charset='utf-8'/>",
    "  <meta name='viewport' content='width=device-width, initial-scale=1' />",
    sprintf("  <title>%s clusters Manhattan</title>", locus_name),
    "  <script src='https://cdn.plot.ly/plotly-2.35.2.min.js'></script>",
    "  <style>",
    "    *{box-sizing:border-box}",
    "    body{margin:0;font-family:'Segoe UI',Arial,sans-serif;background:#f0f4f8;color:#101828}",
    "    .page{display:flex;flex-direction:column;height:100vh;padding:12px 16px;gap:10px}",
    "    .infobar{display:flex;flex-wrap:wrap;gap:8px;background:#fff;",
    "             border:1px solid #d0d5dd;border-radius:10px;padding:10px 14px}",
    "    .card{display:flex;flex-direction:column;min-width:110px;padding:4px 10px;",
    "          border-right:1px solid #e4e7ec}",
    "    .card:last-child{border-right:none}",
    "    .clabel{font-size:10px;font-weight:600;text-transform:uppercase;",
    "            letter-spacing:.05em;color:#667085}",
    "    .cval{font-size:13px;font-weight:700;color:#101828;margin-top:2px}",
    "    .mainrow{display:flex;flex:1;min-height:0;gap:10px}",
    "    .plotwrap{flex:1;min-height:0;border:1px solid #d0d5dd;border-radius:10px;",
    "             background:#fff;overflow:hidden}",
    "    .snpbox{width:350px;max-width:42vw;background:#fff;border:1px solid #d0d5dd;",
    "           border-radius:10px;padding:10px;display:flex;flex-direction:column;gap:8px}",
    "    .snpbox h3{margin:0;font-size:14px}",
    "    .snpgrid{display:grid;grid-template-columns:110px 1fr;gap:6px;font-size:12px}",
    "    .k{color:#667085;font-weight:600}",
    "    .v{font-weight:700;color:#101828;word-break:break-all}",
    "    .snptext{width:100%;height:170px;resize:vertical;border:1px solid #d0d5dd;",
    "            border-radius:8px;padding:8px;font-size:12px;font-family:monospace}",
    "    .btn{align-self:flex-start;border:1px solid #155eef;background:#155eef;color:#fff;",
    "         border-radius:8px;padding:6px 10px;font-size:12px;cursor:pointer}",
    "    .btn:active{transform:translateY(1px)}",
    "    .copymsg{font-size:11px;color:#027a48;min-height:14px}",
    "    @media (max-width: 1100px){.mainrow{flex-direction:column}.snpbox{width:100%;max-width:none}}",
    "    #plot{width:100%;height:100%}",
    "    .hint{font-size:11px;color:#667085;padding:2px 4px}",
    "  </style>",
    "</head>",
    "<body>",
    "  <div class='page'>",
    "    <div class='infobar'>",
    card_html,
    "    </div>",
    "    <div class='mainrow'>",
    "      <div class='plotwrap'><div id='plot'></div></div>",
    "      <div class='snpbox'>",
    "        <h3>Selected SNP Details</h3>",
    "        <div class='snpgrid'>",
    "          <div class='k'>SNP</div><div class='v' id='sel-snp'>Click a point</div>",
    "          <div class='k'>Chromosome</div><div class='v' id='sel-chr'>-</div>",
    "          <div class='k'>Position (bp)</div><div class='v' id='sel-bp'>-</div>",
    "          <div class='k'>P-value</div><div class='v' id='sel-p'>-</div>",
    "          <div class='k'>Cluster</div><div class='v' id='sel-cluster'>-</div>",
    "        </div>",
    "        <textarea id='sel-text' class='snptext' readonly>Click any SNP on the plot to populate copy-ready text.</textarea>",
    "        <button id='copy-btn' class='btn' type='button'>Copy SNP Info</button>",
    "        <div id='copy-msg' class='copymsg'></div>",
    "      </div>",
    "    </div>",
    "    <div class='hint'>&#9432; Click legend items to hide/show clusters &mdash; double-click to isolate one cluster.</div>",
    "  </div>",
    "  <script>",
    sprintf("    const payload = %s;", payload_json),
    sprintf("    const locusName = %s;", locus_name_json),
    "    const fmtInt = (x) => Number.isFinite(Number(x)) ? Number(x).toLocaleString() : 'N/A';",
    "    const fmtSci = (x) => Number.isFinite(Number(x)) ? Number(x).toExponential(3) : 'N/A';",
    "    const el = {",
    "      snp: document.getElementById('sel-snp'),",
    "      chr: document.getElementById('sel-chr'),",
    "      bp: document.getElementById('sel-bp'),",
    "      p: document.getElementById('sel-p'),",
    "      cluster: document.getElementById('sel-cluster'),",
    "      txt: document.getElementById('sel-text'),",
    "      copyBtn: document.getElementById('copy-btn'),",
    "      copyMsg: document.getElementById('copy-msg')",
    "    };",
    "    const updateSnpBox = (point) => {",
    "      const cd = (point && point.customdata) ? point.customdata : [];",
    "      const snp = (cd[0] ?? 'N/A').toString();",
    "      const chr = (cd[1] ?? 'N/A').toString();",
    "      const bp = cd[2];",
    "      const p = cd[3];",
    "      const cluster = (cd[4] ?? 'N/A').toString();",
    "      const chrText = chr === 'N/A' ? 'N/A' : ('chr' + chr);",
    "      el.snp.textContent = snp;",
    "      el.chr.textContent = chrText;",
    "      el.bp.textContent = fmtInt(bp);",
    "      el.p.textContent = fmtSci(p);",
    "      el.cluster.textContent = cluster;",
    "      el.txt.value = [",
    "        'Locus: ' + locusName,",
    "        'SNP: ' + snp,",
    "        'Chromosome: ' + chrText,",
    "        'Position (bp): ' + fmtInt(bp),",
    "        'P-value: ' + fmtSci(p),",
    "        'Cluster: ' + cluster",
    "      ].join('\\n');",
    "      el.copyMsg.textContent = '';",
    "    };",
    "    const copyText = async () => {",
    "      const txt = el.txt.value || '';",
    "      if (!txt.trim()) return;",
    "      try {",
    "        if (navigator.clipboard && navigator.clipboard.writeText) {",
    "          await navigator.clipboard.writeText(txt);",
    "        } else {",
    "          el.txt.focus();",
    "          el.txt.select();",
    "          document.execCommand('copy');",
    "        }",
    "        el.copyMsg.textContent = 'Copied to clipboard.';",
    "      } catch (e) {",
    "        el.copyMsg.textContent = 'Copy failed; select the text manually.';",
    "      }",
    "    };",
    "    el.copyBtn.addEventListener('click', copyText);",
    "    Plotly.newPlot('plot', payload.data, payload.layout, {responsive:true, displaylogo:false}).then(() => {",
    "      const plot = document.getElementById('plot');",
    "      plot.on('plotly_click', (ev) => {",
    "        if (!ev || !ev.points || !ev.points.length) return;",
    "        updateSnpBox(ev.points[0]);",
    "      });",
    "    });",
    "  </script>",
    "</body>",
    "</html>"
  )
  writeLines(html, out_html)
}

dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(opt$`out-dir`, "clusters")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading matched GWAS and LD matrix...")
gwas <- fread(opt$`matched-tsv`, sep = "\t", header = TRUE)
setnames(gwas, tolower(names(gwas)))

needed <- c("snp", "chr", "bp", "p")
miss <- setdiff(needed, names(gwas))
if (length(miss) > 0) {
  stop(sprintf("matched.tsv missing columns: %s", paste(miss, collapse = ",")))
}

gwas[, p := as.numeric(p)]
gwas[, bp := as.integer(bp)]
gwas <- gwas[is.finite(p) & p > 0 & p <= 1 & is.finite(bp)]

ld_raw <- fread(opt$`ld-matrix`, sep = "\t", header = TRUE)
ld_snps <- as.character(ld_raw[[1]])
ld_cols <- names(ld_raw)[-1]
ld_mat <- as.matrix(ld_raw[, -1, with = FALSE])
storage.mode(ld_mat) <- "numeric"
rownames(ld_mat) <- ld_snps
colnames(ld_mat) <- ld_cols

# Read locus boundaries for report columns
loci <- fread(opt$`loci-file`, sep = "\t", header = TRUE)
setnames(loci, tolower(names(loci)))
if (!all(c("loci", "chr", "start", "end") %in% names(loci))) {
  stop("loci file must contain loci, chr, start, end columns")
}

sig <- gwas[p <= opt$`p-threshold`]
sig <- sig[snp %in% rownames(ld_mat) & snp %in% colnames(ld_mat)]

if (nrow(sig) == 0) {
  empty <- data.frame(
    locus = character(0),
    cluster_id = integer(0),
    chr = integer(0),
    cluster_start = integer(0),
    cluster_end = integer(0),
    n_sig_snps = integer(0),
    lead_snp = character(0),
    lead_p = numeric(0),
    stringsAsFactors = FALSE
  )
  fwrite(empty, opt$`out-summary`, sep = "\t")
  base_df <- copy(gwas)
  base_df[, cluster_label := "background"]
  locus_idx_e  <- match(opt$locus, as.character(loci[["loci"]]))
  meta_empty <- list(
    chr         = if (!is.na(locus_idx_e)) loci[["chr"]][locus_idx_e]              else NA,
    locus_start = if (!is.na(locus_idx_e)) as.integer(loci[["start"]][locus_idx_e]) else NA,
    locus_end   = if (!is.na(locus_idx_e)) as.integer(loci[["end"]][locus_idx_e])   else NA,
    n_total     = nrow(base_df),
    n_sig       = 0L,
    n_clusters  = 0L,
    r2_threshold = opt$`r2-threshold`,
    distance_kb  = opt$`distance-kb`,
    lead_snp     = "none",
    lead_p       = NA_real_
  )
  write_interactive_plot(base_df, opt$locus, opt$`p-threshold`, opt$`out-html`, meta = meta_empty)
  diag <- list(
    locus = opt$locus,
    p_threshold = opt$`p-threshold`,
    r2_threshold = opt$`r2-threshold`,
    n_sig_snps = 0,
    n_clusters = 0,
    message = "No SNPs with P <= threshold in locus"
  )
  writeLines(jsonlite::toJSON(diag, auto_unbox = TRUE, pretty = TRUE), opt$`out-diag`)
  writeLines("ok", opt$`done-file`)
  quit(status = 0)
}

sig_snps <- as.character(sig$snp)
r2_sub <- (ld_mat[sig_snps, sig_snps, drop = FALSE])^2

# Build LD graph among significant SNPs using igraph
edge_idx <- which(r2_sub >= opt$`r2-threshold` & upper.tri(r2_sub), arr.ind = TRUE)
if (nrow(edge_idx) > 0) {
  edge_df <- data.frame(
    from = rownames(r2_sub)[edge_idx[, 1]],
    to = colnames(r2_sub)[edge_idx[, 2]],
    stringsAsFactors = FALSE
  )
  g <- graph_from_data_frame(edge_df, directed = FALSE, vertices = sig_snps)
} else {
  g <- make_empty_graph(n = length(sig_snps), directed = FALSE)
  g <- set_vertex_attr(g, "name", value = sig_snps)
}

# Build strict pairwise-LD clusters using maximal cliques.
# This ensures every SNP pair inside a reported cluster passes the r2 threshold.
cluster_members <- list()
if (length(V(g)) > 0) {
  comp <- components(g)
  comp_members <- split(names(comp$membership), comp$membership)

  for (m in comp_members) {
    sub_remaining <- induced_subgraph(g, vids = m)
    while (vcount(sub_remaining) > 0) {
      clq <- maximal.cliques(sub_remaining)
      if (length(clq) == 0) {
        # Isolated vertices: each is its own cluster
        verts <- as_ids(V(sub_remaining))
        cluster_members <- c(cluster_members, as.list(verts))
        break
      }

      # Prefer larger cliques; tie-break by best (lowest) GWAS p among clique SNPs
      clq_sizes <- sapply(clq, length)
      best_size <- max(clq_sizes)
      cand <- which(clq_sizes == best_size)

      if (length(cand) > 1) {
        best_idx <- cand[1]
        best_p <- Inf
        for (ix in cand) {
          s <- as_ids(clq[[ix]])
          pmin_val <- min(sig[snp %in% s]$p, na.rm = TRUE)
          if (is.finite(pmin_val) && pmin_val < best_p) {
            best_p <- pmin_val
            best_idx <- ix
          }
        }
      } else {
        best_idx <- cand[1]
      }

      chosen <- as_ids(clq[[best_idx]])
      cluster_members[[length(cluster_members) + 1]] <- chosen

      # Remove assigned SNPs and continue peeling the component.
      keep_verts <- setdiff(as_ids(V(sub_remaining)), chosen)
      if (length(keep_verts) == 0) {
        break
      }
      sub_remaining <- induced_subgraph(sub_remaining, vids = keep_verts)
    }
  }
}

n_clusters_ld_only <- length(cluster_members)

# Merge nearby clusters if genomic gap is within threshold (kb -> bp).
dist_bp <- as.numeric(opt$`distance-kb`) * 1000
if (!is.finite(dist_bp) || dist_bp < 0) {
  dist_bp <- 50000
}

if (length(cluster_members) > 1) {
  stats <- lapply(seq_along(cluster_members), function(i) {
    snps_i <- cluster_members[[i]]
    cl_i <- sig[snp %in% snps_i]
    data.frame(
      idx = i,
      chr = as.integer(cl_i$chr[1]),
      st = min(cl_i$bp),
      en = max(cl_i$bp),
      stringsAsFactors = FALSE
    )
  })
  st_df <- do.call(rbind, stats)
  st_df <- st_df[order(st_df$chr, st_df$st, st_df$en), , drop = FALSE]

  merged <- list()
  current <- cluster_members[[st_df$idx[1]]]
  cur_chr <- st_df$chr[1]
  cur_st <- st_df$st[1]
  cur_en <- st_df$en[1]

  if (nrow(st_df) > 1) {
    for (k in 2:nrow(st_df)) {
      next_idx <- st_df$idx[k]
      next_chr <- st_df$chr[k]
      next_st <- st_df$st[k]
      next_en <- st_df$en[k]

      # Merge if same chromosome and gap <= distance threshold.
      gap_bp <- next_st - cur_en
      if (next_chr == cur_chr && gap_bp <= dist_bp) {
        current <- unique(c(current, cluster_members[[next_idx]]))
        cur_en <- max(cur_en, next_en)
      } else {
        merged[[length(merged) + 1]] <- current
        current <- cluster_members[[next_idx]]
        cur_chr <- next_chr
        cur_st <- next_st
        cur_en <- next_en
      }
    }
  }
  merged[[length(merged) + 1]] <- current
  cluster_members <- merged
}

make_ld_triangle <- function(r2_matrix) {
  n_s <- nrow(r2_matrix)
  pair_idx <- which(row(r2_matrix) <= col(r2_matrix), arr.ind = TRUE)
  idx_i <- pair_idx[, 1]
  idx_j <- pair_idx[, 2]
  r2_v <- r2_matrix[pair_idx]
  n_p <- length(r2_v)
  xc <- (idx_i + idx_j) / 2
  yc <- (idx_j - idx_i) / 2
  data.frame(
    group = rep(seq_len(n_p), each = 4L),
    r2 = rep(r2_v, each = 4L),
    x = as.vector(rbind(xc - 0.5, xc, xc + 0.5, xc)),
    y = as.vector(rbind(yc, yc - 0.5, yc, yc + 0.5))
  )
}

cluster_rows <- list()
gwas_with_cluster <- copy(gwas)
gwas_with_cluster[, cluster_label := "background"]

for (cid in seq_along(cluster_members)) {
  cl_snps <- cluster_members[[cid]]
  cl <- sig[snp %in% cl_snps]
  setorder(cl, bp)

  lead_idx <- which.min(cl$p)
  lead_snp <- as.character(cl$snp[lead_idx])
  lead_p <- as.numeric(cl$p[lead_idx])
  chr_val <- as.integer(cl$chr[lead_idx])
  st <- min(cl$bp)
  en <- max(cl$bp)

  # Manhattan plot for this cluster against full locus background
  plot_df <- copy(gwas)
  plot_df[, logp := -log10(p)]
  plot_df[, grp := "background"]
  plot_df[snp %in% cl_snps, grp := "cluster"]
  plot_df[snp == lead_snp, grp := "lead"]

  p1 <- ggplot(plot_df, aes(x = bp / 1e6, y = logp, color = grp)) +
    geom_hline(yintercept = -log10(opt$`p-threshold`), linetype = "dashed", color = "red", linewidth = 0.5) +
    geom_point(size = 1.4, alpha = 0.85) +
    scale_color_manual(values = c(background = "grey75", cluster = "#1f78b4", lead = "#e31a1c")) +
    labs(
      title = sprintf("%s | Cluster %d Manhattan", opt$locus, cid),
      subtitle = sprintf("Lead SNP: %s (P=%s)", lead_snp, format(lead_p, scientific = TRUE, digits = 3)),
      x = sprintf("Position on chr%s (Mb)", chr_val),
      y = expression(-log[10](italic(P))),
      color = NULL
    ) +
    theme_bw(base_size = 11)

  ggsave(file.path(plot_dir, sprintf("cluster_%02d_manhattan.png", cid)), p1, width = 10, height = 4.2, dpi = 150)

  # Triangle LD plot for SNPs inside this cluster
  cl_ord <- cl$snp[order(cl$bp)]
  r2_cl <- r2_sub[cl_ord, cl_ord, drop = FALSE]
  poly <- make_ld_triangle(r2_cl)
  max_y <- (nrow(r2_cl) - 1) / 2 + 0.5

  p2 <- ggplot(poly, aes(x = x, y = y, group = group, fill = r2)) +
    geom_polygon(color = NA) +
    scale_fill_gradientn(colours = c("#FFFFFF", "#87CEFA", "#1E90FF", "#0047AB", "#B22222"), limits = c(0, 1)) +
    scale_y_reverse(limits = c(max_y, -0.5), expand = c(0, 0)) +
    labs(
      title = sprintf("%s | Cluster %d LD triangle", opt$locus, cid),
      subtitle = sprintf("%d significant SNPs, r2 >= %.2f edges", length(cl_ord), opt$`r2-threshold`),
      x = "Cluster SNP order by position",
      y = NULL,
      fill = expression(italic(r)^2)
    ) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

  ggsave(file.path(plot_dir, sprintf("cluster_%02d_ld_triangle.png", cid)), p2, width = 8, height = 4.8, dpi = 150)

  cluster_rows[[length(cluster_rows) + 1]] <- data.frame(
    locus = opt$locus,
    cluster_id = cid,
    chr = chr_val,
    cluster_start = st,
    cluster_end = en,
    n_sig_snps = nrow(cl),
    lead_snp = lead_snp,
    lead_p = lead_p,
    stringsAsFactors = FALSE
  )

  gwas_with_cluster[snp %in% cl_snps, cluster_label := sprintf("cluster_%02d", cid)]
}

summary_df <- do.call(rbind, cluster_rows)
summary_df <- summary_df[order(summary_df$chr, summary_df$cluster_start), ]
fwrite(summary_df, opt$`out-summary`, sep = "\t")
locus_idx  <- match(opt$locus, as.character(loci[["loci"]]))
# Pick overall lead SNP (lowest p among all significant SNPs)
lead_row <- sig[which.min(sig$p), ]
meta_final <- list(
  chr          = if (!is.na(locus_idx)) loci[["chr"]][locus_idx]                          else (if (nrow(sig) > 0) sig$chr[1] else NA),
  locus_start  = if (!is.na(locus_idx)) as.integer(loci[["start"]][locus_idx])             else min(gwas$bp),
  locus_end    = if (!is.na(locus_idx)) as.integer(loci[["end"]][locus_idx])               else max(gwas$bp),
  n_total      = nrow(gwas),
  n_sig        = nrow(sig),
  n_clusters   = nrow(summary_df),
  r2_threshold = opt$`r2-threshold`,
  distance_kb  = opt$`distance-kb`,
  lead_snp     = if (nrow(lead_row) > 0) as.character(lead_row$snp[1]) else "N/A",
  lead_p       = if (nrow(lead_row) > 0) as.numeric(lead_row$p[1]) else NA_real_
)
write_interactive_plot(gwas_with_cluster, opt$locus, opt$`p-threshold`, opt$`out-html`, meta = meta_final)

diag <- list(
  locus = opt$locus,
  p_threshold = opt$`p-threshold`,
  r2_threshold = opt$`r2-threshold`,
  distance_kb = opt$`distance-kb`,
  n_sig_snps = nrow(sig),
  n_clusters_ld_only = n_clusters_ld_only,
  n_clusters = nrow(summary_df),
  plot_dir = plot_dir
)
writeLines(jsonlite::toJSON(diag, auto_unbox = TRUE, pretty = TRUE), opt$`out-diag`)
writeLines("ok", opt$`done-file`)
