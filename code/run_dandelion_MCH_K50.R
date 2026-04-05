.libPaths(c('/home/jliucx/R/x86_64-pc-linux-gnu-library/4.4', .libPaths()))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(rtracklayer)
  library(igraph)
  library(scales)
})

# ── 1. Load dfl2 (saved MCH-K50 mediation result) ────────────────────────────
dfl2 <- readRDS("/project/xuanyao/jiaming/paper/output/mediation_analysis/MCH/flames+200kb/k50/result.rds")

# ── 2. ENSG -> symbol (GTF) ───────────────────────────────────────────────────
gtf <- rtracklayer::import(
  "/project/xuanyao/jiaming/Getting_started/data/protein_coding_genes/Homo_sapiens.GRCh38.112.gtf.gz",
  feature.type = "gene")
gtf_df <- as.data.frame(gtf)[, c("gene_id", "gene_name")]
colnames(gtf_df) <- c("ensg", "symbol")
ensg_to_symbol <- setNames(gtf_df$symbol, gtf_df$ensg)

# ── 3. Build rsid -> cluster_id (200kb) ──────────────────────────────────────
gwas_result <- readRDS("/project/xuanyao/jiaming/Getting_started/data/GWAS/MCH/flames+200kb_mapping/result.rds")

rsid_locations <- gwas_result %>%
  dplyr::select(rsid, locations) %>%
  distinct() %>%
  filter(!is.na(locations)) %>%
  mutate(
    chr = sub(":.*", "", locations),
    pos = as.numeric(sub(".*:", "", locations))
  ) %>%
  filter(chr %in% as.character(1:22)) %>%
  arrange(chr, pos) %>%
  group_by(chr) %>%
  mutate(cluster_num = cumsum(c(1, diff(pos) > 200000))) %>%
  ungroup() %>%
  mutate(cluster_id = paste0(chr, "_", cluster_num)) %>%
  dplyr::select(rsid, chr, pos, cluster_id)

symbol_to_rsid <- gwas_result %>%
  filter(!is.na(symbol)) %>%
  dplyr::select(symbol, rsid) %>%
  distinct()

# ── 4. Build merged_df ────────────────────────────────────────────────────────
tbl <- dfl2 %>%
  mutate(ensg   = sub("^gRNA:", "", snp),
         symbol = ensg_to_symbol[ensg]) %>%
  left_join(symbol_to_rsid, by = "symbol", relationship = "many-to-many") %>%
  left_join(rsid_locations %>% dplyr::select(rsid, cluster_id), by = "rsid") %>%
  mutate(cluster_id = ifelse(is.na(cluster_id), paste0("solo_", ensg), cluster_id))

merged_df <- tbl %>%
  group_by(module, cluster_id) %>%
  summarise(
    snp_geneset      = min(snp_geneset,   na.rm = TRUE),
    geneset_trait    = min(geneset_trait, na.rm = TRUE),
    merged_gene_name = paste(sort(unique(symbol[!is.na(symbol)])), collapse = "-"),
    merged_loci      = paste(sort(unique(rsid[!is.na(rsid)])),     collapse = "-"),
    .groups = "drop"
  ) %>%
  mutate(
    summary   = paste0(merged_gene_name, "(", merged_loci, ")"),
    gene_name = merged_gene_name
  ) %>%
  dplyr::select(module, merged_gene_name, merged_loci, summary, gene_name, snp_geneset, geneset_trait)

cat("merged_df dim:", nrow(merged_df), "x", ncol(merged_df), "\n")
print(head(merged_df[, 1:4], 10))

# ── 5. Gene category vectors ──────────────────────────────────────────────────
haemoglobin_genes <- c(
  "KLF1", "NPC1", "SUPT5H", "CCNF", "CAD", "CALR", "CBFA2T3",
  "CCNA2", "CLSPN", "DNAJC24", "NDUFB7", "NUP214", "POLE",
  "ABCB10", "DNAJC13", "MED13L", "MED23", "MTOR", "NCOR1",
  "RPTOR", "SALL4"
)
cell_cycle_genes <- c("HSPA9", "MED14", "MED17", "SUPT5H", "SYMPK")
autophagy_genes  <- c(
  "AP2A1", "ATR", "EIF2B4", "EIF3E", "EIF3F", "EIF3M", "HSPA9",
  "MED13L", "MED14", "MED17", "MED26", "MTOR", "NIPBL", "POLR2G",
  "RNF20", "RPTOR", "RTTN", "SPRED2", "SUPT5H", "SYMPK", "TP73",
  "TRAPPC10", "WDFY3", "CLSPN", "EEF2", "POLE", "BPGM", "UBE2D3",
  "USP10", "XPOT"
)

# ── 6. Plot function ──────────────────────────────────────────────────────────
plot_gene_module_network_MCH <- function(df12,
                                         gene_modules_id = NULL,
                                         top_n_per_module = NULL,
                                         top_n_gene = 80,
                                         pcap = 1e-300,
                                         edge_by = c("none","snp_geneset"),
                                         layout_fun = NULL,
                                         highlight_mod = NULL,
                                         out_pdf = NULL,
                                         width = 12, height = 9,
                                         label_genes = TRUE,
                                         label_modules = TRUE,
                                         size_scale = 0.7,
                                         mod_size_range  = c(15, 40),
                                         gene_size_range = c(4,  20),
                                         highlight_genes = NULL,
                                         sig_GO = NULL,
                                         haemoglobin_genes = NULL,
                                         cell_cycle_genes  = NULL,
                                         autophagy_genes   = NULL,
                                         unclassified_alpha = 0.3) {

  edge_by <- match.arg(edge_by)
  if (is.null(layout_fun)) layout_fun <- igraph::layout_with_fr

  mix_hex_colors <- function(hex_vec) {
    if (length(hex_vec) == 1) return(hex_vec)
    rgb_mat <- grDevices::col2rgb(hex_vec)
    avg <- rowMeans(rgb_mat)
    grDevices::rgb(avg[1], avg[2], avg[3], maxColorValue = 255)
  }
  add_alpha <- function(col, alpha) adjustcolor(col, alpha.f = alpha)

  col_haem  <- "#E74C3C"
  col_cycle <- "#F1C40F"
  col_auto  <- "#3498DB"
  default_gene_col   <- "grey60"
  default_module_col <- "grey10"

  get_gene_color <- function(raw_name) {
    no_snp     <- trimws(sub("\\s*\\(.*", "", raw_name))
    pure_names <- unlist(strsplit(no_snp, "-"))
    cats <- character(0)
    if (!is.null(haemoglobin_genes) && any(pure_names %in% haemoglobin_genes)) cats <- c(cats, col_haem)
    if (!is.null(cell_cycle_genes)  && any(pure_names %in% cell_cycle_genes))  cats <- c(cats, col_cycle)
    if (!is.null(autophagy_genes)   && any(pure_names %in% autophagy_genes))   cats <- c(cats, col_auto)
    if (length(cats) == 0) return(add_alpha(default_gene_col, unclassified_alpha))
    add_alpha(mix_hex_colors(cats), 1.0)
  }

  get_module_color <- function(mod_name) {
    if (is.null(sig_GO)) return(add_alpha(default_module_col, unclassified_alpha))
    rows <- sig_GO[sig_GO$module == mod_name, , drop = FALSE]
    if (nrow(rows) == 0) return(add_alpha(default_module_col, unclassified_alpha))
    summ <- paste(rows$summary, collapse = " ")
    cats <- character(0)
    if (grepl("hemoglobin|haemoglobin", summ, ignore.case = TRUE)) cats <- c(cats, col_haem)
    if (grepl("cell cycle",             summ, ignore.case = TRUE)) cats <- c(cats, col_cycle)
    if (grepl("autophagy",              summ, ignore.case = TRUE)) cats <- c(cats, col_auto)
    if (length(cats) == 0) return(add_alpha(default_module_col, unclassified_alpha))
    mix_hex_colors(cats)
  }

  d <- df12 %>%
    dplyr::transmute(
      module = as.character(module),
      gene   = as.character(gene_name),
      p_snp  = pmax(as.numeric(snp_geneset),  pcap),
      p_mod  = pmax(as.numeric(geneset_trait), pcap)
    ) %>%
    dplyr::filter(!is.na(module), !is.na(gene)) %>%
    dplyr::group_by(module, gene) %>%
    dplyr::summarise(p_snp = min(p_snp, na.rm = TRUE),
                     p_mod = min(p_mod, na.rm = TRUE), .groups = "drop")

  if (!is.null(top_n_per_module)) {
    d <- d %>% dplyr::group_by(module) %>%
      dplyr::arrange(p_snp, .by_group = TRUE) %>%
      dplyr::slice_head(n = top_n_per_module) %>% dplyr::ungroup()
  }
  if (!is.null(top_n_gene)) {
    top_genes <- d %>% dplyr::group_by(gene) %>%
      dplyr::summarise(best_p = min(p_snp, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(best_p) %>% dplyr::slice_head(n = top_n_gene) %>% dplyr::pull(gene)
    d <- d %>% dplyr::filter(gene %in% top_genes)
  }

  mod_tbl <- d %>% dplyr::group_by(module) %>%
    dplyr::summarise(p_mod = min(p_mod, na.rm = TRUE), n_assoc = dplyr::n_distinct(gene), .groups = "drop")
  gene_tbl <- d %>% dplyr::group_by(gene) %>%
    dplyr::summarise(p_snp = min(p_snp, na.rm = TRUE), n_assoc = dplyr::n_distinct(module), .groups = "drop")

  vertices <- dplyr::bind_rows(
    mod_tbl  %>% dplyr::transmute(name = module, type = "module", weight = n_assoc),
    gene_tbl %>% dplyr::transmute(name = gene,   type = "gene",   weight = n_assoc)
  )
  edges <- d %>% dplyr::transmute(from = module, to = gene, p_snp = p_snp)

  g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = vertices)

  igraph::V(g)$shape       <- ifelse(igraph::V(g)$type == "module", "square", "circle")
  igraph::V(g)$label.color <- "white"
  igraph::V(g)$frame.color <- NA
  igraph::V(g)$color       <- ifelse(igraph::V(g)$type == "module",
                                     add_alpha(default_module_col, unclassified_alpha),
                                     add_alpha(default_gene_col,   unclassified_alpha))

  mod_ids  <- which(igraph::V(g)$type == "module")
  gene_ids <- which(igraph::V(g)$type == "gene")

  for (idx in gene_ids) igraph::V(g)$color[idx] <- get_gene_color(igraph::V(g)$name[idx])
  for (idx in mod_ids)  igraph::V(g)$color[idx] <- get_module_color(igraph::V(g)$name[idx])

  if (!is.null(highlight_mod))
    for (idx in mod_ids[igraph::V(g)$name[mod_ids] %in% highlight_mod])
      igraph::V(g)$color[idx] <- "green3"
  if (!is.null(highlight_genes))
    for (idx in gene_ids[igraph::V(g)$name[gene_ids] %in% highlight_genes])
      igraph::V(g)$color[idx] <- "red"

  igraph::V(g)$size <- rep(NA_real_, igraph::vcount(g))
  mod_w  <- igraph::V(g)$weight[mod_ids]
  igraph::V(g)$size[mod_ids]  <- if (length(unique(mod_w))  == 1) rep(mean(mod_size_range)  * size_scale, length(mod_ids))
                                  else scales::rescale(mod_w,  to = mod_size_range  * size_scale)
  gene_w <- igraph::V(g)$weight[gene_ids]
  igraph::V(g)$size[gene_ids] <- if (length(unique(gene_w)) == 1) rep(mean(gene_size_range) * size_scale, length(gene_ids))
                                  else scales::rescale(gene_w, to = gene_size_range * size_scale)

  igraph::E(g)$color <- "grey70"
  igraph::E(g)$width <- if (edge_by == "snp_geneset")
    scales::rescale(-log10(igraph::E(g)$p_snp), to = c(0.5, 4)) else 1

  igraph::V(g)$label     <- NA_character_
  igraph::V(g)$label.cex <- ifelse(igraph::V(g)$type == "module", 0.9, 0.7)
  if (label_modules) for (idx in mod_ids)  igraph::V(g)$label[idx] <- igraph::V(g)$name[idx]
  if (label_genes)   for (idx in gene_ids) igraph::V(g)$label[idx] <- igraph::V(g)$name[idx]

  if (!is.null(out_pdf)) grDevices::pdf(out_pdf, width = width, height = height)

  plot(g, layout = layout_fun(g), vertex.label.dist = 0, edge.curved = 0)

  leg_labels <- c("Haemoglobin","Cell cycle","Autophagy",
                  "Cycle+Haem","Cycle+Auto","Haem+Auto","All three",
                  "Unclassified (gene)","Unclassified (module)")
  leg_cols <- c(col_haem, col_cycle, col_auto,
                mix_hex_colors(c(col_cycle, col_haem)),
                mix_hex_colors(c(col_cycle, col_auto)),
                mix_hex_colors(c(col_haem,  col_auto)),
                mix_hex_colors(c(col_haem, col_cycle, col_auto)),
                add_alpha(default_gene_col,   unclassified_alpha),
                add_alpha(default_module_col, unclassified_alpha))
  graphics::legend("bottomleft", legend = leg_labels, pt.bg = leg_cols, col = leg_cols,
                   pch = c(rep(19, 8), 15), pt.cex = 1.4, cex = 0.75, bty = "n", title = "Category")

  if (!is.null(out_pdf)) grDevices::dev.off()
  invisible(g)
}

# ── 7. Run and save plot ──────────────────────────────────────────────────────
out_png <- "/project/xuanyao/jiaming/paper/output/mediation_analysis/MCH/flames+200kb/k50/dandelion_plot.png"

set.seed(42)
png(out_png, width = 10000, height = 10000, res = 400)
plot_gene_module_network_MCH(
  df12              = merged_df,
  sig_GO            = NULL,          # sig_GO needs gost; skipped here
  haemoglobin_genes = haemoglobin_genes,
  cell_cycle_genes  = cell_cycle_genes,
  autophagy_genes   = autophagy_genes,
  top_n_gene        = 800,
  size_scale        = 0.3,
  edge_by           = "snp_geneset",
  unclassified_alpha = 0.3
)
dev.off()
cat("Plot saved to:", out_png, "\n")
