#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(msigdbr)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
})

proj_dir <- Sys.getenv("OLAPARIB_RESISTANCE_DIR", unset = "")
if (!nzchar(proj_dir)) {
  proj_dir <- "~/bioinfo/projects/olaparib_resistance"
}
proj_dir <- path.expand(proj_dir)
if (!dir.exists(proj_dir)) {
  stop("No existe el directorio del proyecto: ", proj_dir)
}
tables_dir <- file.path(proj_dir, "results", "tables")
out_tab_dir <- file.path(proj_dir, "results", "tables", "gsea_hallmarks")
out_fig_dir <- file.path(proj_dir, "results", "figures", "gsea_hallmarks")

dir.create(out_tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_fig_dir, "lollipop"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_fig_dir, "emap"), recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Parámetros
# ---------------------------
padj_plot_cut <- 0.10   # para mostrar más señal en figura (puedes bajar a 0.05)
show_top <- 12          # pathways por lollipop
minGS <- 15
maxGS <- 500

# ---------------------------
# Hallmarks TERM2GENE
# ---------------------------
msig_h <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::distinct(gs_name, gene_symbol)

TERM2GENE <- msig_h %>% dplyr::select(gs_name, gene_symbol)

# ---------------------------
# Helpers
# ---------------------------

collapse_duplicate_symbols <- function(stats_named) {
  # stats_named: named numeric vector (names = SYMBOL)
  tibble(symbol = names(stats_named), stat = as.numeric(stats_named)) %>%
    filter(!is.na(symbol), symbol != "", !is.na(stat), is.finite(stat)) %>%
    group_by(symbol) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    deframe()
}

map_to_symbol_if_needed <- function(ids) {
  # Si parece Ensembl (ENSG...), mapear a SYMBOL. Si no, asumir ya es SYMBOL.
  ids2 <- gsub("\\..*$", "", ids) # quitar versión ENSGxxxx.1
  is_ens <- grepl("^ENSG", ids2)

  if (!any(is_ens)) return(ids2)

  sym <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = unique(ids2),
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  sym2 <- unname(sym[ids2])
  sym2
}

make_rank_limma <- function(df, id_col = "ID", stat_col = "t_stat") {
  ids <- as.character(df[[id_col]])
  stat <- as.numeric(df[[stat_col]])
  sym <- map_to_symbol_if_needed(ids)
  stats <- stat
  names(stats) <- sym
  stats <- collapse_duplicate_symbols(stats)
  sort(stats, decreasing = TRUE)
}

make_rank_deseq <- function(df, id_col = "gene_id", lfc_col = "log2FoldChange", se_col = "lfcSE") {
  ids <- as.character(df[[id_col]])
  lfc <- as.numeric(df[[lfc_col]])
  se  <- as.numeric(df[[se_col]])
  z <- lfc / se
  sym <- map_to_symbol_if_needed(ids)
  names(z) <- sym
  z <- collapse_duplicate_symbols(z)
  sort(z, decreasing = TRUE)
}

run_gsea_and_plots <- function(ranks, dataset_label, out_prefix) {
  # clusterProfiler GSEA (usa fgsea internamente)
  gsea <- suppressMessages(
    GSEA(
      geneList = ranks,
      TERM2GENE = TERM2GENE,
      minGSSize = minGS,
      maxGSSize = maxGS,
      pvalueCutoff = 1,
      verbose = FALSE,
      eps = 0
    )
  )

  res <- as_tibble(gsea@result)
  keep_cols <- intersect(
    c("ID","Description","setSize","enrichmentScore","NES","pvalue","p.adjust","qvalues"),
    colnames(res)
  )
  res <- res %>%
    dplyr::select(all_of(keep_cols)) %>%
    arrange(p.adjust)

  out_tsv <- file.path(out_tab_dir, paste0(out_prefix, "_GSEA_Hallmarks.tsv"))
  write_tsv(res, out_tsv)

  # Lollipop: top pathways por padj (mostrar también NES negativo/positivo)
  n_half <- ceiling(show_top / 2)
  top_pos <- res %>%
    filter(is.finite(NES), !is.na(p.adjust), NES > 0, p.adjust < padj_plot_cut) %>%
    arrange(p.adjust) %>%
    slice_head(n = n_half)

  top_neg <- res %>%
    filter(is.finite(NES), !is.na(p.adjust), NES < 0, p.adjust < padj_plot_cut) %>%
    arrange(p.adjust) %>%
    slice_head(n = n_half)

  plot_df <- bind_rows(top_neg, top_pos) %>%
    mutate(
      Pathway = str_replace_all(Description, "^HALLMARK_", ""),
      Pathway = str_replace_all(Pathway, "_", " "),
      Pathway = str_to_title(Pathway),
      Pathway = case_when(
        Pathway == "E2f" ~ "E2F",
        Pathway == "Nfkb" ~ "NFKB",
        Pathway == "Uv" ~ "UV",
        TRUE ~ Pathway
      )
    ) %>%
    arrange(NES) %>%
    mutate(Pathway = factor(Pathway, levels = Pathway))

  p_lolli <- ggplot(plot_df, aes(x = NES, y = Pathway, color = NES > 0)) +
    geom_segment(aes(x = 0, xend = NES, yend = Pathway), linewidth = 0.8) +
    geom_point(aes(size = -log10(p.adjust)), alpha = 0.9) +
    scale_color_manual(
      values = c("TRUE" = "#e41a1c", "FALSE" = "#377eb8"),
      labels = c("FALSE" = "Up en parental", "TRUE" = "Up en resistente"),
      name = "Dirección",
      guide = "none"
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = dataset_label,
      x = "NES (Resistant vs Parental)",
      y = NULL,
      size = "-log10(FDR)"
    ) +
    coord_cartesian(xlim = c(-2.5, 2.5))

  out_png_lolli <- file.path(out_fig_dir, "lollipop", paste0(out_prefix, "_lollipop.png"))
  ggsave(out_png_lolli, p_lolli, width = 6.5, height = 4.5, dpi = 300)

  # EnrichmentMap filtrando términos significativos y limpiando nombres
  res_sig <- res %>% filter(!is.na(p.adjust), p.adjust < padj_plot_cut)
  if (nrow(res_sig) < 8) {
    message("Omitiendo EnrichmentMap para ", dataset_label,
            " (solo ", nrow(res_sig), " términos con FDR <", padj_plot_cut, ")")
  } else {
    gsea2 <- gsea
    gsea2@result <- gsea2@result %>%
      as.data.frame() %>%
      filter(p.adjust < padj_plot_cut)

    gsea2@result$Description <- gsub("^HALLMARK_", "", gsea2@result$Description)
    gsea2@result$Description <- gsub("_", " ", gsea2@result$Description)
    gsea2@result$Description <- tools::toTitleCase(tolower(gsea2@result$Description))

    gsea_sim <- pairwise_termsim(gsea2)
    show_n <- min(20, nrow(gsea2@result))
    out_png_emap <- file.path(out_fig_dir, "emap", paste0(out_prefix, "_emap_clean.png"))
    png(out_png_emap, width = 3000, height = 2200, res = 300)
    # Nota: algunas versiones de enrichplot permiten node_label; si no, se ignora.
    print(
      emapplot(
        gsea_sim,
        showCategory = show_n,
        layout = "fr",
        node_label = "none"
      ) + ggtitle(paste0(dataset_label, " — EnrichmentMap (Hallmarks)"))
    )
    dev.off()
  }

  message("✔ ", dataset_label, " -> ", out_tsv)
  invisible(res)
}

# ---------------------------
# Cargar tablas DE y correr GSEA
# ---------------------------

# A2780 (limma)
file_A <- file.path(tables_dir, "GSE153867_A2780_limma_DE_FPKM.tsv")
if (!file.exists(file_A)) stop("No encuentro: ", file_A)
deg_A <- read_tsv(file_A, show_col_types = FALSE)
rA <- make_rank_limma(deg_A, id_col = "ID", stat_col = "t_stat")
run_gsea_and_plots(rA, "GSE153867 — A2780", "A2780")

# PEO1 (limma)
file_P <- file.path(tables_dir, "GSE117765_PEO1_limma_TPM_Res_vs_Par.tsv")
if (!file.exists(file_P)) stop("No encuentro: ", file_P)
deg_P <- read_tsv(file_P, show_col_types = FALSE)
rP <- make_rank_limma(deg_P, id_col = "ID", stat_col = "t_stat")
run_gsea_and_plots(rP, "GSE117765 — PEO1", "PEO1")

# UWB BRCA-def (DESeq2)
file_def <- file.path(tables_dir, "GSE235980_BRCAdef_Res_vs_Par_DESeq2.tsv")
if (!file.exists(file_def)) stop("No encuentro: ", file_def)
deg_def <- read_tsv(file_def, show_col_types = FALSE)
rD <- make_rank_deseq(deg_def, id_col = "gene_id", lfc_col = "log2FoldChange", se_col = "lfcSE")
run_gsea_and_plots(rD, "GSE235980 — UWB BRCA-def", "UWB_BRCAdef")

# UWB BRCA-prof (DESeq2)
file_prof <- file.path(tables_dir, "GSE235980_BRCAprof_Res_vs_Par_DESeq2.tsv")
if (!file.exists(file_prof)) stop("No encuentro: ", file_prof)
deg_prof <- read_tsv(file_prof, show_col_types = FALSE)
rR <- make_rank_deseq(deg_prof, id_col = "gene_id", lfc_col = "log2FoldChange", se_col = "lfcSE")
run_gsea_and_plots(rR, "GSE235980 — UWB BRCA-prof", "UWB_BRCAprof")

message("=== OK: GSEA Hallmarks + Lollipop + EnrichmentMap generados ===")
