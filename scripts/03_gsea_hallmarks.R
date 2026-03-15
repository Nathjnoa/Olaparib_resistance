#!/usr/bin/env Rscript

# ============================================================
# 03_gsea_hallmarks.R
#
# GSEA (Hallmarks) across 3 datasets + integrated NES heatmap:
#
#   Section 1: Run GSEA per dataset → TSV + lollipop + EnrichmentMap
#   Section 2: Load per-dataset results (kept in memory), build
#              integrated NES heatmap (PNG + PDF)
#
# Ranking: limma → t-statistic; DESeq2 → z = log2FC / lfcSE
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(msigdbr)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

source(file.path({f <- commandArgs(trailingOnly=FALSE); f <- grep("--file=",f,value=TRUE); if(length(f)) dirname(normalizePath(sub("--file=","",f[1]))) else path.expand("~/bioinfo/projects/olaparib_resistance/scripts")}, "00_config.R"))

log_handle <- start_log("03_gsea_hallmarks")
set.seed(42)

out_tab_dir <- file.path(tables_dir, "gsea_hallmarks")
out_fig_dir <- file.path(figures_dir, "gsea_hallmarks")

dir.create(out_tab_dir,                          showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_fig_dir, "lollipop"),   showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_fig_dir, "emap"),       showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Hallmarks gene sets
# ============================================================
msig_h   <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::distinct(gs_name, gene_symbol)
TERM2GENE <- msig_h %>% dplyr::select(gs_name, gene_symbol)

# ============================================================
# Ranking helpers
# ============================================================
make_rank_limma <- function(df, id_col = "ID", stat_col = "t_stat") {
  ids  <- as.character(df[[id_col]])
  stat <- as.numeric(df[[stat_col]])
  sym  <- map_to_symbol_if_needed(ids)
  names(stat) <- sym
  stat <- collapse_duplicate_symbols(stat)
  sort(stat, decreasing = TRUE)
}

make_rank_deseq <- function(df, id_col = "gene_id",
                            lfc_col = "log2FoldChange", se_col = "lfcSE") {
  ids <- as.character(df[[id_col]])
  lfc <- as.numeric(df[[lfc_col]])
  se  <- as.numeric(df[[se_col]])
  z   <- lfc / se
  sym <- map_to_symbol_if_needed(ids)
  names(z) <- sym
  z <- collapse_duplicate_symbols(z)
  sort(z, decreasing = TRUE)
}

# ============================================================
# GSEA runner + plots (returns result tibble invisibly)
# ============================================================
run_gsea_and_plots <- function(ranks, dataset_label, out_prefix) {
  gsea <- suppressMessages(
    GSEA(
      geneList     = ranks,
      TERM2GENE    = TERM2GENE,
      minGSSize    = MIN_GS,
      maxGSSize    = MAX_GS,
      pvalueCutoff = 1,
      verbose      = FALSE,
      eps          = 0
    )
  )

  res <- as_tibble(gsea@result) %>%
    dplyr::select(any_of(c("ID","Description","setSize","enrichmentScore",
                            "NES","pvalue","p.adjust","qvalues"))) %>%
    arrange(p.adjust)

  out_tsv <- file.path(out_tab_dir, paste0(out_prefix, "_GSEA_Hallmarks.tsv"))
  readr::write_tsv(res, out_tsv)

  # Lollipop
  n_half  <- ceiling(TOP_N_PLOT / 2)
  top_pos <- res %>%
    filter(is.finite(NES), !is.na(p.adjust), NES > 0, p.adjust < FDR_CUTOFF_VIZ) %>%
    arrange(p.adjust) %>% slice_head(n = n_half)
  top_neg <- res %>%
    filter(is.finite(NES), !is.na(p.adjust), NES < 0, p.adjust < FDR_CUTOFF_VIZ) %>%
    arrange(p.adjust) %>% slice_head(n = n_half)

  plot_df <- bind_rows(top_neg, top_pos) %>%
    mutate(
      Pathway = stringr::str_replace_all(Description, "^HALLMARK_", "") %>%
        stringr::str_replace_all("_", " ") %>%
        stringr::str_to_title()
    ) %>%
    arrange(NES) %>%
    mutate(Pathway = factor(Pathway, levels = Pathway))

  if (nrow(plot_df) > 0) {
    p_lolli <- ggplot(plot_df, aes(x = NES, y = Pathway, color = NES > 0)) +
      geom_segment(aes(x = 0, xend = NES, yend = Pathway), linewidth = 0.8) +
      geom_point(aes(size = -log10(p.adjust)), alpha = 0.9) +
      scale_color_manual(values = c("TRUE" = "#e41a1c", "FALSE" = "#377eb8"),
                         guide = "none") +
      theme_bw(base_size = 12) +
      labs(title = dataset_label, x = "NES (Resistant vs Parental)",
           y = NULL, size = "-log10(FDR)") +
      coord_cartesian(xlim = c(-2.5, 2.5))

    ggsave(file.path(out_fig_dir, "lollipop", paste0(out_prefix, "_lollipop.png")),
           p_lolli, width = 6.5, height = 4.5, dpi = 300)
  }

  # EnrichmentMap
  res_sig <- res %>% filter(!is.na(p.adjust), p.adjust < FDR_CUTOFF_VIZ)
  if (nrow(res_sig) >= 8) {
    gsea2 <- gsea
    gsea2@result <- gsea2@result %>%
      as.data.frame() %>%
      filter(p.adjust < FDR_CUTOFF_VIZ) %>%
      mutate(Description = gsub("^HALLMARK_", "", Description),
             Description = gsub("_", " ", Description),
             Description = tools::toTitleCase(tolower(Description)))
    gsea_sim <- pairwise_termsim(gsea2)
    show_n   <- min(20, nrow(gsea2@result))
    out_emap <- file.path(out_fig_dir, "emap", paste0(out_prefix, "_emap_clean.png"))
    png(out_emap, width = 3000, height = 2200, res = 300)
    print(
      emapplot(gsea_sim, showCategory = show_n, layout = "fr", node_label = "none") +
        ggtitle(paste0(dataset_label, " — EnrichmentMap (Hallmarks)"))
    )
    dev.off()
  }

  message("✔ ", dataset_label, " -> ", out_tsv)
  invisible(res)
}

# ============================================================
# SECTION 1: Run GSEA per dataset (results kept in memory)
# ============================================================

# A2780 (limma)
file_A <- file.path(tables_dir, "de", "GSE153867_A2780_limma_DE_FPKM.tsv")
if (!file.exists(file_A)) stop("No encuentro: ", file_A)
deg_A  <- readr::read_tsv(file_A, show_col_types = FALSE)
rA     <- make_rank_limma(deg_A, id_col = "ID", stat_col = "t_stat")
res_A  <- run_gsea_and_plots(rA, "GSE153867 — A2780", "A2780")

# PEO1 (limma)
file_P <- file.path(tables_dir, "de", "GSE117765_PEO1_limma_TPM_Res_vs_Par.tsv")
if (!file.exists(file_P)) stop("No encuentro: ", file_P)
deg_P  <- readr::read_tsv(file_P, show_col_types = FALSE)
rP     <- make_rank_limma(deg_P, id_col = "ID", stat_col = "t_stat")
res_P  <- run_gsea_and_plots(rP, "GSE117765 — PEO1", "PEO1")

# UWB BRCA-def (DESeq2)
file_def <- file.path(tables_dir, "de", "GSE235980_BRCAdef_Res_vs_Par_DESeq2.tsv")
if (!file.exists(file_def)) stop("No encuentro: ", file_def)
deg_def  <- readr::read_tsv(file_def, show_col_types = FALSE)
rD       <- make_rank_deseq(deg_def, id_col = "gene_id",
                             lfc_col = "log2FoldChange", se_col = "lfcSE")
res_D    <- run_gsea_and_plots(rD, "GSE235980 — UWB BRCA-def", "UWB_BRCAdef")

message("=== GSEA per-dataset COMPLETADO ===")

# ============================================================
# SECTION 2: Integrated NES heatmap (3 datasets)
# ============================================================
message("=== Construyendo heatmap integrado de NES ===")

dataset_levels <- c("A2780", "PEO1", "UWB_BRCAdef")

# Use results already in memory (no re-read needed)
gsea_results <- list(
  A2780       = res_A %>% mutate(dataset = "A2780"),
  PEO1        = res_P %>% mutate(dataset = "PEO1"),
  UWB_BRCAdef = res_D %>% mutate(dataset = "UWB_BRCAdef")
)

tbl <- bind_rows(gsea_results) %>%
  transmute(
    dataset       = factor(dataset, levels = dataset_levels),
    hallmark      = Description,
    NES           = as.numeric(NES),
    padj          = as.numeric(p.adjust)
  ) %>%
  mutate(hallmark_clean = clean_hallmark(hallmark)) %>%
  filter(!is.na(hallmark_clean), hallmark_clean != "")

fdr_cut      <- FDR_CUTOFF_VIZ
min_datasets <- 2
cap_nes      <- 2

recurrent <- tbl %>%
  filter(!is.na(padj), padj < fdr_cut) %>%
  count(hallmark_clean) %>%
  filter(n >= min_datasets) %>%
  pull(hallmark_clean)

if (length(recurrent) == 0) {
  warning("No hay hallmarks recurrentes con FDR<", fdr_cut, " en >=", min_datasets, " datasets.")
} else {
  mat <- tbl %>%
    filter(hallmark_clean %in% recurrent) %>%
    dplyr::select(dataset, hallmark_clean, NES) %>%
    pivot_wider(names_from = dataset, values_from = NES) %>%
    column_to_rownames("hallmark_clean") %>%
    as.matrix()

  mat <- pmax(pmin(mat, cap_nes), -cap_nes)
  row_order <- order(rowMeans(abs(mat), na.rm = TRUE), decreasing = TRUE)
  mat <- mat[row_order, , drop = FALSE]

  col_fun_nes <- colorRamp2(c(-cap_nes, 0, cap_nes), c("#2166ac", "white", "#b2182b"))

  ht <- Heatmap(
    mat,
    name        = "NES",
    col         = col_fun_nes,
    na_col      = "grey90",
    cluster_rows    = TRUE,
    cluster_columns = TRUE,
    show_row_names  = TRUE,
    row_names_side  = "right",
    row_names_gp    = grid::gpar(fontsize = 10),
    show_column_names   = TRUE,
    column_names_rot    = 45,
    column_names_gp     = grid::gpar(fontsize = 11),
    column_title        = paste0("Hallmarks GSEA (NES) — recurrentes en ≥", min_datasets,
                                 " datasets (FDR<", fdr_cut, ")"),
    column_title_gp     = grid::gpar(fontsize = 16, fontface = "bold"),
    heatmap_legend_param = list(title = "NES")
  )

  out_png <- file.path(out_fig_dir, "Hallmarks_NES_heatmap.png")
  out_pdf <- sub("\\.png$", ".pdf", out_png)

  png(out_png, width = 2800, height = 2100, res = 300)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
  dev.off()

  pdf(out_pdf, width = 10, height = 8, useDingbats = FALSE)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
  dev.off()

  message("✔ NES heatmap -> ", out_png)
}

stop_log(log_handle)
