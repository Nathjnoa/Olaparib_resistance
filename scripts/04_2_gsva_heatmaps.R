#!/usr/bin/env Rscript

# ============================================================
# 04_2_gsva_heatmaps.R
#
# Heatmaps from GSVA results (3 datasets only):
#
#   Section 1: Per-sample Z-score heatmap (Fig 3A)
#              Input: GSVA_ALL_merged.tsv + GSVA_metadata.tsv
#   Section 2: Integrated logFC heatmap (Fig 3B)
#              Input: GSVA_DE_ALL_datasets.tsv
#
# Pathways shown: recurrent across >= 2 datasets (FDR < 0.10)
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

source(file.path({f <- commandArgs(trailingOnly=FALSE); f <- grep("--file=",f,value=TRUE); if(length(f)) dirname(normalizePath(sub("--file=","",f[1]))) else path.expand("~/bioinfo/projects/olaparib_resistance/scripts")}, "00_config.R"))

log_handle <- start_log("04_2_gsva_heatmaps")
set.seed(42)

gsva_dir   <- file.path(tables_dir, "gsva")
gsvaDE_dir <- file.path(tables_dir, "gsva_DE")
out_dir    <- file.path(figures_dir, "gsva_heatmaps")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dataset_levels <- c("A2780", "PEO1", "UWB_BRCAdef")

fdr_cut      <- FDR_CUTOFF_VIZ
min_datasets <- 2

# ============================================================
# SECTION 1: Per-sample Z-score heatmap (Fig 3A)
# ============================================================
message("=== Fig 3A: per-sample heatmap ===")

in_gsva   <- file.path(gsva_dir,   "GSVA_ALL_merged.tsv")
in_meta   <- file.path(gsva_dir,   "GSVA_metadata.tsv")
in_gsvaDE <- file.path(gsvaDE_dir, "GSVA_DE_ALL_datasets.tsv")

missing_files <- c(in_gsva, in_meta, in_gsvaDE)[!file.exists(c(in_gsva, in_meta, in_gsvaDE))]
if (length(missing_files) > 0) {
  stop("Archivos de entrada no encontrados:\n", paste(missing_files, collapse = "\n"))
}

gsva_tbl <- readr::read_tsv(in_gsva, show_col_types = FALSE)
gsva_mat  <- as.matrix(gsva_tbl[, -1, drop = FALSE])
rownames(gsva_mat) <- gsva_tbl$Pathway

meta <- readr::read_tsv(in_meta, show_col_types = FALSE) %>%
  mutate(
    dataset   = factor(dataset,   levels = dataset_levels),
    condition = factor(condition, levels = c("Parental", "Resistant"))
  )

# Keep only 3-dataset samples (exclude any residual BRCAprof if present)
meta <- meta %>% filter(dataset %in% dataset_levels)

stopifnot(all(meta$sample_id %in% colnames(gsva_mat)))

gsvaDE <- readr::read_tsv(in_gsvaDE, show_col_types = FALSE) %>%
  filter(Dataset %in% dataset_levels) %>%
  mutate(Dataset = factor(Dataset, levels = dataset_levels))

# Select recurrent pathways
core_tbl <- gsvaDE %>%
  filter(!is.na(adj.P.Val), adj.P.Val < fdr_cut) %>%
  group_by(Pathway) %>%
  summarise(
    n_datasets     = n_distinct(Dataset),
    mean_abs_logFC = mean(abs(logFC), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_datasets >= min_datasets) %>%
  arrange(desc(mean_abs_logFC))

core_paths <- core_tbl$Pathway
if (length(core_paths) == 0) {
  warning("No pathways recurrentes con FDR<", fdr_cut, " en >= ", min_datasets, " datasets.")
} else {
  top_max <- 20
  if (length(core_paths) > top_max) core_paths <- core_paths[1:top_max]

  mat_3 <- gsva_mat[intersect(core_paths, rownames(gsva_mat)), meta$sample_id, drop = FALSE]

  meta_ord <- meta %>% arrange(dataset, condition, sample_id)
  mat_3    <- mat_3[, meta_ord$sample_id, drop = FALSE]

  message("Pathways recurrentes seleccionados: ", nrow(mat_3))
  print(table(meta_ord$dataset, meta_ord$condition))

  mat_z <- t(scale(t(mat_3)))
  mat_z[is.na(mat_z)] <- 0
  cap <- 2
  mat_z <- pmax(pmin(mat_z, cap), -cap)
  col_fun <- colorRamp2(c(-cap, 0, cap), c("#2166ac", "white", "#b2182b"))

  ha <- HeatmapAnnotation(
    Dataset   = meta_ord$dataset,
    Condition = meta_ord$condition,
    annotation_name_gp = grid::gpar(fontsize = 10)
  )

  ht_A <- Heatmap(
    mat_z,
    name         = "Z(GSVA)",
    col          = col_fun,
    na_col       = "grey90",
    cluster_rows    = TRUE,
    cluster_columns = FALSE,
    column_split    = meta_ord$dataset,
    column_gap      = unit(2, "mm"),
    show_column_names = TRUE,
    show_row_names    = TRUE,
    row_names_side    = "right",
    row_names_gp      = grid::gpar(fontsize = 9),
    column_names_gp   = grid::gpar(fontsize = 9),
    column_names_rot  = 45,
    top_annotation    = ha,
    column_title      = paste0("GSVA Hallmarks (por muestra) — recurrentes (FDR<", fdr_cut,
                               ", ≥", min_datasets, " datasets)"),
    column_title_gp   = grid::gpar(fontsize = 14, fontface = "bold"),
    heatmap_legend_param = list(title = "Z(GSVA)")
  )

  out_3A_png <- file.path(out_dir, "Fig3A_GSVA_scores_per_sample_heatmap.png")
  out_3A_pdf <- sub("\\.png$", ".pdf", out_3A_png)

  png(out_3A_png, width = 3200, height = 2200, res = 300)
  ComplexHeatmap::draw(ht_A, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  pdf(out_3A_pdf, width = 12, height = 8, useDingbats = FALSE)
  ComplexHeatmap::draw(ht_A, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  message("✔ Fig 3A -> ", out_3A_png)
}

# ============================================================
# SECTION 2: Integrated logFC heatmap (Fig 3B)
# ============================================================
message("=== Fig 3B: integrated logFC heatmap ===")

tbl <- readr::read_tsv(in_gsvaDE, show_col_types = FALSE) %>%
  filter(Dataset %in% dataset_levels) %>%
  mutate(Dataset = factor(Dataset, levels = dataset_levels))

stopifnot(all(c("Pathway","Dataset","logFC","adj.P.Val") %in% colnames(tbl)))

core <- tbl %>%
  filter(!is.na(adj.P.Val), adj.P.Val < fdr_cut) %>%
  group_by(Pathway) %>%
  summarise(
    n_datasets     = n_distinct(Dataset),
    mean_abs_logFC = mean(abs(logFC), na.rm = TRUE),
    sign_consistent = {
      s <- sign(logFC[!is.na(adj.P.Val) & adj.P.Val < fdr_cut])
      if (length(s) == 0) NA else (length(unique(s)) == 1)
    },
    .groups = "drop"
  ) %>%
  filter(n_datasets >= min_datasets) %>%
  arrange(desc(mean_abs_logFC))

if (nrow(core) == 0) {
  warning("No hay pathways recurrentes para Fig 3B.")
} else {
  core_paths_B <- core$Pathway

  mat_B <- tbl %>%
    filter(Pathway %in% core_paths_B) %>%
    dplyr::select(Pathway, Dataset, logFC) %>%
    pivot_wider(names_from = Dataset, values_from = logFC) %>%
    tibble::column_to_rownames("Pathway") %>%
    as.matrix()

  row_order <- order(rowMeans(abs(mat_B), na.rm = TRUE), decreasing = TRUE)
  mat_B <- mat_B[row_order, dataset_levels, drop = FALSE]

  out_tsv_B <- file.path(gsvaDE_dir, "GSVA_core_pathways_summary.tsv")
  readr::write_tsv(
    core %>% mutate(Pathway = factor(Pathway, levels = rownames(mat_B))) %>% arrange(Pathway),
    out_tsv_B
  )

  cap_B <- 1
  mat_B2 <- pmax(pmin(mat_B, cap_B), -cap_B)
  col_fun_B <- colorRamp2(c(-cap_B, 0, cap_B), c("#2166ac", "white", "#b2182b"))

  ht_B <- Heatmap(
    mat_B2,
    name         = "logFC",
    col          = col_fun_B,
    na_col       = "grey90",
    cluster_rows    = TRUE,
    cluster_columns = TRUE,
    show_row_names    = TRUE,
    row_names_side    = "right",
    row_names_gp      = grid::gpar(fontsize = 10),
    column_names_gp   = grid::gpar(fontsize = 11),
    column_names_rot  = 45,
    column_title      = paste0("ΔGSVA (Resistant − Parental) — Hallmarks recurrentes (FDR<", fdr_cut,
                               ", ≥", min_datasets, " datasets)"),
    column_title_gp   = grid::gpar(fontsize = 14, fontface = "bold"),
    heatmap_legend_param = list(title = "logFC (Δscore)")
  )

  out_3B_png <- file.path(out_dir, "Fig3B_GSVA_logFC_heatmap_integrated.png")
  out_3B_pdf <- sub("\\.png$", ".pdf", out_3B_png)

  png(out_3B_png, width = 2800, height = 2200, res = 300)
  ComplexHeatmap::draw(ht_B, heatmap_legend_side = "right")
  dev.off()

  pdf(out_3B_pdf, width = 10, height = 8, useDingbats = FALSE)
  ComplexHeatmap::draw(ht_B, heatmap_legend_side = "right")
  dev.off()

  message("✔ Fig 3B -> ", out_3B_png)
  message("✔ Core summary -> ", out_tsv_B)
}

stop_log(log_handle)
