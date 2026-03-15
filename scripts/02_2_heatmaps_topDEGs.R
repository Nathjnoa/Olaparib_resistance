#!/usr/bin/env Rscript

# ============================================================
# 02_2_heatmaps_topDEGs.R
#
# ComplexHeatmap of top DEGs per dataset (Z-score of log2-expression):
#   - GSE153867: A2780 (FPKM)
#   - GSE235980: UWB BRCA-deficient (VST counts)
#   - GSE117765: PEO1 (TPM)
#
# Top 20 up + 20 down (FDR < 0.05) per dataset.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(ComplexHeatmap)
  library(circlize)
  library(R.utils)
})

source(file.path({f <- commandArgs(trailingOnly=FALSE); f <- grep("--file=",f,value=TRUE); if(length(f)) dirname(normalizePath(sub("--file=","",f[1]))) else path.expand("~/bioinfo/projects/olaparib_resistance/scripts")}, "00_config.R"))

log_handle <- start_log("02_2_heatmaps_topDEGs")

out_dir <- file.path(figures_dir, "heatmaps")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Shared color scale
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b"))

require_file <- function(path) {
  if (!file.exists(path)) stop("No encuentro: ", path)
}

# ============================================================
# Helpers
# ============================================================
get_top_genes <- function(df, id_col, fdr_col, logfc_col,
                          n_up = TOP_N_GENES, n_down = TOP_N_GENES,
                          fdr_cut = FDR_CUTOFF) {
  df_sig <- df %>% filter(!is.na(.data[[fdr_col]]), .data[[fdr_col]] < fdr_cut)

  up <- df_sig %>%
    filter(.data[[logfc_col]] > 0) %>%
    arrange(desc(.data[[logfc_col]])) %>%
    head(n_up)

  down <- df_sig %>%
    filter(.data[[logfc_col]] < 0) %>%
    arrange(.data[[logfc_col]]) %>%
    head(n_down)

  genes <- c(up[[id_col]], down[[id_col]])
  dir   <- c(rep("Up", nrow(up)), rep("Down", nrow(down)))
  names(dir) <- genes
  genes <- unique(genes)
  dir   <- dir[genes]
  list(genes = genes, direction = dir)
}

make_heatmap <- function(expr_mat, meta, top_info, dataset_name, outfile_png) {
  meta     <- meta %>% arrange(condition)
  expr_mat <- expr_mat[, meta$sample_id, drop = FALSE]

  genes <- intersect(top_info$genes, rownames(expr_mat))
  if (length(genes) < 5) stop("Muy pocos genes para heatmap: ", dataset_name)

  mat_z <- t(scale(t(expr_mat[genes, , drop = FALSE])))
  mat_z[is.na(mat_z)] <- 0

  gene_dir <- factor(top_info$direction[genes], levels = c("Up", "Down"))

  ha_col <- HeatmapAnnotation(
    Condition = meta$condition,
    col = list(Condition = c(Parental = "#377eb8", Resistant = "#e41a1c")),
    annotation_name_gp = gpar(fontsize = 10)
  )

  ht <- Heatmap(
    mat_z,
    name = "Z-score",
    col  = col_fun,
    show_row_names = TRUE,
    row_names_side = "right",
    row_names_gp   = gpar(fontsize = 6),
    show_column_names = TRUE,
    column_names_rot  = 45,
    column_names_gp   = gpar(fontsize = 8),
    cluster_rows      = TRUE,
    cluster_columns   = FALSE,
    clustering_distance_rows = "pearson",
    row_split       = gene_dir,
    top_annotation  = ha_col,
    column_title    = dataset_name,
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(title = "Z-score",
                                title_gp  = gpar(fontsize = 10),
                                labels_gp = gpar(fontsize = 8))
  )

  png(outfile_png, width = 2200, height = 2600, res = 300)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  message("Heatmap -> ", outfile_png)
}

# ============================================================
# 1) A2780 (FPKM)
# ============================================================
message("=== A2780 ===")
fpkm_file <- file.path(raw_dir, "GSE153867", "GSE153867_fpkm.txt")
require_file(fpkm_file)

fpkm_tbl <- readr::read_tsv(fpkm_file, show_col_types = FALSE)
gene_ids  <- fpkm_tbl[[1]]
fpkm_mat  <- as.matrix(fpkm_tbl[, -1])
rownames(fpkm_mat) <- gene_ids

samples_A  <- c(A2780_PAR, A2780_RES)
expr_A_log <- log2(fpkm_mat[, samples_A, drop = FALSE] + 1)

meta_A <- tibble(
  sample_id = samples_A,
  condition = factor(ifelse(grepl("^O-", samples_A), "Parental", "Resistant"),
                     levels = c("Parental", "Resistant"))
)

deg_A_file <- file.path(tables_dir, "de", "GSE153867_A2780_limma_DE_FPKM.tsv")
require_file(deg_A_file)
deg_A <- readr::read_tsv(deg_A_file, show_col_types = FALSE)

# A2780 uses "ID" column (gene symbol from limma)
top_A <- get_top_genes(deg_A, id_col = "ID", fdr_col = "adj.P.Val", logfc_col = "log2FC")

make_heatmap(expr_A_log, meta_A, top_A,
             "GSE153867 — A2780",
             file.path(out_dir, "heatmap_GSE153867_A2780_topDEGs.png"))

# ============================================================
# 2) UWB BRCA-deficient (VST counts)
# ============================================================
message("=== UWB BRCA-def ===")
cts_gz  <- file.path(raw_dir, "GSE235980", "GSE235980_CountReads.txt.gz")
cts_txt <- file.path(raw_dir, "GSE235980", "GSE235980_CountReads.txt")

if (!file.exists(cts_txt) && file.exists(cts_gz)) R.utils::gunzip(cts_gz, overwrite = TRUE)
require_file(cts_txt)

cts_tbl <- readr::read_tsv(cts_txt, show_col_types = FALSE)
gene_ids <- cts_tbl[[1]]
cts_mat  <- as.matrix(cts_tbl[, -1])
rownames(cts_mat) <- gene_ids

samples_def <- c(UWB_DEF_PAR, UWB_DEF_RES)
meta_def <- tibble(
  sample_id = samples_def,
  condition = factor(c(rep("Parental", length(UWB_DEF_PAR)), rep("Resistant", length(UWB_DEF_RES))),
                     levels = c("Parental", "Resistant"))
)

dds_def  <- DESeqDataSetFromMatrix(
  countData = round(cts_mat[, samples_def, drop = FALSE]),
  colData   = as.data.frame(meta_def),
  design    = ~ condition
)
vsd_def  <- vst(dds_def, blind = TRUE)
expr_def_vst <- assay(vsd_def)

deg_def_file <- file.path(tables_dir, "de", "GSE235980_BRCAdef_Res_vs_Par_DESeq2.tsv")
require_file(deg_def_file)
deg_def <- readr::read_tsv(deg_def_file, show_col_types = FALSE)

top_def <- get_top_genes(deg_def, id_col = "gene_id", fdr_col = "padj", logfc_col = "log2FoldChange")

make_heatmap(expr_def_vst, meta_def, top_def,
             "GSE235980 — UWB BRCA-def",
             file.path(out_dir, "heatmap_GSE235980_BRCAdef_topDEGs.png"))

# ============================================================
# 3) PEO1 (TPM)
# ============================================================
message("=== PEO1 ===")
tpm_gz  <- file.path(raw_dir, "GSE117765", "GSE117765_matrix.txt.gz")
tpm_txt <- file.path(raw_dir, "GSE117765", "GSE117765_matrix.txt")

if (!file.exists(tpm_txt) && file.exists(tpm_gz)) R.utils::gunzip(tpm_gz, overwrite = TRUE)
require_file(tpm_txt)

tpm_tbl <- readr::read_tsv(tpm_txt, show_col_types = FALSE)
gene_ids <- tpm_tbl[[1]]
tpm_mat  <- as.matrix(tpm_tbl[, -1])
rownames(tpm_mat) <- gene_ids

sample_ids <- colnames(tpm_mat)
meta_P <- tibble(
  sample_id = sample_ids,
  condition = factor(
    case_when(
      grepl("Adherent", sample_id, ignore.case = TRUE) ~ "Parental",
      grepl("Clone",    sample_id, ignore.case = TRUE) ~ "Resistant",
      TRUE ~ NA_character_
    ),
    levels = c("Parental", "Resistant")
  )
)

expr_P_log <- log2(tpm_mat[, meta_P$sample_id, drop = FALSE] + 1)

deg_P_file <- file.path(tables_dir, "de", "GSE117765_PEO1_limma_TPM_Res_vs_Par.tsv")
require_file(deg_P_file)
deg_P <- readr::read_tsv(deg_P_file, show_col_types = FALSE)

top_P <- get_top_genes(deg_P, id_col = "ID", fdr_col = "adj.P.Val", logfc_col = "log2FC")

make_heatmap(expr_P_log, meta_P, top_P,
             "GSE117765 — PEO1",
             file.path(out_dir, "heatmap_GSE117765_PEO1_topDEGs.png"))

message(">> Heatmaps ComplexHeatmap generados para los 3 datasets.")
stop_log(log_handle)
