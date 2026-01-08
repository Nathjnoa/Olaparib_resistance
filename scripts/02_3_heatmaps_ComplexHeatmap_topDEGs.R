#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(ComplexHeatmap)
  library(circlize)
})

# ------------------------------------------------------------
# Configuración general
# ------------------------------------------------------------
proj_dir <- Sys.getenv("OLAPARIB_RESISTANCE_DIR", unset = "")
if (!nzchar(proj_dir)) {
  proj_dir <- "~/bioinfo/projects/olaparib_resistance"
}
proj_dir <- path.expand(proj_dir)
if (!dir.exists(proj_dir)) {
  stop("No existe el directorio del proyecto: ", proj_dir)
}
raw_dir     <- file.path(proj_dir, "data", "raw")
tables_dir  <- file.path(proj_dir, "results", "tables")
figures_dir <- file.path(proj_dir, "results", "figures", "heatmaps")

dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

fdr_cutoff   <- 0.05
top_n_up     <- 20
top_n_down   <- 20

# paleta común a todos los heatmaps
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b"))

require_file <- function(path) {
  if (!file.exists(path)) stop("No encuentro archivo: ", path)
}

# ------------------------------------------------------------
# Helper: obtener top genes (limma o DESeq2)
# ------------------------------------------------------------
get_top_genes <- function(df,
                          id_col,
                          fdr_col,
                          logfc_col,
                          n_up   = 20,
                          n_down = 20,
                          fdr_cut = 0.05) {

  df_sig <- df %>%
    filter(!is.na(.data[[fdr_col]]),
           .data[[fdr_col]] < fdr_cut)

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

# ------------------------------------------------------------
# Helper: construir heatmap con ComplexHeatmap
# ------------------------------------------------------------
make_heatmap <- function(expr_mat,   # matriz log2-normalizada: genes x muestras
                         meta,       # data.frame con 'sample_id' y 'condition'
                         top_info,   # lista con $genes, $direction
                         dataset_name,
                         outfile_png) {

  # Orden fijo de columnas: Parental -> Resistant
  meta <- meta %>%
    arrange(condition)
  expr_mat <- expr_mat[, meta$sample_id, drop = FALSE]

  genes <- intersect(top_info$genes, rownames(expr_mat))
  if (length(genes) < 5) {
    stop("Muy pocos genes en la intersección para ", dataset_name)
  }

  mat <- expr_mat[genes, , drop = FALSE]

  # Z-score por gen
  mat_z <- t(scale(t(mat)))
  mat_z[is.na(mat_z)] <- 0

  # Ordenar filas por dirección (Up / Down)
  gene_dir <- top_info$direction[genes]
  gene_dir <- factor(gene_dir, levels = c("Up", "Down"))

  # Anotación de columnas (condición)
  ha_col <- HeatmapAnnotation(
    Condition = meta$condition,
    col = list(
      Condition = c(
        Parental  = "#377eb8",
        Resistant = "#e41a1c"
      )
    ),
    annotation_name_gp = gpar(fontsize = 10)
  )

  ht <- Heatmap(
    mat_z,
    name = "Z-score",
    col = col_fun,
    show_row_names = TRUE,
    row_names_side = "right",
    row_names_gp = gpar(fontsize = 6),
    show_column_names = TRUE,
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 8),
    cluster_rows = TRUE,
    cluster_columns = FALSE,  # mantiene orden Parental -> Resistant
    clustering_distance_rows = "pearson",
    row_split = gene_dir,
    top_annotation = ha_col,
    column_title = dataset_name,
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(
      title = "Z-score",
      title_gp = gpar(fontsize = 10),
      labels_gp = gpar(fontsize = 8)
    )
  )

  png(outfile_png, width = 2200, height = 2600, res = 300)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  message("Heatmap guardado en: ", outfile_png)
}

# ============================================================
# 1) GSE153867 - A2780 (FPKM + limma)
# ============================================================
message("=== GSE153867 - A2780 ===")

fpkm_file <- file.path(raw_dir, "GSE153867", "GSE153867_fpkm.txt")
require_file(fpkm_file)
fpkm_tbl  <- readr::read_tsv(fpkm_file, show_col_types = FALSE)

gene_ids <- fpkm_tbl[[1]]
fpkm_mat <- as.matrix(fpkm_tbl[, -1])
rownames(fpkm_mat) <- gene_ids
colnames(fpkm_mat) <- colnames(fpkm_tbl)[-1]

# seleccionar sólo las muestras A2780 (C-* y O-*)
par_ids <- paste0("O-", 1:8)
res_ids <- paste0("C-", 1:8)
samples_A <- c(par_ids, res_ids)  # Parental primero

missing_A <- setdiff(samples_A, colnames(fpkm_mat))
if (length(missing_A) > 0) {
  stop("A2780: faltan columnas en FPKM: ", paste(missing_A, collapse = ", "))
}

expr_A <- fpkm_mat[, samples_A, drop = FALSE]
expr_A_log <- log2(expr_A + 1)

meta_A <- tibble(
  sample_id = samples_A,
  condition = factor(
    ifelse(grepl("^O-", sample_id), "Parental", "Resistant"),
    levels = c("Parental", "Resistant")
  )
)

deg_A_file <- file.path(tables_dir, "GSE153867_A2780_limma_DE_FPKM.tsv")
require_file(deg_A_file)
deg_A <- readr::read_tsv(deg_A_file, show_col_types = FALSE)

top_A <- get_top_genes(
  df        = deg_A,
  id_col    = "ID",          # símbolos de gen
  fdr_col   = "adj.P.Val",
  logfc_col = "log2FC",
  n_up      = top_n_up,
  n_down    = top_n_down,
  fdr_cut   = fdr_cutoff
)

make_heatmap(
  expr_mat    = expr_A_log,
  meta        = meta_A,
  top_info    = top_A,
  dataset_name = "GSE153867 — A2780",
  outfile_png  = file.path(figures_dir, "heatmap_GSE153867_A2780_topDEGs.png")
)

# ============================================================
# 2) GSE235980 - UWB1.289 BRCA1-deficient
# ============================================================
message("=== GSE235980 - UWB BRCA-def ===")

cts_gz  <- file.path(raw_dir, "GSE235980", "GSE235980_CountReads.txt.gz")
cts_txt <- file.path(raw_dir, "GSE235980", "GSE235980_CountReads.txt")

if (!file.exists(cts_txt) && file.exists(cts_gz)) {
  R.utils::gunzip(cts_gz, overwrite = TRUE)
}
require_file(cts_txt)
cts_tbl <- readr::read_tsv(cts_txt, show_col_types = FALSE)

gene_ids <- cts_tbl[[1]]
cts_mat  <- as.matrix(cts_tbl[, -1])
rownames(cts_mat) <- gene_ids
colnames(cts_mat) <- colnames(cts_tbl)[-1]

# Parental primero, luego Resistant
par_def <- c("AN14028320", "AN14028322")
res_def <- c("AN14028316", "AN14028318")
samples_def <- c(par_def, res_def)

missing_def <- setdiff(samples_def, colnames(cts_mat))
if (length(missing_def) > 0) {
  stop("UWB BRCA-def: faltan columnas en counts: ", paste(missing_def, collapse = ", "))
}

meta_def <- tibble(
  sample_id = samples_def,
  condition = factor(
    c(rep("Parental", length(par_def)),
      rep("Resistant", length(res_def))),
    levels = c("Parental", "Resistant")
  )
)

dds_def <- DESeqDataSetFromMatrix(
  countData = round(cts_mat[, samples_def, drop = FALSE]),
  colData   = as.data.frame(meta_def),
  design    = ~ condition
)

vsd_def <- vst(dds_def, blind = TRUE)
expr_def_vst <- assay(vsd_def)

deg_def_file <- file.path(tables_dir, "GSE235980_BRCAdef_Res_vs_Par_DESeq2.tsv")
require_file(deg_def_file)
deg_def <- readr::read_tsv(deg_def_file, show_col_types = FALSE)

top_def <- get_top_genes(
  df        = deg_def,
  id_col    = "gene_id",
  fdr_col   = "padj",
  logfc_col = "log2FoldChange",
  n_up      = top_n_up,
  n_down    = top_n_down,
  fdr_cut   = fdr_cutoff
)

make_heatmap(
  expr_mat     = expr_def_vst,
  meta         = meta_def,
  top_info     = top_def,
  dataset_name = "GSE235980 — UWB BRCA-def",
  outfile_png  = file.path(figures_dir, "heatmap_GSE235980_BRCAdef_topDEGs.png")
)

# ============================================================
# 3) GSE235980 - UWB1.289 BRCA1-proficient
# ============================================================
message("=== GSE235980 - UWB BRCA-prof ===")

par_prof <- c("AN14028328", "AN14028330")
res_prof <- c("AN14028324", "AN14028326")
samples_prof <- c(par_prof, res_prof)

missing_prof <- setdiff(samples_prof, colnames(cts_mat))
if (length(missing_prof) > 0) {
  stop("UWB BRCA-prof: faltan columnas en counts: ", paste(missing_prof, collapse = ", "))
}

meta_prof <- tibble(
  sample_id = samples_prof,
  condition = factor(
    c(rep("Parental", length(par_prof)),
      rep("Resistant", length(res_prof))),
    levels = c("Parental", "Resistant")
  )
)

dds_prof <- DESeqDataSetFromMatrix(
  countData = round(cts_mat[, samples_prof, drop = FALSE]),
  colData   = as.data.frame(meta_prof),
  design    = ~ condition
)

vsd_prof <- vst(dds_prof, blind = TRUE)
expr_prof_vst <- assay(vsd_prof)

deg_prof_file <- file.path(tables_dir, "GSE235980_BRCAprof_Res_vs_Par_DESeq2.tsv")
require_file(deg_prof_file)
deg_prof <- readr::read_tsv(deg_prof_file, show_col_types = FALSE)

top_prof <- get_top_genes(
  df        = deg_prof,
  id_col    = "gene_id",
  fdr_col   = "padj",
  logfc_col = "log2FoldChange",
  n_up      = top_n_up,
  n_down    = top_n_down,
  fdr_cut   = fdr_cutoff
)

make_heatmap(
  expr_mat     = expr_prof_vst,
  meta         = meta_prof,
  top_info     = top_prof,
  dataset_name = "GSE235980 — UWB BRCA-prof",
  outfile_png  = file.path(figures_dir, "heatmap_GSE235980_BRCAprof_topDEGs.png")
)

# ============================================================
# 4) GSE117765 - PEO1 (TPM + limma)
# ============================================================
message("=== GSE117765 - PEO1 ===")

tpm_gz  <- file.path(raw_dir, "GSE117765", "GSE117765_matrix.txt.gz")
tpm_txt <- file.path(raw_dir, "GSE117765", "GSE117765_matrix.txt")

if (!file.exists(tpm_txt) && file.exists(tpm_gz)) {
  R.utils::gunzip(tpm_gz, overwrite = TRUE)
}
require_file(tpm_txt)
tpm_tbl <- readr::read_tsv(tpm_txt, show_col_types = FALSE)

gene_ids <- tpm_tbl[[1]]
tpm_mat  <- as.matrix(tpm_tbl[, -1])
rownames(tpm_mat) <- gene_ids
colnames(tpm_mat) <- colnames(tpm_tbl)[-1]

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
) %>%
  arrange(condition)

expr_P <- tpm_mat[, meta_P$sample_id, drop = FALSE]
expr_P_log <- log2(expr_P + 1)

deg_P_file <- file.path(tables_dir, "GSE117765_PEO1_limma_TPM_Res_vs_Par.tsv")
require_file(deg_P_file)
deg_P <- readr::read_tsv(deg_P_file, show_col_types = FALSE)

top_P <- get_top_genes(
  df        = deg_P,
  id_col    = "ID",          # símbolos
  fdr_col   = "adj.P.Val",
  logfc_col = "log2FC",
  n_up      = top_n_up,
  n_down    = top_n_down,
  fdr_cut   = fdr_cutoff
)

make_heatmap(
  expr_mat     = expr_P_log,
  meta         = meta_P,
  top_info     = top_P,
  dataset_name = "GSE117765 — PEO1",
  outfile_png  = file.path(figures_dir, "heatmap_GSE117765_PEO1_topDEGs.png")
)

message(">> Heatmaps de ComplexHeatmap generados para los cuatro datasets.")
