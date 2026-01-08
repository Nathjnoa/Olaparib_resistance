#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(GSVA)
  library(msigdbr)
  library(DESeq2)
  library(R.utils)
})

# =====================================================
# CONFIG paths
# =====================================================
proj_dir <- Sys.getenv("OLAPARIB_RESISTANCE_DIR", unset = "")
if (!nzchar(proj_dir)) {
  proj_dir <- "~/bioinfo/projects/olaparib_resistance"
}
proj_dir <- path.expand(proj_dir)
if (!dir.exists(proj_dir)) {
  stop("No existe el directorio del proyecto: ", proj_dir)
}

raw_dir       <- file.path(proj_dir, "data", "raw")
processed_dir <- file.path(proj_dir, "data", "processed")

gsva_tab_dir  <- file.path(proj_dir, "results", "tables", "gsva")
dir.create(gsva_tab_dir, showWarnings = FALSE, recursive = TRUE)

# =====================================================
# Helpers
# =====================================================
load_expr_matrix <- function(path, gene_col = 1) {
  if (!file.exists(path)) stop("No existe: ", path)
  message(">> Cargando: ", path)
  df <- readr::read_tsv(path, show_col_types = FALSE)

  gene_col_name <- gene_col
  if (is.numeric(gene_col)) gene_col_name <- colnames(df)[gene_col]

  if (is.na(gene_col_name) || !(gene_col_name %in% colnames(df))) {
    stop("No se encontró columna de gen en: ", path)
  }

  colnames(df)[colnames(df) == gene_col_name] <- "Gene"
  df <- df %>% distinct(Gene, .keep_all = TRUE)

  mat <- as.matrix(df[, -1, drop = FALSE])
  rownames(mat) <- df$Gene
  mat
}

safe_gunzip <- function(gz_path, out_path) {
  if (!file.exists(out_path) && file.exists(gz_path)) {
    message(">> Descomprimiendo: ", gz_path)
    R.utils::gunzip(gz_path, overwrite = TRUE, remove = FALSE)
  }
  if (!file.exists(out_path)) stop("No existe: ", out_path)
}

# Filtro robusto para GSVA: finitos + varianza > 0
clean_for_gsva <- function(mat, dataset_name) {
  keep_finite <- apply(mat, 1, function(x) all(is.finite(x)))
  keep_var    <- apply(mat[keep_finite, , drop = FALSE], 1, sd, na.rm = TRUE) > 0

  n_drop_finite <- sum(!keep_finite)
  n_drop_var    <- sum(keep_finite) - sum(keep_var)

  if (n_drop_finite > 0) message(dataset_name, ": genes excluidos por NA/Inf = ", n_drop_finite)
  if (n_drop_var > 0)    message(dataset_name, ": genes excluidos por varianza 0 = ", n_drop_var)

  mat2 <- mat[keep_finite, , drop = FALSE]
  mat2 <- mat2[keep_var, , drop = FALSE]
  mat2
}

run_gsva_hallmarks <- function(mat, pathways, dataset_name) {
  mat2 <- clean_for_gsva(mat, dataset_name)

  param <- gsvaParam(
    exprData = mat2,
    geneSets = pathways,
    kcdf     = "Gaussian"
  )

  gsva(param)
}

# =====================================================
# 1) A2780 (usar processed: 16 muestras O-* y C-*)
# =====================================================
file_A_proc <- file.path(processed_dir, "GSE153867_A2780_FPKM.tsv")
mat_A_raw <- load_expr_matrix(file_A_proc, gene_col = "gene_id")

# log2(FPKM+1)
mat_A <- log2(mat_A_raw + 1)

col_par_A <- grep("^O-", colnames(mat_A), value = TRUE)
col_res_A <- grep("^C-", colnames(mat_A), value = TRUE)

if (length(col_par_A) == 0 || length(col_res_A) == 0) {
  stop("A2780: no detecté columnas O-* (Parental) y/o C-* (Resistant). Revisa colnames().")
}

mat_A <- mat_A[, c(col_par_A, col_res_A), drop = FALSE]
message("A2780 dims: ", nrow(mat_A), " x ", ncol(mat_A))

meta_A <- tibble(
  sample_id = colnames(mat_A),
  dataset   = "A2780",
  condition = factor(ifelse(grepl("^O-", sample_id), "Parental", "Resistant"),
                     levels = c("Parental", "Resistant"))
)

# =====================================================
# 2) PEO1 (TPM -> log2)
# =====================================================
file_P_gz  <- file.path(raw_dir, "GSE117765", "GSE117765_matrix.txt.gz")
file_P_txt <- file.path(raw_dir, "GSE117765", "GSE117765_matrix.txt")
safe_gunzip(file_P_gz, file_P_txt)

mat_P_raw <- load_expr_matrix(file_P_txt, gene_col = 1)
mat_P <- log2(mat_P_raw + 1)

col_par_P <- grep("Adherent", colnames(mat_P), value = TRUE)
col_res_P <- grep("Clone",    colnames(mat_P), value = TRUE)

if (length(col_par_P) == 0 || length(col_res_P) == 0) {
  stop("PEO1: no detecté columnas Adherent (Parental) y/o Clone (Resistant).")
}

mat_P <- mat_P[, c(col_par_P, col_res_P), drop = FALSE]
message("PEO1 dims: ", nrow(mat_P), " x ", ncol(mat_P))

meta_P <- tibble(
  sample_id = colnames(mat_P),
  dataset   = "PEO1",
  condition = factor(ifelse(grepl("Adherent", sample_id, ignore.case = TRUE),
                            "Parental", "Resistant"),
                     levels = c("Parental", "Resistant"))
)

# =====================================================
# 3) UWB (counts -> VST por subset)
# =====================================================
file_U_gz  <- file.path(raw_dir, "GSE235980", "GSE235980_CountReads.txt.gz")
file_U_txt <- file.path(raw_dir, "GSE235980", "GSE235980_CountReads.txt")
safe_gunzip(file_U_gz, file_U_txt)

uwb_raw <- readr::read_tsv(file_U_txt, show_col_types = FALSE)
colnames(uwb_raw)[1] <- "Gene"
uwb_raw <- uwb_raw %>% distinct(Gene, .keep_all = TRUE)

count_mat <- as.matrix(uwb_raw[, -1, drop = FALSE])
rownames(count_mat) <- uwb_raw$Gene

col_def  <- c("AN14028316", "AN14028318", "AN14028320", "AN14028322")
col_prof <- c("AN14028324", "AN14028326", "AN14028328", "AN14028330")

missing_def  <- setdiff(col_def,  colnames(count_mat))
missing_prof <- setdiff(col_prof, colnames(count_mat))
if (length(missing_def) > 0 || length(missing_prof) > 0) {
  stop("UWB: faltan columnas.\nDef: ", paste(missing_def, collapse = ", "),
       "\nProf: ", paste(missing_prof, collapse = ", "))
}

# condición real para metadata (Parental/Resistant)
cond_def  <- c("Resistant", "Resistant", "Parental", "Parental")
cond_prof <- c("Resistant", "Resistant", "Parental", "Parental")

meta_def <- tibble(
  sample_id = col_def,
  dataset   = "UWB_BRCAdef",
  BRCA      = "deficient",
  condition = factor(cond_def, levels = c("Parental", "Resistant"))
)

meta_prof <- tibble(
  sample_id = col_prof,
  dataset   = "UWB_BRCAprof",
  BRCA      = "proficient",
  condition = factor(cond_prof, levels = c("Parental", "Resistant"))
)

make_vst <- function(counts_sub, meta_sub, dataset_name) {
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_sub),
    colData   = as.data.frame(meta_sub %>% column_to_rownames("sample_id")),
    design    = ~ 1
  )
  vsd <- vst(dds, blind = TRUE)
  assay(vsd)
}

vst_def  <- make_vst(count_mat[, col_def,  drop = FALSE], meta_def,  "UWB_BRCAdef")
vst_prof <- make_vst(count_mat[, col_prof, drop = FALSE], meta_prof, "UWB_BRCAprof")

# reordenar columnas Parental -> Resistant en VST
vst_def  <- vst_def[,  meta_def %>% arrange(condition)  %>% pull(sample_id), drop = FALSE]
vst_prof <- vst_prof[, meta_prof %>% arrange(condition) %>% pull(sample_id), drop = FALSE]

message("UWB BRCAdef dims:  ", nrow(vst_def),  " x ", ncol(vst_def))
message("UWB BRCAprof dims: ", nrow(vst_prof), " x ", ncol(vst_prof))

# =====================================================
# 4) Genes comunes (para integración consistente)
# =====================================================
genes_common <- Reduce(intersect, list(
  rownames(mat_A),
  rownames(mat_P),
  rownames(vst_def),
  rownames(vst_prof)
))

message("Genes comunes entre datasets: ", length(genes_common))

mat_A    <- mat_A[genes_common, , drop = FALSE]
mat_P    <- mat_P[genes_common, , drop = FALSE]
vst_def  <- vst_def[genes_common, , drop = FALSE]
vst_prof <- vst_prof[genes_common, , drop = FALSE]

# =====================================================
# 5) MSigDB Hallmarks (50)
# =====================================================
msig <- msigdbr(species = "Homo sapiens", collection = "H")
df_h <- msig %>% dplyr::select(gs_name, gene_symbol) %>% distinct()
pathways <- split(df_h$gene_symbol, df_h$gs_name)
message("Número de pathways (Hallmarks): ", length(pathways))

# =====================================================
# 6) Ejecutar GSVA
# =====================================================
message(">> Ejecutando GSVA A2780 ...")
gsva_A <- run_gsva_hallmarks(mat_A, pathways, "A2780")

message(">> Ejecutando GSVA PEO1 ...")
gsva_P <- run_gsva_hallmarks(mat_P, pathways, "PEO1")

message(">> Ejecutando GSVA UWB BRCA-def ...")
gsva_def <- run_gsva_hallmarks(vst_def, pathways, "UWB_BRCAdef")

message(">> Ejecutando GSVA UWB BRCA-prof ...")
gsva_prof <- run_gsva_hallmarks(vst_prof, pathways, "UWB_BRCAprof")

# =====================================================
# 7) Guardar resultados + metadata
# =====================================================
readr::write_tsv(as.data.frame(gsva_A)    %>% rownames_to_column("Pathway"),
                 file.path(gsva_tab_dir, "GSVA_A2780.tsv"))
readr::write_tsv(as.data.frame(gsva_P)    %>% rownames_to_column("Pathway"),
                 file.path(gsva_tab_dir, "GSVA_PEO1.tsv"))
readr::write_tsv(as.data.frame(gsva_def)  %>% rownames_to_column("Pathway"),
                 file.path(gsva_tab_dir, "GSVA_UWB_BRCAdef.tsv"))
readr::write_tsv(as.data.frame(gsva_prof) %>% rownames_to_column("Pathway"),
                 file.path(gsva_tab_dir, "GSVA_UWB_BRCAprof.tsv"))

gsva_merged <- cbind(gsva_A, gsva_P, gsva_def, gsva_prof)
readr::write_tsv(as.data.frame(gsva_merged) %>% rownames_to_column("Pathway"),
                 file.path(gsva_tab_dir, "GSVA_ALL_merged.tsv"))

meta_all <- bind_rows(meta_A, meta_P, meta_def, meta_prof)
readr::write_tsv(meta_all, file.path(gsva_tab_dir, "GSVA_metadata.tsv"))

message(">> GSVA COMPLETADO. Outputs en: ", gsva_tab_dir)
