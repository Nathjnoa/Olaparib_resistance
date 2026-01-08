#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(org.Hs.eg.db)
})

# ============================================================
# Paths
# ============================================================
proj_dir <- Sys.getenv("OLAPARIB_RESISTANCE_DIR", unset = "")
if (!nzchar(proj_dir)) {
  proj_dir <- "~/bioinfo/projects/olaparib_resistance"
}
proj_dir <- path.expand(proj_dir)
if (!dir.exists(proj_dir)) {
  stop("No existe el directorio del proyecto: ", proj_dir)
}
tables_dir <- file.path(proj_dir, "results", "tables")
out_dir    <- file.path(tables_dir, "integrated")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_long <- file.path(out_dir, "meta_DE_input_long.tsv")
out_wide <- file.path(out_dir, "meta_DE_input_wide.tsv")

# ============================================================
# CONFIG: opcional, por si algún dataset tiene el contraste invertido
# (Resistant - Parental) debe ser POSITIVO cuando "sube en resistente".
# ============================================================
flip_log2FC <- c(
  A2780        = FALSE,
  PEO1         = FALSE,
  UWB_BRCAdef  = FALSE,
  UWB_BRCAprof = FALSE
)

# ============================================================
# Helpers
# ============================================================

.pick_first <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

# ---------- limma loader (flexible) ----------
load_limma_deg <- function(file, dataset_id) {
  message(">> Cargando limma: ", dataset_id, " | ", file)
  df <- readr::read_tsv(file, show_col_types = FALSE)

  # Columnas posibles
  gene_col <- intersect(c("ID", "gene", "Gene", "symbol", "SYMBOL"), colnames(df)) |> .pick_first()
  logfc_col <- intersect(c("log2FC", "logFC"), colnames(df)) |> .pick_first()
  t_col <- intersect(c("t_stat", "t", "tstat"), colnames(df)) |> .pick_first()
  p_col <- intersect(c("P.Value", "pvalue", "pval"), colnames(df)) |> .pick_first()
  padj_col <- intersect(c("adj.P.Val", "padj", "FDR"), colnames(df)) |> .pick_first()

  needed <- c(gene_col, logfc_col, t_col, p_col, padj_col)
  if (any(is.na(needed))) {
    stop(
      "Faltan columnas esperadas en limma para ", dataset_id, ".\n",
      "Encontradas: ", paste(colnames(df), collapse = ", "), "\n",
      "Necesarias (alguna variante): gene/logFC/t/p/adj"
    )
  }

  df2 <- df %>%
    transmute(
      gene_symbol = as.character(.data[[gene_col]]),
      log2FC      = as.numeric(.data[[logfc_col]]),
      t_stat      = as.numeric(.data[[t_col]]),
      pval        = as.numeric(.data[[p_col]]),
      padj        = as.numeric(.data[[padj_col]]),
      dataset     = dataset_id
    ) %>%
    filter(!is.na(gene_symbol), gene_symbol != "")

  # SE aproximado desde t (si no tienes stdev.unscaled + s2.post)
  df2 <- df2 %>%
    mutate(
      SE = if_else(is.na(t_stat) | t_stat == 0, NA_real_, abs(log2FC / t_stat))
    ) %>%
    dplyr::select(gene_symbol, log2FC, SE, pval, padj, dataset)

  # Resolver símbolos duplicados: quedarnos con menor padj
  df2 <- df2 %>%
    group_by(gene_symbol) %>%
    slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
    ungroup()

  df2
}

# ---------- DESeq2 loader (flexible + ENSEMBL->SYMBOL) ----------
load_deseq_deg <- function(file, dataset_id) {
  message(">> Cargando DESeq2: ", dataset_id, " | ", file)
  df <- readr::read_tsv(file, show_col_types = FALSE)

  gene_col <- intersect(c("gene_id", "GeneID", "gene", "Gene", "id"), colnames(df)) |> .pick_first()
  lfc_col  <- intersect(c("log2FoldChange", "logFC", "log2FC"), colnames(df)) |> .pick_first()
  se_col   <- intersect(c("lfcSE", "SE"), colnames(df)) |> .pick_first()
  p_col    <- intersect(c("pvalue", "pval", "P.Value"), colnames(df)) |> .pick_first()
  padj_col <- intersect(c("padj", "adj.P.Val", "FDR"), colnames(df)) |> .pick_first()

  needed <- c(gene_col, lfc_col, se_col, p_col, padj_col)
  if (any(is.na(needed))) {
    stop(
      "Faltan columnas esperadas en DESeq2 para ", dataset_id, ".\n",
      "Encontradas: ", paste(colnames(df), collapse = ", "), "\n",
      "Necesarias (alguna variante): gene/log2FC/lfcSE/p/padj"
    )
  }

  gene_vec <- as.character(df[[gene_col]])
  is_ensg  <- all(grepl("^ENS", gene_vec))

  if (is_ensg) {
    df <- df %>% mutate(ensg = sub("\\.\\d+$", "", .data[[gene_col]]))

    map <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys    = unique(df$ensg),
      keytype = "ENSEMBL",
      columns = "SYMBOL"
    )

    # Si hay 1->many, nos quedamos con el primer símbolo por ENSEMBL
    map <- map %>%
      arrange(ENSEMBL, SYMBOL) %>%
      group_by(ENSEMBL) %>%
      summarise(SYMBOL = .pick_first(SYMBOL), .groups = "drop")

    df2 <- df %>%
      left_join(map, by = c("ensg" = "ENSEMBL")) %>%
      rename(gene_symbol = SYMBOL)
  } else {
    df2 <- df %>% mutate(gene_symbol = .data[[gene_col]])
  }

  df2 <- df2 %>%
    transmute(
      gene_symbol = as.character(gene_symbol),
      log2FC      = as.numeric(.data[[lfc_col]]),
      SE          = as.numeric(.data[[se_col]]),
      pval        = as.numeric(.data[[p_col]]),
      padj        = as.numeric(.data[[padj_col]]),
      dataset     = dataset_id
    ) %>%
    filter(!is.na(gene_symbol), gene_symbol != "")

  # Resolver símbolos duplicados: menor padj
  df2 <- df2 %>%
    group_by(gene_symbol) %>%
    slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
    ungroup()

  df2
}

# ============================================================
# Load the 4 datasets (mismos paths que tu script original)
# ============================================================
file_A2780   <- file.path(tables_dir, "GSE153867_A2780_limma_DE_FPKM.tsv")
file_PEO1    <- file.path(tables_dir, "GSE117765_PEO1_limma_TPM_Res_vs_Par.tsv")
file_UWBdef  <- file.path(tables_dir, "GSE235980_BRCAdef_Res_vs_Par_DESeq2.tsv")
file_UWBprof <- file.path(tables_dir, "GSE235980_BRCAprof_Res_vs_Par_DESeq2.tsv")

stopifnot(file.exists(file_A2780),
          file.exists(file_PEO1),
          file.exists(file_UWBdef),
          file.exists(file_UWBprof))

deg_A2780   <- load_limma_deg(file_A2780,  "A2780")
deg_PEO1    <- load_limma_deg(file_PEO1,   "PEO1")
deg_UWBdef  <- load_deseq_deg(file_UWBdef, "UWB_BRCAdef")
deg_UWBprof <- load_deseq_deg(file_UWBprof,"UWB_BRCAprof")

deg_long <- bind_rows(deg_A2780, deg_PEO1, deg_UWBdef, deg_UWBprof) %>%
  mutate(
    # aplicar flips si se requiere
    log2FC = if_else(flip_log2FC[dataset], -log2FC, log2FC),
    contrast = "Resistant_minus_Parental"
  ) %>%
  arrange(gene_symbol, dataset)

# ============================================================
# QC rápido
# ============================================================
qc <- deg_long %>%
  group_by(dataset) %>%
  summarise(
    n_genes = n(),
    n_ok_se = sum(!is.na(SE) & SE > 0),
    n_fdr_01 = sum(!is.na(padj) & padj < 0.10),
    up_fdr_01 = sum(!is.na(padj) & padj < 0.10 & log2FC > 0),
    down_fdr_01 = sum(!is.na(padj) & padj < 0.10 & log2FC < 0),
    median_abs_log2FC = median(abs(log2FC), na.rm = TRUE),
    .groups = "drop"
  )

message("---- QC resumen ----")
print(qc)

message("Total filas en meta_input_long: ", nrow(deg_long))
message("Total genes únicos: ", n_distinct(deg_long$gene_symbol))

readr::write_tsv(deg_long, out_long)
message(">> Guardado: ", out_long)

# ============================================================
# Formato WIDE
# ============================================================
deg_wide <- deg_long %>%
  dplyr::select(gene_symbol, dataset, log2FC, SE, pval, padj) %>%
  tidyr::pivot_longer(
    cols = c(log2FC, SE, pval, padj),
    names_to = "stat",
    values_to = "value"
  ) %>%
  tidyr::unite("colname", stat, dataset) %>%
  tidyr::pivot_wider(names_from = colname, values_from = value)

readr::write_tsv(deg_wide, out_wide)
message(">> Guardado: ", out_wide)

message(">> Construcción de tablas de entrada para metaanálisis COMPLETADA.")
