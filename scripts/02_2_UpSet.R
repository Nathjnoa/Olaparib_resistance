#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(UpSetR)  # install.packages("UpSetR") si hace falta
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
figures_dir <- file.path(proj_dir, "results", "figures")
dir.create(file.path(figures_dir, "upset"), showWarnings = FALSE, recursive = TRUE)

fdr_cutoff   <- 0.05
logfc_cutoff <- 0.584  # ~1.5x

# ---------------------------
# 1) Cargar tablas de DE
# ---------------------------

file_A      <- file.path(tables_dir, "GSE153867_A2780_limma_DE_FPKM.tsv")
file_U_def  <- file.path(tables_dir, "GSE235980_BRCAdef_Res_vs_Par_DESeq2.tsv")
file_U_prof <- file.path(tables_dir, "GSE235980_BRCAprof_Res_vs_Par_DESeq2.tsv")
file_P      <- file.path(tables_dir, "GSE117765_PEO1_limma_TPM_Res_vs_Par.tsv")

missing_files <- c(file_A, file_U_def, file_U_prof, file_P)
missing_files <- missing_files[!file.exists(missing_files)]
if (length(missing_files) > 0) {
  stop("No se encontraron archivos de entrada:\n", paste(missing_files, collapse = "\n"))
}

A     <- readr::read_tsv(file_A)      %>% mutate(gene_id = as.character(gene_id))
U_def <- readr::read_tsv(file_U_def)  %>% mutate(gene_id = as.character(gene_id))
U_prof<- readr::read_tsv(file_U_prof) %>% mutate(gene_id = as.character(gene_id))
P     <- readr::read_tsv(file_P)      %>% mutate(gene_id = as.character(gene_id))

# ---------------------------
# 2) Listas de DEGs (UP / DOWN)
# ---------------------------

# A2780 (limma)
A_up <- A %>%
  filter(!is.na(adj.P.Val), !is.na(log2FC),
         adj.P.Val < fdr_cutoff, log2FC >  logfc_cutoff) %>%
  pull(gene_id) %>% unique()

A_down <- A %>%
  filter(!is.na(adj.P.Val), !is.na(log2FC),
         adj.P.Val < fdr_cutoff, log2FC < -logfc_cutoff) %>%
  pull(gene_id) %>% unique()

# UWB (DESeq2) - unión def + prof
U_def_up <- U_def %>%
  filter(!is.na(padj), !is.na(log2FoldChange),
         padj < fdr_cutoff, log2FoldChange >  logfc_cutoff) %>%
  pull(gene_id) %>% unique()

U_def_down <- U_def %>%
  filter(!is.na(padj), !is.na(log2FoldChange),
         padj < fdr_cutoff, log2FoldChange < -logfc_cutoff) %>%
  pull(gene_id) %>% unique()

U_prof_up <- U_prof %>%
  filter(!is.na(padj), !is.na(log2FoldChange),
         padj < fdr_cutoff, log2FoldChange >  logfc_cutoff) %>%
  pull(gene_id) %>% unique()

U_prof_down <- U_prof %>%
  filter(!is.na(padj), !is.na(log2FoldChange),
         padj < fdr_cutoff, log2FoldChange < -logfc_cutoff) %>%
  pull(gene_id) %>% unique()

U_up   <- union(U_def_up,   U_prof_up)
U_down <- union(U_def_down, U_prof_down)

# PEO1 (limma)
P_up <- P %>%
  filter(!is.na(adj.P.Val), !is.na(log2FC),
         adj.P.Val < fdr_cutoff, log2FC >  logfc_cutoff) %>%
  pull(gene_id) %>% unique()

P_down <- P %>%
  filter(!is.na(adj.P.Val), !is.na(log2FC),
         adj.P.Val < fdr_cutoff, log2FC < -logfc_cutoff) %>%
  pull(gene_id) %>% unique()

# ---------------------------
# 3) UpSet usando fromList()
# ---------------------------

# Lista de genes UP por dataset
up_list <- list(
  A2780 = A_up,
  UWB   = U_up,
  PEO1  = P_up
)

# Lista de genes DOWN por dataset
down_list <- list(
  A2780 = A_down,
  UWB   = U_down,
  PEO1  = P_down
)

# Convertir listas a matriz de presencia/ausencia que entiende UpSetR
up_mat   <- UpSetR::fromList(up_list)
down_mat <- UpSetR::fromList(down_list)

# UP
png(file.path(figures_dir, "upset/DEGs_UP_UpSet_A2780_UWB_PEO1.png"),
    width = 2000, height = 1500, res = 300)
upset(up_mat,
      nsets = 3,
      nintersects = NA,
      sets = c("A2780", "UWB", "PEO1"),
      order.by = "freq",
      main.bar.color = "grey30",
      sets.bar.color = "grey40",
      mainbar.y.label = "Genes en intersección",
      sets.x.label = "Genes por dataset",
      text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.2))
dev.off()

# DOWN
png(file.path(figures_dir, "upset/DEGs_DOWN_UpSet_A2780_UWB_PEO1.png"),
    width = 2000, height = 1500, res = 300)
upset(down_mat,
      nsets = 3,
      nintersects = NA,
      sets = c("A2780", "UWB", "PEO1"),
      order.by = "freq",
      main.bar.color = "grey30",
      sets.bar.color = "grey40",
      mainbar.y.label = "Genes en intersección",
      sets.x.label = "Genes por dataset",
      text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.2))
dev.off()

message(">> UpSet plots generados correctamente.")
