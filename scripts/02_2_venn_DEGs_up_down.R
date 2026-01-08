#!/usr/bin/env Rscript

# ============================================================
# Venn diagrams de DEGs entre A2780, UWB (def+prof) y PEO1
# Usa:
#  - GSE153867_A2780_limma_DE_FPKM.tsv
#  - GSE235980_BRCAdef_Res_vs_Par_DESeq2.tsv
#  - GSE235980_BRCAprof_Res_vs_Par_DESeq2.tsv
#  - GSE117765_PEO1_limma_TPM_Res_vs_Par.tsv
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggVennDiagram)  # install.packages("ggVennDiagram") si hace falta
  library(ggplot2)
})

proj_dir <- Sys.getenv("OLAPARIB_RESISTANCE_DIR", unset = "")
if (!nzchar(proj_dir)) {
  proj_dir <- "~/bioinfo/projects/olaparib_resistance"
}
proj_dir <- path.expand(proj_dir)
if (!dir.exists(proj_dir)) {
  stop("No existe el directorio del proyecto: ", proj_dir)
}
tables_dir  <- file.path(proj_dir, "results", "tables")
figures_dir <- file.path(proj_dir, "results", "figures")

dir.create(file.path(figures_dir, "venn"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(tables_dir,  "venn"), showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Parámetros globales
# ------------------------------------------------------------
fdr_cutoff   <- 0.05
logfc_cutoff <- 0.584  # FC >= 1.5

venn_palette <- c(
  "A2780" = "#1f78b4",  # azul
  "UWB"   = "#33a02c",  # verde
  "PEO1"  = "#e31a1c"   # rojo
)

# ============================================================
# 1) Cargar resultados de cada dataset
# ============================================================

# A2780 - GSE153867 (limma sobre FPKM)
file_A <- file.path(tables_dir, "GSE153867_A2780_limma_DE_FPKM.tsv")
if (!file.exists(file_A)) stop("No encuentro: ", file_A)
A <- readr::read_tsv(file_A) %>%
  mutate(gene_id = as.character(gene_id))

# UWB - GSE235980 (DESeq2)
file_U_def  <- file.path(tables_dir, "GSE235980_BRCAdef_Res_vs_Par_DESeq2.tsv")
file_U_prof <- file.path(tables_dir, "GSE235980_BRCAprof_Res_vs_Par_DESeq2.tsv")
if (!file.exists(file_U_def)) stop("No encuentro: ", file_U_def)
if (!file.exists(file_U_prof)) stop("No encuentro: ", file_U_prof)
U_def  <- readr::read_tsv(file_U_def)  %>% mutate(gene_id = as.character(gene_id))
U_prof <- readr::read_tsv(file_U_prof) %>% mutate(gene_id = as.character(gene_id))

# PEO1 - GSE117765 (limma sobre TPM)
file_P <- file.path(tables_dir, "GSE117765_PEO1_limma_TPM_Res_vs_Par.tsv")
if (!file.exists(file_P)) stop("No encuentro: ", file_P)
P <- readr::read_tsv(file_P) %>%
  mutate(gene_id = as.character(gene_id))

# ============================================================
# 2) Definir listas de DEGs (up / down) por dataset
# ============================================================

# A2780 (limma: adj.P.Val, log2FC)
A_up <- A %>%
  filter(!is.na(adj.P.Val), !is.na(log2FC),
         adj.P.Val < fdr_cutoff, log2FC >  logfc_cutoff) %>%
  pull(gene_id) %>%
  unique()

A_down <- A %>%
  filter(!is.na(adj.P.Val), !is.na(log2FC),
         adj.P.Val < fdr_cutoff, log2FC < -logfc_cutoff) %>%
  pull(gene_id) %>%
  unique()

# UWB (DESeq2: padj, log2FoldChange)
U_def_up <- U_def %>%
  filter(!is.na(padj), !is.na(log2FoldChange),
         padj < fdr_cutoff, log2FoldChange >  logfc_cutoff) %>%
  pull(gene_id) %>%
  unique()

U_def_down <- U_def %>%
  filter(!is.na(padj), !is.na(log2FoldChange),
         padj < fdr_cutoff, log2FoldChange < -logfc_cutoff) %>%
  pull(gene_id) %>%
  unique()

U_prof_up <- U_prof %>%
  filter(!is.na(padj), !is.na(log2FoldChange),
         padj < fdr_cutoff, log2FoldChange >  logfc_cutoff) %>%
  pull(gene_id) %>%
  unique()

U_prof_down <- U_prof %>%
  filter(!is.na(padj), !is.na(log2FoldChange),
         padj < fdr_cutoff, log2FoldChange < -logfc_cutoff) %>%
  pull(gene_id) %>%
  unique()

# UWB (unión de genes DE en BRCA-deficient y BRCA-proficient)
U_up   <- union(U_def_up,   U_prof_up)
U_down <- union(U_def_down, U_prof_down)

# PEO1 (limma: adj.P.Val, log2FC)
P_up <- P %>%
  filter(!is.na(adj.P.Val), !is.na(log2FC),
         adj.P.Val < fdr_cutoff, log2FC >  logfc_cutoff) %>%
  pull(gene_id) %>%
  unique()

P_down <- P %>%
  filter(!is.na(adj.P.Val), !is.na(log2FC),
         adj.P.Val < fdr_cutoff, log2FC < -logfc_cutoff) %>%
  pull(gene_id) %>%
  unique()

# Guardar tablas con conteos básicos
summary_list <- tibble(
  dataset = c("A2780_up", "A2780_down",
              "UWB_up", "UWB_down",
              "PEO1_up", "PEO1_down"),
  n_genes = c(length(A_up), length(A_down),
              length(U_up), length(U_down),
              length(P_up), length(P_down))
)

readr::write_tsv(summary_list,
                 file.path(tables_dir, "venn/DEG_counts_for_venn.tsv"))

# ============================================================
# 3) Venn de genes UP
# ============================================================

venn_up_list <- list(
  A2780 = A_up,
  UWB   = U_up,
  PEO1  = P_up
)

p_venn_up <- ggVennDiagram(
  venn_up_list,
  set_color   = venn_palette[names(venn_up_list)],
  label_alpha = 0,
  label = "count"
) +
  scale_fill_gradient(low = "white", high = "grey80") +
  scale_color_manual(values = venn_palette, guide = "none") +
  ggtitle(paste0("DEGs UP en resistentes\nFDR<",
                 fdr_cutoff, ", |log2FC|>", logfc_cutoff)) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold", color = "#111111"),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

out_venn_up <- file.path(figures_dir, "venn/DEGs_UP_Venn_A2780_UWB_PEO1.png")
ggsave(out_venn_up, p_venn_up, width = 7.5, height = 7.5, dpi = 300)

# Extra: guardar intersecciones explícitas
up_triple <- Reduce(intersect, venn_up_list)
readr::write_tsv(
  tibble(gene_id = up_triple),
  file.path(tables_dir, "venn/DEGs_UP_triple_intersection_A2780_UWB_PEO1.tsv")
)

# ============================================================
# 4) Venn de genes DOWN
# ============================================================

venn_down_list <- list(
  A2780 = A_down,
  UWB   = U_down,
  PEO1  = P_down
)

p_venn_down <- ggVennDiagram(
  venn_down_list,
  set_color   = venn_palette[names(venn_down_list)],
  label_alpha = 0,
  label = "count"
) +
  scale_fill_gradient(low = "white", high = "grey80") +
  scale_color_manual(values = venn_palette, guide = "none") +
  ggtitle(paste0("DEGs DOWN en resistentes (UP en parentales)\nFDR<",
                 fdr_cutoff, ", |log2FC|>", logfc_cutoff)) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold", color = "#111111"),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

out_venn_down <- file.path(figures_dir, "venn/DEGs_DOWN_Venn_A2780_UWB_PEO1.png")
ggsave(out_venn_down, p_venn_down, width = 7.5, height = 7.5, dpi = 300)

down_triple <- Reduce(intersect, venn_down_list)
readr::write_tsv(
  tibble(gene_id = down_triple),
  file.path(tables_dir, "venn/DEGs_DOWN_triple_intersection_A2780_UWB_PEO1.tsv")
)

# ============================================================
# 5) Tablas de membresía (binario) para genes UP y DOWN
# ============================================================
all_up_genes <- unique(c(A_up, U_up, P_up))

up_membership <- tibble(
  gene_id  = all_up_genes,
  A2780_up = gene_id %in% A_up,
  UWB_up   = gene_id %in% U_up,
  PEO1_up  = gene_id %in% P_up
)

readr::write_tsv(
  up_membership,
  file.path(tables_dir, "venn/DEGs_UP_membership_A2780_UWB_PEO1.tsv")
)

all_down_genes <- unique(c(A_down, U_down, P_down))

down_membership <- tibble(
  gene_id    = all_down_genes,
  A2780_down = gene_id %in% A_down,
  UWB_down   = gene_id %in% U_down,
  PEO1_down  = gene_id %in% P_down
)

readr::write_tsv(
  down_membership,
  file.path(tables_dir, "venn/DEGs_DOWN_membership_A2780_UWB_PEO1.tsv")
)

message(">> Venns generados correctamente.")
