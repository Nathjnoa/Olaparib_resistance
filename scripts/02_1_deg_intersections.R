#!/usr/bin/env Rscript

# ============================================================
# 02_1_deg_intersections.R
#
# DEG intersections across 3 datasets (A2780, UWB BRCA-def, PEO1):
#   Section 1: Venn diagrams (ggVennDiagram)
#   Section 2: UpSet plots (UpSetR)
#
# Input:  DE tables from results/tables/
# Output: Venn PNGs, UpSet PNGs, membership TSVs
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggVennDiagram)
  library(UpSetR)
})

source(file.path({f <- commandArgs(trailingOnly=FALSE); f <- grep("--file=",f,value=TRUE); if(length(f)) dirname(normalizePath(sub("--file=","",f[1]))) else path.expand("~/bioinfo/projects/olaparib_resistance/scripts")}, "00_config.R"))

log_handle <- start_log("02_1_deg_intersections")

dir.create(file.path(figures_dir, "venn"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(figures_dir, "upset"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(tables_dir,  "venn"),  showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1) Load DE tables (single shared loading block)
# ============================================================
file_A   <- file.path(tables_dir, "de", "GSE153867_A2780_limma_DE_FPKM.tsv")
file_def <- file.path(tables_dir, "de", "GSE235980_BRCAdef_Res_vs_Par_DESeq2.tsv")
file_P   <- file.path(tables_dir, "de", "GSE117765_PEO1_limma_TPM_Res_vs_Par.tsv")

missing <- c(file_A, file_def, file_P)[!file.exists(c(file_A, file_def, file_P))]
if (length(missing) > 0) stop("Archivos de entrada no encontrados:\n", paste(missing, collapse = "\n"))

# limma topTable uses "ID"; DESeq2 uses "gene_id" — normalize to gene_id on load
A    <- readr::read_tsv(file_A,   show_col_types = FALSE) %>% dplyr::rename(gene_id = ID)
Udef <- readr::read_tsv(file_def, show_col_types = FALSE)
P    <- readr::read_tsv(file_P,   show_col_types = FALSE) %>% dplyr::rename(gene_id = ID)

# ============================================================
# 2) Build DEG lists (UP / DOWN per dataset)
# ============================================================

# A2780 (limma: adj.P.Val, log2FC, ID column)
A_up <- A %>%
  filter(!is.na(adj.P.Val), !is.na(log2FC),
         adj.P.Val < FDR_CUTOFF, log2FC >  LOGFC_CUTOFF) %>%
  pull(gene_id) %>% unique()

A_down <- A %>%
  filter(!is.na(adj.P.Val), !is.na(log2FC),
         adj.P.Val < FDR_CUTOFF, log2FC < -LOGFC_CUTOFF) %>%
  pull(gene_id) %>% unique()

# UWB BRCA-def (DESeq2: padj, log2FoldChange)
U_up <- Udef %>%
  filter(!is.na(padj), !is.na(log2FoldChange),
         padj < FDR_CUTOFF, log2FoldChange >  LOGFC_CUTOFF) %>%
  pull(gene_id) %>% unique()

U_down <- Udef %>%
  filter(!is.na(padj), !is.na(log2FoldChange),
         padj < FDR_CUTOFF, log2FoldChange < -LOGFC_CUTOFF) %>%
  pull(gene_id) %>% unique()

# PEO1 (limma: adj.P.Val, log2FC)
P_up <- P %>%
  filter(!is.na(adj.P.Val), !is.na(log2FC),
         adj.P.Val < FDR_CUTOFF, log2FC >  LOGFC_CUTOFF) %>%
  pull(gene_id) %>% unique()

P_down <- P %>%
  filter(!is.na(adj.P.Val), !is.na(log2FC),
         adj.P.Val < FDR_CUTOFF, log2FC < -LOGFC_CUTOFF) %>%
  pull(gene_id) %>% unique()

# Counts summary
summary_counts <- tibble(
  dataset = c("A2780_up", "A2780_down", "UWB_BRCAdef_up", "UWB_BRCAdef_down",
              "PEO1_up", "PEO1_down"),
  n_genes = c(length(A_up), length(A_down), length(U_up), length(U_down),
              length(P_up), length(P_down))
)
message("DEG counts:")
print(summary_counts)
readr::write_tsv(summary_counts, file.path(tables_dir, "venn/DEG_counts_for_venn.tsv"))

# ============================================================
# SECTION 1: Venn diagrams
# ============================================================

# UP genes
venn_up_list <- list(A2780 = A_up, UWB = U_up, PEO1 = P_up)

p_venn_up <- ggVennDiagram(
  venn_up_list,
  set_color   = venn_palette[names(venn_up_list)],
  label_alpha = 0,
  label = "count"
) +
  scale_fill_gradient(low = "white", high = "grey80") +
  scale_color_manual(values = venn_palette, guide = "none") +
  ggtitle(paste0("DEGs UP en resistentes\nFDR<", FDR_CUTOFF, ", |log2FC|>", LOGFC_CUTOFF)) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold", color = "#111111"),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

ggsave(file.path(figures_dir, "venn/DEGs_UP_Venn_A2780_UWB_PEO1.png"),
       p_venn_up, width = 7.5, height = 7.5, dpi = 300)

# Triple intersection UP
up_triple <- Reduce(intersect, venn_up_list)
readr::write_tsv(tibble(gene_id = up_triple),
                 file.path(tables_dir, "venn/DEGs_UP_triple_intersection_A2780_UWB_PEO1.tsv"))

# DOWN genes
venn_down_list <- list(A2780 = A_down, UWB = U_down, PEO1 = P_down)

p_venn_down <- ggVennDiagram(
  venn_down_list,
  set_color   = venn_palette[names(venn_down_list)],
  label_alpha = 0,
  label = "count"
) +
  scale_fill_gradient(low = "white", high = "grey80") +
  scale_color_manual(values = venn_palette, guide = "none") +
  ggtitle(paste0("DEGs DOWN en resistentes\nFDR<", FDR_CUTOFF, ", |log2FC|>", LOGFC_CUTOFF)) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold", color = "#111111"),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

ggsave(file.path(figures_dir, "venn/DEGs_DOWN_Venn_A2780_UWB_PEO1.png"),
       p_venn_down, width = 7.5, height = 7.5, dpi = 300)

down_triple <- Reduce(intersect, venn_down_list)
readr::write_tsv(tibble(gene_id = down_triple),
                 file.path(tables_dir, "venn/DEGs_DOWN_triple_intersection_A2780_UWB_PEO1.tsv"))

# Membership tables
all_up_genes <- unique(c(A_up, U_up, P_up))
up_membership <- tibble(
  gene_id       = all_up_genes,
  A2780_up      = gene_id %in% A_up,
  UWB_BRCAdef_up = gene_id %in% U_up,
  PEO1_up       = gene_id %in% P_up
)
readr::write_tsv(up_membership,
                 file.path(tables_dir, "venn/DEGs_UP_membership_A2780_UWB_PEO1.tsv"))

all_down_genes <- unique(c(A_down, U_down, P_down))
down_membership <- tibble(
  gene_id         = all_down_genes,
  A2780_down      = gene_id %in% A_down,
  UWB_BRCAdef_down = gene_id %in% U_down,
  PEO1_down       = gene_id %in% P_down
)
readr::write_tsv(down_membership,
                 file.path(tables_dir, "venn/DEGs_DOWN_membership_A2780_UWB_PEO1.tsv"))

message("✔ Venn diagrams generados")

# ============================================================
# SECTION 2: UpSet plots
# ============================================================

up_list   <- list(A2780 = A_up,   UWB = U_up,   PEO1 = P_up)
down_list <- list(A2780 = A_down, UWB = U_down, PEO1 = P_down)

up_mat   <- UpSetR::fromList(up_list)
down_mat <- UpSetR::fromList(down_list)

upset_params <- list(
  nsets = 3,
  nintersects = NA,
  sets = c("A2780", "UWB", "PEO1"),
  order.by = "freq",
  main.bar.color = "grey30",
  sets.bar.color = "grey40",
  mainbar.y.label = "Genes en intersección",
  sets.x.label = "Genes por dataset",
  text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.2)
)

png(file.path(figures_dir, "upset/DEGs_UP_UpSet_A2780_UWB_PEO1.png"),
    width = 2000, height = 1500, res = 300)
do.call(upset, c(list(data = up_mat), upset_params))
dev.off()

png(file.path(figures_dir, "upset/DEGs_DOWN_UpSet_A2780_UWB_PEO1.png"),
    width = 2000, height = 1500, res = 300)
do.call(upset, c(list(data = down_mat), upset_params))
dev.off()

message("✔ UpSet plots generados")

stop_log(log_handle)
