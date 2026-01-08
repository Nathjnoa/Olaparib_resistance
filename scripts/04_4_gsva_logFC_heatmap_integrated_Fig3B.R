#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

proj_dir <- Sys.getenv("OLAPARIB_RESISTANCE_DIR", unset = "")
if (!nzchar(proj_dir)) {
  proj_dir <- "~/bioinfo/projects/olaparib_resistance"
}
proj_dir <- path.expand(proj_dir)
if (!dir.exists(proj_dir)) {
  stop("No existe el directorio del proyecto: ", proj_dir)
}

in_gsvaDE <- file.path(proj_dir, "results", "tables", "gsva_DE", "GSVA_DE_ALL_datasets.tsv")
out_dir   <- file.path(proj_dir, "results", "figures", "gsva_heatmaps")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_png <- file.path(out_dir, "Fig3B_GSVA_logFC_heatmap_integrated.png")
out_pdf <- sub("\\.png$", ".pdf", out_png)
out_tsv <- file.path(proj_dir, "results", "tables", "gsva_DE", "GSVA_core_pathways_summary.tsv")

if (!file.exists(in_gsvaDE)) {
  stop("No encuentro archivo de entrada: ", in_gsvaDE)
}

# -----------------------------
# Parámetros
# -----------------------------
fdr_cut <- 0.10
min_datasets <- 2

dataset_levels <- c("A2780", "PEO1", "UWB_BRCAdef", "UWB_BRCAprof")

# -----------------------------
# Cargar tabla integrada GSVA-DE
# -----------------------------
tbl <- readr::read_tsv(in_gsvaDE, show_col_types = FALSE) %>%
  mutate(Dataset = factor(Dataset, levels = dataset_levels))

stopifnot(all(c("Pathway","Dataset","logFC","adj.P.Val") %in% colnames(tbl)))

# -----------------------------
# Seleccionar pathways recurrentes
# -----------------------------
core <- tbl %>%
  filter(!is.na(adj.P.Val), adj.P.Val < fdr_cut) %>%
  group_by(Pathway) %>%
  summarise(
    n_datasets = n_distinct(Dataset),
    mean_abs_logFC = mean(abs(logFC), na.rm = TRUE),
    # “concordancia” simple: si en los significativos todos tienen mismo signo
    sign_consistent = {
      s <- sign(logFC[!is.na(adj.P.Val) & adj.P.Val < fdr_cut])
      if (length(s) == 0) NA else (length(unique(s)) == 1)
    },
    .groups = "drop"
  ) %>%
  filter(n_datasets >= min_datasets) %>%
  arrange(desc(mean_abs_logFC))

if (nrow(core) == 0) {
  stop("No hay pathways recurrentes con FDR<", fdr_cut, " en >= ", min_datasets, " datasets.")
}

core_paths <- core$Pathway

# matriz logFC
mat <- tbl %>%
  filter(Pathway %in% core_paths) %>%
  select(Pathway, Dataset, logFC) %>%
  pivot_wider(names_from = Dataset, values_from = logFC) %>%
  column_to_rownames("Pathway") %>%
  as.matrix()

# ordenar filas por señal promedio
row_order <- order(rowMeans(abs(mat), na.rm = TRUE), decreasing = TRUE)
mat <- mat[row_order, , drop = FALSE]
mat <- mat[, dataset_levels, drop = FALSE]

# guardar resumen
core_out <- core %>%
  mutate(Pathway = factor(Pathway, levels = rownames(mat))) %>%
  arrange(Pathway)
readr::write_tsv(core_out, out_tsv)

# -----------------------------
# Color: rango fijo y simétrico
# -----------------------------
cap <- 1
mat2 <- pmax(pmin(mat, cap), -cap)
col_fun <- circlize::colorRamp2(c(-cap, 0, cap), c("#2166ac", "white", "#b2182b"))

# -----------------------------
# Heatmap
# -----------------------------
ht <- Heatmap(
  mat2,
  name = "logFC",
  col = col_fun,
  na_col = "grey90",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = grid::gpar(fontsize = 10),
  column_names_gp = grid::gpar(fontsize = 11),
  column_names_rot = 45,
  column_title = paste0("ΔGSVA (Resistant − Parental) — Hallmarks recurrentes (FDR<", fdr_cut, ", ≥", min_datasets, " datasets)"),
  column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(title = "logFC (Δscore)")
)

png(out_png, width = 2800, height = 2200, res = 300)
draw(ht, heatmap_legend_side = "right")
dev.off()

pdf(out_pdf, width = 10, height = 8, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right")
dev.off()

message("✔ Fig3B guardada: ", out_png)
message("✔ PDF guardado: ", out_pdf)
message("✔ Resumen core TSV: ", out_tsv)
