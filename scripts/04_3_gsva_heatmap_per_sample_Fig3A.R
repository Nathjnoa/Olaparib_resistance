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

gsva_dir   <- file.path(proj_dir, "results", "tables", "gsva")
gsvaDE_dir <- file.path(proj_dir, "results", "tables", "gsva_DE")
out_dir    <- file.path(proj_dir, "results", "figures", "gsva_heatmaps")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

in_gsva   <- file.path(gsva_dir, "GSVA_ALL_merged.tsv")
in_meta   <- file.path(gsva_dir, "GSVA_metadata.tsv")
in_gsvaDE <- file.path(gsvaDE_dir, "GSVA_DE_ALL_datasets.tsv")

out_png <- file.path(out_dir, "Fig3A_GSVA_scores_per_sample_heatmap.png")
out_pdf <- sub("\\.png$", ".pdf", out_png)

missing_files <- c(in_gsva, in_meta, in_gsvaDE)
missing_files <- missing_files[!file.exists(missing_files)]
if (length(missing_files) > 0) {
  stop("No encuentro archivos de entrada:\n", paste(missing_files, collapse = "\n"))
}

# -----------------------------
# Parámetros
# -----------------------------
fdr_cut <- 0.10        # para seleccionar pathways “recurrentes”
min_datasets <- 2      # recurrentes en >=2 datasets
top_max <- 20          # límite (si hay más recurrentes, toma top por |logFC| promedio)

dataset_levels <- c("A2780", "PEO1", "UWB_BRCAdef", "UWB_BRCAprof")

# -----------------------------
# Cargar datos
# -----------------------------
gsva_tbl <- readr::read_tsv(in_gsva, show_col_types = FALSE)
stopifnot("Pathway" %in% colnames(gsva_tbl))

gsva_mat <- as.matrix(gsva_tbl[, -1, drop = FALSE])
rownames(gsva_mat) <- gsva_tbl$Pathway

meta <- readr::read_tsv(in_meta, show_col_types = FALSE) %>%
  mutate(
    dataset = factor(dataset, levels = dataset_levels),
    condition = factor(condition, levels = c("Parental", "Resistant"))
  )

stopifnot(all(meta$sample_id %in% colnames(gsva_mat)))
stopifnot(all(colnames(gsva_mat) %in% meta$sample_id))

gsvaDE <- readr::read_tsv(in_gsvaDE, show_col_types = FALSE) %>%
  mutate(
    Dataset = factor(Dataset, levels = dataset_levels)
  )

# -----------------------------
# Selección de pathways “recurrentes” por GSVA-DE
# -----------------------------
core_tbl <- gsvaDE %>%
  filter(!is.na(adj.P.Val), adj.P.Val < fdr_cut) %>%
  group_by(Pathway) %>%
  summarise(
    n_datasets = n_distinct(Dataset),
    mean_abs_logFC = mean(abs(logFC), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_datasets >= min_datasets) %>%
  arrange(desc(mean_abs_logFC))

core_paths <- core_tbl$Pathway
if (length(core_paths) == 0) {
  stop("No hay pathways recurrentes con FDR<", fdr_cut, " en >= ", min_datasets, " datasets.")
}

# limitar a top_max si hay demasiados
if (length(core_paths) > top_max) core_paths <- core_paths[1:top_max]

# -----------------------------
# Subset + orden de columnas
# -----------------------------
mat <- gsva_mat[intersect(core_paths, rownames(gsva_mat)), , drop = FALSE]
if (nrow(mat) < 2) stop("Muy pocos pathways tras el subset para heatmap.")

meta_ord <- meta %>%
  arrange(dataset, condition, sample_id)

mat <- mat[, meta_ord$sample_id, drop = FALSE]

message("Conteo de muestras por dataset/condición:")
print(table(meta_ord$dataset, meta_ord$condition))
message("Primeras 10 columnas tras el ordenamiento: ",
        paste(head(colnames(mat), 10), collapse = ", "))

# -----------------------------
# Escalado (Z-score por pathway)
# -----------------------------
mat_z <- t(scale(t(mat)))
mat_z[is.na(mat_z)] <- 0

# Rango fijo para color
cap <- 2
mat_z <- pmax(pmin(mat_z, cap), -cap)

col_fun <- circlize::colorRamp2(c(-cap, 0, cap), c("#2166ac", "white", "#b2182b"))

# -----------------------------
# Anotaciones de columnas
# -----------------------------
ha <- HeatmapAnnotation(
  Dataset = meta_ord$dataset,
  Condition = meta_ord$condition,
  annotation_name_gp = grid::gpar(fontsize = 10)
)

# -----------------------------
# Heatmap
# -----------------------------
ht <- Heatmap(
  mat_z,
  name = "Z(GSVA)",
  col = col_fun,
  na_col = "grey90",
  cluster_rows = TRUE,
  cluster_columns = FALSE, # columnas ya ordenadas por dataset/condición
  column_split = meta_ord$dataset,
  column_gap = unit(2, "mm"),
  show_column_names = TRUE,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = grid::gpar(fontsize = 9),
  column_names_gp = grid::gpar(fontsize = 9),
  column_names_rot = 45,
  top_annotation = ha,
  column_title = paste0("GSVA Hallmarks (por muestra) — recurrentes (FDR<", fdr_cut, ", ≥", min_datasets, " datasets)"),
  column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(title = "Z(GSVA)")
)

png(out_png, width = 3200, height = 2200, res = 300)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

pdf(out_pdf, width = 12, height = 8, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

message("✔ Fig3A guardada: ", out_png)
message("✔ PDF guardado: ", out_pdf)
