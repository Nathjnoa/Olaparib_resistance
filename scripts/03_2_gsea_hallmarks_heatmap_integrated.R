#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)   # <- clave para usar grid::gpar explícito
})

proj_dir <- Sys.getenv("OLAPARIB_RESISTANCE_DIR", unset = "")
if (!nzchar(proj_dir)) {
  proj_dir <- "~/bioinfo/projects/olaparib_resistance"
}
proj_dir <- path.expand(proj_dir)
if (!dir.exists(proj_dir)) {
  stop("No existe el directorio del proyecto: ", proj_dir)
}
in_dir   <- file.path(proj_dir, "results", "tables", "gsea_hallmarks")
out_png  <- file.path(proj_dir, "results", "figures", "gsea_hallmarks", "Hallmarks_NES_heatmap.png")
out_pdf  <- sub("\\.png$", ".pdf", out_png)

# Parámetros
fdr_cut      <- 0.10
min_datasets <- 2
cap_nes      <- 2

dataset_levels <- c("A2780", "PEO1", "UWB_BRCAdef", "UWB_BRCAprof")

files <- tibble(
  dataset = dataset_levels,
  path = file.path(in_dir, c(
    "A2780_GSEA_Hallmarks.tsv",
    "PEO1_GSEA_Hallmarks.tsv",
    "UWB_BRCAdef_GSEA_Hallmarks.tsv",
    "UWB_BRCAprof_GSEA_Hallmarks.tsv"
  ))
)

missing_files <- files$path[!file.exists(files$path)]
if (length(missing_files) > 0) {
  stop("No encuentro archivos de GSEA Hallmarks:\n", paste(missing_files, collapse = "\n"))
}

clean_hallmark <- function(x) {
  y <- x %>%
    str_replace_all("^HALLMARK_", "") %>%
    str_replace_all("_", " ") %>%
    str_to_lower() %>%
    str_replace_all("\\s+", " ") %>%
    str_trim() %>%
    str_to_title()

  y <- y %>%
    str_replace_all("\\bE2f\\b", "E2F") %>%
    str_replace_all("\\bNfkb\\b", "NFKB") %>%
    str_replace_all("\\bTnfa\\b", "TNFA") %>%
    str_replace_all("\\bIl6 Jak Stat3\\b", "IL6 JAK STAT3") %>%
    str_replace_all("\\bMtorc1\\b", "MTORC1") %>%
    str_replace_all("\\bUv\\b", "UV") %>%
    str_replace_all("\\bG2m\\b", "G2M") %>%
    str_replace_all("\\bMyc\\b", "MYC")
  y
}

tbl <- files %>%
  mutate(df = map(path, ~ readr::read_tsv(.x, show_col_types = FALSE))) %>%
  unnest(df) %>%
  transmute(
    dataset = factor(dataset, levels = dataset_levels),
    hallmark = Description,
    NES = as.numeric(NES),
    padj = as.numeric(p.adjust)
  ) %>%
  mutate(hallmark_clean = clean_hallmark(hallmark)) %>%
  filter(!is.na(hallmark_clean), hallmark_clean != "")

recurrent <- tbl %>%
  filter(!is.na(padj), padj < fdr_cut) %>%
  count(hallmark_clean) %>%
  filter(n >= min_datasets) %>%
  pull(hallmark_clean)

tbl2 <- tbl %>% filter(hallmark_clean %in% recurrent)

mat <- tbl2 %>%
  select(dataset, hallmark_clean, NES) %>%
  pivot_wider(names_from = dataset, values_from = NES) %>%
  column_to_rownames("hallmark_clean") %>%
  as.matrix()

# Cap NES para rango fijo comparable
mat <- pmax(pmin(mat, cap_nes), -cap_nes)

# Ordenar filas por señal promedio
row_order <- order(rowMeans(abs(mat), na.rm = TRUE), decreasing = TRUE)
mat <- mat[row_order, , drop = FALSE]

col_fun <- circlize::colorRamp2(c(-cap_nes, 0, cap_nes), c("#2166ac", "white", "#b2182b"))

ht <- Heatmap(
  mat,
  name = "NES",
  col = col_fun,
  na_col = "grey90",
  cluster_rows = TRUE,
  cluster_columns = TRUE,

  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = grid::gpar(fontsize = 10),

  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp = grid::gpar(fontsize = 11),

  column_title = paste0("Hallmarks GSEA (NES), recurrentes en ≥", min_datasets, " datasets (FDR<", fdr_cut, ")"),
  column_title_gp = grid::gpar(fontsize = 16, fontface = "bold"),

  heatmap_legend_param = list(title = "NES")
)

dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)

png(out_png, width = 2800, height = 2100, res = 300)
draw(ht, heatmap_legend_side = "right")
dev.off()

pdf(out_pdf, width = 10, height = 8, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right")
dev.off()

message("✔ Heatmap integrador guardado en: ", out_png)
message("✔ PDF guardado en: ", out_pdf)
