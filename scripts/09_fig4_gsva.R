#!/usr/bin/env Rscript

# ============================================================
# 09_fig4_gsva.R
#
# Figure 4: Sample-level and integrated pathway activity in olaparib resistance
#
# Layout (double_col, 200 × 150 mm):
#   [A: per-sample Z(GSVA) heatmap (26 cols) | B: ΔGSVA logFC (3 datasets)]
#
# Question: Is the differential pathway activation from GSEA (Fig 3) consistent
#           at the individual sample level, and confirmed by the independent
#           GSVA + limma statistical framework?
#
# Prerequisites: 04_1_gsva_analysis.R
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

source(file.path({
  f <- commandArgs(trailingOnly = FALSE)
  f <- grep("--file=", f, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("--file=", "", f[1])))
  else path.expand("~/bioinfo/projects/olaparib_resistance/scripts")
}, "00_config.R"))

log_handle <- start_log("09_fig4_gsva")
set.seed(42)

# ============================================================
# 0. Validate inputs
# ============================================================
req <- c(
  gsva_mat  = file.path(tables_dir, "gsva", "GSVA_ALL_merged.tsv"),
  gsva_meta = file.path(tables_dir, "gsva", "GSVA_metadata.tsv"),
  gsva_de   = file.path(tables_dir, "gsva_DE", "GSVA_DE_ALL_datasets.tsv")
)
missing_req <- req[!file.exists(req)]
if (length(missing_req) > 0) {
  stop("Archivos faltantes (re-ejecutar 04_1_gsva_analysis.R):\n",
       paste0("  ", names(missing_req), ": ", missing_req, collapse = "\n"))
}

fig4_dir <- file.path(figures_dir, "fig4")
dir.create(fig4_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Style constants — double_col preset
# ============================================================
FIG_W_MM  <- 200
FIG_H_MM  <- 150
BASE_SIZE <- 8.5

PAL_DATASET <- c("A2780"    = "#1f78b4",
                 "PEO1"     = "#e31a1c",
                 "UWB1.289" = "#33a02c")
PAL_COND    <- c("Parental"  = "#56B4E9",
                 "Resistant" = "#E69F00")
PAL_BRCA    <- c("BRCA-WT"   = "#a6cee3",
                 "BRCA1-mut" = "#fb9a99",
                 "BRCA1-def" = "#b2df8a")

BRCA_MAP     <- c("A2780"       = "BRCA-WT",
                  "PEO1"        = "BRCA1-mut",
                  "UWB_BRCAdef" = "BRCA1-def")
DATASET_DISP <- c("A2780"       = "A2780",
                  "PEO1"        = "PEO1",
                  "UWB_BRCAdef" = "UWB1.289")

DATASET_LEVELS_DATA <- c("A2780", "PEO1", "UWB_BRCAdef")
DATASET_LEVELS_DISP <- c("A2780", "PEO1", "UWB1.289")

# Abbreviate long Hallmark names for compact heatmap labels
abbrev_pathway <- function(x) {
  x %>%
    str_replace_all("Epithelial Mesenchymal Transition", "EMT") %>%
    str_replace_all("Interferon Gamma Response",         "IFN-\u03b3 Response") %>%
    str_replace_all("Interferon Alpha Response",         "IFN-\u03b1 Response") %>%
    str_replace_all("TNFA Signaling Via NFKB",           "TNFA via NFKB") %>%
    str_replace_all("Il6 Jak Stat3 Signaling",           "IL6-JAK-STAT3") %>%
    str_replace_all("Oxidative Phosphorylation",         "Oxidative Phosphoryl.")
}

clean_path <- function(x) abbrev_pathway(clean_hallmark(x))

# ============================================================
# 2. Load data
# ============================================================
gsva_tbl <- read_tsv(req["gsva_mat"], show_col_types = FALSE)
gsva_mat  <- as.matrix(gsva_tbl[, -1, drop = FALSE])
rownames(gsva_mat) <- gsva_tbl$Pathway

meta <- read_tsv(req["gsva_meta"], show_col_types = FALSE) %>%
  filter(dataset %in% DATASET_LEVELS_DATA) %>%
  mutate(
    dataset_disp = DATASET_DISP[dataset],
    brca         = BRCA_MAP[dataset],
    dataset_disp = factor(dataset_disp, levels = DATASET_LEVELS_DISP),
    condition    = factor(condition, levels = c("Parental", "Resistant"))
  )

gsva_de <- read_tsv(req["gsva_de"], show_col_types = FALSE) %>%
  filter(Dataset %in% DATASET_LEVELS_DATA) %>%
  mutate(dataset_disp = factor(DATASET_DISP[Dataset], levels = DATASET_LEVELS_DISP))

# ============================================================
# 3. Select core pathways: GSVA+limma FDR < 0.10, recurrent in >= 2 datasets
# ============================================================
core_paths <- gsva_de %>%
  filter(!is.na(adj.P.Val), adj.P.Val < FDR_CUTOFF_VIZ) %>%
  group_by(Pathway) %>%
  summarise(
    n_datasets     = n_distinct(Dataset),
    mean_abs_logFC = mean(abs(logFC), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_datasets >= 2) %>%
  arrange(desc(mean_abs_logFC)) %>%
  pull(Pathway)

stopifnot("No core pathways found — check 04_1_gsva_analysis.R" = length(core_paths) > 0)

message("Core pathways selected (n=", length(core_paths), "):")
walk(core_paths, ~ message("  ", clean_path(.x)))

# ============================================================
# 4. Panel A: per-sample Z(GSVA) heatmap
# ============================================================
meta_ord <- meta %>% arrange(dataset_disp, condition, sample_id)

mat_raw <- gsva_mat[core_paths, meta_ord$sample_id, drop = FALSE]

# Z-score within each dataset so FPKM/TPM/count scales are comparable
mat_z <- do.call(cbind, lapply(DATASET_LEVELS_DATA, function(ds) {
  sids <- meta_ord$sample_id[meta_ord$dataset == ds]
  t(scale(t(mat_raw[, sids, drop = FALSE])))
}))
mat_z[is.na(mat_z)] <- 0
stopifnot(identical(colnames(mat_z), meta_ord$sample_id))

cap_z   <- 2
mat_z   <- pmax(pmin(mat_z, cap_z), -cap_z)
rownames(mat_z) <- clean_path(core_paths)

# Pre-compute row clustering so Panel B shows the same order
row_clust <- hclust(dist(mat_z), method = "complete")
row_ord   <- row_clust$order

col_z <- colorRamp2(c(-cap_z, 0, cap_z), c("#2166ac", "white", "#b2182b"))

ha_top <- HeatmapAnnotation(
  Dataset   = meta_ord$dataset_disp,
  Condition = meta_ord$condition,
  BRCA      = meta_ord$brca,
  col = list(
    Dataset   = PAL_DATASET,
    Condition = PAL_COND,
    BRCA      = PAL_BRCA
  ),
  annotation_name_gp   = gpar(fontsize = 6.5),
  annotation_height    = unit(rep(2.5, 3), "mm"),
  show_annotation_name = TRUE,
  annotation_name_side = "right",
  annotation_legend_param = list(
    Dataset   = list(labels_gp = gpar(fontsize = 6),
                     title_gp  = gpar(fontsize = 6.5, fontface = "bold")),
    Condition = list(labels_gp = gpar(fontsize = 6),
                     title_gp  = gpar(fontsize = 6.5, fontface = "bold")),
    BRCA      = list(labels_gp = gpar(fontsize = 6),
                     title_gp  = gpar(fontsize = 6.5, fontface = "bold"))
  )
)

ht_A <- Heatmap(
  mat_z,
  name              = "Z(GSVA)",
  col               = col_z,
  na_col            = "grey90",
  width             = unit(75, "mm"),
  cluster_rows      = row_clust,
  cluster_columns   = FALSE,
  column_split      = meta_ord$dataset_disp,
  column_gap        = unit(1.5, "mm"),
  show_column_names = FALSE,
  show_row_names    = FALSE,
  column_title_gp   = gpar(fontsize = 7.5, fontface = "bold"),
  top_annotation    = ha_top,
  heatmap_legend_param = list(
    title         = "Z(GSVA)",
    title_gp      = gpar(fontsize = 7, fontface = "bold"),
    labels_gp     = gpar(fontsize = 6.5),
    legend_height = unit(18, "mm"),
    at            = c(-2, -1, 0, 1, 2)
  )
)

# ============================================================
# 5. Panel B: integrated ΔGSVA (logFC, limma), same row order as Panel A
# ============================================================
mat_b <- gsva_de %>%
  filter(Pathway %in% core_paths) %>%
  select(Pathway, dataset_disp, logFC) %>%
  pivot_wider(names_from = dataset_disp, values_from = logFC) %>%
  column_to_rownames("Pathway") %>%
  as.matrix()
mat_b <- mat_b[core_paths, DATASET_LEVELS_DISP, drop = FALSE]
rownames(mat_b) <- clean_path(core_paths)

sig_mat <- gsva_de %>%
  filter(Pathway %in% core_paths) %>%
  select(Pathway, dataset_disp, adj.P.Val) %>%
  pivot_wider(names_from = dataset_disp, values_from = adj.P.Val) %>%
  column_to_rownames("Pathway") %>%
  as.matrix()
sig_mat <- sig_mat[core_paths, DATASET_LEVELS_DISP, drop = FALSE]
rownames(sig_mat) <- clean_path(core_paths)

cap_b  <- 1
mat_b2 <- pmax(pmin(mat_b, cap_b), -cap_b)
col_b  <- colorRamp2(c(-cap_b, 0, cap_b), c("#2166ac", "white", "#b2182b"))

ht_B <- Heatmap(
  mat_b2,
  name              = "\u0394GSVA",
  col               = col_b,
  na_col            = "grey90",
  width             = unit(26, "mm"),
  cluster_rows      = FALSE,
  row_order         = row_ord,
  cluster_columns   = FALSE,
  show_column_names = TRUE,
  show_row_names    = TRUE,
  row_names_side    = "right",
  row_names_gp      = gpar(fontsize = 7),
  column_names_gp   = gpar(fontsize = 7.5, fontface = "bold"),
  column_names_rot  = 45,
  cell_fun = function(j, i, x, y, width, height, fill) {
    pv  <- sig_mat[i, j]
    lbl <- if (!is.na(pv) && pv < 0.05)  "*"
           else if (!is.na(pv) && pv < 0.10) "\u00b7"
           else ""
    if (nzchar(lbl)) {
      grid.text(lbl, x, y, gp = gpar(fontsize = 9, col = "black"))
    }
  },
  heatmap_legend_param = list(
    title         = "\u0394GSVA",
    title_gp      = gpar(fontsize = 7, fontface = "bold"),
    labels_gp     = gpar(fontsize = 6.5),
    legend_height = unit(18, "mm"),
    at            = c(-1, -0.5, 0, 0.5, 1)
  )
)

# ============================================================
# 6. Compose and export
# ============================================================
out_pdf <- file.path(fig4_dir, "Fig4.pdf")
out_png <- file.path(fig4_dir, "Fig4.png")

FIG_W_IN <- FIG_W_MM / 25.4
FIG_H_IN <- FIG_H_MM / 25.4

draw_fig4 <- function() {
  grid.newpage()
  lt <- grid.layout(
    nrow   = 1,
    ncol   = 2,
    widths = unit(c(FIG_W_MM * 0.65, FIG_W_MM * 0.35), "mm")
  )
  pushViewport(viewport(layout = lt))

  # Panel A — per-sample heatmap
  pushViewport(viewport(layout.pos.col = 1))
  draw(ht_A, newpage = FALSE,
       heatmap_legend_side    = "right",
       annotation_legend_side = "right",
       padding = unit(c(4, 3, 4, 5), "mm"))
  grid.text("A", x = 0.01, y = 0.99, just = c("left", "top"),
            gp = gpar(fontsize = 9, fontface = "bold"))
  popViewport()

  # Panel B — integrated logFC heatmap
  pushViewport(viewport(layout.pos.col = 2))
  draw(ht_B, newpage = FALSE,
       heatmap_legend_side    = "right",
       annotation_legend_side = "right",
       padding = unit(c(4, 5, 4, 3), "mm"))
  grid.text("B", x = 0.01, y = 0.99, just = c("left", "top"),
            gp = gpar(fontsize = 9, fontface = "bold"))
  popViewport()

  popViewport()
}

pdf(out_pdf, width = FIG_W_IN, height = FIG_H_IN)
draw_fig4()
invisible(dev.off())

png(out_png, width = FIG_W_MM, height = FIG_H_MM, units = "mm", res = 600)
draw_fig4()
invisible(dev.off())

message("Figure 4 saved:")
message("  PDF: ", out_pdf)
message("  PNG: ", out_png)

stop_log(log_handle)
