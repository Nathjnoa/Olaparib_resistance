#!/usr/bin/env Rscript

# ============================================================
# 07_fig2_overlap.R
#
# Figure 2: Shared transcriptomic signature across olaparib-resistant cell lines
#
# Layout (double_col, 200 × 230 mm):
#   Row 1 (42%): [A: UpSet — upregulated DEGs | B: UpSet — downregulated DEGs]
#   Row 2 (58%): [C: ComplexHeatmap — triple-intersection genes, full width]
#
# Prerequisites: 01_differential_expression.R + 02_1_deg_intersections.R
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(edgeR)
  library(grid)
})

source(file.path({
  f <- commandArgs(trailingOnly = FALSE)
  f <- grep("--file=", f, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("--file=", "", f[1])))
  else path.expand("~/bioinfo/projects/olaparib_resistance/scripts")
}, "00_config.R"))

log_handle <- start_log("07_fig2_overlap")
set.seed(42)

# ============================================================
# 0. Validate required inputs
# ============================================================
req <- c(
  memb_up    = file.path(tables_dir, "venn", "DEGs_UP_membership_A2780_UWB_PEO1.tsv"),
  memb_down  = file.path(tables_dir, "venn", "DEGs_DOWN_membership_A2780_UWB_PEO1.tsv"),
  trip_up    = file.path(tables_dir, "venn", "DEGs_UP_triple_intersection_A2780_UWB_PEO1.tsv"),
  trip_down  = file.path(tables_dir, "venn", "DEGs_DOWN_triple_intersection_A2780_UWB_PEO1.tsv"),
  expr_a2780 = file.path(processed_dir, "GSE153867_A2780_FPKM.tsv"),
  expr_peo1  = file.path(raw_dir, "GSE117765", "GSE117765_matrix.txt"),
  expr_uwb   = file.path(raw_dir, "GSE235980", "GSE235980_CountReads.txt")
)
missing_req <- req[!file.exists(req)]
if (length(missing_req) > 0) {
  stop("Missing files (re-run prerequisites):\n",
       paste0("  ", names(missing_req), ": ", missing_req, collapse = "\n"))
}

fig2_dir <- file.path(figures_dir, "fig2")
dir.create(fig2_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Style constants — double_col preset, consistent with Fig 1
# ============================================================
FIG_W_MM  <- 200
FIG_H_MM  <- 260

PAL_COND    <- c("Parental"   = "#56B4E9", "Resistant"  = "#E69F00")
PAL_DATASET <- c("A2780"      = "#1f78b4", "PEO1"       = "#e31a1c", "UWB" = "#33a02c")
PAL_DIR     <- c("UP"         = "#D81B60", "DOWN"       = "#1E88E5")
PAL_BRCA    <- c("BRCA-WT"    = "#a6cee3", "BRCA1-mut"  = "#fb9a99", "BRCA1-def" = "#b2df8a")

# ============================================================
# 2. Load membership tables → gene lists for UpSet
# ============================================================
memb_up   <- read_tsv(req["memb_up"],   show_col_types = FALSE)
memb_down <- read_tsv(req["memb_down"], show_col_types = FALSE)
trip_up   <- read_tsv(req["trip_up"],   show_col_types = FALSE)$gene_id
trip_down <- read_tsv(req["trip_down"], show_col_types = FALSE)$gene_id

message("Triple intersection: ", length(trip_up), " UP / ", length(trip_down), " DOWN genes")

# Detect boolean columns dynamically (order: A2780, UWB, PEO1)
bool_cols_up   <- setdiff(colnames(memb_up),   "gene_id")
bool_cols_down <- setdiff(colnames(memb_down), "gene_id")

lists_up <- setNames(
  lapply(bool_cols_up,   function(col) memb_up$gene_id[memb_up[[col]]]),
  c("A2780", "UWB", "PEO1")
)
lists_down <- setNames(
  lapply(bool_cols_down, function(col) memb_down$gene_id[memb_down[[col]]]),
  c("A2780", "UWB", "PEO1")
)

# ============================================================
# 3. UpSet plots — Panels A and B
# ============================================================
make_upset <- function(gene_lists, title) {
  m <- make_comb_mat(gene_lists, mode = "distinct")
  UpSet(
    m,
    set_order  = c("A2780", "UWB", "PEO1"),
    comb_order = order(comb_size(m), decreasing = TRUE),
    top_annotation = upset_top_annotation(
      m,
      add_numbers        = TRUE,
      numbers_gp         = gpar(fontsize = 5.5),
      annotation_name_gp = gpar(fontsize = 7),
      height             = unit(2.5, "cm")
    ),
    right_annotation = upset_right_annotation(
      m,
      add_numbers        = TRUE,
      annotation_name_gp = gpar(fontsize = 7),
      width              = unit(2, "cm")
    ),
    row_names_gp    = gpar(fontsize = 7.5, fontface = "bold",
                           col = PAL_DATASET[c("A2780", "UWB", "PEO1")]),
    pt_size         = unit(3.5, "mm"),
    lwd             = 2,
    column_title    = title,
    column_title_gp = gpar(fontsize = 8, fontface = "bold")
  )
}

ht_upset_up   <- make_upset(lists_up,   "Upregulated in resistance")
ht_upset_down <- make_upset(lists_down, "Downregulated in resistance")

# ============================================================
# 4. Load expression matrices for Panel C
# ============================================================

# --- A2780: log2(FPKM+1), parental = O-1..8, resistant = C-1..8 ---
a2780_raw  <- read_tsv(req["expr_a2780"], show_col_types = FALSE)
par_a2780  <- paste0("O-", 1:8)
res_a2780  <- paste0("C-", 1:8)

expr_a2780 <- a2780_raw %>%
  filter(gene_id %in% c(trip_up, trip_down)) %>%
  column_to_rownames("gene_id") %>%
  select(all_of(c(par_a2780, res_a2780))) %>%
  mutate(across(everything(), ~ log2(. + 1))) %>%
  as.matrix()

# --- PEO1: log2(TPM+1), parental = Adherent, resistant = Clone ---
peo1_raw  <- read_tsv(req["expr_peo1"], show_col_types = FALSE) %>%
  rename(gene_id = 1)
par_peo1  <- grep("Adherent", colnames(peo1_raw), value = TRUE)
res_peo1  <- grep("Clone",    colnames(peo1_raw), value = TRUE)

expr_peo1 <- peo1_raw %>%
  filter(gene_id %in% c(trip_up, trip_down)) %>%
  column_to_rownames("gene_id") %>%
  select(all_of(c(par_peo1, res_peo1))) %>%
  mutate(across(everything(), ~ log2(. + 1))) %>%
  as.matrix()

# --- UWB: log2(CPM+1), BRCAdef only (AN14028316/18 = resistant, AN14028320/22 = parental) ---
# Load full count matrix for proper library-size normalization, then filter to BRCAdef
uwb_full <- read_tsv(req["expr_uwb"], show_col_types = FALSE) %>%
  rename(gene_id = 1) %>%
  filter(!startsWith(gene_id, "__"))  # remove HTSeq summary lines

uwb_keep <- c(UWB_DEF_PAR, UWB_DEF_RES)
uwb_counts_full <- uwb_full %>%
  column_to_rownames("gene_id") %>%
  select(any_of(uwb_keep)) %>%
  as.matrix()

uwb_lcpm_full <- cpm(uwb_counts_full, log = TRUE, prior.count = 1)

par_uwb  <- UWB_DEF_PAR[UWB_DEF_PAR %in% colnames(uwb_lcpm_full)]
res_uwb  <- UWB_DEF_RES[UWB_DEF_RES %in% colnames(uwb_lcpm_full)]
uwb_genes <- intersect(c(trip_up, trip_down), rownames(uwb_lcpm_full))
expr_uwb  <- uwb_lcpm_full[uwb_genes, c(par_uwb, res_uwb), drop = FALSE]

# ============================================================
# 5. Within-dataset Z-score and combined matrix
# ============================================================
z_score_rows <- function(mat) {
  row_sd <- apply(mat, 1, sd, na.rm = TRUE)
  row_sd <- pmax(row_sd, 1e-6)  # guard against zero-variance genes
  sweep(sweep(mat, 1, rowMeans(mat, na.rm = TRUE), "-"), 1, row_sd, "/")
}

zmat_a2780 <- z_score_rows(expr_a2780)
zmat_peo1  <- z_score_rows(expr_peo1)
zmat_uwb   <- z_score_rows(expr_uwb)

# Genes present in all 3 datasets
genes_avail <- Reduce(intersect,
                      list(rownames(zmat_a2780), rownames(zmat_peo1), rownames(zmat_uwb)))
missing_g   <- setdiff(c(trip_up, trip_down), genes_avail)
if (length(missing_g) > 0) {
  message("Genes absent from >=1 dataset (excluded): ", paste(missing_g, collapse = ", "))
}

# Ordered: UP genes first, DOWN genes second
gene_order_up   <- intersect(trip_up,   genes_avail)
gene_order_down <- intersect(trip_down, genes_avail)
gene_order      <- c(gene_order_up, gene_order_down)

mat_combined <- cbind(
  zmat_a2780[gene_order, , drop = FALSE],
  zmat_peo1[gene_order,  , drop = FALSE],
  zmat_uwb[gene_order,   , drop = FALSE]
)

message("Heatmap matrix: ", nrow(mat_combined), " genes × ",
        ncol(mat_combined), " samples")

# ============================================================
# 6. ComplexHeatmap — Panel C
# ============================================================

# Column annotation
anno_df <- data.frame(
  Dataset   = factor(c(rep("A2780", ncol(expr_a2780)),
                       rep("PEO1",  ncol(expr_peo1)),
                       rep("UWB",   ncol(expr_uwb))),
                     levels = c("A2780", "PEO1", "UWB")),
  Condition = factor(c(rep("Parental",  length(par_a2780)),
                       rep("Resistant", length(res_a2780)),
                       rep("Parental",  length(par_peo1)),
                       rep("Resistant", length(res_peo1)),
                       rep("Parental",  length(par_uwb)),
                       rep("Resistant", length(res_uwb))),
                     levels = c("Parental", "Resistant")),
  BRCA      = factor(c(rep("BRCA-WT",   ncol(expr_a2780)),
                       rep("BRCA1-mut", ncol(expr_peo1)),
                       rep("BRCA1-def", ncol(expr_uwb))),
                     levels = c("BRCA-WT", "BRCA1-mut", "BRCA1-def")),
  row.names = colnames(mat_combined),
  stringsAsFactors = FALSE
)

# Row annotation
row_df <- data.frame(
  Direction = factor(c(rep("UP",   length(gene_order_up)),
                       rep("DOWN", length(gene_order_down))),
                     levels = c("UP", "DOWN")),
  row.names = gene_order
)

legend_param <- function(title, ...) {
  list(title = title,
       title_gp  = gpar(fontsize = 7, fontface = "bold"),
       labels_gp = gpar(fontsize = 6),
       grid_height = unit(3, "mm"),
       grid_width  = unit(3, "mm"),
       ...)
}

col_ha <- HeatmapAnnotation(
  Dataset   = anno_df$Dataset,
  Condition = anno_df$Condition,
  BRCA      = anno_df$BRCA,
  col = list(
    Dataset   = PAL_DATASET,
    Condition = PAL_COND,
    BRCA      = PAL_BRCA
  ),
  annotation_height  = unit(c(1.5, 1.5, 1.5), "mm"),
  annotation_name_gp = gpar(fontsize = 6.5),
  show_legend        = c(TRUE, TRUE, TRUE),
  annotation_legend_param = list(
    Dataset   = legend_param("Dataset"),
    Condition = legend_param("Condition"),
    BRCA      = legend_param("BRCA status")
  )
)

row_ha <- rowAnnotation(
  Direction = row_df$Direction,
  col       = list(Direction = PAL_DIR),
  width     = unit(3, "mm"),
  annotation_name_side = "top",
  annotation_name_gp   = gpar(fontsize = 6.5),
  show_legend = FALSE
)

ht_heatmap <- Heatmap(
  mat_combined,
  name = "Z-score",
  col  = colorRamp2(c(-2, 0, 2), c("#1E88E5", "white", "#D81B60")),
  top_annotation        = col_ha,
  left_annotation       = row_ha,
  column_split          = anno_df$Dataset,
  row_split             = row_df$Direction,
  cluster_row_slices    = FALSE,   # UP block on top, DOWN on bottom
  cluster_column_slices = FALSE,   # dataset order preserved
  clustering_distance_rows    = "pearson",
  clustering_method_rows      = "ward.D2",
  clustering_distance_columns = "pearson",
  clustering_method_columns   = "ward.D2",
  show_column_names = FALSE,
  show_row_names    = TRUE,
  row_names_gp      = gpar(fontsize = 6.5),
  column_title_gp   = gpar(fontsize = 7.5, fontface = "bold"),
  row_title_gp      = gpar(fontsize = 8,   fontface = "bold"),
  row_title_rot     = 0,
  column_gap        = unit(2, "mm"),
  row_gap           = unit(1, "mm"),
  border            = FALSE,
  heatmap_legend_param = list(
    title        = "Z-score",
    title_gp     = gpar(fontsize = 7, fontface = "bold"),
    labels_gp    = gpar(fontsize = 6),
    grid_width   = unit(3, "mm"),
    direction    = "horizontal",
    legend_width = unit(3, "cm")
  )
)

# ============================================================
# 7. Compose panels with grid layout and export
# ============================================================
out_pdf <- file.path(fig2_dir, "Fig2.pdf")
out_png <- file.path(fig2_dir, "Fig2.png")

FIG_W_IN <- FIG_W_MM / 25.4
FIG_H_IN <- FIG_H_MM / 25.4

draw_fig2 <- function() {
  grid.newpage()
  lt <- grid.layout(
    nrow    = 2,
    ncol    = 2,
    heights = unit(c(FIG_H_MM * 0.38, FIG_H_MM * 0.62), "mm"),
    widths  = unit(c(FIG_W_MM * 0.5,  FIG_W_MM * 0.5),  "mm")
  )
  pushViewport(viewport(layout = lt))

  # Panel A — UpSet, upregulated (top-left)
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  draw(ht_upset_up, newpage = FALSE, padding = unit(c(4, 3, 2, 3), "mm"))
  grid.text("A", x = 0.02, y = 0.97, just = c("left", "top"),
            gp = gpar(fontsize = 9, fontface = "bold"))
  popViewport()

  # Panel B — UpSet, downregulated (top-right)
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  draw(ht_upset_down, newpage = FALSE, padding = unit(c(4, 3, 2, 3), "mm"))
  grid.text("B", x = 0.02, y = 0.97, just = c("left", "top"),
            gp = gpar(fontsize = 9, fontface = "bold"))
  popViewport()

  # Panel C — heatmap, full width (bottom)
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:2))
  draw(ht_heatmap, newpage = FALSE,
       heatmap_legend_side    = "right",
       annotation_legend_side = "right",
       padding = unit(c(2, 4, 3, 3), "mm"))
  grid.text("C", x = 0.01, y = 0.99, just = c("left", "top"),
            gp = gpar(fontsize = 9, fontface = "bold"))
  popViewport()

  popViewport()  # root layout
}

pdf(out_pdf, width = FIG_W_IN, height = FIG_H_IN)
draw_fig2()
invisible(dev.off())

png(out_png, width = FIG_W_MM, height = FIG_H_MM, units = "mm", res = 600)
draw_fig2()
invisible(dev.off())

message("Figure 2 saved:")
message("  PDF: ", out_pdf)
message("  PNG: ", out_png)

stop_log(log_handle)
