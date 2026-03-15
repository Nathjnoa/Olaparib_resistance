#!/usr/bin/env Rscript

# ============================================================
# 08_fig3_hallmarks.R
#
# Figure 3: Hallmark pathway enrichment in olaparib resistance
#
# Layout (double_col, 200 × 220 mm):
#   Row 1 (45%): [A: lollipop A2780 | B: lollipop PEO1 | C: lollipop UWB]
#   Row 2 (55%): [D: integrated NES heatmap, full width]
#
# Question: Do transcriptomic changes in olaparib resistance converge on
#           the same biological programs across all 3 cell line models?
#
# Prerequisites: 03_gsea_hallmarks.R
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(patchwork)
  library(grid)
})

source(file.path({
  f <- commandArgs(trailingOnly = FALSE)
  f <- grep("--file=", f, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("--file=", "", f[1])))
  else path.expand("~/bioinfo/projects/olaparib_resistance/scripts")
}, "00_config.R"))

log_handle <- start_log("08_fig3_hallmarks")
set.seed(42)

# ============================================================
# 0. Validate inputs
# ============================================================
req <- c(
  gsea_a2780 = file.path(tables_dir, "gsea_hallmarks", "A2780_GSEA_Hallmarks.tsv"),
  gsea_peo1  = file.path(tables_dir, "gsea_hallmarks", "PEO1_GSEA_Hallmarks.tsv"),
  gsea_uwb   = file.path(tables_dir, "gsea_hallmarks", "UWB_BRCAdef_GSEA_Hallmarks.tsv")
)
missing_req <- req[!file.exists(req)]
if (length(missing_req) > 0) {
  stop("Missing files (re-run 03_gsea_hallmarks.R):\n",
       paste0("  ", names(missing_req), ": ", missing_req, collapse = "\n"))
}

fig3_dir <- file.path(figures_dir, "fig3")
dir.create(fig3_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Style constants — double_col preset
# ============================================================
FIG_W_MM  <- 200
FIG_H_MM  <- 190
BASE_SIZE <- 8.5

PAL_DIR  <- c("UP" = "#D81B60", "DOWN" = "#1E88E5")
PAL_BRCA <- c("BRCA-WT" = "#a6cee3", "BRCA1-mut" = "#fb9a99", "BRCA1-def" = "#b2df8a")

# Abbreviate the longest Hallmark names for lollipop labels
abbrev_pathway <- function(x) {
  x %>%
    str_replace_all("Epithelial Mesenchymal Transition", "EMT") %>%
    str_replace_all("Interferon Gamma Response",         "IFN-\u03b3 Response") %>%
    str_replace_all("Interferon Alpha Response",         "IFN-\u03b1 Response") %>%
    str_replace_all("TNFA Signaling Via NFKB",           "TNFA via NFKB") %>%
    str_replace_all("IL6 JAK STAT3 Signaling",           "IL6-JAK-STAT3") %>%
    str_replace_all("Oxidative Phosphorylation",         "Oxidative Phosphoryl.") %>%
    str_replace_all("Estrogen Response Early",           "Estrogen Resp. Early") %>%
    str_replace_all("Estrogen Response Late",            "Estrogen Resp. Late")
}

theme_pub <- function() {
  theme_bw(base_size = BASE_SIZE) +
    theme(
      plot.title        = element_text(size = 9, face = "bold", hjust = 0.5,
                                       margin = margin(b = 2)),
      axis.title        = element_text(size = 8),
      axis.text         = element_text(size = 7),
      legend.title      = element_text(size = 7, face = "bold"),
      legend.text       = element_text(size = 7),
      legend.key.size   = unit(3, "mm"),
      legend.background = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin       = margin(3, 3, 3, 3, "mm"),
      plot.tag          = element_text(size = 9, face = "bold")
    )
}

# ============================================================
# 2. Load and prepare GSEA tables
# ============================================================
load_gsea <- function(path, dataset_label) {
  read_tsv(path, show_col_types = FALSE) %>%
    mutate(
      dataset  = dataset_label,
      pathway  = abbrev_pathway(clean_hallmark(ID)),
      log10fdr = -log10(p.adjust),
      dir      = if_else(NES > 0, "UP", "DOWN")
    ) %>%
    select(dataset, pathway, NES, fdr = p.adjust, log10fdr, dir)
}

gsea_list <- list(
  A2780 = load_gsea(req["gsea_a2780"], "A2780"),
  PEO1  = load_gsea(req["gsea_peo1"],  "PEO1"),
  UWB   = load_gsea(req["gsea_uwb"],   "UWB1.289")
)

# Global NES limits for consistent x-axis across lollipops
sig_any <- bind_rows(gsea_list) %>% filter(fdr < FDR_CUTOFF)
nes_lim <- max(abs(sig_any$NES), na.rm = TRUE)
nes_lim <- ceiling(nes_lim * 10) / 10 + 0.2  # round up + small margin

message("Global NES range for lollipops: ±", round(nes_lim, 2))

# ============================================================
# 3. Lollipop function — Panels A, B, C
# ============================================================
make_lollipop <- function(df, title, show_legend = FALSE) {
  # Select top N per direction (FDR < 0.05), sorted by NES
  plot_df <- df %>%
    filter(fdr < FDR_CUTOFF) %>%
    group_by(dir) %>%
    slice_max(order_by = abs(NES), n = 5) %>%
    ungroup() %>%
    arrange(NES) %>%
    mutate(pathway = factor(pathway, levels = pathway))

  if (nrow(plot_df) == 0) {
    # Fallback: relax to FDR < 0.1 if nothing passes strict cutoff
    plot_df <- df %>%
      filter(fdr < FDR_CUTOFF_VIZ) %>%
      group_by(dir) %>%
      slice_max(order_by = abs(NES), n = 5) %>%
      ungroup() %>%
      arrange(NES) %>%
      mutate(pathway = factor(pathway, levels = pathway))
  }

  ggplot(plot_df, aes(x = NES, y = pathway, color = dir)) +
    geom_segment(aes(x = 0, xend = NES, y = pathway, yend = pathway),
                 linewidth = 0.6) +
    geom_point(aes(size = log10fdr)) +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey40") +
    scale_color_manual(
      values = PAL_DIR,
      name   = "Direction",
      guide  = if (show_legend) guide_legend(order = 1, override.aes = list(size = 3))
               else "none"
    ) +
    scale_size_continuous(
      name   = expression(-log[10](FDR)),
      range  = c(1.5, 5),
      breaks = c(1, 2, 3),
      guide  = if (show_legend) guide_legend(order = 2) else "none"
    ) +
    scale_x_continuous(
      limits = c(-nes_lim, nes_lim),
      expand = expansion(mult = 0.02)
    ) +
    labs(
      title = title,
      x     = "NES (Resistant vs Parental)",
      y     = NULL
    ) +
    theme_pub() +
    theme(
      legend.position  = if (show_legend) "right" else "none",
      axis.text.y      = element_text(size = 7),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)
    )
}

p_A <- make_lollipop(gsea_list$A2780, "A2780",    show_legend = FALSE)
p_B <- make_lollipop(gsea_list$PEO1,  "PEO1",     show_legend = FALSE)
p_C <- make_lollipop(gsea_list$UWB,   "UWB1.289", show_legend = TRUE)

# Shared legend for dot size (extracted from p_C and placed via patchwork)
row1 <- (p_A | p_B | p_C) +
  plot_annotation(tag_levels = list(c("A", "B", "C"))) &
  theme(plot.tag = element_text(size = 9, face = "bold"))

# ============================================================
# 4. Integrated NES heatmap — Panel D
# ============================================================

# Build NES and FDR matrices (rows = pathways, cols = datasets)
all_pathways <- gsea_list$A2780$pathway  # same 50 for all (all Hallmarks tested)

nes_mat <- sapply(gsea_list, function(df) {
  v <- setNames(df$NES, df$pathway)
  v[all_pathways]
})
rownames(nes_mat) <- all_pathways

fdr_mat <- sapply(gsea_list, function(df) {
  v <- setNames(df$fdr, df$pathway)
  v[all_pathways]
})
rownames(fdr_mat) <- all_pathways

# Select pathways significant (FDR < 0.1) in ≥2 datasets
n_sig <- rowSums(fdr_mat < FDR_CUTOFF_VIZ, na.rm = TRUE)
keep  <- names(n_sig)[n_sig >= 2]
message("Pathways in integrated heatmap (FDR<0.1 in ≥2 datasets): ", length(keep))

nes_hm  <- nes_mat[keep, , drop = FALSE]
fdr_hm  <- fdr_mat[keep, , drop = FALSE]

# Column names: clean dataset labels
colnames(nes_hm) <- c("A2780", "PEO1", "UWB1.289")
colnames(fdr_hm) <- c("A2780", "PEO1", "UWB1.289")

# Significance annotation function (asterisks per cell)
sig_cell_fun <- function(j, i, x, y, width, height, fill) {
  fdr_val <- fdr_hm[i, j]
  if (!is.na(fdr_val) && fdr_val < FDR_CUTOFF) {
    grid.text("*", x, y,
              gp = gpar(fontsize = 9, col = "black", fontface = "bold"))
  } else if (!is.na(fdr_val) && fdr_val < FDR_CUTOFF_VIZ) {
    grid.text("\u00b7", x, y,              # middle dot for 0.05 < FDR < 0.1
              gp = gpar(fontsize = 11, col = "grey40"))
  }
}

# NES color scale: symmetric, centered at 0
nes_max  <- ceiling(max(abs(nes_hm), na.rm = TRUE) * 10) / 10
col_nes  <- colorRamp2(c(-nes_max, 0, nes_max),
                       c("#1E88E5", "white", "#D81B60"))

# Transpose: datasets as rows (3), pathways as columns (17)
nes_hm_t <- t(nes_hm)
fdr_hm_t <- t(fdr_hm)

# cell_fun for transposed matrix (i = dataset row, j = pathway col)
sig_cell_fun_t <- function(j, i, x, y, width, height, fill) {
  fdr_val <- fdr_hm_t[i, j]
  if (!is.na(fdr_val) && fdr_val < FDR_CUTOFF) {
    grid.text("*", x, y,
              gp = gpar(fontsize = 9, col = "black", fontface = "bold"))
  } else if (!is.na(fdr_val) && fdr_val < FDR_CUTOFF_VIZ) {
    grid.text("\u00b7", x, y,
              gp = gpar(fontsize = 11, col = "grey40"))
  }
}

# Left row annotation: BRCA status (one per dataset row)
row_left_ha <- rowAnnotation(
  "BRCA" = c("BRCA-WT", "BRCA1-mut", "BRCA1-def"),
  col    = list("BRCA" = PAL_BRCA),
  width  = unit(3, "mm"),
  annotation_name_gp   = gpar(fontsize = 6.5),
  annotation_name_side = "top",
  annotation_legend_param = list(
    "BRCA" = list(
      title     = "BRCA status",
      title_gp  = gpar(fontsize = 7, fontface = "bold"),
      labels_gp = gpar(fontsize = 6),
      grid_height = unit(3, "mm"), grid_width = unit(3, "mm")
    )
  )
)

ht_nes <- Heatmap(
  nes_hm_t,
  name             = "NES",
  col              = col_nes,
  left_annotation  = row_left_ha,
  cluster_rows     = FALSE,           # keep A2780 → PEO1 → UWB1.289 order
  cluster_columns  = TRUE,
  clustering_distance_columns = "euclidean",
  clustering_method_columns   = "ward.D2",
  column_dend_side  = "top",
  row_names_side    = "left",
  column_names_side = "bottom",
  column_names_rot  = 45,
  row_names_gp      = gpar(fontsize = 8, fontface = "bold"),
  column_names_gp   = gpar(fontsize = 7),
  cell_fun          = sig_cell_fun_t,
  column_title      = NULL,
  rect_gp           = gpar(col = "white", lwd = 0.8),
  height            = unit(30, "mm"),   # 3 rows × 10 mm
  heatmap_legend_param = list(
    title      = "NES",
    title_gp   = gpar(fontsize = 7, fontface = "bold"),
    labels_gp  = gpar(fontsize = 6),
    grid_width = unit(3, "mm"),
    direction  = "vertical",
    at         = c(-round(nes_max), -1, 0, 1, round(nes_max))
  )
)

# ============================================================
# 5. Compose figure with grid viewports and export
# ============================================================
out_pdf <- file.path(fig3_dir, "Fig3.pdf")
out_png <- file.path(fig3_dir, "Fig3.png")

FIG_W_IN <- FIG_W_MM / 25.4
FIG_H_IN <- FIG_H_MM / 25.4

draw_fig3 <- function() {
  grid.newpage()
  lt <- grid.layout(
    nrow    = 2,
    ncol    = 1,
    heights = unit(c(FIG_H_MM * 0.54, FIG_H_MM * 0.46), "mm")
  )
  pushViewport(viewport(layout = lt))

  # Row 1: lollipop panels A, B, C (patchwork)
  pushViewport(viewport(layout.pos.row = 1))
  print(row1, newpage = FALSE)
  popViewport()

  # Row 2: NES heatmap panel D
  pushViewport(viewport(layout.pos.row = 2))
  draw(ht_nes, newpage = FALSE,
       heatmap_legend_side    = "right",
       annotation_legend_side = "right",
       padding = unit(c(3, 4, 3, 4), "mm"))
  grid.text("D", x = 0.01, y = 0.99, just = c("left", "top"),
            gp = gpar(fontsize = 9, fontface = "bold"))
  popViewport()

  popViewport()
}

pdf(out_pdf, width = FIG_W_IN, height = FIG_H_IN)
draw_fig3()
invisible(dev.off())

png(out_png, width = FIG_W_MM, height = FIG_H_MM, units = "mm", res = 600)
draw_fig3()
invisible(dev.off())

message("Figure 3 saved:")
message("  PDF: ", out_pdf)
message("  PNG: ", out_png)

stop_log(log_handle)
