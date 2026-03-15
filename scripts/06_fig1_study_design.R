#!/usr/bin/env Rscript

# ============================================================
# 06_fig1_study_design.R
#
# Figure 1: Study overview and per-dataset transcriptomic validation
#
# Layout:
#   Row 1: [A: Study design diagram (60%)] [D: DEG counts bar (40%)]
#   Row 2: [B1: A2780 PCA] [B2: PEO1 PCA] [B3: UWB PCA]
#   Row 3: [C1: A2780 Volcano] [C2: PEO1 Volcano] [C3: UWB Volcano]
#
# Preset: double_col  — 200 mm × 195 mm, base_size = 8.5 pt, 600 dpi
#
# Prerequisites: run 01_differential_expression.R first to generate
#   data/processed/pca_coords_{A2780,PEO1,UWB}.tsv
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})

source(file.path({
  f <- commandArgs(trailingOnly = FALSE)
  f <- grep("--file=", f, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("--file=", "", f[1])))
  else path.expand("~/bioinfo/projects/olaparib_resistance/scripts")
}, "00_config.R"))

log_handle <- start_log("06_fig1_study_design")
set.seed(42)

# ============================================================
# 0. Validate required inputs
# ============================================================
req <- c(
  pca_a2780  = file.path(processed_dir, "pca_coords_A2780.tsv"),
  pca_peo1   = file.path(processed_dir, "pca_coords_PEO1.tsv"),
  pca_uwb    = file.path(processed_dir, "pca_coords_UWB.tsv"),
  de_a2780   = file.path(tables_dir, "de", "GSE153867_A2780_limma_DE_FPKM.tsv"),
  de_peo1    = file.path(tables_dir, "de", "GSE117765_PEO1_limma_TPM_Res_vs_Par.tsv"),
  de_uwb     = file.path(tables_dir, "de", "GSE235980_BRCAdef_Res_vs_Par_DESeq2.tsv"),
  deg_counts = file.path(tables_dir, "venn", "DEG_counts_for_venn.tsv")
)
missing_req <- req[!file.exists(req)]
if (length(missing_req) > 0) {
  stop("Archivos faltantes (re-ejecutar 01_differential_expression.R):\n",
       paste0("  ", names(missing_req), ": ", missing_req, collapse = "\n"))
}

# ============================================================
# 1. Theme y paletas — double_col preset
# ============================================================
FIG_W_MM  <- 200
FIG_H_MM  <- 195
BASE_SIZE <- 8.5

theme_pub <- function() {
  theme_bw(base_size = BASE_SIZE) +
    theme(
      plot.title        = element_text(size = 9, face = "bold", hjust = 0.5,
                                       margin = margin(b = 2)),
      plot.subtitle     = element_text(size = 7, hjust = 0.5, color = "grey40",
                                       margin = margin(t = 0, b = 2)),
      axis.title        = element_text(size = 8),
      axis.text         = element_text(size = 7),
      legend.title      = element_text(size = 7, face = "bold"),
      legend.text       = element_text(size = 7),
      legend.key.size   = unit(3, "mm"),
      legend.background = element_blank(),
      panel.grid.minor  = element_blank(),
      plot.margin       = margin(3, 3, 3, 3, "mm"),
      plot.tag          = element_text(size = 9, face = "bold")
    )
}

# Colorblind-safe (Okabe-Ito): Parental/Resistant
PAL_COND  <- c("Parental"  = "#56B4E9", "Resistant" = "#E69F00")
SHP_COND  <- c("Parental"  = 16,        "Resistant" = 17)

# ============================================================
# 2. Panel A — Diagrama de diseño del estudio
# ============================================================
# Coordinate space: x [0, 10], y [0, 3.5]
# 3 filas de datos (y=0.5, 1.5, 2.5) + encabezado (y=3.1)

ds_rows <- tibble(
  y       = c(2.5,       1.5,           0.5),
  ds_name = c("A2780",   "PEO1",        "UWB1.289"),
  gse     = c("GSE153867", "GSE117765", "GSE235980"),
  brca    = c("BRCA-WT", "BRCA1-mut",   "BRCA1-def"),
  n_par   = c(8L, 2L, 2L),
  n_res   = c(8L, 4L, 2L),
  bg      = c("#D6E8F5", "#F5D6D6", "#D6F0D6")
)

# posiciones x de cada columna (sin Method — redistribuidas)
X_DS   <- 1.0
X_BRCA <- 3.0
X_PAR  <- 5.2
X_ARR  <- 6.5    # centro de la flecha
X_RES  <- 7.8

panel_A <-
  ggplot() +
  # ── fondos de fila ──
  geom_rect(data = ds_rows,
            aes(xmin = 0, xmax = 10, ymin = y - 0.44, ymax = y + 0.44, fill = bg),
            alpha = 0.55, show.legend = FALSE) +
  scale_fill_identity() +

  # ── separadores horizontales ──
  geom_hline(yintercept = c(2.95, -0.01), linewidth = 0.4, color = "grey50") +
  geom_hline(yintercept = c(2.06, 1.06),  linewidth = 0.2, color = "grey75",
             linetype = "dashed") +

  # ── encabezados de columna ──
  annotate("text", x = X_DS,   y = 3.2, label = "Dataset",     size = 2.4,
           fontface = "bold", hjust = 0.5) +
  annotate("text", x = X_BRCA, y = 3.2, label = "BRCA status", size = 2.4,
           fontface = "bold", hjust = 0.5) +
  annotate("text", x = (X_PAR + X_RES) / 2, y = 3.2,
           label = "Design", size = 2.4, fontface = "bold", hjust = 0.5) +

  # ── nombre dataset + GSE ──
  geom_text(data = ds_rows,
            aes(x = X_DS, y = y + 0.12, label = ds_name),
            size = 2.6, fontface = "bold", hjust = 0.5) +
  geom_text(data = ds_rows,
            aes(x = X_DS, y = y - 0.15, label = gse),
            size = 1.9, color = "grey45", hjust = 0.5) +

  # ── badge BRCA ──
  geom_tile(data = ds_rows,
            aes(x = X_BRCA, y = y, width = 1.6, height = 0.42, fill = bg),
            color = "grey40", linewidth = 0.35, show.legend = FALSE) +
  geom_text(data = ds_rows,
            aes(x = X_BRCA, y = y, label = brca),
            size = 2.2, fontface = "italic", hjust = 0.5) +

  # ── círculos Parental ──
  geom_point(data = ds_rows,
             aes(x = X_PAR, y = y),
             shape = 21, size = 9, fill = PAL_COND["Parental"],
             color = "#2B7CBE", stroke = 0.5) +
  geom_text(data = ds_rows,
            aes(x = X_PAR, y = y + 0.30, label = "Parental"),
            size = 1.9, hjust = 0.5) +
  geom_text(data = ds_rows,
            aes(x = X_PAR, y = y - 0.31, label = paste0("n = ", n_par)),
            size = 1.9, color = "grey35", hjust = 0.5) +

  # ── flecha ──
  geom_segment(data = ds_rows,
               aes(x = X_PAR + 0.42, xend = X_RES - 0.42, y = y, yend = y),
               arrow = arrow(length = unit(1.8, "mm"), type = "closed"),
               linewidth = 0.5, color = "grey25") +

  # ── círculos Resistant ──
  geom_point(data = ds_rows,
             aes(x = X_RES, y = y),
             shape = 21, size = 9, fill = PAL_COND["Resistant"],
             color = "#A86800", stroke = 0.5) +
  geom_text(data = ds_rows,
            aes(x = X_RES, y = y + 0.30, label = "Resistant"),
            size = 1.9, hjust = 0.5) +
  geom_text(data = ds_rows,
            aes(x = X_RES, y = y - 0.31, label = paste0("n = ", n_res)),
            size = 1.9, color = "grey35", hjust = 0.5) +

  scale_x_continuous(limits = c(0, 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.02, 3.42), expand = c(0, 0)) +
  labs(tag = "A") +
  theme_void() +
  theme(
    plot.tag    = element_text(size = 9, face = "bold"),
    plot.margin = margin(3, 2, 3, 3, "mm")
  )

# ============================================================
# 3. Panel D — DEG counts (barras espejo UP/DOWN)
# ============================================================
deg_raw <- readr::read_tsv(req["deg_counts"], show_col_types = FALSE)

deg_long <- deg_raw %>%
  tidyr::separate(dataset,
                  into = c("ds_raw", "direction"),
                  sep  = "_(?=[^_]+$)") %>%
  mutate(
    direction = str_to_title(direction),
    ds_label  = case_when(
      ds_raw == "A2780"      ~ "A2780",
      ds_raw == "UWB_BRCAdef" ~ "UWB1.289",
      ds_raw == "PEO1"       ~ "PEO1",
      TRUE                   ~ ds_raw
    ),
    ds_label = factor(ds_label, levels = c("A2780", "PEO1", "UWB1.289")),
    y_val    = ifelse(direction == "Up", n_genes, -n_genes),
    label_y  = ifelse(direction == "Up",
                      n_genes + 180,
                      -(n_genes + 180))
  )

y_max <- max(deg_long$n_genes) * 1.20

panel_D <-
  ggplot(deg_long, aes(x = ds_label, y = y_val, fill = direction)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "grey20") +
  geom_text(aes(y = label_y, label = comma(n_genes)),
            size = 2, hjust = 0.5) +
  scale_fill_manual(
    values = c("Up" = "#D81B60", "Down" = "#1E88E5"),
    name   = "Direction",
    breaks = c("Up", "Down")
  ) +
  scale_y_continuous(
    limits = c(-y_max, y_max),
    labels = function(x) comma(abs(x)),
    breaks = pretty(c(-y_max, y_max), n = 5)
  ) +
  labs(
    title = "DEGs per dataset",
    x     = NULL,
    y     = "Number of DEGs\n(FDR < 0.05, |log\u2082FC| \u2265 0.58)",
    tag   = "D"
  ) +
  theme_pub() +
  theme(legend.position = "right")

# ============================================================
# 4. Paneles B — PCA por dataset
# ============================================================
pca_a2780 <- readr::read_tsv(req["pca_a2780"], show_col_types = FALSE)
pca_peo1  <- readr::read_tsv(req["pca_peo1"],  show_col_types = FALSE)
pca_uwb   <- readr::read_tsv(req["pca_uwb"],   show_col_types = FALSE)

make_pca <- function(pca_df, title, show_legend = FALSE) {
  pc1_pct <- unique(pca_df$pct_PC1)
  pc2_pct <- unique(pca_df$pct_PC2)

  p <- ggplot(pca_df,
              aes(PC1, PC2, color = condition, shape = condition)) +
    geom_point(size = 2.0, alpha = 0.9, stroke = 0.3) +
    scale_color_manual(values = PAL_COND, name = "Condition",
                       breaks = c("Parental", "Resistant")) +
    scale_shape_manual(values = SHP_COND,  name = "Condition",
                       breaks = c("Parental", "Resistant")) +
    labs(
      title = title,
      x     = paste0("PC1 (", pc1_pct, "%)"),
      y     = paste0("PC2 (", pc2_pct, "%)")
    ) +
    theme_pub()

  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

panel_B1 <- make_pca(pca_a2780, "A2780",          show_legend = FALSE) + labs(tag = "B")
panel_B2 <- make_pca(pca_peo1,  "PEO1",            show_legend = FALSE)
panel_B3 <- make_pca(pca_uwb,   "UWB1.289 BRCA-def", show_legend = TRUE) +
  theme(legend.position = "right")

# ============================================================
# 5. Paneles C — Volcano por dataset
# ============================================================
# Carga y estandariza tablas DE
load_de <- function(path, logfc_col, pval_col, dataset_label) {
  readr::read_tsv(path, show_col_types = FALSE) %>%
    dplyr::rename(log2FC = !!logfc_col, FDR = !!pval_col) %>%
    dplyr::select(log2FC, FDR) %>%
    filter(!is.na(log2FC), !is.na(FDR), is.finite(log2FC)) %>%
    mutate(
      # Limitar y-axis a percentil 99.5 para no comprimir la nube
      neg_log10_FDR = pmin(-log10(FDR), quantile(-log10(FDR), 0.995, na.rm = TRUE)),
      sig = case_when(
        FDR < FDR_CUTOFF & log2FC >  LOGFC_CUTOFF ~ "Up (Resistant)",
        FDR < FDR_CUTOFF & log2FC < -LOGFC_CUTOFF ~ "Down (Resistant)",
        TRUE ~ "NS"
      ) %>% factor(levels = c("Up (Resistant)", "Down (Resistant)", "NS")),
      dataset = dataset_label
    )
}

de_a2780 <- load_de(req["de_a2780"], "log2FC",          "adj.P.Val", "A2780")
de_peo1  <- load_de(req["de_peo1"],  "log2FC",          "adj.P.Val", "PEO1")
de_uwb   <- load_de(req["de_uwb"],   "log2FoldChange",  "padj",      "UWB1.289")

# Subtítulo con recuento de DEGs
make_subtitle <- function(de_df) {
  n_up   <- sum(de_df$sig == "Up (Resistant)")
  n_down <- sum(de_df$sig == "Down (Resistant)")
  paste0("\u2191 ", comma(n_up), "  \u2193 ", comma(n_down))
}

make_volcano <- function(de_df, title, show_legend = FALSE) {
  # Límites simétricos del eje x (percentil 99 ± 10%)
  x_lim <- max(abs(quantile(de_df$log2FC, c(0.005, 0.995), na.rm = TRUE))) * 1.1

  p <- ggplot(de_df, aes(log2FC, neg_log10_FDR, color = sig, shape = sig)) +
    geom_point(alpha = 0.55, size = 0.7, stroke = 0) +
    geom_vline(xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF),
               linetype = "dashed", linewidth = 0.35, color = "grey50") +
    geom_hline(yintercept = -log10(FDR_CUTOFF),
               linetype = "dashed", linewidth = 0.35, color = "grey50") +
    scale_color_manual(values = volcano_pal,    name = "Significance",
                       drop = FALSE) +
    scale_shape_manual(values = volcano_shapes, name = "Significance",
                       drop = FALSE) +
    scale_x_continuous(limits = c(-x_lim, x_lim), oob = squish) +
    labs(
      title    = title,
      subtitle = make_subtitle(de_df),
      x        = "log\u2082FC (Resistant / Parental)",
      y        = "-log\u2081\u2080(FDR)"
    ) +
    theme_pub()

  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

panel_C1 <- make_volcano(de_a2780, "A2780",             show_legend = FALSE) + labs(tag = "C")
panel_C2 <- make_volcano(de_peo1,  "PEO1",              show_legend = FALSE)
panel_C3 <- make_volcano(de_uwb,   "UWB1.289 BRCA-def", show_legend = TRUE) +
  theme(legend.position = "right")

# ============================================================
# 6. Ensamblaje con patchwork
# ============================================================
row1 <- (panel_A | panel_D)     + plot_layout(widths = c(3, 2))
row2 <- (panel_B1 | panel_B2 | panel_B3)
row3 <- (panel_C1 | panel_C2 | panel_C3)

fig1 <- (row1 / row2 / row3) +
  plot_layout(heights = c(1.3, 1, 1))

# ============================================================
# 7. Export
# ============================================================
fig1_dir <- file.path(figures_dir, "fig1")
dir.create(fig1_dir, showWarnings = FALSE, recursive = TRUE)

out_pdf <- file.path(fig1_dir, "Fig1.pdf")
out_png <- file.path(fig1_dir, "Fig1.png")

# PDF vectorial (cairo para fuentes embebidas)
pdf_dev <- if (capabilities("cairo")) cairo_pdf else pdf
ggsave(out_pdf,
       plot    = fig1,
       device  = pdf_dev,
       width   = FIG_W_MM, height = FIG_H_MM,
       units   = "mm",
       limitsize = FALSE)

# PNG de alta resolución
ggsave(out_png,
       plot    = fig1,
       width   = FIG_W_MM, height = FIG_H_MM,
       units   = "mm",
       dpi     = 600,
       limitsize = FALSE)

message("Fig1 exportada:")
message("  PDF -> ", out_pdf)
message("  PNG -> ", out_png)

# ============================================================
# 8. QC: verificar que los archivos existen y tienen tamaño razonable
# ============================================================
for (f in c(out_pdf, out_png)) {
  sz <- file.size(f)
  if (is.na(sz) || sz < 50000) {
    warning("Archivo sospechosamente pequeño o inexistente: ", f, " (", sz, " bytes)")
  } else {
    message("QC OK: ", basename(f), " (", round(sz / 1024), " KB)")
  }
}

stop_log(log_handle)
