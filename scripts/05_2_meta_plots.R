#!/usr/bin/env Rscript

# ============================================================
# 05_2_meta_plots.R
#
# Meta-analysis visualizations (3 datasets):
#
#   Section 1: Meta-volcano plot (beta_meta vs -log10(FDR))
#   Section 2: Lollipop of top UP/DOWN genes from meta
#   Section 3: Forest plots for top UP + top DOWN genes
#
# CLI: Rscript 05_2_meta_plots.R [fdr_thr] [top_each_dir] [beta_thr] [label_n]
#   fdr_thr      — FDR threshold (default: 0.05)
#   top_each_dir — top genes per direction (default: 10)
#   beta_thr     — |beta| guide lines in volcano (default: 0.5)
#   label_n      — genes labeled in volcano (default: 12)
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

source(file.path({f <- commandArgs(trailingOnly=FALSE); f <- grep("--file=",f,value=TRUE); if(length(f)) dirname(normalizePath(sub("--file=","",f[1]))) else path.expand("~/bioinfo/projects/olaparib_resistance/scripts")}, "00_config.R"))

log_handle <- start_log("05_2_meta_plots")

int_dir <- file.path(tables_dir, "integrated")
fig_dir <- file.path(figures_dir, "meta_analysis")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# CLI args
# ============================================================
fdr_thr      <- 0.05
beta_thr     <- 0.5
top_each_dir <- 10
label_n      <- 12
use_highconf <- TRUE

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) fdr_thr      <- as.numeric(args[1])
if (length(args) >= 2) top_each_dir <- as.integer(args[2])
if (length(args) >= 3) beta_thr     <- as.numeric(args[3])
if (length(args) >= 4) label_n      <- as.integer(args[4])

message("Params: fdr_thr=", fdr_thr,
        " | top_each_dir=", top_each_dir,
        " | beta_thr=", beta_thr,
        " | label_n=", label_n)

# ============================================================
# Load + sanitize meta-analysis results
# ============================================================
meta_file <- file.path(int_dir, "meta_DE_random_effects.tsv")
long_file <- file.path(int_dir, "meta_DE_input_long.tsv")

missing_files <- c(meta_file, long_file)[!file.exists(c(meta_file, long_file))]
if (length(missing_files) > 0) {
  stop("Archivos de entrada no encontrados:\n", paste(missing_files, collapse = "\n"))
}

meta <- readr::read_tsv(meta_file, show_col_types = FALSE)

num_cols <- intersect(
  c("beta_meta","se_meta","z_meta","p_meta","FDR_meta","I2","k",
    "beta_meta_DL","se_meta_DL","p_meta_DL","tau2_DL"),
  colnames(meta)
)
meta <- meta %>%
  mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.x))))

if (!("p_meta" %in% names(meta))) stop("Columna p_meta no encontrada en ", meta_file)

meta <- meta %>%
  mutate(
    p_meta2   = pmax(p_meta, 1e-300),
    FDR_meta2 = p.adjust(p_meta2, method = "BH")
  )

if (!("gene_symbol" %in% names(meta))) stop("Columna gene_symbol no encontrada.")
if (!("beta_meta"   %in% names(meta))) stop("Columna beta_meta no encontrada.")

# ============================================================
# SECTION 1: Meta-volcano plot
# ============================================================
message("=== Section 1: Meta-volcano ===")

pal_dir2 <- c("UP" = "#D81B60", "DOWN" = "#1E88E5", "NS" = "#9E9E9E")

volc <- meta %>%
  mutate(
    neglog10FDR     = -log10(FDR_meta2),
    sig = case_when(
      !is.na(FDR_meta2) & FDR_meta2 < fdr_thr & beta_meta >  0 ~ "UP",
      !is.na(FDR_meta2) & FDR_meta2 < fdr_thr & beta_meta <  0 ~ "DOWN",
      TRUE ~ "NS"
    ),
    neglog10FDR_cap = pmin(neglog10FDR, 100)
  )

lab_df <- volc %>%
  filter(FDR_meta2 < fdr_thr, is.finite(beta_meta), is.finite(neglog10FDR)) %>%
  arrange(FDR_meta2, desc(abs(beta_meta))) %>%
  slice_head(n = label_n)

p_volcano <- ggplot(volc, aes(x = beta_meta, y = neglog10FDR_cap)) +
  geom_point(aes(color = sig, shape = sig), alpha = 0.7, size = 1.6) +
  scale_color_manual(values = pal_dir2, name = "Status", drop = FALSE) +
  scale_shape_manual(values = c(UP = 15, DOWN = 16, NS = 17), name = "Status", drop = FALSE) +
  geom_vline(xintercept = c(-beta_thr, beta_thr), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(fdr_thr), linetype = "dashed", color = "grey50") +
  labs(
    title    = "Meta-analysis (random effects): gene-level effects",
    subtitle = paste0(
      "Volcano: beta_meta vs -log10(FDR).\n",
      "Threshold: FDR<", fdr_thr, " | guide lines at |beta|=", beta_thr
    ),
    x     = "beta_meta (log2FC, Resistant - Parental)",
    y     = "-log10(FDR) (capped at 100)",
    shape = "Status",
    color = "Status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(face = "bold"),
    plot.subtitle     = element_text(hjust = 0.5, size = 10, lineheight = 1.1, margin = margin(b = 6)),
    legend.position   = "right"
  )

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_volcano <- p_volcano +
    ggrepel::geom_text_repel(
      data        = lab_df,
      aes(label   = gene_symbol),
      size        = 3,
      max.overlaps = Inf,
      box.padding  = 0.3,
      point.padding = 0.2
    )
} else {
  p_volcano <- p_volcano +
    geom_text(data = lab_df, aes(label = gene_symbol), size = 2.8, vjust = -0.8)
}

out_vol_png <- file.path(fig_dir, "MetaVolcano_random_effects.png")
out_vol_pdf <- sub("\\.png$", ".pdf", out_vol_png)
ggsave(out_vol_png, p_volcano, width = 7.2, height = 5.8, dpi = 300)
ggsave(out_vol_pdf, p_volcano, width = 7.2, height = 5.8)
message("✔ Volcano -> ", out_vol_png)

# ============================================================
# SECTION 2: Lollipop — top UP / top DOWN genes
# ============================================================
message("=== Section 2: Lollipop top genes ===")

core_file      <- file.path(int_dir, "meta_core_highconf.tsv")
df_for_lollipop <- meta

if (use_highconf && file.exists(core_file)) {
  message("Usando meta_core_highconf.tsv para lollipop.")
  core <- readr::read_tsv(core_file, show_col_types = FALSE)
  if (!("gene_symbol" %in% names(core))) stop("meta_core_highconf.tsv no tiene gene_symbol.")
  df_for_lollipop <- meta %>% semi_join(core %>% dplyr::select(gene_symbol), by = "gene_symbol")
} else {
  message("Usando meta_DE_random_effects.tsv completo para lollipop.")
}

cap_val <- 100
sig_df <- df_for_lollipop %>%
  filter(!is.na(FDR_meta2), FDR_meta2 < fdr_thr,
         !is.na(beta_meta), !is.na(gene_symbol)) %>%
  mutate(
    direction  = if_else(beta_meta > 0, "UP", "DOWN"),
    mlog10     = -log10(FDR_meta2),
    mlog10_cap = pmin(mlog10, cap_val)
  )

top_up <- sig_df %>%
  filter(direction == "UP") %>%
  arrange(FDR_meta2, desc(abs(beta_meta))) %>%
  slice_head(n = top_each_dir)

top_down <- sig_df %>%
  filter(direction == "DOWN") %>%
  arrange(FDR_meta2, desc(abs(beta_meta))) %>%
  slice_head(n = top_each_dir)

plot_df <- bind_rows(top_down, top_up) %>%
  mutate(
    gene_plot = factor(gene_symbol, levels = c(top_down$gene_symbol, rev(top_up$gene_symbol)))
  )

if (nrow(plot_df) == 0) {
  stop("No hay genes significativos con FDR<", fdr_thr, " para lollipop.")
}

out_top_tsv <- file.path(int_dir, "meta_topgenes_for_lollipop.tsv")
readr::write_tsv(plot_df, out_top_tsv)
message("✔ Tabla top genes lollipop: ", out_top_tsv)

p_sep <- if (nrow(top_down) > 0 && nrow(top_up) > 0) nrow(top_down) + 0.5 else NA_real_

pal_dir_plot <- c("UP" = "#D81B60", "DOWN" = "#1E88E5")

p_lolli <- ggplot(plot_df, aes(x = beta_meta, y = gene_plot, color = direction)) +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.4, alpha = 0.6) +
  geom_segment(aes(x = 0, xend = beta_meta, yend = gene_plot), linewidth = 0.7, alpha = 0.8) +
  geom_point(aes(size = mlog10_cap), alpha = 0.9) +
  scale_color_manual(values = pal_dir_plot, name = "Direction") +
  labs(
    title    = "Meta-analysis: top significant genes (UP/DOWN)",
    subtitle = paste0(
      "Filter: FDR<", fdr_thr,
      " | Select: top ", top_each_dir, " UP + top ", top_each_dir,
      " DOWN\nOrder: FDR then |beta| | Size capped at ", cap_val
    ),
    x    = "beta_meta (log2FC, Resistant - Parental)",
    y    = NULL,
    size = "-log10(FDR) (capped)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y   = element_text(size = 9),
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, lineheight = 1.1, margin = margin(b = 6))
  )

if (!is.na(p_sep)) {
  p_lolli <- p_lolli + geom_hline(yintercept = p_sep, linetype = "dashed", color = "grey50")
}

out_lol_png <- file.path(fig_dir, "MetaLollipop_topGenes.png")
out_lol_pdf <- sub("\\.png$", ".pdf", out_lol_png)
ggsave(out_lol_png, p_lolli, width = 8.5, height = max(6, 0.22 * nrow(plot_df) + 2.5), dpi = 300)
ggsave(out_lol_pdf, p_lolli, width = 8.5, height = max(6, 0.22 * nrow(plot_df) + 2.5))
message("✔ Lollipop -> ", out_lol_png)

# ============================================================
# SECTION 3: Forest plots — top UP + top DOWN genes
# ============================================================
message("=== Section 3: Forest plots ===")

deg_long <- readr::read_tsv(long_file, show_col_types = FALSE) %>%
  dplyr::select(gene_symbol, dataset, log2FC, SE) %>%
  filter(!is.na(log2FC), !is.na(SE), SE > 0)

meta_forest <- readr::read_tsv(meta_file, show_col_types = FALSE) %>%
  dplyr::select(
    gene_symbol, k, beta_meta, se_meta, p_meta, FDR_meta, ci_lb, ci_ub, I2, direction
  ) %>%
  filter(!is.na(beta_meta), !is.na(ci_lb), !is.na(ci_ub), !is.na(FDR_meta))

n_each <- 3

meta_sig <- meta_forest %>%
  filter(FDR_meta < fdr_thr, !is.na(direction))

top_up_f <- meta_sig %>%
  filter(direction == "UP") %>%
  arrange(FDR_meta, desc(abs(beta_meta))) %>%
  slice_head(n = n_each)

top_down_f <- meta_sig %>%
  filter(direction == "DOWN") %>%
  arrange(FDR_meta, desc(abs(beta_meta))) %>%
  slice_head(n = n_each)

genes_plot <- bind_rows(top_up_f, top_down_f) %>% pull(gene_symbol) %>% unique()

if (length(genes_plot) == 0) {
  warning("No hay genes para forest plots con FDR_meta < ", fdr_thr, ". Omitiendo.")
} else {
  meta_plot_f <- meta_forest %>%
    filter(gene_symbol %in% genes_plot) %>%
    mutate(
      FDR_fmt = if_else(FDR_meta <= 1e-300, "FDR < 1e-300",
                        format(FDR_meta, scientific = TRUE, digits = 2)),
      gene_lab = paste0(
        gene_symbol,
        "  (k=", k, ", I\u00b2=", sprintf("%.1f", I2), "%, ", FDR_fmt, ")"
      )
    )

  df_studies <- deg_long %>%
    filter(gene_symbol %in% genes_plot) %>%
    mutate(
      ci_lb = log2FC - 1.96 * SE,
      ci_ub = log2FC + 1.96 * SE,
      type  = "Study",
      label = dataset
    ) %>%
    dplyr::select(gene_symbol, dataset, label, log2FC, ci_lb, ci_ub, type)

  df_meta_f <- meta_plot_f %>%
    transmute(
      gene_symbol,
      dataset = "META",
      label   = "META (RE)",
      log2FC  = beta_meta,
      ci_lb   = ci_lb,
      ci_ub   = ci_ub,
      type    = "Meta"
    )

  # Fix: only 3 datasets — no UWB_BRCAprof in factor levels
  study_levels <- rev(c("A2780", "PEO1", "UWB_BRCAdef", "META (RE)"))

  df_forest <- bind_rows(df_studies, df_meta_f) %>%
    left_join(meta_plot_f %>% dplyr::select(gene_symbol, gene_lab, direction), by = "gene_symbol") %>%
    mutate(label = factor(label, levels = study_levels))

  gene_order_f <- meta_plot_f %>%
    arrange(factor(direction, levels = c("UP","DOWN")), FDR_meta, desc(abs(beta_meta))) %>%
    pull(gene_lab)

  df_forest <- df_forest %>%
    mutate(
      gene_lab = stringr::str_wrap(gene_lab, width = 30),
      gene_lab = factor(gene_lab, levels = stringr::str_wrap(gene_order_f, width = 30))
    )

  abs_ci <- abs(c(df_forest$ci_lb, df_forest$ci_ub))
  xmax   <- quantile(abs_ci[is.finite(abs_ci)], 0.99, na.rm = TRUE)
  if (!is.finite(xmax) || xmax <= 0) xmax <- max(abs_ci[is.finite(abs_ci)], na.rm = TRUE)
  if (!is.finite(xmax) || xmax <= 0) xmax <- 1
  xmax <- xmax * 1.05

  df_forest_plot <- df_forest %>%
    filter(!is.na(log2FC), !is.na(ci_lb), !is.na(ci_ub))

  p_forest <- ggplot(df_forest_plot, aes(x = label, y = log2FC)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey40") +
    geom_linerange(aes(ymin = ci_lb, ymax = ci_ub), linewidth = 0.8, color = "black") +
    geom_point(
      aes(shape = type),
      size = 2.8, fill = "white", color = "black", stroke = 0.7
    ) +
    geom_point(
      data = df_forest_plot %>% filter(type == "Meta"),
      aes(x = label, y = log2FC),
      size = 4, shape = 23, fill = "black", color = "black"
    ) +
    coord_flip() +
    facet_wrap(~ gene_lab, ncol = 2, scales = "fixed") +
    scale_shape_manual(values = c(Study = 21, Meta = 23), guide = "none") +
    scale_y_continuous(limits = c(-xmax, xmax)) +
    labs(
      title    = "Meta-analysis (random effects): per-gene consistency",
      subtitle = paste0(
        "Study effects (log2FC Res-Par) + meta effect (RE) with 95% CI.\n",
        "Genes: top ", n_each, " UP + top ", n_each, " DOWN (FDR_meta<", fdr_thr, ")"
      ),
      x = NULL,
      y = "log2FC (Resistant - Parental)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title        = element_text(face = "bold", size = 17),
      plot.subtitle     = element_text(size = 12),
      strip.text        = element_text(face = "bold", size = 12),
      axis.text.y       = element_text(size = 10),
      axis.text.x       = element_text(size = 10),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor  = element_blank()
    )

  out_forest_png <- file.path(fig_dir, "meta_forest_top_genes.png")
  out_forest_pdf <- sub("\\.png$", ".pdf", out_forest_png)
  ggsave(out_forest_png, p_forest, width = 12.5, height = 7.5, dpi = 300)
  ggsave(out_forest_pdf, p_forest, width = 12.5, height = 7.5)
  message("✔ Forest plots -> ", out_forest_png)
}

stop_log(log_handle)
