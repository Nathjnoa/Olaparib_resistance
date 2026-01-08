#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

proj_dir <- Sys.getenv("OLAPARIB_RESISTANCE_DIR", unset = "")
if (!nzchar(proj_dir)) {
  proj_dir <- "~/bioinfo/projects/olaparib_resistance"
}
proj_dir <- path.expand(proj_dir)
if (!dir.exists(proj_dir)) {
  stop("No existe el directorio del proyecto: ", proj_dir)
}
int_dir  <- file.path(proj_dir, "results", "tables", "integrated")
fig_dir  <- file.path(proj_dir, "results", "figures", "meta")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

in_long <- file.path(int_dir, "meta_DE_input_long.tsv")
in_meta <- file.path(int_dir, "meta_DE_random_effects.tsv")

out_png <- file.path(fig_dir, "meta_forest_top_genes.png")
out_pdf <- sub("\\.png$", ".pdf", out_png)

missing_files <- c(in_long, in_meta)
missing_files <- missing_files[!file.exists(missing_files)]
if (length(missing_files) > 0) {
  stop("No encuentro archivos de entrada:\n", paste(missing_files, collapse = "\n"))
}

# -------------------------------
# Load
# -------------------------------
deg_long <- readr::read_tsv(in_long, show_col_types = FALSE) %>%
  dplyr::select(gene_symbol, dataset, log2FC, SE) %>%
  filter(!is.na(log2FC), !is.na(SE), SE > 0)

meta <- readr::read_tsv(in_meta, show_col_types = FALSE) %>%
  dplyr::select(
    gene_symbol, k, beta_meta, se_meta, p_meta, FDR_meta, ci_lb, ci_ub, I2, direction
  ) %>%
  filter(!is.na(beta_meta), !is.na(ci_lb), !is.na(ci_ub), !is.na(FDR_meta))

# -------------------------------
# Select genes to plot
# -------------------------------
thr_fdr <- 0.05
n_each  <- 3          # top UP + top DOWN
max_I2  <- 100        # cambia a 75 si quieres evitar heterogeneidad alta

meta_sig <- meta %>%
  filter(FDR_meta < thr_fdr, !is.na(direction), I2 <= max_I2)

top_up <- meta_sig %>%
  filter(direction == "UP") %>%
  arrange(FDR_meta, desc(abs(beta_meta))) %>%
  slice_head(n = n_each)

top_down <- meta_sig %>%
  filter(direction == "DOWN") %>%
  arrange(FDR_meta, desc(abs(beta_meta))) %>%
  slice_head(n = n_each)

genes_plot <- bind_rows(top_up, top_down) %>% pull(gene_symbol) %>% unique()

if (length(genes_plot) == 0) {
  stop("No hay genes para graficar con FDR_meta < ", thr_fdr, ".")
}

meta_plot <- meta %>%
  filter(gene_symbol %in% genes_plot) %>%
  mutate(
    FDR_fmt = if_else(FDR_meta <= 1e-300, "FDR < 1e-300",
                      format(FDR_meta, scientific = TRUE, digits = 2)),
    gene_lab = paste0(
      gene_symbol,
      "  (k=", k, ", I²=", sprintf("%.1f", I2), "%, ", FDR_fmt, ")"
    )
  )

# -------------------------------
# Build forest dataframe
# -------------------------------
df_studies <- deg_long %>%
  filter(gene_symbol %in% genes_plot) %>%
  mutate(
    ci_lb = log2FC - 1.96 * SE,
    ci_ub = log2FC + 1.96 * SE,
    type  = "Study",
    label = dataset
  ) %>%
  select(gene_symbol, dataset, label, log2FC, ci_lb, ci_ub, type)

df_meta <- meta_plot %>%
  transmute(
    gene_symbol,
    dataset = "META",
    label   = "META (RE)",
    log2FC  = beta_meta,
    ci_lb   = ci_lb,
    ci_ub   = ci_ub,
    type    = "Meta"
  )

df_forest <- bind_rows(df_studies, df_meta) %>%
  left_join(meta_plot %>% select(gene_symbol, gene_lab, direction), by = "gene_symbol") %>%
  mutate(
    label = factor(label, levels = rev(c("A2780", "PEO1", "UWB_BRCAdef", "UWB_BRCAprof", "META (RE)")))
  )

# Para ordenar genes (arriba: UP, abajo: DOWN; dentro por FDR)
gene_order <- meta_plot %>%
  arrange(factor(direction, levels = c("UP","DOWN")), FDR_meta, desc(abs(beta_meta))) %>%
  pull(gene_lab)

df_forest <- df_forest %>%
  mutate(
    gene_lab = str_wrap(gene_lab, width = 30),
    gene_lab = factor(gene_lab, levels = str_wrap(gene_order, width = 30))
  )

# xlim global (log2FC) basado en CI (estudios + meta), usando percentil 99% para evitar outliers extremos
abs_ci <- abs(c(df_forest$ci_lb, df_forest$ci_ub))
xmax <- quantile(abs_ci[is.finite(abs_ci)], 0.99, na.rm = TRUE)
if (!is.finite(xmax) || xmax <= 0) {
  xmax <- max(abs_ci[is.finite(abs_ci)], na.rm = TRUE)
}
if (!is.finite(xmax) || xmax <= 0) xmax <- 1
xmax <- xmax * 1.05

# -------------------------------
# Plot (ggplot forest)
# -------------------------------
df_forest_plot <- df_forest %>%
  filter(!is.na(log2FC), !is.na(ci_lb), !is.na(ci_ub))

p <- ggplot(df_forest_plot, aes(x = label, y = log2FC)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey40") +
  geom_linerange(aes(ymin = ci_lb, ymax = ci_ub), linewidth = 0.8, color = "black") +
  geom_point(
    aes(shape = type),
    size = 2.8,
    fill = "white",
    color = "black",
    stroke = 0.7
  ) +
  # resaltar META
  geom_point(
    data = df_forest %>% filter(type == "Meta"),
    aes(x = label, y = log2FC),
    size = 4,
    shape = 23, fill = "black", color = "black"
  ) +
  coord_flip() +
  facet_wrap(~ gene_lab, ncol = 2, scales = "fixed") +
  scale_shape_manual(values = c(Study = 21, Meta = 23), guide = "none") +
  scale_y_continuous(limits = c(-xmax, xmax)) +
  labs(
    title = "Meta-analysis (random effects): consistencia por gen",
    subtitle = paste0("Efectos por dataset (log2FC Res-Par) y efecto meta (RE) con IC95%.  Genes: top ", n_each, " UP + top ", n_each, " DOWN (FDR_meta<", thr_fdr, ")"),
    x = NULL,
    y = "log2FC (Resistant - Parental)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 12),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_blank()
  )

ggsave(out_png, p, width = 12.5, height = 7.5, dpi = 300)
ggsave(out_pdf, p, width = 12.5, height = 7.5)

message("✔ Forest plot guardado en:\n  ", out_png, "\n  ", out_pdf)
