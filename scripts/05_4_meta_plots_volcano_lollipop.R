#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# Paletas
pal_dir <- c("UP" = "#D81B60", "DOWN" = "#1E88E5", "NS" = "#9E9E9E")
pal_dir2 <- c("UP" = "#D81B60", "DOWN" = "#1E88E5")

# ============================================================
# Config
# ============================================================
proj_dir <- Sys.getenv("OLAPARIB_RESISTANCE_DIR", unset = "")
if (!nzchar(proj_dir)) {
  proj_dir <- "~/bioinfo/projects/olaparib_resistance"
}
proj_dir <- path.expand(proj_dir)
if (!dir.exists(proj_dir)) {
  stop("No existe el directorio del proyecto: ", proj_dir)
}
int_dir  <- file.path(proj_dir, "results", "tables", "integrated")
fig_dir  <- file.path(proj_dir, "results", "figures", "meta_analysis")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

meta_file <- file.path(int_dir, "meta_DE_random_effects.tsv")
stopifnot(file.exists(meta_file))

# Defaults (ajustables por CLI)
fdr_thr      <- 0.05
beta_thr     <- 0.5     # solo para líneas verticales del volcano (no filtra)
top_each_dir <- 10      # top UP y top DOWN para lollipop
label_n      <- 12      # genes etiquetados en volcano
use_highconf_if_exists <- TRUE

# CLI: Rscript script.R [fdr_thr] [top_each_dir] [beta_thr] [label_n]
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
# Load + sanitize
# ============================================================
meta <- readr::read_tsv(meta_file, show_col_types = FALSE)

# Asegurar tipos numéricos en columnas clave (por si hay parsing raro)
num_cols <- intersect(
  c("beta_meta","se_meta","z_meta","p_meta","FDR_meta","I2","k",
    "beta_meta_DL","se_meta_DL","p_meta_DL","tau2_DL"),
  colnames(meta)
)
meta <- meta %>%
  mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.x))))

# FDR estable (si p=0 por underflow)
if (!("p_meta" %in% names(meta))) stop("No encuentro columna p_meta en meta.")
meta <- meta %>%
  mutate(
    p_meta2 = pmax(p_meta, 1e-300),
    FDR_meta2 = if ("FDR_meta" %in% names(meta)) {
      # recalcular por seguridad desde p_meta2 (y conservar FDR_meta original)
      p.adjust(p_meta2, method = "BH")
    } else {
      p.adjust(p_meta2, method = "BH")
    }
  )

if (!("gene_symbol" %in% names(meta))) stop("No encuentro gene_symbol.")
if (!("beta_meta" %in% names(meta))) stop("No encuentro beta_meta.")

# ============================================================
# (A) Volcano plot
# ============================================================
volc <- meta %>%
  mutate(
    neglog10FDR = -log10(FDR_meta2),
    sig = case_when(
      !is.na(FDR_meta2) & FDR_meta2 < fdr_thr & beta_meta >  0 ~ "UP",
      !is.na(FDR_meta2) & FDR_meta2 < fdr_thr & beta_meta <  0 ~ "DOWN",
      TRUE ~ "NS"
    ),
    neglog10FDR_cap = pmin(neglog10FDR, 100) # cap estético
  )

# genes a etiquetar: significativos, ordenar por FDR y magnitud
lab_df <- volc %>%
  filter(FDR_meta2 < fdr_thr, is.finite(beta_meta), is.finite(neglog10FDR)) %>%
  arrange(FDR_meta2, desc(abs(beta_meta))) %>%
  slice_head(n = label_n)

p_volcano <- ggplot(volc, aes(x = beta_meta, y = neglog10FDR_cap)) +
  geom_point(aes(color = sig, shape = sig), alpha = 0.7, size = 1.6) +
  scale_color_manual(values = pal_dir, name = "Status", drop = FALSE) +
  scale_shape_manual(values = c(UP = 15, DOWN = 16, NS = 17), name = "Status", drop = FALSE) +
  geom_vline(xintercept = c(-beta_thr, beta_thr), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(fdr_thr), linetype = "dashed", color = "grey50") +
  labs(
    title = "Meta-analysis (random effects): gene-level effects",
    subtitle = paste0(
      "Volcano: beta_meta vs -log10(FDR).\n",
      "Threshold: FDR<", fdr_thr, " | guide lines at |beta|=", beta_thr
    ),
    x = "beta_meta (log2FC, Resistant - Parental)",
    y = "-log10(FDR) (capped at 100)",
    shape = "Status",
    color = "Status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, lineheight = 1.1, margin = margin(b = 6)),
    legend.position = "right"
  )

# Etiquetado (ggrepel si existe, si no fallback)
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_volcano <- p_volcano +
    ggrepel::geom_text_repel(
      data = lab_df,
      aes(label = gene_symbol),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.3,
      point.padding = 0.2
    )
} else {
  p_volcano <- p_volcano +
    geom_text(
      data = lab_df,
      aes(label = gene_symbol),
      size = 2.8,
      vjust = -0.8
    )
}

out_vol_png <- file.path(fig_dir, "MetaVolcano_random_effects.png")
out_vol_pdf <- file.path(fig_dir, "MetaVolcano_random_effects.pdf")
ggsave(out_vol_png, p_volcano, width = 7.2, height = 5.8, dpi = 300)
ggsave(out_vol_pdf, p_volcano, width = 7.2, height = 5.8)
message("✔ Volcano: ", out_vol_png)
message("✔ Volcano: ", out_vol_pdf)

# ============================================================
# (B) Lollipop plot: top UP / top DOWN genes from meta
# Preferir meta_core_highconf si existe
# ============================================================
core_file <- file.path(int_dir, "meta_core_highconf.tsv")

df_for_lollipop <- meta
if (use_highconf_if_exists && file.exists(core_file)) {
  message("Usando meta_core_highconf.tsv para lollipop.")
  core <- readr::read_tsv(core_file, show_col_types = FALSE)
  if (!("gene_symbol" %in% names(core))) stop("meta_core_highconf.tsv no tiene gene_symbol.")
  df_for_lollipop <- meta %>% semi_join(core %>% select(gene_symbol), by = "gene_symbol")
} else {
  message("Usando meta_DE_random_effects.tsv completo para lollipop.")
}

# Selección "paper-ready": filtra por FDR, separa UP/DOWN, ordena y cap de tamaño
cap_val <- 100
top_n   <- top_each_dir

sig_df <- df_for_lollipop %>%
  filter(!is.na(FDR_meta2), FDR_meta2 < fdr_thr,
         !is.na(beta_meta), !is.na(gene_symbol)) %>%
  mutate(
    direction = if_else(beta_meta > 0, "UP", "DOWN"),
    mlog10 = -log10(FDR_meta2),
    mlog10_cap = pmin(mlog10, cap_val)
  )

top_up <- sig_df %>%
  filter(direction == "UP") %>%
  arrange(FDR_meta2, desc(abs(beta_meta))) %>%
  slice_head(n = top_n)

top_down <- sig_df %>%
  filter(direction == "DOWN") %>%
  arrange(FDR_meta2, desc(abs(beta_meta))) %>%
  slice_head(n = top_n)

# Orden: DOWN al fondo, UP arriba
plot_df <- bind_rows(top_down, top_up) %>%
  mutate(
    gene_plot = factor(gene_symbol, levels = c(top_down$gene_symbol, rev(top_up$gene_symbol)))
  )

if (nrow(plot_df) == 0) {
  stop("No hay genes significativos con el umbral FDR<", fdr_thr,
       ". Prueba fdr_thr=0.10 o revisa meta_DE_random_effects.tsv")
}

out_top_tsv <- file.path(int_dir, "meta_topgenes_for_lollipop.tsv")
readr::write_tsv(plot_df, out_top_tsv)
message("✔ Tabla top genes lollipop: ", out_top_tsv)

p_sep <- if (nrow(top_down) > 0 && nrow(top_up) > 0) nrow(top_down) + 0.5 else NA_real_

p_lolli <- ggplot(plot_df, aes(x = beta_meta, y = gene_plot, color = direction)) +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.4, alpha = 0.6) +
  geom_segment(aes(x = 0, xend = beta_meta, yend = gene_plot), linewidth = 0.7, alpha = 0.8) +
  geom_point(aes(size = mlog10_cap), alpha = 0.9) +
  scale_color_manual(values = pal_dir2, name = "Direction") +
  labs(
    title = "Meta-analysis: top significant genes (UP/DOWN)",
    subtitle = paste0(
      "Filter: FDR<", fdr_thr,
      " | Select: top ", top_n, " UP + top ", top_n,
      " DOWN\nOrder: FDR then |beta| | Size capped at ", cap_val
    ),
    x = "beta_meta (log2FC, Resistant - Parental)",
    y = NULL,
    size = "-log10(FDR) (capped)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, lineheight = 1.1, margin = margin(b = 6))
  )

if (!is.na(p_sep)) {
  p_lolli <- p_lolli + geom_hline(yintercept = p_sep, linetype = "dashed", color = "grey50")
}

out_lol_png <- file.path(fig_dir, "MetaLollipop_topGenes.png")
out_lol_pdf <- file.path(fig_dir, "MetaLollipop_topGenes.pdf")
ggsave(out_lol_png, p_lolli, width = 8.5, height = max(6, 0.22 * nrow(plot_df) + 2.5), dpi = 300)
ggsave(out_lol_pdf, p_lolli, width = 8.5, height = max(6, 0.22 * nrow(plot_df) + 2.5))
message("✔ Lollipop: ", out_lol_png)
message("✔ Lollipop: ", out_lol_pdf)

message(">> DONE: meta volcano + lollipop.")
