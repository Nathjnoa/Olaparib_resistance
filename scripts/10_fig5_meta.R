#!/usr/bin/env Rscript

# ============================================================
# 10_fig5_meta.R
#
# Figure 5: Integrated transcriptomic signature from random-effects meta-analysis
#
# Layout (double_col, 200 × 220 mm):
#   Row 1 (57%): [A: Meta-volcano (45%) | B: Top genes lollipop — high-confidence (55%)]
#   Row 2 (43%): [C: Meta-GSEA Hallmarks (beta_meta ranking) — full width]
#
# Question: What is the integrated transcriptomic signature of olaparib resistance
#           across all 3 models, and what biological programs does it implicate?
#
# Prerequisites: 05_1_meta_analysis.R, 05_3_meta_gsea.R
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
  library(grid)
})

source(file.path({
  f <- commandArgs(trailingOnly = FALSE)
  f <- grep("--file=", f, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("--file=", "", f[1])))
  else path.expand("~/bioinfo/projects/olaparib_resistance/scripts")
}, "00_config.R"))

log_handle <- start_log("10_fig5_meta")
set.seed(42)

# ============================================================
# 0. Validate inputs
# ============================================================
int_dir      <- file.path(tables_dir, "integrated")
meta_gsea_dir <- file.path(tables_dir, "meta_gsea")

req <- c(
  meta_results  = file.path(int_dir,       "meta_DE_random_effects.tsv"),
  top_genes     = file.path(int_dir,       "meta_topgenes_for_lollipop.tsv"),
  highconf      = file.path(int_dir,       "meta_core_highconf.tsv"),
  gsea_hallmarks = file.path(meta_gsea_dir, "metaGSEA_HALLMARKS_BETA.tsv")
)
missing_req <- req[!file.exists(req)]
if (length(missing_req) > 0) {
  stop("Archivos faltantes:\n",
       paste0("  ", names(missing_req), ": ", missing_req, collapse = "\n"))
}

fig5_dir <- file.path(figures_dir, "fig5")
dir.create(fig5_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Style constants — double_col preset
# ============================================================
FIG_W_MM  <- 200
FIG_H_MM  <- 220
BASE_SIZE <- 8.5

PAL_DIR <- c("UP" = "#D81B60", "DOWN" = "#1E88E5", "NS" = "#9E9E9E")

theme_pub <- function() {
  theme_bw(base_size = BASE_SIZE) +
    theme(
      plot.title        = element_text(size = 9, face = "bold", hjust = 0.5,
                                       margin = margin(b = 2)),
      plot.subtitle     = element_text(size = 7, hjust = 0.5, color = "grey40",
                                       margin = margin(b = 3)),
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

# Abbreviate long Hallmark names
abbrev_hallmark <- function(x) {
  x %>%
    str_replace_all("Epithelial Mesenchymal Transition", "EMT") %>%
    str_replace_all("Interferon Gamma Response",         "IFN-\u03b3 Response") %>%
    str_replace_all("Interferon Alpha Response",         "IFN-\u03b1 Response") %>%
    str_replace_all("TNFA Signaling Via NFKB",           "TNFA via NFKB") %>%
    str_replace_all("Il2 Stat5 Signaling",               "IL2-STAT5 Signaling") %>%
    str_replace_all("Kras Signaling Dn",                 "KRAS Signaling DN") %>%
    str_replace_all("Kras Signaling Up",                 "KRAS Signaling UP") %>%
    str_replace_all("E2f Targets",                       "E2F Targets") %>%
    str_replace_all("G2m Checkpoint",                    "G2M Checkpoint") %>%
    str_replace_all("Mtorc1 Signaling",                  "MTORC1 Signaling") %>%
    str_replace_all("Myc Targets V1",                    "MYC Targets V1") %>%
    str_replace_all("Myc Targets V2",                    "MYC Targets V2") %>%
    str_replace_all("Uv Response Dn",                    "UV Response DN") %>%
    str_replace_all("Uv Response Up",                    "UV Response UP")
}

clean_gsea_path <- function(x) {
  abbrev_hallmark(
    str_to_title(str_replace_all(str_replace_all(x, "^HALLMARK_", ""), "_", " "))
  )
}

# ============================================================
# 2. Load data
# ============================================================
meta <- read_tsv(req["meta_results"], show_col_types = FALSE) %>%
  mutate(across(c(beta_meta, se_meta, z_meta, p_meta, FDR_meta, I2, k),
                ~ suppressWarnings(as.numeric(.x)))) %>%
  mutate(
    p_meta2    = pmax(p_meta, 1e-300),
    FDR_meta2  = p.adjust(p_meta2, method = "BH"),
    mlog10_fdr = pmin(-log10(FDR_meta2), 100),
    status     = case_when(
      !is.na(FDR_meta2) & FDR_meta2 < FDR_CUTOFF & beta_meta >  0 ~ "UP",
      !is.na(FDR_meta2) & FDR_meta2 < FDR_CUTOFF & beta_meta <  0 ~ "DOWN",
      TRUE ~ "NS"
    )
  )

top_genes <- read_tsv(req["top_genes"], show_col_types = FALSE)

highconf <- read_tsv(req["highconf"], show_col_types = FALSE)

gsea_h <- read_tsv(req["gsea_hallmarks"], show_col_types = FALSE) %>%
  mutate(
    pathway_clean = clean_gsea_path(pathway),
    direction     = if_else(NES > 0, "UP", "DOWN"),
    mlog10_fdr    = -log10(padj)
  ) %>%
  filter(!is.na(padj), padj < FDR_CUTOFF) %>%
  arrange(NES) %>%
  mutate(pathway_clean = factor(pathway_clean, levels = pathway_clean))

stopifnot("No significant meta-GSEA Hallmarks (beta)" = nrow(gsea_h) > 0)
message("Meta-GSEA Hallmarks (beta, FDR<0.05): n=", nrow(gsea_h))

# ============================================================
# 3. Panel A — Meta-volcano
# ============================================================

# Label genes: the same genes in the lollipop (Panel B) to link both panels
lollipop_genes <- top_genes$gene_symbol
lab_df <- meta %>%
  filter(gene_symbol %in% lollipop_genes,
         is.finite(beta_meta), is.finite(mlog10_fdr))

# Count sig genes for subtitle
n_up   <- sum(meta$status == "UP",   na.rm = TRUE)
n_down <- sum(meta$status == "DOWN", na.rm = TRUE)

p_vol <- ggplot(meta, aes(x = beta_meta, y = mlog10_fdr)) +
  geom_point(aes(color = status, shape = status),
             alpha = 0.55, size = 1.2, stroke = 0) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_text_repel(
    data          = lab_df,
    aes(label     = gene_symbol, color = status),
    size          = 2.2,
    max.overlaps  = Inf,
    box.padding   = 0.25,
    point.padding = 0.15,
    segment.size  = 0.25,
    segment.color = "grey60",
    show.legend   = FALSE,
    min.segment.length = 0.1
  ) +
  scale_color_manual(values = PAL_DIR, name = "Status",
                     guide  = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  scale_shape_manual(values = c("UP" = 15, "DOWN" = 16, "NS" = 17),
                     name   = "Status") +
  scale_x_continuous(breaks = seq(-5, 5, 1)) +
  scale_y_continuous(limits = c(0, 102), breaks = c(0, 25, 50, 75, 100)) +
  labs(
    x = expression(beta[meta]~(log[2]*FC*","~Resistant-Parental)),
    y = expression(-log[10](FDR)~~(capped~at~100))
  ) +
  theme_pub() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey92"),
    legend.position    = c(0.84, 0.82),
    legend.box.background = element_rect(color = "grey80", linewidth = 0.3)
  )

# ============================================================
# 4. Panel B — Top genes lollipop (high-confidence)
# ============================================================

# Build plot dataframe: sort DOWN then UP by beta_meta for layout
lollipop_df <- top_genes %>%
  arrange(beta_meta) %>%
  mutate(gene_symbol = factor(gene_symbol, levels = gene_symbol),
         mlog10_cap  = pmin(-log10(FDR_meta2), 100))

# Add a divider row between DOWN and UP
divider_y <- sum(lollipop_df$direction == "DOWN") + 0.5

p_lollipop <- ggplot(lollipop_df,
                     aes(x = beta_meta, y = gene_symbol, color = direction)) +
  geom_segment(aes(x = 0, xend = beta_meta,
                   y = gene_symbol, yend = gene_symbol),
               linewidth = 0.55) +
  geom_point(aes(size = mlog10_cap)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey40") +
  geom_hline(yintercept = divider_y, linetype = "dashed",
             color = "grey55", linewidth = 0.4) +
  scale_color_manual(
    values = c("UP" = "#D81B60", "DOWN" = "#1E88E5"),
    name   = "Direction",
    guide  = guide_legend(order = 1, override.aes = list(size = 2.5))
  ) +
  scale_size_continuous(
    name   = expression(-log[10](FDR)),
    range  = c(1.5, 5),
    breaks = c(40, 70, 100),
    labels = c("40", "70", "\u2265100"),
    guide  = guide_legend(order = 2,
                          override.aes = list(color = "grey40"))
  ) +
  labs(
    x = expression(beta[meta]~(log[2]*FC)),
    y = NULL
  ) +
  theme_pub() +
  theme(
    axis.text.y      = element_text(size = 7, face = "italic"),
    panel.grid.major.x = element_line(color = "grey92"),
    legend.position  = "right"
  )

# ============================================================
# 5. Panel C — Meta-GSEA Hallmarks (beta_meta ranking)
# ============================================================

nes_lim <- max(abs(gsea_h$NES), na.rm = TRUE)
nes_lim <- ceiling(nes_lim * 5) / 5 + 0.2

p_gsea <- ggplot(gsea_h,
                 aes(x = NES, y = pathway_clean, color = direction)) +
  geom_segment(aes(x = 0, xend = NES,
                   y = pathway_clean, yend = pathway_clean),
               linewidth = 0.6) +
  geom_point(aes(size = mlog10_fdr)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey40") +
  scale_color_manual(
    values = c("UP" = "#D81B60", "DOWN" = "#1E88E5"),
    name   = "Direction in\nresistant cells",
    guide  = guide_legend(order = 1, override.aes = list(size = 3))
  ) +
  scale_size_continuous(
    name   = expression(-log[10](FDR)),
    range  = c(2, 5),
    limits = c(1.0, 2.5),
    breaks = c(1.1, 1.3, 1.5),
    labels = c("1.1", "1.3", "1.5"),
    guide  = guide_legend(order = 2, override.aes = list(color = "grey40"))
  ) +
  scale_x_continuous(limits = c(-nes_lim, nes_lim),
                     breaks = seq(-2, 2, 0.5)) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  theme_pub() +
  theme(
    axis.text.y       = element_text(size = 7.5),
    panel.grid.major.x = element_line(color = "grey92"),
    legend.position   = "right"
  )

# ============================================================
# 6. Compose and export
# ============================================================
out_pdf <- file.path(fig5_dir, "Fig5.pdf")
out_png <- file.path(fig5_dir, "Fig5.png")

FIG_W_IN <- FIG_W_MM / 25.4
FIG_H_IN <- FIG_H_MM / 25.4

ROW1_FRAC <- 0.57   # fraction of height for row 1
ROW2_FRAC <- 1 - ROW1_FRAC

draw_fig5 <- function() {
  grid.newpage()
  lt <- grid.layout(
    nrow    = 2,
    ncol    = 1,
    heights = unit(c(FIG_H_MM * ROW1_FRAC, FIG_H_MM * ROW2_FRAC), "mm")
  )
  pushViewport(viewport(layout = lt))

  # Row 1: volcano + lollipop side by side
  pushViewport(viewport(layout.pos.row = 1))
  row1 <- (p_vol | p_lollipop) +
    plot_layout(widths = c(0.54, 0.46)) +
    plot_annotation(tag_levels = list(c("A", "B"))) &
    theme(plot.tag = element_text(size = 9, face = "bold"))
  print(row1, newpage = FALSE)
  popViewport()

  # Row 2: meta-GSEA full width
  pushViewport(viewport(layout.pos.row = 2))
  row2 <- p_gsea +
    plot_annotation(tag_levels = list(c("C"))) &
    theme(plot.tag = element_text(size = 9, face = "bold"))
  print(row2, newpage = FALSE)
  popViewport()

  popViewport()
}

pdf(out_pdf, width = FIG_W_IN, height = FIG_H_IN)
draw_fig5()
invisible(dev.off())

png(out_png, width = FIG_W_MM, height = FIG_H_MM, units = "mm", res = 600)
draw_fig5()
invisible(dev.off())

message("Figure 5 saved:")
message("  PDF: ", out_pdf)
message("  PNG: ", out_png)

stop_log(log_handle)
