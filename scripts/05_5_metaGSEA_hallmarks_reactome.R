#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(msigdbr)
  library(fgsea)
})

# ============================================================
# Config paths
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

rank_z_file    <- file.path(int_dir, "meta_DE_meta_rank_z.tsv")
rank_beta_file <- file.path(int_dir, "meta_DE_meta_rank_beta.tsv")
meta_file      <- file.path(int_dir, "meta_DE_random_effects.tsv")

out_tab_dir <- file.path(proj_dir, "results", "tables", "meta_gsea")
out_fig_dir <- file.path(proj_dir, "results", "figures", "meta_gsea")
dir.create(out_tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Parameters (CLI)
# Uso:
#   Rscript scripts/11_metaGSEA_hallmarks_reactome.R [stat=z|beta] [padj_thr] [top_n_plot] [make_plots=0|1]
# ============================================================
stat_choice <- "z"
padj_thr    <- 0.10
top_n_plot  <- 10
make_plots  <- 1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) stat_choice <- tolower(args[1])
if (length(args) >= 2) padj_thr    <- as.numeric(args[2])
if (length(args) >= 3) top_n_plot  <- as.integer(args[3])
if (length(args) >= 4) make_plots  <- as.integer(args[4])

if (!stat_choice %in% c("z", "beta")) {
  stop("stat_choice debe ser 'z' o 'beta'. Recibido: ", stat_choice)
}

message("Params: stat=", stat_choice, " | padj_thr=", padj_thr,
        " | top_n_plot=", top_n_plot, " | make_plots=", make_plots)

# ============================================================
# Helpers
# ============================================================
make_stats <- function(df, stat_col, gene_col = "gene_symbol") {
  if (!gene_col %in% names(df)) stop("No encuentro columna ", gene_col)
  if (!stat_col %in% names(df)) stop("No encuentro columna ", stat_col)

  tmp <- df %>%
    transmute(
      gene = as.character(.data[[gene_col]]),
      stat = suppressWarnings(as.numeric(.data[[stat_col]]))
    ) %>%
    filter(!is.na(gene), gene != "", !is.na(stat), is.finite(stat))

  # Resolver duplicados por gen: quedarnos con el mayor |stat|
  tmp2 <- tmp %>%
    group_by(gene) %>%
    summarise(stat = stat[which.max(abs(stat))], .groups = "drop")

  stats <- tmp2$stat
  names(stats) <- tmp2$gene

  # Jitter determinístico minúsculo para evitar empates exactos
  stats <- stats + seq_along(stats) * 1e-12
  sort(stats, decreasing = TRUE)
}

get_pathways_msigdbr <- function(collection, subcollection = NULL) {
  msig <- msigdbr(
    species = "Homo sapiens",
    collection = collection,
    subcollection = subcollection
  )

  pathways <- msig %>%
    select(gs_name, gene_symbol) %>%
    mutate(gene_symbol = as.character(gene_symbol)) %>%
    filter(!is.na(gene_symbol), gene_symbol != "") %>%
    split(x = .$gene_symbol, f = .$gs_name)

  lapply(pathways, unique)
}

run_fgsea <- function(pathways, stats, minSize = 15, maxSize = 500) {
  fgseaMultilevel(
    pathways = pathways,
    stats    = stats,
    minSize  = minSize,
    maxSize  = maxSize
  ) %>%
    as_tibble() %>%
    arrange(padj, desc(abs(NES)))
}

# ---- Etiquetas "bonitas" ----
clean_pathway_label <- function(x, wrap_width = 40, truncate_chars = 120) {
  lab <- x %>%
    # Hallmarks / Reactome prefijos
    str_replace("^HALLMARK_", "") %>%
    str_replace("^REACTOME_", "") %>%
    str_replace("^REACTOME\\s+", "") %>%
    str_replace("^Reactome\\s+", "") %>%
    # separadores
    str_replace_all("_", " ") %>%
    str_squish() %>%
    str_to_title()

  # acrónimos comunes
  lab <- lab %>%
    str_replace_all("\\bE2f\\b", "E2F") %>%
    str_replace_all("\\bNfkb\\b", "NFKB") %>%
    str_replace_all("\\bIl\\s*\\b", "IL") %>%
    str_replace_all("\\bJak\\b", "JAK") %>%
    str_replace_all("\\bStat\\b", "STAT") %>%
    str_replace_all("\\bDna\\b", "DNA") %>%
    str_replace_all("\\bRna\\b", "RNA") %>%
    str_replace_all("\\bTnf\\b", "TNF") %>%
    str_replace_all("\\bUv\\b", "UV")

  # truncado suave si algún nombre es ridículamente largo
  if (!is.null(truncate_chars) && is.finite(truncate_chars)) {
    lab <- ifelse(nchar(lab) > truncate_chars,
                  paste0(substr(lab, 1, truncate_chars - 1), "…"),
                  lab)
  }

  str_wrap(lab, width = wrap_width)
}

# ---- Plot lollipop: top UP + top DOWN entre padj<thr ----
plot_lollipop_pathways <- function(res,
                                   title,
                                   out_png,
                                   out_pdf,
                                   top_n_each = 10,
                                   padj_thr = 0.10,
                                   wrap_width = 40,
                                   truncate_chars = 120,
                                   width = 10,
                                   height = 6,
                                   mlog10_cap_max = 12) {

  df <- res %>%
    filter(!is.na(padj), !is.na(NES)) %>%
    mutate(
      direction = if_else(NES > 0, "UP", "DOWN"),
      mlog10 = -log10(pmax(padj, 1e-300)),
      mlog10_cap = pmin(mlog10, mlog10_cap_max)
    )

  df_sig <- df %>% filter(padj < padj_thr)

  n_up_sig   <- sum(df_sig$direction == "UP")
  n_down_sig <- sum(df_sig$direction == "DOWN")

  # helper: top por padj y luego |NES|
  pick_top <- function(dat, n_show) {
    dat %>%
      arrange(padj, desc(abs(NES))) %>%
      slice_head(n = n_show)
  }

  top_up   <- df_sig %>% filter(direction == "UP")   %>% pick_top(top_n_each)
  top_down <- df_sig %>% filter(direction == "DOWN") %>% pick_top(top_n_each)

  plot_df <- bind_rows(top_down, top_up) %>%
    distinct(pathway, .keep_all = TRUE) %>%
    mutate(pathway_label = clean_pathway_label(pathway,
                                               wrap_width = wrap_width,
                                               truncate_chars = truncate_chars))

  if (nrow(plot_df) == 0) {
    message("Sin términos para plotear en ", title, " (padj<", padj_thr, ")")
    return(invisible(NULL))
  }

  # Orden en Y: DOWN arriba (más negativo) y UP abajo (más positivo), manteniendo selección
  levels_raw <- c(rev(top_down$pathway), rev(top_up$pathway))
  if (length(levels_raw) == 0) levels_raw <- plot_df$pathway

  levels_clean <- clean_pathway_label(levels_raw,
                                      wrap_width = wrap_width,
                                      truncate_chars = truncate_chars)

  plot_df <- plot_df %>%
    mutate(pathway_label = factor(pathway_label, levels = levels_clean))

  xlim <- max(abs(plot_df$NES), na.rm = TRUE)
  if (!is.finite(xlim) || xlim <= 0) xlim <- 1
  xlim <- xlim * 1.12

  subtitle_txt <- paste0(
    "Top ", top_n_each, " UP + Top ", top_n_each, " DOWN (padj<", padj_thr, "). ",
    "Disponibles: UP=", n_up_sig, ", DOWN=", n_down_sig
  )

  p <- ggplot(plot_df, aes(x = NES, y = pathway_label, color = direction)) +
    geom_vline(xintercept = 0, linewidth = 0.6, color = "grey55") +
    geom_segment(aes(x = 0, xend = NES, yend = pathway_label),
                 linewidth = 1.0, color = "grey25") +
    geom_point(aes(size = mlog10_cap), alpha = 0.95) +
    coord_cartesian(xlim = c(-xlim, xlim)) +
    scale_color_manual(values = c("UP" = "#D81B60", "DOWN" = "#1E88E5"),
                       name = "Direction") +
    scale_size_continuous(name = expression(-log[10]*"(FDR)"),
                          range = c(2.2, 7)) +
    labs(title = title, subtitle = subtitle_txt, x = "NES", y = NULL) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
      plot.subtitle = element_text(size = 11.5, hjust = 0.5),
      axis.title.x = element_text(size = 14),
      axis.text.y = element_text(size = 10.5),
      axis.text.x = element_text(size = 11.5),
      legend.title = element_text(size = 12.5),
      legend.text = element_text(size = 11.5),
      legend.position = "right",
      legend.box = "vertical",
      plot.margin = margin(10, 18, 10, 12)
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(size = 4)),
      size  = guide_legend(order = 2)
    )

  ggsave(out_png, p, width = width, height = height, dpi = 300)
  ggsave(out_pdf, p, width = width, height = height)
  invisible(p)
}

# ============================================================
# Load ranking (prefer export)
# ============================================================
rank_df <- NULL
stat_col <- NULL

if (stat_choice == "z" && file.exists(rank_z_file)) {
  message("Usando ranking Z: ", rank_z_file)
  rank_df <- readr::read_tsv(rank_z_file, show_col_types = FALSE)

  if ("z_meta" %in% names(rank_df)) stat_col <- "z_meta"
  else if ("z" %in% names(rank_df)) stat_col <- "z"
  else if ("stat" %in% names(rank_df)) stat_col <- "stat"
  else if (all(c("beta_meta", "se_meta") %in% names(rank_df))) {
    rank_df <- rank_df %>% mutate(z_meta = beta_meta / se_meta)
    stat_col <- "z_meta"
  } else stop("No encuentro columna z en ", rank_z_file)

} else if (stat_choice == "beta" && file.exists(rank_beta_file)) {
  message("Usando ranking BETA: ", rank_beta_file)
  rank_df <- readr::read_tsv(rank_beta_file, show_col_types = FALSE)

  if ("beta_meta" %in% names(rank_df)) stat_col <- "beta_meta"
  else if ("beta" %in% names(rank_df)) stat_col <- "beta"
  else stop("No encuentro columna beta en ", rank_beta_file)

} else {
  message("No encontré ranking exportado para stat_choice=", stat_choice, ". Usando meta_DE_random_effects.tsv")
  stopifnot(file.exists(meta_file))
  rank_df <- readr::read_tsv(meta_file, show_col_types = FALSE)

  if (stat_choice == "z") {
    if ("z_meta" %in% names(rank_df)) {
      stat_col <- "z_meta"
    } else if (all(c("beta_meta","se_meta") %in% names(rank_df))) {
      rank_df <- rank_df %>% mutate(z_meta = beta_meta / se_meta)
      stat_col <- "z_meta"
    } else stop("No puedo construir z_meta (faltan z_meta o beta_meta+se_meta).")
  } else {
    if (!"beta_meta" %in% names(rank_df)) stop("No encuentro beta_meta en meta_DE_random_effects.tsv")
    stat_col <- "beta_meta"
  }
}

stats <- make_stats(rank_df, stat_col = stat_col, gene_col = "gene_symbol")
message("Stats vector: n_genes=", length(stats),
        " | stat_col=", stat_col,
        " | min=", signif(min(stats), 3),
        " | max=", signif(max(stats), 3))

# ============================================================
# Build pathways (Hallmarks + Reactome)
# ============================================================
message("Descargando MSigDB (msigdbr) para Hallmarks y Reactome...")
pathways_H <- get_pathways_msigdbr(collection = "H")
pathways_R <- get_pathways_msigdbr(collection = "C2", subcollection = "CP:REACTOME")

message("Pathways Hallmarks: ", length(pathways_H))
message("Pathways Reactome:  ", length(pathways_R))

# ============================================================
# Run meta-GSEA
# ============================================================
set.seed(1)

message(">> Corriendo fgseaMultilevel: Hallmarks")
res_H <- run_fgsea(pathways_H, stats, minSize = 15, maxSize = 500) %>%
  mutate(collection = "HALLMARKS")

message(">> Corriendo fgseaMultilevel: Reactome")
res_R <- run_fgsea(pathways_R, stats, minSize = 15, maxSize = 500) %>%
  mutate(collection = "REACTOME")

# ============================================================
# Export tables
# ============================================================
out_H <- file.path(out_tab_dir, paste0("metaGSEA_HALLMARKS_", toupper(stat_choice), ".tsv"))
out_R <- file.path(out_tab_dir, paste0("metaGSEA_REACTOME_",  toupper(stat_choice), ".tsv"))

readr::write_tsv(res_H, out_H)
readr::write_tsv(res_R, out_R)

message("✔ Hallmarks TSV: ", out_H)
message("✔ Reactome TSV:  ", out_R)

res_ALL <- bind_rows(res_H, res_R) %>%
  select(collection, pathway, NES, pval, padj, size, leadingEdge)
out_ALL <- file.path(out_tab_dir, paste0("metaGSEA_ALL_", toupper(stat_choice), ".tsv"))
readr::write_tsv(res_ALL, out_ALL)
message("✔ Combined TSV:  ", out_ALL)

# ============================================================
# Optional: lollipop plots
# ============================================================
if (make_plots == 1) {
  message(">> Generando lollipop plots (top pathways)")

  plot_lollipop_pathways(
    res_H,
    title   = paste0("Meta-GSEA Hallmarks (", stat_col, ")"),
    out_png = file.path(out_fig_dir, paste0("meta_gsea_hallmarks_", stat_choice, "_lollipop.png")),
    out_pdf = file.path(out_fig_dir, paste0("meta_gsea_hallmarks_", stat_choice, "_lollipop.pdf")),
    top_n_each = top_n_plot,
    padj_thr = padj_thr,
    wrap_width = 36,
    truncate_chars = 110,
    width = 10,
    height = 6,
    mlog10_cap_max = 12
  )

  # Reactome: más ancho + wrap más agresivo para que no se “coma” el panel
  plot_lollipop_pathways(
    res_R,
    title   = paste0("Meta-GSEA Reactome (", stat_col, ")"),
    out_png = file.path(out_fig_dir, paste0("meta_gsea_reactome_", stat_choice, "_lollipop.png")),
    out_pdf = file.path(out_fig_dir, paste0("meta_gsea_reactome_", stat_choice, "_lollipop.pdf")),
    top_n_each = top_n_plot,
    padj_thr = padj_thr,
    wrap_width = 30,
    truncate_chars = 120,
    width = 13,
    height = 7,
    mlog10_cap_max = 20
  )

  message("✔ Lollipops guardados en: ", out_fig_dir)
}

message(">> Meta-GSEA COMPLETADO.")
