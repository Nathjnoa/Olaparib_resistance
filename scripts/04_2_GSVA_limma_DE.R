#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(patchwork)
})

# ============================================================
# Paths
# ============================================================
proj_dir <- Sys.getenv("OLAPARIB_RESISTANCE_DIR", unset = "")
if (!nzchar(proj_dir)) {
  proj_dir <- "~/bioinfo/projects/olaparib_resistance"
}
proj_dir <- path.expand(proj_dir)
if (!dir.exists(proj_dir)) {
  stop("No existe el directorio del proyecto: ", proj_dir)
}
gsva_tab_dir <- file.path(proj_dir, "results", "tables", "gsva")

meta_file    <- file.path(gsva_tab_dir, "GSVA_metadata.tsv")

out_tab_dir  <- file.path(proj_dir, "results", "tables", "gsva_DE")
out_fig_dir  <- file.path(proj_dir, "results", "figures", "gsva_DE")

dir.create(out_tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Helper: load GSVA matrix
# ============================================================
load_gsva <- function(path) {
  if (!file.exists(path)) stop("No existe: ", path)
  df <- readr::read_tsv(path, show_col_types = FALSE)
  stopifnot("Pathway" %in% colnames(df))
  mat <- as.matrix(df[, -1, drop = FALSE])
  rownames(mat) <- df$Pathway
  mat
}

# ============================================================
# Helper: run limma on GSVA matrix with explicit contrast
# ============================================================
run_limma_gsva <- function(mat, meta, dataset_label) {
  # meta must have: sample_id, condition (Parental/Resistant)
  meta <- meta %>%
    mutate(
      condition = factor(condition, levels = c("Parental", "Resistant"))
    ) %>%
    arrange(condition, sample_id)

  # align matrix columns
  missing_cols <- setdiff(meta$sample_id, colnames(mat))
  if (length(missing_cols) > 0) {
    stop(dataset_label, ": faltan columnas en GSVA matrix: ",
         paste(missing_cols, collapse = ", "))
  }
  mat <- mat[, meta$sample_id, drop = FALSE]

  # design and contrast
  design <- model.matrix(~ 0 + condition, data = meta)
  colnames(design) <- levels(meta$condition)

  fit <- lmFit(mat, design)
  contr <- makeContrasts(Res_vs_Par = Resistant - Parental, levels = design)
  fit2 <- contrasts.fit(fit, contr)
  fit2 <- eBayes(fit2)

  res <- topTable(fit2, coef = "Res_vs_Par", number = Inf, sort.by = "P") %>%
    rownames_to_column("Pathway") %>%
    as_tibble() %>%
    rename(
      logFC = logFC,
      P.Value = P.Value,
      adj.P.Val = adj.P.Val
    ) %>%
    mutate(
      Dataset = dataset_label,
      direction = case_when(
        !is.na(adj.P.Val) & adj.P.Val < 0.05 & logFC > 0  ~ "UP",
        !is.na(adj.P.Val) & adj.P.Val < 0.05 & logFC < 0  ~ "DOWN",
        TRUE ~ "NS"
      )
    )

  list(res = res, design = design)
}

# ============================================================
# Helpers: lollipop plots
# ============================================================
clean_pathway <- function(x) {
  x %>%
    str_replace_all("^HALLMARK_", "") %>%
    str_replace_all("_", " ") %>%
    str_to_title() %>%
    str_replace("^E2f$", "E2F") %>%
    str_replace("^Nfkb$", "NFKB") %>%
    str_replace("^Il6 Jak Stat3$", "IL6 JAK STAT3") %>%
    str_replace("^Tnfa", "TNFA") %>%
    str_replace("^Uv$", "UV")
}

select_paths <- function(res, dataset_label, fdr_cut = 0.10, top_up = 8, top_down = 8, rank_mode = "fdr") {
  required <- c("Pathway", "logFC", "adj.P.Val")
  miss <- setdiff(required, colnames(res))
  if (length(miss) > 0) {
    stop("Resultados de ", dataset_label, " carecen de columnas: ", paste(miss, collapse = ", "))
  }

  res_sig <- res %>% filter(!is.na(adj.P.Val), adj.P.Val < fdr_cut)
  res_up_all <- res_sig %>% filter(logFC > 0)
  res_down_all <- res_sig %>% filter(logFC < 0)

  if (rank_mode == "effect") {
    res_up <- res_up_all %>% arrange(desc(abs(logFC)), adj.P.Val) %>% slice_head(n = top_up)
    res_down <- res_down_all %>% arrange(desc(abs(logFC)), adj.P.Val) %>% slice_head(n = top_down)
  } else {
    res_up <- res_up_all %>% arrange(adj.P.Val) %>% slice_head(n = top_up)
    res_down <- res_down_all %>% arrange(adj.P.Val) %>% slice_head(n = top_down)
  }

  res_plot <- bind_rows(res_up, res_down) %>%
    mutate(
      direction = if_else(logFC > 0, "UP", "DOWN"),
      Pathway_clean = clean_pathway(Pathway),
      neglog10FDR = -log10(adj.P.Val),
      dataset = dataset_label
    )

  list(res_plot = res_plot, n_sig = nrow(res_sig), n_up = nrow(res_up), n_down = nrow(res_down))
}

plot_lollipop <- function(df, dataset_label, n_sig, n_up, n_down, rank_mode, x_limits) {
  pal_dir <- c("UP" = "#D81B60", "DOWN" = "#1E88E5")
  df <- df %>% mutate(Pathway = forcats::fct_reorder(Pathway_clean, logFC))
  subtitle_txt <- paste0("FDR<0.1: n=", n_sig, " | mostrado UP=", n_up, ", DOWN=", n_down, " | rank=", rank_mode)

  ggplot(df, aes(x = logFC, y = Pathway, color = direction)) +
    geom_segment(aes(x = 0, xend = logFC, yend = Pathway), linewidth = 0.7, alpha = 0.8) +
    geom_point(aes(size = -log10(adj.P.Val)), alpha = 0.9) +
    scale_color_manual(values = pal_dir, name = "Direction") +
    scale_size_continuous(name = "-log10(FDR)", range = c(2, 6)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    labs(
      title = dataset_label,
      subtitle = subtitle_txt,
      x = "logFC (GSVA score: Resistant - Parental)",
      y = NULL
    ) +
    coord_cartesian(xlim = c(-x_limits, x_limits))
}

# ============================================================
# Load metadata (generated in GSVA compute script)
# ============================================================
if (!file.exists(meta_file)) {
  stop("No encuentro GSVA_metadata.tsv en: ", meta_file,
       "\nAsegúrate de correr primero el script de GSVA (08_gsva_hallmarks_compute.R).")
}

meta_all <- readr::read_tsv(meta_file, show_col_types = FALSE)

required_cols <- c("sample_id", "dataset", "condition")
if (!all(required_cols %in% colnames(meta_all))) {
  stop("GSVA_metadata.tsv debe contener columnas: ",
       paste(required_cols, collapse = ", "))
}

# ============================================================
# Load GSVA matrices
# ============================================================
mat_A    <- load_gsva(file.path(gsva_tab_dir, "GSVA_A2780.tsv"))
mat_P    <- load_gsva(file.path(gsva_tab_dir, "GSVA_PEO1.tsv"))
mat_def  <- load_gsva(file.path(gsva_tab_dir, "GSVA_UWB_BRCAdef.tsv"))
mat_prof <- load_gsva(file.path(gsva_tab_dir, "GSVA_UWB_BRCAprof.tsv"))

# split metadata by dataset
meta_A    <- meta_all %>% filter(dataset == "A2780")       %>% select(sample_id, condition)
meta_P    <- meta_all %>% filter(dataset == "PEO1")        %>% select(sample_id, condition)
meta_def  <- meta_all %>% filter(dataset == "UWB_BRCAdef") %>% select(sample_id, condition)
meta_prof <- meta_all %>% filter(dataset == "UWB_BRCAprof")%>% select(sample_id, condition)

# ============================================================
# Run limma
# ============================================================
cat(">> Ejecutando limma sobre scores GSVA...\n")

outA    <- run_limma_gsva(mat_A,    meta_A,    "A2780")
outP    <- run_limma_gsva(mat_P,    meta_P,    "PEO1")
outDef  <- run_limma_gsva(mat_def,  meta_def,  "UWB_BRCAdef")
outProf <- run_limma_gsva(mat_prof, meta_prof, "UWB_BRCAprof")

res_A    <- outA$res
res_P    <- outP$res
res_def  <- outDef$res
res_prof <- outProf$res

# ============================================================
# Export results
# ============================================================
readr::write_tsv(res_A,    file.path(out_tab_dir, "GSVA_DE_A2780.tsv"))
readr::write_tsv(res_P,    file.path(out_tab_dir, "GSVA_DE_PEO1.tsv"))
readr::write_tsv(res_def,  file.path(out_tab_dir, "GSVA_DE_UWB_BRCAdef.tsv"))
readr::write_tsv(res_prof, file.path(out_tab_dir, "GSVA_DE_UWB_BRCAprof.tsv"))

# Integrated table
res_all <- bind_rows(res_A, res_P, res_def, res_prof)
readr::write_tsv(res_all, file.path(out_tab_dir, "GSVA_DE_ALL_datasets.tsv"))

# ============================================================
# Lollipop plots (FDR<0.1, top_up/down)
# ============================================================
fdr_cut_lolli <- 0.10
top_up_lolli   <- 8
top_down_lolli <- 8
rank_mode <- "fdr"  # opciones: "fdr" o "effect"

dataset_results <- list(
  list(label = "A2780",        res = res_A,    file_stub = "GSVA_DE_A2780"),
  list(label = "PEO1",         res = res_P,    file_stub = "GSVA_DE_PEO1"),
  list(label = "UWB BRCAdef",  res = res_def,  file_stub = "GSVA_DE_UWB_BRCAdef"),
  list(label = "UWB BRCAprof", res = res_prof, file_stub = "GSVA_DE_UWB_BRCAprof")
)

plot_entries <- list()
caption_lines <- character()
all_logfc <- c()

for (ds in dataset_results) {
  top_info <- select_paths(ds$res, ds$label,
                           fdr_cut = fdr_cut_lolli,
                           top_up = top_up_lolli,
                           top_down = top_down_lolli,
                           rank_mode = rank_mode)
  res_plot <- top_info$res_plot
  n_sig    <- top_info$n_sig
  if (n_sig == 0 || nrow(res_plot) == 0) {
    caption_lines <- c(caption_lines,
                       paste0(ds$label, ": 0 pathways (FDR<", fdr_cut_lolli, ")"))
    next
  }
  caption_lines <- c(caption_lines,
                     paste0(ds$label, ": FDR<", fdr_cut_lolli,
                            " n_sig=", n_sig,
                            " | mostrado UP=", top_info$n_up,
                            ", DOWN=", top_info$n_down,
                            " | rank=", rank_mode))
  all_logfc <- c(all_logfc, res_plot$logFC)
  plot_entries[[ds$label]] <- list(
    plot_df = res_plot,
    n_sig   = n_sig,
    n_up    = top_info$n_up,
    n_down  = top_info$n_down,
    file_stub = ds$file_stub
  )
}

x_lim <- if (length(all_logfc) > 0) max(1, max(abs(all_logfc), na.rm = TRUE)) else 1

plot_list <- list()
if (length(plot_entries) > 0) {
  for (nm in names(plot_entries)) {
    entry <- plot_entries[[nm]]
    p <- plot_lollipop(entry$plot_df, nm, entry$n_sig, entry$n_up, entry$n_down, rank_mode, x_limits = x_lim)
    plot_list[[nm]] <- p
    out_stub <- entry$file_stub
    ggsave(file.path(out_fig_dir, paste0(out_stub, "_lollipop_top10.png")), p, width = 7, height = 5, dpi = 300)
    ggsave(file.path(out_fig_dir, paste0(out_stub, "_lollipop_top10.pdf")), p, width = 7, height = 5)
  }
}

combined_caption <- paste(caption_lines, collapse = " | ")
if (length(plot_list) > 0) {
  combined <- wrap_plots(plot_list, ncol = 2) +
    plot_annotation(
      title = "GSVA differential pathways (Resistant - Parental): Top significant (FDR<0.1)",
      caption = combined_caption
    )

  ggsave(file.path(out_fig_dir, "GSVA_DE_lollipop_top10_all.png"), combined, width = 12, height = 10, dpi = 300)
  ggsave(file.path(out_fig_dir, "GSVA_DE_lollipop_top10_all.pdf"), combined, width = 12, height = 10)
} else {
  message("No hay pathways significativos en ningún dataset (FDR<", fdr_cut_lolli, ").")
}

cat(">> GSVA differential analysis COMPLETADO.\n")
cat("Tablas: ", out_tab_dir, "\n")
cat("Figuras:", out_fig_dir, "\n")
