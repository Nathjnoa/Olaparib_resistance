#!/usr/bin/env Rscript

# ============================================================
# 04_1_gsva_analysis.R
#
# GSVA scoring + limma differential analysis on pathway scores:
#
#   Section 1: Compute GSVA scores (Hallmarks) per dataset
#   Section 2: limma on GSVA scores (Resistant vs Parental)
#   Section 3: Per-dataset lollipop plots + combined figure
#
# Datasets: A2780 (FPKM), PEO1 (TPM), UWB BRCA-def (VST counts)
# kcdf = "Gaussian" (continuous expression data)
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(GSVA)
  library(msigdbr)
  library(DESeq2)
  library(limma)
  library(patchwork)
  library(R.utils)
})

source(file.path({f <- commandArgs(trailingOnly=FALSE); f <- grep("--file=",f,value=TRUE); if(length(f)) dirname(normalizePath(sub("--file=","",f[1]))) else path.expand("~/bioinfo/projects/olaparib_resistance/scripts")}, "00_config.R"))

log_handle <- start_log("04_1_gsva_analysis")

gsva_tab_dir <- file.path(tables_dir, "gsva")
gsva_DE_dir  <- file.path(tables_dir, "gsva_DE")
gsva_fig_dir <- file.path(figures_dir, "gsva_DE")

dir.create(gsva_tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(gsva_DE_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(gsva_fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Helpers
# ============================================================
load_expr_matrix <- function(path, gene_col = 1) {
  if (!file.exists(path)) stop("No existe: ", path)
  df   <- readr::read_tsv(path, show_col_types = FALSE)
  gcol <- if (is.numeric(gene_col)) colnames(df)[gene_col] else gene_col
  if (!(gcol %in% colnames(df))) stop("No se encontró columna gen en: ", path)
  colnames(df)[colnames(df) == gcol] <- "Gene"
  df  <- df %>% dplyr::distinct(Gene, .keep_all = TRUE)
  mat <- as.matrix(df[, -1, drop = FALSE])
  rownames(mat) <- df$Gene
  mat
}

safe_gunzip <- function(gz_path, out_path) {
  if (!file.exists(out_path) && file.exists(gz_path)) {
    message(">> Descomprimiendo: ", gz_path)
    R.utils::gunzip(gz_path, overwrite = TRUE, remove = FALSE)
  }
  if (!file.exists(out_path)) stop("No existe: ", out_path)
}

clean_for_gsva <- function(mat, dataset_name) {
  keep_finite <- apply(mat, 1, function(x) all(is.finite(x)))
  keep_var    <- apply(mat[keep_finite, , drop = FALSE], 1, sd, na.rm = TRUE) > 0
  n_drop_finite <- sum(!keep_finite)
  n_drop_var    <- sum(keep_finite) - sum(keep_var)
  if (n_drop_finite > 0) message(dataset_name, ": excluidos NA/Inf = ", n_drop_finite)
  if (n_drop_var > 0)    message(dataset_name, ": excluidos var=0 = ", n_drop_var)
  mat2 <- mat[keep_finite, , drop = FALSE]
  mat2[keep_var, , drop = FALSE]
}

run_gsva_hallmarks <- function(mat, pathways, dataset_name) {
  mat2  <- clean_for_gsva(mat, dataset_name)
  param <- gsvaParam(exprData = mat2, geneSets = pathways, kcdf = "Gaussian")
  gsva(param)
}

# ============================================================
# SECTION 1: Compute GSVA scores
# ============================================================

# A2780 (log2 FPKM)
file_A_proc <- file.path(processed_dir, "GSE153867_A2780_FPKM.tsv")
mat_A_raw   <- load_expr_matrix(file_A_proc, gene_col = "gene_id")
mat_A       <- log2(mat_A_raw + 1)

col_par_A <- grep("^O-", colnames(mat_A), value = TRUE)
col_res_A <- grep("^C-", colnames(mat_A), value = TRUE)
if (length(col_par_A) == 0 || length(col_res_A) == 0) {
  stop("A2780: no detecté columnas O-* / C-*")
}
mat_A <- mat_A[, c(col_par_A, col_res_A), drop = FALSE]
message("A2780 dims: ", nrow(mat_A), " x ", ncol(mat_A))

meta_A <- tibble(
  sample_id = colnames(mat_A),
  dataset   = "A2780",
  condition = factor(ifelse(grepl("^O-", colnames(mat_A)), "Parental", "Resistant"),
                     levels = c("Parental", "Resistant"))
)

# PEO1 (log2 TPM)
file_P_gz  <- file.path(raw_dir, "GSE117765", "GSE117765_matrix.txt.gz")
file_P_txt <- file.path(raw_dir, "GSE117765", "GSE117765_matrix.txt")
safe_gunzip(file_P_gz, file_P_txt)

mat_P_raw <- load_expr_matrix(file_P_txt, gene_col = 1)
mat_P     <- log2(mat_P_raw + 1)

col_par_P <- grep("Adherent", colnames(mat_P), value = TRUE)
col_res_P <- grep("Clone",    colnames(mat_P), value = TRUE)
if (length(col_par_P) == 0 || length(col_res_P) == 0) {
  stop("PEO1: no detecté columnas Adherent / Clone")
}
mat_P <- mat_P[, c(col_par_P, col_res_P), drop = FALSE]
message("PEO1 dims: ", nrow(mat_P), " x ", ncol(mat_P))

meta_P <- tibble(
  sample_id = colnames(mat_P),
  dataset   = "PEO1",
  condition = factor(ifelse(grepl("Adherent", colnames(mat_P), ignore.case = TRUE),
                            "Parental", "Resistant"),
                     levels = c("Parental", "Resistant"))
)

# UWB BRCA-def (VST counts)
file_U_gz  <- file.path(raw_dir, "GSE235980", "GSE235980_CountReads.txt.gz")
file_U_txt <- file.path(raw_dir, "GSE235980", "GSE235980_CountReads.txt")
safe_gunzip(file_U_gz, file_U_txt)

uwb_raw <- readr::read_tsv(file_U_txt, show_col_types = FALSE)
colnames(uwb_raw)[1] <- "Gene"
uwb_raw  <- uwb_raw %>% dplyr::distinct(Gene, .keep_all = TRUE)
count_mat <- as.matrix(uwb_raw[, -1, drop = FALSE])
rownames(count_mat) <- uwb_raw$Gene

col_def <- c(UWB_DEF_PAR, UWB_DEF_RES)
missing <- setdiff(col_def, colnames(count_mat))
if (length(missing) > 0) stop("UWB: faltan columnas: ", paste(missing, collapse = ", "))

meta_def <- tibble(
  sample_id = col_def,
  dataset   = "UWB_BRCAdef",
  condition = factor(c(rep("Parental", length(UWB_DEF_PAR)), rep("Resistant", length(UWB_DEF_RES))),
                     levels = c("Parental", "Resistant"))
)

make_vst <- function(counts_sub, meta_sub) {
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_sub),
    colData   = as.data.frame(meta_sub %>% tibble::column_to_rownames("sample_id")),
    design    = ~ 1
  )
  vsd <- vst(dds, blind = TRUE)
  assay(vsd)
}

vst_def <- make_vst(count_mat[, col_def, drop = FALSE], meta_def)
vst_def <- vst_def[, meta_def %>% arrange(condition) %>% pull(sample_id), drop = FALSE]
message("UWB BRCAdef dims: ", nrow(vst_def), " x ", ncol(vst_def))

# Common genes (for consistent meta-analysis downstream)
genes_common <- Reduce(intersect, list(rownames(mat_A), rownames(mat_P), rownames(vst_def)))
message("Genes comunes: ", length(genes_common))

mat_A   <- mat_A[genes_common, , drop = FALSE]
mat_P   <- mat_P[genes_common, , drop = FALSE]
vst_def <- vst_def[genes_common, , drop = FALSE]

# Hallmarks gene sets
msig      <- msigdbr(species = "Homo sapiens", collection = "H")
df_h      <- msig %>% dplyr::select(gs_name, gene_symbol) %>% dplyr::distinct()
pathways  <- split(df_h$gene_symbol, df_h$gs_name)
message("Hallmarks: ", length(pathways))

# Run GSVA
message(">> GSVA A2780 ...")
gsva_A <- run_gsva_hallmarks(mat_A, pathways, "A2780")

message(">> GSVA PEO1 ...")
gsva_P <- run_gsva_hallmarks(mat_P, pathways, "PEO1")

message(">> GSVA UWB BRCA-def ...")
gsva_def <- run_gsva_hallmarks(vst_def, pathways, "UWB_BRCAdef")

# Save GSVA score matrices
readr::write_tsv(as.data.frame(gsva_A)   %>% tibble::rownames_to_column("Pathway"),
                 file.path(gsva_tab_dir, "GSVA_A2780.tsv"))
readr::write_tsv(as.data.frame(gsva_P)   %>% tibble::rownames_to_column("Pathway"),
                 file.path(gsva_tab_dir, "GSVA_PEO1.tsv"))
readr::write_tsv(as.data.frame(gsva_def) %>% tibble::rownames_to_column("Pathway"),
                 file.path(gsva_tab_dir, "GSVA_UWB_BRCAdef.tsv"))

gsva_merged <- cbind(gsva_A, gsva_P, gsva_def)
readr::write_tsv(as.data.frame(gsva_merged) %>% tibble::rownames_to_column("Pathway"),
                 file.path(gsva_tab_dir, "GSVA_ALL_merged.tsv"))

meta_all <- bind_rows(meta_A, meta_P, meta_def)
readr::write_tsv(meta_all, file.path(gsva_tab_dir, "GSVA_metadata.tsv"))

message(">> GSVA scores guardados en: ", gsva_tab_dir)

# ============================================================
# SECTION 2: limma on GSVA scores
# ============================================================
run_limma_gsva <- function(mat, meta, dataset_label) {
  meta <- meta %>%
    mutate(condition = factor(condition, levels = c("Parental", "Resistant"))) %>%
    arrange(condition, sample_id)

  missing_cols <- setdiff(meta$sample_id, colnames(mat))
  if (length(missing_cols) > 0) {
    stop(dataset_label, ": faltan columnas: ", paste(missing_cols, collapse = ", "))
  }
  mat    <- mat[, meta$sample_id, drop = FALSE]
  design <- model.matrix(~ 0 + condition, data = meta)
  colnames(design) <- levels(meta$condition)

  fit  <- lmFit(mat, design)
  contr <- makeContrasts(Res_vs_Par = Resistant - Parental, levels = design)
  fit2 <- contrasts.fit(fit, contr)
  fit2 <- eBayes(fit2)

  topTable(fit2, coef = "Res_vs_Par", number = Inf, sort.by = "P") %>%
    tibble::rownames_to_column("Pathway") %>%
    as_tibble() %>%
    mutate(
      Dataset   = dataset_label,
      direction = case_when(
        !is.na(adj.P.Val) & adj.P.Val < FDR_CUTOFF & logFC > 0 ~ "UP",
        !is.na(adj.P.Val) & adj.P.Val < FDR_CUTOFF & logFC < 0 ~ "DOWN",
        TRUE ~ "NS"
      )
    )
}

cat(">> limma sobre scores GSVA...\n")
res_A    <- run_limma_gsva(gsva_A,   meta_A,   "A2780")
res_P    <- run_limma_gsva(gsva_P,   meta_P,   "PEO1")
res_def  <- run_limma_gsva(gsva_def, meta_def, "UWB_BRCAdef")

readr::write_tsv(res_A,   file.path(gsva_DE_dir, "GSVA_DE_A2780.tsv"))
readr::write_tsv(res_P,   file.path(gsva_DE_dir, "GSVA_DE_PEO1.tsv"))
readr::write_tsv(res_def, file.path(gsva_DE_dir, "GSVA_DE_UWB_BRCAdef.tsv"))

res_all <- bind_rows(res_A, res_P, res_def)
readr::write_tsv(res_all, file.path(gsva_DE_dir, "GSVA_DE_ALL_datasets.tsv"))
message(">> GSVA-DE tablas guardadas en: ", gsva_DE_dir)

# ============================================================
# SECTION 3: Lollipop plots
# ============================================================
clean_pathway <- function(x) {
  x %>%
    stringr::str_replace_all("^HALLMARK_", "") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_to_title() %>%
    stringr::str_replace("^E2f$", "E2F") %>%
    stringr::str_replace("^Nfkb$", "NFKB") %>%
    stringr::str_replace("^Il6 Jak Stat3$", "IL6 JAK STAT3") %>%
    stringr::str_replace("^Tnfa", "TNFA") %>%
    stringr::str_replace("^Uv$", "UV")
}

select_paths <- function(res, dataset_label, fdr_cut = FDR_CUTOFF_VIZ,
                         top_up = 8, top_down = 8) {
  res_sig   <- res %>% filter(!is.na(adj.P.Val), adj.P.Val < fdr_cut)
  res_up    <- res_sig %>% filter(logFC > 0) %>% arrange(adj.P.Val) %>% slice_head(n = top_up)
  res_down  <- res_sig %>% filter(logFC < 0) %>% arrange(adj.P.Val) %>% slice_head(n = top_down)
  res_plot  <- bind_rows(res_up, res_down) %>%
    mutate(direction = if_else(logFC > 0, "UP", "DOWN"),
           Pathway_clean = clean_pathway(Pathway),
           dataset = dataset_label)
  list(res_plot = res_plot, n_sig = nrow(res_sig),
       n_up = nrow(res_up), n_down = nrow(res_down))
}

plot_lollipop <- function(df, dataset_label, n_sig, n_up, n_down, x_limits) {
  df <- df %>% mutate(Pathway = forcats::fct_reorder(Pathway_clean, logFC))
  ggplot(df, aes(x = logFC, y = Pathway, color = direction)) +
    geom_segment(aes(x = 0, xend = logFC, yend = Pathway), linewidth = 0.7, alpha = 0.8) +
    geom_point(aes(size = -log10(adj.P.Val)), alpha = 0.9) +
    scale_color_manual(values = pal_dir, name = "Direction") +
    scale_size_continuous(name = "-log10(FDR)", range = c(2, 6)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5)) +
    labs(title = dataset_label,
         subtitle = paste0("FDR<0.1: n=", n_sig, " | UP=", n_up, ", DOWN=", n_down),
         x = "logFC (GSVA score: Resistant - Parental)", y = NULL) +
    coord_cartesian(xlim = c(-x_limits, x_limits))
}

dataset_results <- list(
  list(label = "A2780",       res = res_A,   file_stub = "GSVA_DE_A2780"),
  list(label = "PEO1",        res = res_P,   file_stub = "GSVA_DE_PEO1"),
  list(label = "UWB BRCAdef", res = res_def, file_stub = "GSVA_DE_UWB_BRCAdef")
)

plot_list  <- list()
all_logfc  <- c()
for (ds in dataset_results) {
  top_info  <- select_paths(ds$res, ds$label)
  res_plot  <- top_info$res_plot
  if (nrow(res_plot) == 0) next
  all_logfc <- c(all_logfc, res_plot$logFC)
  plot_list[[ds$label]] <- list(
    plot_df   = res_plot,
    n_sig     = top_info$n_sig,
    n_up      = top_info$n_up,
    n_down    = top_info$n_down,
    file_stub = ds$file_stub
  )
}

x_lim <- if (length(all_logfc) > 0) max(1, max(abs(all_logfc), na.rm = TRUE)) else 1

plots_out <- list()
for (nm in names(plot_list)) {
  entry  <- plot_list[[nm]]
  p      <- plot_lollipop(entry$plot_df, nm, entry$n_sig, entry$n_up, entry$n_down, x_lim)
  plots_out[[nm]] <- p
  ggsave(file.path(gsva_fig_dir, paste0(entry$file_stub, "_lollipop_top10.png")),
         p, width = 7, height = 5, dpi = 300)
  ggsave(file.path(gsva_fig_dir, paste0(entry$file_stub, "_lollipop_top10.pdf")),
         p, width = 7, height = 5)
}

if (length(plots_out) > 0) {
  combined <- patchwork::wrap_plots(plots_out, ncol = 2) +
    patchwork::plot_annotation(
      title = "GSVA differential pathways (Resistant - Parental): Top significant (FDR<0.1)"
    )
  ggsave(file.path(gsva_fig_dir, "GSVA_DE_lollipop_top10_all.png"),
         combined, width = 12, height = 10, dpi = 300)
  ggsave(file.path(gsva_fig_dir, "GSVA_DE_lollipop_top10_all.pdf"),
         combined, width = 12, height = 10)
}

message(">> GSVA analysis COMPLETADO.")
stop_log(log_handle)
