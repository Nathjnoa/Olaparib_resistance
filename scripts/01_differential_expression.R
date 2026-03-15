#!/usr/bin/env Rscript

# ============================================================
# 01_differential_expression.R
#
# Differential expression analysis for 3 olaparib-resistance datasets:
#   - GSE153867: A2780 (FPKM, limma)
#   - GSE117765: PEO1  (TPM,  limma)
#   - GSE235980: UWB1.289 BRCA-deficient only (counts, DESeq2 ~ resistance)
#
# NOTE on method choice:
#   limma (without voom) is used for FPKM/TPM data. voom is designed
#   for raw count data and exploits the mean-variance relationship of
#   Poisson/NB counts. FPKM and TPM are pre-normalized continuous
#   estimates; standard limma on log2(x+1) is the recommended approach
#   (Ritchie et al. 2015 NAR; Law et al. 2014 Genome Biol).
#
# Outputs: DE tables (TSV), volcano plots (PNG), PCA plots (PNG).
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(DESeq2)
  library(R.utils)
  library(matrixStats)
})

source(file.path({f <- commandArgs(trailingOnly=FALSE); f <- grep("--file=",f,value=TRUE); if(length(f)) dirname(normalizePath(sub("--file=","",f[1]))) else path.expand("~/bioinfo/projects/olaparib_resistance/scripts")}, "00_config.R"))

log_handle <- start_log("01_differential_expression")

de_fig_dir <- file.path(figures_dir, "de")
de_tab_dir <- file.path(tables_dir,  "de")
dir.create(de_fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(de_tab_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# GSE153867 — A2780 (FPKM, limma)
# ============================================================
run_a2780 <- function() {
  raw_dir_ds <- file.path(raw_dir, "GSE153867")
  fpkm_file  <- file.path(raw_dir_ds, "GSE153867_fpkm.txt")

  if (!file.exists(fpkm_file)) stop("No encuentro: ", fpkm_file)

  message(">> A2780: leyendo FPKM")
  fpkm_raw  <- readr::read_tsv(fpkm_file, show_col_types = FALSE)
  gene_ids  <- fpkm_raw[[1]]
  expr_mat  <- as.matrix(fpkm_raw[, -1])
  rownames(expr_mat) <- gene_ids

  samples_keep <- c(A2780_RES, A2780_PAR)
  missing_cols <- setdiff(samples_keep, colnames(expr_mat))
  if (length(missing_cols) > 0) stop("A2780: faltan columnas: ", paste(missing_cols, collapse = ", "))

  expr_a2780 <- expr_mat[, samples_keep]

  # Save processed FPKM matrix
  readr::write_tsv(
    as_tibble(expr_a2780, rownames = "gene_id"),
    file.path(processed_dir, "GSE153867_A2780_FPKM.tsv")
  )

  log_expr <- log2(expr_a2780 + 1)
  gene_sd  <- matrixStats::rowSds(log_expr)
  log_expr <- log_expr[is.finite(gene_sd) & gene_sd > 0, , drop = FALSE]

  group <- factor(
    ifelse(colnames(log_expr) %in% A2780_RES, "Resistant", "Parental"),
    levels = c("Parental", "Resistant")
  )

  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  fit  <- lmFit(log_expr, design)
  contr <- makeContrasts(Resistant_vs_Parental = Resistant - Parental, levels = design)
  fit2 <- contrasts.fit(fit, contr)
  fit2 <- eBayes(fit2)

  res <- topTable(fit2, coef = "Resistant_vs_Parental", number = Inf, sort.by = "P") %>%
    dplyr::rename(log2FC = logFC, t_stat = t)

  out_de <- file.path(de_tab_dir, "GSE153867_A2780_limma_DE_FPKM.tsv")
  readr::write_tsv(res, out_de)
  message("A2780 DE -> ", out_de)

  # Volcano
  volc_data <- res %>%
    drop_na(log2FC, adj.P.Val) %>%
    mutate(
      neg_log10_FDR = -log10(adj.P.Val),
      sig = case_when(
        adj.P.Val < FDR_CUTOFF & log2FC >  1 ~ "Up (Resistant)",
        adj.P.Val < FDR_CUTOFF & log2FC < -1 ~ "Down (Resistant)",
        TRUE ~ "NS"
      ) %>% factor(levels = c("Up (Resistant)", "Down (Resistant)", "NS"))
    )

  p_vol <- ggplot(volc_data, aes(x = log2FC, y = neg_log10_FDR, color = sig, shape = sig)) +
    geom_point(alpha = 0.7, size = 1.7) +
    scale_color_manual(values = volcano_pal, drop = FALSE) +
    scale_shape_manual(values = volcano_shapes, drop = FALSE) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed") +
    labs(title = "GSE153867 — A2780 (FPKM)\nResistant vs Parental",
         x = "log2FC (Resistant / Parental)", y = "-log10(FDR)") +
    theme_bw(base_size = 12)

  ggsave(file.path(de_fig_dir, "GSE153867_A2780_limma_FPKM_volcano.png"),
         p_vol, width = 6, height = 5, dpi = 300)

  # PCA
  pca    <- prcomp(t(log_expr), center = TRUE, scale. = TRUE)
  pca_df <- as_tibble(pca$x[, 1:2]) %>%
    mutate(sample_id = rownames(pca$x), condition = group)

  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample_id)) +
    geom_point(size = 3) +
    geom_text(vjust = -1, size = 3) +
    theme_bw(base_size = 12) +
    labs(title = "GSE153867 — A2780 PCA log2(FPKM+1)", color = "Condition")

  ggsave(file.path(de_fig_dir, "GSE153867_A2780_limma_FPKM_PCA.png"),
         p_pca, width = 6, height = 5, dpi = 300)

  message(">> A2780 COMPLETADO")
}

# ============================================================
# GSE235980 — UWB1.289 BRCA-deficient only (counts, DESeq2)
# ============================================================
run_uwb_brcadef <- function() {
  raw_dir_ds <- file.path(raw_dir, "GSE235980")
  counts_gz  <- file.path(raw_dir_ds, "GSE235980_CountReads.txt.gz")
  counts_txt <- file.path(raw_dir_ds, "GSE235980_CountReads.txt")

  if (!file.exists(counts_txt) && file.exists(counts_gz)) {
    message(">> Descomprimiendo GSE235980_CountReads.txt.gz ...")
    R.utils::gunzip(counts_gz, overwrite = TRUE)
  }
  if (!file.exists(counts_txt)) stop("No encuentro counts UWB: ", counts_txt)

  message(">> UWB BRCA-def: leyendo counts")
  cts_raw  <- readr::read_tsv(counts_txt, show_col_types = FALSE)
  gene_ids <- cts_raw[[1]]
  cts_mat  <- as.matrix(cts_raw[, -1])
  rownames(cts_mat) <- gene_ids

  # Only BRCA-deficient samples (4 samples total)
  samples_def <- c(UWB_DEF_PAR, UWB_DEF_RES)
  missing <- setdiff(samples_def, colnames(cts_mat))
  if (length(missing) > 0) stop("UWB: faltan columnas: ", paste(missing, collapse = ", "))

  cts_def <- cts_mat[, samples_def, drop = FALSE]

  coldata <- data.frame(
    sample_id  = samples_def,
    resistance = factor(
      c(rep("Parental", length(UWB_DEF_PAR)), rep("Resistant", length(UWB_DEF_RES))),
      levels = c("Parental", "Resistant")
    ),
    row.names = samples_def
  )

  message("UWB coldata:")
  print(coldata)

  # DESeq2 design: ~ resistance (only BRCA-def samples; BRCA_status is constant)
  dds <- DESeqDataSetFromMatrix(
    countData = round(cts_def),
    colData   = coldata,
    design    = ~ resistance
  )

  keep <- rowSums(counts(dds) >= 10) >= 2
  dds  <- dds[keep, ]
  message("UWB genes tras filtrado: ", nrow(dds))
  dds <- DESeq(dds)

  res_def <- results(dds, contrast = c("resistance", "Resistant", "Parental"))

  res_def_shrunk <- lfcShrink(
    dds,
    coef = "resistance_Resistant_vs_Parental",
    type = "apeglm"
  )

  res_tbl <- as_tibble(res_def_shrunk, rownames = "gene_id") %>% arrange(padj)

  out_def <- file.path(de_tab_dir, "GSE235980_BRCAdef_Res_vs_Par_DESeq2.tsv")
  readr::write_tsv(res_tbl, out_def)
  message("UWB BRCA-def DE -> ", out_def)

  # Volcano
  volc_def <- res_tbl %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      neg_log10_padj = -log10(padj),
      sig = case_when(
        padj < FDR_CUTOFF & log2FoldChange >  1 ~ "Up (Resistant)",
        padj < FDR_CUTOFF & log2FoldChange < -1 ~ "Down (Resistant)",
        TRUE ~ "NS"
      ) %>% factor(levels = c("Up (Resistant)", "Down (Resistant)", "NS"))
    )

  p_vol_def <- ggplot(volc_def, aes(x = log2FoldChange, y = neg_log10_padj, color = sig, shape = sig)) +
    geom_point(alpha = 0.7, size = 1.7) +
    scale_color_manual(values = volcano_pal, drop = FALSE) +
    scale_shape_manual(values = volcano_shapes, drop = FALSE) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed") +
    labs(title = "GSE235980 — UWB BRCA-deficient\nResistant vs Parental",
         x = "log2FC (Resistant / Parental)", y = "-log10(FDR)") +
    theme_bw(base_size = 12)

  ggsave(file.path(de_fig_dir, "GSE235980_BRCAdef_volcano.png"),
         p_vol_def, width = 6, height = 5, dpi = 300)

  # PCA (4 samples, colored by resistance only)
  vsd    <- vst(dds, blind = TRUE)
  pca    <- prcomp(t(assay(vsd)), center = TRUE, scale. = TRUE)
  pca_df <- as_tibble(pca$x[, 1:2]) %>%
    mutate(sample_id = rownames(pca$x)) %>%
    left_join(as_tibble(coldata) %>% dplyr::select(sample_id, resistance), by = "sample_id")

  p_pca_def <- ggplot(pca_df, aes(x = PC1, y = PC2, color = resistance, label = sample_id)) +
    geom_point(size = 3) +
    geom_text(vjust = -1, size = 3) +
    theme_bw(base_size = 12) +
    labs(title = "GSE235980 — UWB BRCA-def PCA (VST)",
         color = "Resistance")

  ggsave(file.path(de_fig_dir, "GSE235980_BRCAdef_PCA.png"),
         p_pca_def, width = 6, height = 5, dpi = 300)

  message(">> UWB BRCA-def COMPLETADO")
}

# ============================================================
# GSE117765 — PEO1 (TPM, limma)
# ============================================================
run_peo1 <- function() {
  raw_dir_ds <- file.path(raw_dir, "GSE117765")
  mat_gz  <- file.path(raw_dir_ds, "GSE117765_matrix.txt.gz")
  mat_txt <- file.path(raw_dir_ds, "GSE117765_matrix.txt")

  if (!file.exists(mat_txt) && file.exists(mat_gz)) {
    message(">> Descomprimiendo GSE117765_matrix.txt.gz ...")
    R.utils::gunzip(mat_gz, overwrite = TRUE)
  }
  if (!file.exists(mat_txt)) stop("No encuentro PEO1 matrix: ", mat_txt)

  message(">> PEO1: leyendo TPM")
  mat_raw  <- readr::read_tsv(mat_txt, show_col_types = FALSE)
  gene_ids <- mat_raw[[1]]
  expr_mat <- as.matrix(mat_raw[, -1])
  rownames(expr_mat) <- gene_ids

  sample_ids <- colnames(expr_mat)
  meta <- tibble(
    sample_id = sample_ids,
    condition = case_when(
      grepl("Adherent", sample_id, ignore.case = TRUE) ~ "Parental",
      grepl("Clone",    sample_id, ignore.case = TRUE) ~ "Resistant",
      TRUE ~ NA_character_
    )
  )

  if (any(is.na(meta$condition))) stop("PEO1: muestras sin condicion asignada.")

  meta$condition <- factor(meta$condition, levels = c("Parental", "Resistant"))

  readr::write_tsv(meta, file.path(processed_dir, "GSE117765_PEO1_metadata.tsv"))

  log_expr <- log2(expr_mat + 1)
  gene_sd  <- matrixStats::rowSds(log_expr)
  log_expr <- log_expr[is.finite(gene_sd) & gene_sd > 0, , drop = FALSE]

  group  <- meta$condition
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  fit  <- lmFit(log_expr, design)
  contr <- makeContrasts(Resistant_vs_Parental = Resistant - Parental, levels = design)
  fit2 <- contrasts.fit(fit, contr)
  fit2 <- eBayes(fit2)

  res <- topTable(fit2, coef = "Resistant_vs_Parental", number = Inf, sort.by = "P") %>%
    dplyr::rename(log2FC = logFC, t_stat = t)

  out_peo1 <- file.path(de_tab_dir, "GSE117765_PEO1_limma_TPM_Res_vs_Par.tsv")
  readr::write_tsv(res, out_peo1)
  message("PEO1 DE -> ", out_peo1)

  # Volcano
  volc_p <- res %>%
    drop_na(log2FC, adj.P.Val) %>%
    mutate(
      neg_log10_FDR = -log10(adj.P.Val),
      sig = case_when(
        adj.P.Val < FDR_CUTOFF & log2FC >  1 ~ "Up (Resistant)",
        adj.P.Val < FDR_CUTOFF & log2FC < -1 ~ "Down (Resistant)",
        TRUE ~ "NS"
      ) %>% factor(levels = c("Up (Resistant)", "Down (Resistant)", "NS"))
    )

  p_vol_p <- ggplot(volc_p, aes(x = log2FC, y = neg_log10_FDR, color = sig, shape = sig)) +
    geom_point(alpha = 0.7, size = 1.7) +
    scale_color_manual(values = volcano_pal, drop = FALSE) +
    scale_shape_manual(values = volcano_shapes, drop = FALSE) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed") +
    labs(title = "GSE117765 — PEO1 (TPM)\nResistant vs Parental",
         x = "log2FC (Resistant / Parental)", y = "-log10(FDR)") +
    theme_bw(base_size = 12)

  ggsave(file.path(de_fig_dir, "GSE117765_PEO1_limma_TPM_volcano.png"),
         p_vol_p, width = 6, height = 5, dpi = 300)

  # PCA
  pca    <- prcomp(t(log_expr), center = TRUE, scale. = TRUE)
  pca_df <- as_tibble(pca$x[, 1:2]) %>%
    mutate(sample_id = rownames(pca$x)) %>%
    left_join(meta, by = "sample_id")

  p_pca_p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample_id)) +
    geom_point(size = 3) +
    geom_text(vjust = -1, size = 3) +
    theme_bw(base_size = 12) +
    labs(title = "GSE117765 — PEO1 PCA log2(TPM+1)", color = "Condition")

  ggsave(file.path(de_fig_dir, "GSE117765_PEO1_limma_TPM_PCA.png"),
         p_pca_p, width = 6, height = 5, dpi = 300)

  message(">> PEO1 COMPLETADO")
}

# ============================================================
# Run all
# ============================================================
run_a2780()
run_uwb_brcadef()
run_peo1()

stop_log(log_handle)
