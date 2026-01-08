#!/usr/bin/env Rscript

# ============================================================
# Olaparib resistance - unified runner
# Reúne:
#  - GSE153867 (A2780, FPKM, limma)
#  - GSE235980 (UWB1.289, counts, DESeq2)
#  - GSE117765 (PEO1, TPM, limma)
# Mantiene los mismos outputs TSV/PNG que los scripts originales.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(DESeq2)
  library(R.utils)
  library(matrixStats)
})

# ------------ RUTAS DEL PROYECTO ----------------------------
proj_dir <- Sys.getenv("OLAPARIB_RESISTANCE_DIR", unset = "")
if (!nzchar(proj_dir)) {
  proj_dir <- "~/bioinfo/projects/olaparib_resistance"
}
proj_dir <- path.expand(proj_dir)
if (!dir.exists(proj_dir)) {
  stop("No existe el directorio del proyecto: ", proj_dir)
}

raw_dir        <- file.path(proj_dir, "data", "raw")
processed_dir  <- file.path(proj_dir, "data", "processed")
tables_dir     <- file.path(proj_dir, "results", "tables")
figures_dir    <- file.path(proj_dir, "results", "figures")

dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir,   showWarnings = FALSE, recursive = TRUE)

volcano_pal <- c(
  "Up (Resistant)"   = "#D81B60",
  "Down (Resistant)" = "#1E88E5",
  "NS"               = "#9E9E9E"
)
volcano_shapes <- c(
  "Up (Resistant)"   = 15,
  "Down (Resistant)" = 16,
  "NS"               = 17
)

# ============================================================
# GSE153867 - A2780 parental vs Olaparib-resistant (FPKM, limma)
# ============================================================
run_gse153867_limma_fpkm <- function() {
  raw_dir_ds <- file.path(raw_dir, "GSE153867")
  fpkm_file  <- file.path(raw_dir_ds, "GSE153867_fpkm.txt")

  if (!file.exists(fpkm_file)) {
    stop("No encuentro GSE153867_fpkm.txt en: ", raw_dir_ds)
  }

  message(">> Leyendo matriz de FPKM: ", fpkm_file)
  fpkm_raw <- readr::read_tsv(fpkm_file)

  gene_ids <- fpkm_raw[[1]]
  expr_mat <- as.matrix(fpkm_raw[, -1])
  rownames(expr_mat) <- gene_ids

  resistant_ids <- paste0("C-", 1:8)
  parental_ids  <- paste0("O-", 1:8)
  samples_keep  <- c(resistant_ids, parental_ids)

  missing_cols <- setdiff(samples_keep, colnames(expr_mat))
  if (length(missing_cols) > 0) {
    stop("Las siguientes columnas esperadas no están en la matriz FPKM:\n",
         paste(missing_cols, collapse = ", "),
         "\nRevisa los nombres de columnas en GSE153867_fpkm.txt")
  }

  expr_a2780 <- expr_mat[, samples_keep]
  message("Dimensiones submatriz A2780: ",
          nrow(expr_a2780), " genes x ", ncol(expr_a2780), " muestras")

  fpkm_a2780_tsv <- file.path(processed_dir, "GSE153867_A2780_FPKM.tsv")
  fpkm_a2780_tbl <- as_tibble(expr_a2780, rownames = "gene_id")
  readr::write_tsv(fpkm_a2780_tbl, fpkm_a2780_tsv)
  message("Matriz FPKM A2780 guardada en: ", fpkm_a2780_tsv)

  message(">> Transformando a log2(FPKM + 1) ...")
  log_expr <- log2(expr_a2780 + 1)

  # Filtrar genes con varianza > 0 antes de limma
  gene_sd <- matrixStats::rowSds(log_expr)
  keep_genes <- is.finite(gene_sd) & gene_sd > 0
  log_expr <- log_expr[keep_genes, , drop = FALSE]

  group <- factor(
    ifelse(colnames(log_expr) %in% resistant_ids, "Resistant", "Parental"),
    levels = c("Parental", "Resistant")
  )

  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  message("Diseño limma:")
  print(design)

  fit <- lmFit(log_expr, design)
  contrast_matrix <- makeContrasts(
    Resistant_vs_Parental = Resistant - Parental,
    levels = design
  )

  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  res_limma <- topTable(
    fit2,
    coef     = "Resistant_vs_Parental",
    number   = Inf,
    sort.by  = "P"
  )

  res_limma <- as.data.frame(res_limma)
  available_cols <- colnames(res_limma)
  logfc_candidates <- c("logFC",
                        "logFC.Resistant_vs_Parental",
                        "Resistant_vs_Parental",
                        "coef",
                        "Estimate")
  logfc_col <- intersect(available_cols, logfc_candidates)

  if (length(logfc_col) == 0 &&
      "Resistant_vs_Parental" %in% colnames(fit2$coefficients)) {
    # Fallback: tomar coeficiente directamente del objeto fit2
    res_limma$logFC <- fit2$coefficients[, "Resistant_vs_Parental"]
    available_cols <- colnames(res_limma)
    logfc_col <- "logFC"
  }

  if (length(logfc_col) == 0) {
    stop("No se encontró la columna logFC en topTable. Columnas disponibles: ",
         paste(available_cols, collapse = ", "))
  }

  res_limma <- res_limma %>% rownames_to_column("gene_id")
  names(res_limma)[names(res_limma) == logfc_col[1]] <- "log2FC"
  if ("AveExpr" %in% names(res_limma)) names(res_limma)[names(res_limma) == "AveExpr"] <- "AveExpr"
  if ("t" %in% names(res_limma)) names(res_limma)[names(res_limma) == "t"] <- "t_stat"
  if ("P.Value" %in% names(res_limma)) names(res_limma)[names(res_limma) == "P.Value"] <- "P.Value"
  if ("adj.P.Val" %in% names(res_limma)) names(res_limma)[names(res_limma) == "adj.P.Val"] <- "adj.P.Val"
  if ("B" %in% names(res_limma)) names(res_limma)[names(res_limma) == "B"] <- "B_stat"

  out_res_file <- file.path(tables_dir, "GSE153867_A2780_limma_DE_FPKM.tsv")
  readr::write_tsv(res_limma, out_res_file)
  message("Tabla DE limma guardada en: ", out_res_file)

  message(">> Generando volcano plot...")
  volcano_data <- res_limma %>%
    drop_na(log2FC, adj.P.Val) %>%
    mutate(
      neg_log10_FDR = -log10(adj.P.Val),
      sig = case_when(
        adj.P.Val < 0.05 & log2FC > 1  ~ "Up (Resistant)",
        adj.P.Val < 0.05 & log2FC < -1 ~ "Down (Resistant)",
        TRUE                           ~ "NS"
      ) %>% factor(levels = c("Up (Resistant)", "Down (Resistant)", "NS"))
    )

  p_volcano <- ggplot(volcano_data,
                      aes(x = log2FC, y = neg_log10_FDR,
                          color = sig, shape = sig)) +
    geom_point(alpha = 0.7, size = 1.7) +
    scale_color_manual(values = volcano_pal, drop = FALSE) +
    scale_shape_manual(values = volcano_shapes, drop = FALSE) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(
      title = "GSE153867 - A2780 (FPKM)\nResistant (C-1..8) vs Parental (O-1..8)",
      x = "log2FC (Resistant / Parental)",
      y = "-log10(FDR)"
    ) +
    theme_bw(base_size = 12)

  out_volcano_file <- file.path(figures_dir, "GSE153867_A2780_limma_FPKM_volcano.png")
  ggsave(out_volcano_file, p_volcano, width = 6, height = 5, dpi = 300)
  message("Volcano plot guardado en: ", out_volcano_file)

  message(">> Generando PCA...")
  pca <- prcomp(t(log_expr), center = TRUE, scale. = TRUE)
  pca_df <- as_tibble(pca$x[, 1:2]) %>%
    mutate(sample_id = rownames(pca$x),
           condition = group)

  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample_id)) +
    geom_point(size = 3) +
    geom_text(vjust = -1, size = 3) +
    theme_bw(base_size = 12) +
    labs(title = "GSE153867 - A2780 (FPKM)\nPCA log2(FPKM+1)", color = "Condition")

  out_pca_file <- file.path(figures_dir, "GSE153867_A2780_limma_FPKM_PCA.png")
  ggsave(out_pca_file, p_pca, width = 6, height = 5, dpi = 300)
  message("PCA plot guardado en: ", out_pca_file)
  message(">> Análisis GSE153867 A2780 (FPKM + limma) COMPLETADO.")
}

# ============================================================
# GSE235980 - UWB1.289 ± BRCA1 (counts, DESeq2)
# ============================================================
run_gse235980_deseq2 <- function() {
  raw_dir_ds <- file.path(raw_dir, "GSE235980")
  counts_gz  <- file.path(raw_dir_ds, "GSE235980_CountReads.txt.gz")
  counts_txt <- file.path(raw_dir_ds, "GSE235980_CountReads.txt")

  if (!file.exists(counts_txt) && file.exists(counts_gz)) {
    message(">> Descomprimiendo GSE235980_CountReads.txt.gz ...")
    R.utils::gunzip(counts_gz, overwrite = TRUE)
  }

  if (!file.exists(counts_txt)) {
    stop("No encuentro ni GSE235980_CountReads.txt.gz ni GSE235980_CountReads.txt en: ", raw_dir_ds)
  }

  message(">> Leyendo matriz de counts: ", counts_txt)
  cts_raw <- readr::read_tsv(counts_txt)

  gene_ids <- cts_raw[[1]]
  cts_mat  <- as.matrix(cts_raw[, -1])
  rownames(cts_mat) <- gene_ids

  message("Dimensiones de la matriz de counts: ",
          nrow(cts_mat), " genes x ", ncol(cts_mat), " muestras")

  sample_ids <- colnames(cts_mat)
  print(sample_ids)

  BRCA_status <- c("deficient", "deficient",
                  "deficient", "deficient",
                  "proficient", "proficient",
                  "proficient", "proficient")

  resistance  <- c("Resistant", "Resistant",
                  "Parental",  "Parental",
                  "Resistant", "Resistant",
                  "Parental",  "Parental")

  coldata <- tibble(
    sample_id   = sample_ids,
    BRCA_status = factor(BRCA_status, levels = c("deficient", "proficient")),
    resistance  = factor(resistance,  levels = c("Parental", "Resistant"))
  )

  message(">> colData:")
  print(coldata)

  dds <- DESeqDataSetFromMatrix(
    countData = round(cts_mat),
    colData   = as.data.frame(coldata),
    design    = ~ BRCA_status + resistance + BRCA_status:resistance
  )

  keep <- rowSums(counts(dds) >= 10) >= 2
  dds  <- dds[keep, ]
  message("Genes retenidos tras filtrado: ", nrow(dds))
  dds <- DESeq(dds)

  group <- factor(
    paste(coldata$BRCA_status, coldata$resistance, sep = "_"),
    levels = c("deficient_Parental",
               "deficient_Resistant",
               "proficient_Parental",
               "proficient_Resistant")
  )

  coldata_group <- coldata %>%
    mutate(group = group)

  dds_group <- DESeqDataSetFromMatrix(
    countData = round(cts_mat[rownames(dds), ]),
    colData   = as.data.frame(coldata_group),
    design    = ~ group
  )

  dds_group <- DESeq(dds_group)

  res_def <- results(dds_group,
                     contrast = c("group",
                                  "deficient_Resistant",
                                  "deficient_Parental"))

  res_def_shrunk <- lfcShrink(dds_group,
                              coef = "group_deficient_Resistant_vs_deficient_Parental",
                              type = "apeglm")

  res_def_tbl <- as_tibble(res_def_shrunk, rownames = "gene_id") %>%
    arrange(padj)

  out_def_file <- file.path(tables_dir, "GSE235980_BRCAdef_Res_vs_Par_DESeq2.tsv")
  readr::write_tsv(res_def_tbl, out_def_file)
  message("Resultados BRCA1-deficient guardados en: ", out_def_file)

  volcano_def <- res_def_tbl %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      neg_log10_padj = -log10(padj),
      sig = case_when(
        padj < 0.05 & log2FoldChange > 1  ~ "Up (Resistant)",
        padj < 0.05 & log2FoldChange < -1 ~ "Down (Resistant)",
        TRUE                              ~ "NS"
      ) %>% factor(levels = c("Up (Resistant)", "Down (Resistant)", "NS"))
    )

  p_volcano_def <- ggplot(volcano_def,
                          aes(x = log2FoldChange, y = neg_log10_padj,
                              color = sig, shape = sig)) +
    geom_point(alpha = 0.7, size = 1.7) +
    scale_color_manual(values = volcano_pal, drop = FALSE) +
    scale_shape_manual(values = volcano_shapes, drop = FALSE) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(
      title = "GSE235980 - BRCA1-deficient\nResistant vs Parental",
      x = "log2FC (Resistant / Parental)",
      y = "-log10(FDR)"
    ) +
    theme_bw(base_size = 12)

  out_volcano_def <- file.path(figures_dir, "GSE235980_BRCAdef_volcano.png")
  ggsave(out_volcano_def, p_volcano_def, width = 6, height = 5, dpi = 300)
  message("Volcano BRCA-def guardado en: ", out_volcano_def)

  # Re-nivelar para obtener coeficiente directo proficiente Resistant vs Parental con apeglm
  dds_group_prof <- dds_group
  dds_group_prof$group <- relevel(dds_group_prof$group, "proficient_Parental")
  dds_group_prof <- DESeq(dds_group_prof)

  res_prof <- results(dds_group_prof,
                      name = "group_proficient_Resistant_vs_proficient_Parental")

  res_prof_shrunk <- lfcShrink(dds_group_prof,
                               coef = "group_proficient_Resistant_vs_proficient_Parental",
                               type = "apeglm")

  res_prof_tbl <- as_tibble(res_prof_shrunk, rownames = "gene_id") %>%
    arrange(padj)

  out_prof_file <- file.path(tables_dir, "GSE235980_BRCAprof_Res_vs_Par_DESeq2.tsv")
  readr::write_tsv(res_prof_tbl, out_prof_file)
  message("Resultados BRCA1-proficient guardados en: ", out_prof_file)

  volcano_prof <- res_prof_tbl %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      neg_log10_padj = -log10(padj),
      sig = case_when(
        padj < 0.05 & log2FoldChange > 1  ~ "Up (Resistant)",
        padj < 0.05 & log2FoldChange < -1 ~ "Down (Resistant)",
        TRUE                              ~ "NS"
      ) %>% factor(levels = c("Up (Resistant)", "Down (Resistant)", "NS"))
    )

  p_volcano_prof <- ggplot(volcano_prof,
                           aes(x = log2FoldChange, y = neg_log10_padj,
                               color = sig, shape = sig)) +
    geom_point(alpha = 0.7, size = 1.7) +
    scale_color_manual(values = volcano_pal, drop = FALSE) +
    scale_shape_manual(values = volcano_shapes, drop = FALSE) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(
      title = "GSE235980 - BRCA1-proficient\nResistant vs Parental",
      x = "log2FC (Resistant / Parental)",
      y = "-log10(FDR)"
    ) +
    theme_bw(base_size = 12)

  out_volcano_prof <- file.path(figures_dir, "GSE235980_BRCAprof_volcano.png")
  ggsave(out_volcano_prof, p_volcano_prof, width = 6, height = 5, dpi = 300)
  message("Volcano BRCA-prof guardado en: ", out_volcano_prof)

  message(">> Generando PCA global...")
  vsd <- vst(dds_group, blind = TRUE)

  pca <- prcomp(t(assay(vsd)), center = TRUE, scale. = TRUE)
  pca_df <- as_tibble(pca$x[, 1:2]) %>%
    mutate(sample_id = rownames(pca$x)) %>%
    left_join(coldata_group, by = "sample_id")

  p_pca <- ggplot(pca_df,
                  aes(x = PC1, y = PC2,
                      color = BRCA_status,
                      shape = resistance,
                      label = sample_id)) +
    geom_point(size = 3) +
    geom_text(vjust = -1, size = 3) +
    theme_bw(base_size = 12) +
    labs(title = "GSE235980 - PCA VST",
         color = "BRCA1 status",
         shape = "Resistance")

  out_pca <- file.path(figures_dir, "GSE235980_PCA.png")
  ggsave(out_pca, p_pca, width = 6, height = 5, dpi = 300)
  message("PCA global guardado en: ", out_pca)
  message(">> Análisis GSE235980 (DESeq2) COMPLETADO.")
}

# ============================================================
# GSE117765 - PEO1 parental vs olaparib-resistant clones (TPM, limma)
# ============================================================
run_gse117765_limma_tpm <- function() {
  raw_dir_ds <- file.path(raw_dir, "GSE117765")
  mat_gz  <- file.path(raw_dir_ds, "GSE117765_matrix.txt.gz")
  mat_txt <- file.path(raw_dir_ds, "GSE117765_matrix.txt")

  if (!file.exists(mat_txt) && file.exists(mat_gz)) {
    message(">> Descomprimiendo GSE117765_matrix.txt.gz ...")
    R.utils::gunzip(mat_gz, overwrite = TRUE)
  }

  if (!file.exists(mat_txt)) {
    stop("No encuentro ni GSE117765_matrix.txt.gz ni GSE117765_matrix.txt en: ", raw_dir_ds)
  }

  message(">> Leyendo matriz: ", mat_txt)
  mat_raw <- readr::read_tsv(mat_txt)

  gene_ids <- mat_raw[[1]]
  expr_mat <- as.matrix(mat_raw[, -1])
  rownames(expr_mat) <- gene_ids

  message("Dimensiones de la matriz: ",
          nrow(expr_mat), " genes x ", ncol(expr_mat), " muestras")

  sample_ids <- colnames(expr_mat)
  message("IDs de columnas:")
  print(sample_ids)

  meta <- tibble(
    sample_id = sample_ids,
    condition = case_when(
      grepl("Adherent", sample_id, ignore.case = TRUE) ~ "Parental",
      grepl("Clone",    sample_id, ignore.case = TRUE) ~ "Resistant",
      TRUE ~ NA_character_
    )
  )

  if (any(is.na(meta$condition))) {
    stop("Hay muestras sin condition asignada. Revisa los nombres de columnas.")
  }

  meta$condition <- factor(meta$condition, levels = c("Parental", "Resistant"))
  message("colData:")
  print(meta)

  meta_file <- file.path(processed_dir, "GSE117765_PEO1_metadata.tsv")
  readr::write_tsv(meta, meta_file)
  message("Metadata de muestras guardada en: ", meta_file)

  message(">> Transformando a log2(TPM + 1) ...")
  log_expr <- log2(expr_mat + 1)

  # Filtrar genes con varianza > 0 antes de limma
  gene_sd <- matrixStats::rowSds(log_expr)
  keep_genes <- is.finite(gene_sd) & gene_sd > 0
  log_expr <- log_expr[keep_genes, , drop = FALSE]

  group <- meta$condition
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  message("Diseño limma:")
  print(design)

  fit <- lmFit(log_expr, design)
  contrast_matrix <- makeContrasts(
    Resistant_vs_Parental = Resistant - Parental,
    levels = design
  )

  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  res_limma <- topTable(
    fit2,
    coef     = "Resistant_vs_Parental",
    number   = Inf,
    sort.by  = "P"
  )

  res_limma <- as.data.frame(res_limma)
  available_cols <- colnames(res_limma)
  logfc_candidates <- c("logFC",
                        "logFC.Resistant_vs_Parental",
                        "Resistant_vs_Parental",
                        "coef",
                        "Estimate")
  logfc_col <- intersect(available_cols, logfc_candidates)

  if (length(logfc_col) == 0 &&
      "Resistant_vs_Parental" %in% colnames(fit2$coefficients)) {
    res_limma$logFC <- fit2$coefficients[, "Resistant_vs_Parental"]
    available_cols <- colnames(res_limma)
    logfc_col <- "logFC"
  }

  if (length(logfc_col) == 0) {
    stop("No se encontró la columna logFC en topTable. Columnas disponibles: ",
         paste(available_cols, collapse = ", "))
  }

  res_limma <- res_limma %>% rownames_to_column("gene_id")
  names(res_limma)[names(res_limma) == logfc_col[1]] <- "log2FC"
  if ("AveExpr" %in% names(res_limma)) names(res_limma)[names(res_limma) == "AveExpr"] <- "AveExpr"
  if ("t" %in% names(res_limma)) names(res_limma)[names(res_limma) == "t"] <- "t_stat"
  if ("P.Value" %in% names(res_limma)) names(res_limma)[names(res_limma) == "P.Value"] <- "P.Value"
  if ("adj.P.Val" %in% names(res_limma)) names(res_limma)[names(res_limma) == "adj.P.Val"] <- "adj.P.Val"
  if ("B" %in% names(res_limma)) names(res_limma)[names(res_limma) == "B"] <- "B_stat"

  out_res_file <- file.path(tables_dir, "GSE117765_PEO1_limma_TPM_Res_vs_Par.tsv")
  readr::write_tsv(res_limma, out_res_file)
  message("Tabla limma guardada en: ", out_res_file)

  message(">> Generando volcano plot...")
  volcano <- res_limma %>%
    drop_na(log2FC, adj.P.Val) %>%
    mutate(
      neg_log10_FDR = -log10(adj.P.Val),
      sig = case_when(
        adj.P.Val < 0.05 & log2FC > 1  ~ "Up (Resistant)",
        adj.P.Val < 0.05 & log2FC < -1 ~ "Down (Resistant)",
        TRUE                           ~ "NS"
      ) %>% factor(levels = c("Up (Resistant)", "Down (Resistant)", "NS"))
    )

  p_volc <- ggplot(volcano,
                   aes(x = log2FC, y = neg_log10_FDR,
                       color = sig, shape = sig)) +
    geom_point(alpha = 0.7, size = 1.7) +
    scale_color_manual(values = volcano_pal, drop = FALSE) +
    scale_shape_manual(values = volcano_shapes, drop = FALSE) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(
      title = "GSE117765 - PEO1 (TPM, limma)\nResistant vs Parental",
      x = "log2FC (Resistant / Parental)",
      y = "-log10(FDR)"
    ) +
    theme_bw(base_size = 12)

  out_volc <- file.path(figures_dir, "GSE117765_PEO1_limma_TPM_volcano.png")
  ggsave(out_volc, p_volc, width = 6, height = 5, dpi = 300)
  message("Volcano guardado en: ", out_volc)

  message(">> Generando PCA...")
  pca <- prcomp(t(log_expr), center = TRUE, scale. = TRUE)
  pca_df <- as_tibble(pca$x[, 1:2]) %>%
    mutate(sample_id = rownames(pca$x)) %>%
    left_join(meta, by = "sample_id")

  p_pca <- ggplot(pca_df,
                  aes(x = PC1, y = PC2,
                      color = condition,
                      label = sample_id)) +
    geom_point(size = 3) +
    geom_text(vjust = -1, size = 3) +
    theme_bw(base_size = 12) +
    labs(title = "GSE117765 - PEO1 (TPM, limma)\nPCA log2(TPM+1)",
         color = "Condition")

  out_pca <- file.path(figures_dir, "GSE117765_PEO1_limma_TPM_PCA.png")
  ggsave(out_pca, p_pca, width = 6, height = 5, dpi = 300)
  message("PCA guardado en: ", out_pca)
  message(">> Análisis GSE117765 (TPM + limma) COMPLETADO.")
}

# ============================================================
# CLI sencillo
# ============================================================
run_selected <- function(dataset) {
  dataset <- tolower(dataset)
  if (dataset %in% c("all")) {
    run_gse153867_limma_fpkm()
    run_gse235980_deseq2()
    run_gse117765_limma_tpm()
  } else if (dataset %in% c("153867", "gse153867")) {
    run_gse153867_limma_fpkm()
  } else if (dataset %in% c("235980", "gse235980")) {
    run_gse235980_deseq2()
  } else if (dataset %in% c("117765", "gse117765")) {
    run_gse117765_limma_tpm()
  } else {
    stop("Dataset no reconocido: ", dataset,
         " (usa all | 153867 | 235980 | 117765)")
  }
}

args <- commandArgs(trailingOnly = TRUE)
sel <- if (length(args) >= 1) args[1] else "all"
run_selected(sel)
