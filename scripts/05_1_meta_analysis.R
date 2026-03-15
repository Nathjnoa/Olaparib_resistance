#!/usr/bin/env Rscript

# ============================================================
# 05_1_meta_analysis.R
#
# Random-effects meta-analysis across 3 datasets (A2780, PEO1,
# UWB BRCA-deficient). REML estimator via metafor::rma.uni().
#
#   Section 1: Build long/wide input tables
#   Section 2: Run REML random-effects model per gene
#   Section 3: Export signatures (UP/DOWN), rankings, I2, concordance
#
# SE derivation: SE = |logFC| / |t|  (exact for limma t-statistic)
#                SE = lfcSE           (exact from DESeq2)
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(metafor)
  library(org.Hs.eg.db)
})

source(file.path({f <- commandArgs(trailingOnly=FALSE); f <- grep("--file=",f,value=TRUE); if(length(f)) dirname(normalizePath(sub("--file=","",f[1]))) else path.expand("~/bioinfo/projects/olaparib_resistance/scripts")}, "00_config.R"))

log_handle <- start_log("05_1_meta_analysis")

int_dir <- file.path(tables_dir, "integrated")
dir.create(int_dir, showWarnings = FALSE, recursive = TRUE)

out_long <- file.path(int_dir, "meta_DE_input_long.tsv")
out_wide <- file.path(int_dir, "meta_DE_input_wide.tsv")
out_re   <- file.path(int_dir, "meta_DE_random_effects.tsv")

# ============================================================
# SECTION 1: Build input tables
# ============================================================

# flip_log2FC: set TRUE if Resistant - Parental is inverted for a dataset
flip_log2FC <- c(A2780 = FALSE, PEO1 = FALSE, UWB_BRCAdef = FALSE)

.pick_first <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

load_limma_deg <- function(file, dataset_id) {
  message(">> Cargando limma: ", dataset_id)
  df <- readr::read_tsv(file, show_col_types = FALSE)

  gene_col  <- intersect(c("ID", "gene", "Gene", "symbol", "SYMBOL"), colnames(df)) |> .pick_first()
  logfc_col <- intersect(c("log2FC", "logFC"),          colnames(df)) |> .pick_first()
  t_col     <- intersect(c("t_stat", "t", "tstat"),     colnames(df)) |> .pick_first()
  p_col     <- intersect(c("P.Value", "pvalue", "pval"), colnames(df)) |> .pick_first()
  padj_col  <- intersect(c("adj.P.Val", "padj", "FDR"), colnames(df)) |> .pick_first()

  needed <- c(gene_col, logfc_col, t_col, p_col, padj_col)
  if (any(is.na(needed))) stop("Faltan columnas en limma para ", dataset_id)

  df2 <- df %>%
    transmute(
      gene_symbol = as.character(.data[[gene_col]]),
      log2FC      = as.numeric(.data[[logfc_col]]),
      t_stat      = as.numeric(.data[[t_col]]),
      pval        = as.numeric(.data[[p_col]]),
      padj        = as.numeric(.data[[padj_col]]),
      dataset     = dataset_id
    ) %>%
    filter(!is.na(gene_symbol), gene_symbol != "") %>%
    mutate(SE = if_else(is.na(t_stat) | t_stat == 0, NA_real_, abs(log2FC / t_stat))) %>%
    dplyr::select(gene_symbol, log2FC, SE, pval, padj, dataset) %>%
    group_by(gene_symbol) %>%
    slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
    ungroup()

  df2
}

load_deseq_deg <- function(file, dataset_id) {
  message(">> Cargando DESeq2: ", dataset_id)
  df <- readr::read_tsv(file, show_col_types = FALSE)

  gene_col <- intersect(c("gene_id","GeneID","gene","Gene","id"), colnames(df)) |> .pick_first()
  lfc_col  <- intersect(c("log2FoldChange","logFC","log2FC"),     colnames(df)) |> .pick_first()
  se_col   <- intersect(c("lfcSE","SE"),                          colnames(df)) |> .pick_first()
  p_col    <- intersect(c("pvalue","pval","P.Value"),             colnames(df)) |> .pick_first()
  padj_col <- intersect(c("padj","adj.P.Val","FDR"),             colnames(df)) |> .pick_first()

  needed <- c(gene_col, lfc_col, se_col, p_col, padj_col)
  if (any(is.na(needed))) stop("Faltan columnas en DESeq2 para ", dataset_id)

  gene_vec <- as.character(df[[gene_col]])
  is_ensg  <- all(grepl("^ENS", gene_vec))

  if (is_ensg) {
    df <- df %>% mutate(ensg = sub("\\.\\d+$", "", .data[[gene_col]]))
    map <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys    = unique(df$ensg),
      keytype = "ENSEMBL",
      columns = "SYMBOL"
    ) %>%
      arrange(ENSEMBL, SYMBOL) %>%
      group_by(ENSEMBL) %>%
      summarise(SYMBOL = .pick_first(SYMBOL), .groups = "drop")
    df2 <- df %>%
      left_join(map, by = c("ensg" = "ENSEMBL")) %>%
      rename(gene_symbol = SYMBOL)
  } else {
    df2 <- df %>% mutate(gene_symbol = .data[[gene_col]])
  }

  df2 <- df2 %>%
    transmute(
      gene_symbol = as.character(gene_symbol),
      log2FC      = as.numeric(.data[[lfc_col]]),
      SE          = as.numeric(.data[[se_col]]),
      pval        = as.numeric(.data[[p_col]]),
      padj        = as.numeric(.data[[padj_col]]),
      dataset     = dataset_id
    ) %>%
    filter(!is.na(gene_symbol), gene_symbol != "") %>%
    group_by(gene_symbol) %>%
    slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
    ungroup()

  df2
}

# Load DE tables
file_A2780  <- file.path(tables_dir, "de", "GSE153867_A2780_limma_DE_FPKM.tsv")
file_PEO1   <- file.path(tables_dir, "de", "GSE117765_PEO1_limma_TPM_Res_vs_Par.tsv")
file_UWBdef <- file.path(tables_dir, "de", "GSE235980_BRCAdef_Res_vs_Par_DESeq2.tsv")

stopifnot(file.exists(file_A2780), file.exists(file_PEO1), file.exists(file_UWBdef))

deg_A2780  <- load_limma_deg(file_A2780,  "A2780")
deg_PEO1   <- load_limma_deg(file_PEO1,   "PEO1")
deg_UWBdef <- load_deseq_deg(file_UWBdef, "UWB_BRCAdef")

deg_long <- bind_rows(deg_A2780, deg_PEO1, deg_UWBdef) %>%
  mutate(
    log2FC   = if_else(flip_log2FC[dataset], -log2FC, log2FC),
    contrast = "Resistant_minus_Parental"
  ) %>%
  arrange(gene_symbol, dataset)

# QC
qc <- deg_long %>%
  group_by(dataset) %>%
  summarise(
    n_genes          = n(),
    n_ok_se          = sum(!is.na(SE) & SE > 0),
    n_fdr_01         = sum(!is.na(padj) & padj < 0.10),
    up_fdr_01        = sum(!is.na(padj) & padj < 0.10 & log2FC > 0),
    down_fdr_01      = sum(!is.na(padj) & padj < 0.10 & log2FC < 0),
    median_abs_log2FC = median(abs(log2FC), na.rm = TRUE),
    .groups = "drop"
  )
message("---- QC resumen ----")
print(qc)
message("Total filas: ", nrow(deg_long), " | Genes únicos: ", n_distinct(deg_long$gene_symbol))

readr::write_tsv(deg_long, out_long)
message(">> Guardado: ", out_long)

# Wide format
deg_wide <- deg_long %>%
  dplyr::select(gene_symbol, dataset, log2FC, SE, pval, padj) %>%
  tidyr::pivot_longer(cols = c(log2FC, SE, pval, padj), names_to = "stat", values_to = "value") %>%
  tidyr::unite("colname", stat, dataset) %>%
  tidyr::pivot_wider(names_from = colname, values_from = value)

readr::write_tsv(deg_wide, out_wide)
message(">> Guardado: ", out_wide)

# ============================================================
# SECTION 2: Random-effects meta-analysis (REML)
# ============================================================
message(">> Iniciando meta-análisis por gen...")

run_meta_for_gene <- function(df_gene) {
  df_ok <- df_gene %>%
    filter(!is.na(log2FC), is.finite(log2FC), !is.na(SE), is.finite(SE), SE > 0)

  k       <- nrow(df_ok)
  ds_list <- paste(sort(unique(df_ok$dataset)), collapse = ",")

  empty <- tibble(
    k = k, datasets = ds_list, reml_warning = FALSE,
    beta_meta = NA_real_, se_meta = NA_real_, z_meta = NA_real_, p_meta = NA_real_,
    ci_lb = NA_real_, ci_ub = NA_real_,
    tau2 = NA_real_, Q = NA_real_, Q_pval = NA_real_, I2 = NA_real_,
    beta_fixed = NA_real_, se_fixed = NA_real_, z_fixed = NA_real_, p_fixed = NA_real_,
    beta_meta_DL = NA_real_, se_meta_DL = NA_real_, p_meta_DL = NA_real_, tau2_DL = NA_real_
  )

  if (k < 2) return(empty)

  reml_warn_flag <- FALSE
  fit_re <- tryCatch(
    withCallingHandlers(
      rma.uni(yi = df_ok$log2FC, sei = df_ok$SE, method = "REML"),
      warning = function(w) {
        if (grepl("Fisher scoring algorithm may have gotten stuck", conditionMessage(w))) {
          reml_warn_flag <<- TRUE
          invokeRestart("muffleWarning")
        }
      }
    ),
    error = function(e) NULL
  )

  fit_dl <- NULL
  if (reml_warn_flag) {
    fit_dl <- tryCatch(rma.uni(yi = df_ok$log2FC, sei = df_ok$SE, method = "DL"),
                       error = function(e) NULL)
  }

  fit_fe <- tryCatch(rma.uni(yi = df_ok$log2FC, sei = df_ok$SE, method = "FE"),
                     error = function(e) NULL)

  tibble(
    k = k, datasets = ds_list, reml_warning = reml_warn_flag,
    beta_meta = if (is.null(fit_re)) NA_real_ else as.numeric(fit_re$b),
    se_meta   = if (is.null(fit_re)) NA_real_ else fit_re$se,
    z_meta    = if (is.null(fit_re)) NA_real_ else fit_re$zval,
    p_meta    = if (is.null(fit_re)) NA_real_ else fit_re$pval,
    ci_lb     = if (is.null(fit_re)) NA_real_ else fit_re$ci.lb,
    ci_ub     = if (is.null(fit_re)) NA_real_ else fit_re$ci.ub,
    tau2      = if (is.null(fit_re)) NA_real_ else fit_re$tau2,
    Q         = if (is.null(fit_re)) NA_real_ else fit_re$QE,
    Q_pval    = if (is.null(fit_re)) NA_real_ else fit_re$QEp,
    I2        = if (is.null(fit_re)) NA_real_ else fit_re$I2,
    beta_fixed = if (is.null(fit_fe)) NA_real_ else as.numeric(fit_fe$b),
    se_fixed   = if (is.null(fit_fe)) NA_real_ else fit_fe$se,
    z_fixed    = if (is.null(fit_fe)) NA_real_ else fit_fe$zval,
    p_fixed    = if (is.null(fit_fe)) NA_real_ else fit_fe$pval,
    beta_meta_DL = if (is.null(fit_dl)) NA_real_ else as.numeric(fit_dl$b),
    se_meta_DL   = if (is.null(fit_dl)) NA_real_ else as.numeric(fit_dl$se),
    p_meta_DL    = if (is.null(fit_dl)) NA_real_ else as.numeric(fit_dl$pval),
    tau2_DL      = if (is.null(fit_dl)) NA_real_ else as.numeric(fit_dl$tau2)
  )
}

meta_results <- deg_long %>%
  group_by(gene_symbol) %>%
  group_modify(~ run_meta_for_gene(.x)) %>%
  ungroup() %>%
  mutate(
    reml_warning = as.logical(reml_warning),
    across(c(beta_meta_DL, se_meta_DL, p_meta_DL, tau2_DL), as.numeric)
  ) %>%
  mutate(
    FDR_meta  = p.adjust(p_meta, method = "BH"),
    direction = case_when(
      !is.na(beta_meta) & beta_meta > 0 ~ "UP",
      !is.na(beta_meta) & beta_meta < 0 ~ "DOWN",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(FDR_meta, p_meta)

readr::write_tsv(meta_results, out_re)
message(">> Meta-análisis completado: ", out_re)

warn_genes <- meta_results %>% filter(isTRUE(reml_warning))
if (nrow(warn_genes) > 0) {
  out_warn <- file.path(int_dir, "meta_DE_REML_warnings_genes.tsv")
  readr::write_tsv(warn_genes, out_warn)
  message(">> Genes con warning REML: ", out_warn)
}

message("\nTop 10 por FDR_meta:")
print(meta_results %>% dplyr::select(gene_symbol, k, beta_meta, FDR_meta, I2, datasets) %>% head(10))

# ============================================================
# SECTION 3: Export signatures
# ============================================================
message(">> Exportando firmas y rankings...")

meta <- meta_results %>%
  mutate(
    p_meta2  = if_else(!is.na(p_meta) & p_meta == 0, 1e-300, p_meta),
    FDR_meta = p.adjust(p_meta2, method = "BH"),
    z_meta   = if_else(!is.na(beta_meta) & !is.na(se_meta) & se_meta > 0,
                       beta_meta / se_meta, z_meta)
  )

# Concordance
concord <- deg_long %>%
  filter(!is.na(log2FC)) %>%
  group_by(gene_symbol) %>%
  summarise(
    k_eff       = n(),
    n_pos       = sum(log2FC > 0, na.rm = TRUE),
    n_neg       = sum(log2FC < 0, na.rm = TRUE),
    concordance = if_else(k_eff > 0, pmax(n_pos, n_neg) / k_eff, NA_real_),
    .groups = "drop"
  )

meta <- meta %>% left_join(concord, by = "gene_symbol")

# UP / DOWN lists
meta_up <- meta %>%
  filter(!is.na(FDR_meta), FDR_meta < FDR_CUTOFF, !is.na(beta_meta), beta_meta > 0) %>%
  arrange(FDR_meta)
meta_down <- meta %>%
  filter(!is.na(FDR_meta), FDR_meta < FDR_CUTOFF, !is.na(beta_meta), beta_meta < 0) %>%
  arrange(FDR_meta)

readr::write_tsv(meta_up,   file.path(int_dir, "meta_DE_meta_up.tsv"))
readr::write_tsv(meta_down, file.path(int_dir, "meta_DE_meta_down.tsv"))
message("Genes UP (FDR<", FDR_CUTOFF, "): ",   nrow(meta_up))
message("Genes DOWN (FDR<", FDR_CUTOFF, "): ", nrow(meta_down))

# Rankings
meta_rank_z    <- meta %>% filter(!is.na(z_meta))   %>% arrange(desc(z_meta))
meta_rank_beta <- meta %>% filter(!is.na(beta_meta)) %>% arrange(desc(beta_meta))
readr::write_tsv(meta_rank_z,    file.path(int_dir, "meta_DE_meta_rank_z.tsv"))
readr::write_tsv(meta_rank_beta, file.path(int_dir, "meta_DE_meta_rank_beta.tsv"))

# I2 summary
I2_summary <- meta %>%
  filter(!is.na(I2)) %>%
  summarise(
    n_genes    = n(), mean_I2 = mean(I2), median_I2 = median(I2),
    I2_q25 = quantile(I2, 0.25), I2_q75 = quantile(I2, 0.75),
    n_I2_lt25  = sum(I2 < 25),
    n_I2_25_75 = sum(I2 >= 25 & I2 <= 75),
    n_I2_gt75  = sum(I2 > 75)
  )
readr::write_tsv(I2_summary, file.path(int_dir, "meta_DE_I2_summary.tsv"))
message("Mediana I2: ", round(I2_summary$median_I2, 1), "%")

# k distribution
readr::write_tsv(meta %>% count(k, name = "n_genes"),
                 file.path(int_dir, "meta_DE_k_distribution.tsv"))

# Concordance and high-confidence signatures
readr::write_tsv(concord, file.path(int_dir, "meta_gene_concordance.tsv"))

meta_core_highconf <- meta %>%
  filter(!is.na(FDR_meta), FDR_meta < FDR_CUTOFF,
         !is.na(k), k >= 3,
         !is.na(concordance), concordance == 1,
         !is.na(I2), I2 <= 50)
meta_sig_hetero <- meta %>%
  filter(!is.na(FDR_meta), FDR_meta < FDR_CUTOFF, !is.na(I2), I2 > 75)

readr::write_tsv(meta_core_highconf, file.path(int_dir, "meta_core_highconf.tsv"))
readr::write_tsv(meta_sig_hetero,    file.path(int_dir, "meta_sig_heterogeneous.tsv"))

message(">> Meta-análisis y exportación COMPLETADOS.")
stop_log(log_handle)
