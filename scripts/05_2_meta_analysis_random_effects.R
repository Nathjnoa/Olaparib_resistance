#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(metafor)
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
in_file  <- file.path(int_dir, "meta_DE_input_long.tsv")
out_file <- file.path(int_dir, "meta_DE_random_effects.tsv")

if (!file.exists(in_file)) stop("No se encontró archivo de entrada: ", in_file)

meta_long <- readr::read_tsv(in_file, show_col_types = FALSE)

message("Filas en meta_long: ", nrow(meta_long))
message("Genes únicos: ", dplyr::n_distinct(meta_long$gene_symbol))
message("Datasets encontrados: ", paste(sort(unique(meta_long$dataset)), collapse = ", "))

# ------------------------------------------------------------
# Metaanálisis por gen
# ------------------------------------------------------------
run_meta_for_gene <- function(df_gene) {
  df_ok <- df_gene %>%
    filter(
      !is.na(log2FC),
      is.finite(log2FC),
      !is.na(SE),
      is.finite(SE),
      SE > 0
    )

  k <- nrow(df_ok)
  ds_list <- paste(sort(unique(df_ok$dataset)), collapse = ",")

  # salida vacía si k<2
  if (k < 2) {
    return(tibble(
      k = k,
      datasets = ds_list,
      reml_warning = FALSE,
      beta_meta = NA_real_, se_meta = NA_real_, z_meta = NA_real_, p_meta = NA_real_,
      ci_lb = NA_real_, ci_ub = NA_real_,
      tau2 = NA_real_, Q = NA_real_, Q_pval = NA_real_, I2 = NA_real_,
      beta_fixed = NA_real_, se_fixed = NA_real_, z_fixed = NA_real_, p_fixed = NA_real_,
      beta_meta_DL = NA_real_, se_meta_DL = NA_real_, p_meta_DL = NA_real_, tau2_DL = NA_real_
    ))
  }

  reml_warn_flag <- FALSE

  # Random effects REML con captura de warning Fisher scoring
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
    fit_dl <- tryCatch(
      rma.uni(yi = df_ok$log2FC, sei = df_ok$SE, method = "DL"),
      error = function(e) NULL
    )
  }

  # Fixed effects (sensibilidad)
  fit_fe <- tryCatch(
    rma.uni(yi = df_ok$log2FC, sei = df_ok$SE, method = "FE"),
    error = function(e) NULL
  )

  tibble(
    k = k,
    datasets = ds_list,

    reml_warning = reml_warn_flag,

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

message(">> Iniciando metaanálisis por gen...")

meta_results <- meta_long %>%
  group_by(gene_symbol) %>%
  group_modify(~ run_meta_for_gene(.x)) %>%
  ungroup() %>%
  mutate(
    reml_warning = as.logical(reml_warning),
    across(c(beta_meta_DL, se_meta_DL, p_meta_DL, tau2_DL), as.numeric)
  ) %>%
  mutate(
    FDR_meta = p.adjust(p_meta, method = "BH"),
    direction = case_when(
      !is.na(beta_meta) & beta_meta > 0 ~ "UP",
      !is.na(beta_meta) & beta_meta < 0 ~ "DOWN",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(FDR_meta, p_meta)

readr::write_tsv(meta_results, out_file)
message(">> Metaanálisis completado. Resultados en: ", out_file)

# Exportar genes con warning en REML
warn_genes <- meta_results %>% filter(isTRUE(reml_warning))
if (nrow(warn_genes) > 0) {
  out_warn <- file.path(int_dir, "meta_DE_REML_warnings_genes.tsv")
  readr::write_tsv(warn_genes, out_warn)
  message(">> Genes con warning REML guardados en: ", out_warn)
} else {
  message(">> No hubo warnings REML.")
}

# QC extra: top hits
message("\nTop 10 por FDR_meta:")
print(meta_results %>% select(gene_symbol, k, beta_meta, FDR_meta, I2, datasets) %>% head(10))

warnings()
