#!/usr/bin/env Rscript

# ============================================================
# 220_export_meta_results.R
#
# Exporta tablas derivadas del metaanálisis:
#  - Lista de genes meta-UP y meta-DOWN (FDR_meta < 0.05)
#  - Ranking completo por p_meta
#  - Resumen de I2
#  - Distribución de k (nº datasets por gen)
# ============================================================

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
in_file  <- file.path(int_dir, "meta_DE_random_effects.tsv")
in_long  <- file.path(int_dir, "meta_DE_input_long.tsv")

if (!file.exists(in_file)) {
  stop("No se encontró el archivo de metaanálisis: ", in_file)
}

if (!file.exists(in_long)) {
  stop("No se encontró el archivo de entrada long: ", in_long)
}

meta <- readr::read_tsv(in_file, show_col_types = FALSE)
meta_long <- readr::read_tsv(in_long, show_col_types = FALSE)

message("Filas en meta: ", nrow(meta))
message("Genes únicos: ", dplyr::n_distinct(meta$gene_symbol))

# ------------------------------------------------------------
# 1) Listas UP / DOWN (FDR_meta < 0.05)
# ------------------------------------------------------------
thr_FDR <- 0.05
thr_log2FC <- 0 # si luego quieres exigir magnitud, cámbialo aquí

# Manejo de p_meta == 0 para BH estable y cálculo de z_meta si falta
meta <- meta %>%
  mutate(
    p_meta2 = if_else(!is.na(p_meta) & p_meta == 0, 1e-300, p_meta),
    FDR_meta = p.adjust(p_meta2, method = "BH"),
    z_meta = if ("z_meta" %in% colnames(meta)) z_meta else ifelse(is.na(beta_meta) | is.na(se_meta) | se_meta == 0, NA_real_, beta_meta / se_meta)
  )

# Concordancia por gen desde meta_long
concord <- meta_long %>%
  filter(!is.na(log2FC)) %>%
  group_by(gene_symbol) %>%
  summarise(
    k_eff = n(),
    n_pos = sum(log2FC > 0, na.rm = TRUE),
    n_neg = sum(log2FC < 0, na.rm = TRUE),
    concordance = if_else(k_eff > 0, pmax(n_pos, n_neg) / k_eff, NA_real_),
    .groups = "drop"
  )

meta <- meta %>% left_join(concord, by = "gene_symbol")

meta_up <- meta %>%
  filter(
    !is.na(FDR_meta),
    FDR_meta < thr_FDR,
    !is.na(beta_meta),
    beta_meta > thr_log2FC
  ) %>%
  arrange(FDR_meta)

meta_down <- meta %>%
  filter(
    !is.na(FDR_meta),
    FDR_meta < thr_FDR,
    !is.na(beta_meta),
    beta_meta < -thr_log2FC
  ) %>%
  arrange(FDR_meta)

out_up   <- file.path(int_dir, "meta_DE_meta_up.tsv")
out_down <- file.path(int_dir, "meta_DE_meta_down.tsv")

readr::write_tsv(meta_up,   out_up)
readr::write_tsv(meta_down, out_down)

message("Genes UP (FDR<", thr_FDR, "):   ", nrow(meta_up))
message("Genes DOWN (FDR<", thr_FDR, "): ", nrow(meta_down))
message(">> Listas guardadas en:")
message("   ", out_up)
message("   ", out_down)

# ------------------------------------------------------------
# 2) Rankings globales
# ------------------------------------------------------------

meta_rank_z <- meta %>%
  filter(!is.na(z_meta)) %>%
  arrange(desc(z_meta))

meta_rank_beta <- meta %>%
  filter(!is.na(beta_meta)) %>%
  arrange(desc(beta_meta))

out_rank_z <- file.path(int_dir, "meta_DE_meta_rank_z.tsv")
out_rank_beta <- file.path(int_dir, "meta_DE_meta_rank_beta.tsv")
readr::write_tsv(meta_rank_z, out_rank_z)
readr::write_tsv(meta_rank_beta, out_rank_beta)
message(">> Rankings guardados en: ", out_rank_z, " y ", out_rank_beta)

# ------------------------------------------------------------
# 3) Resumen de I2
# ------------------------------------------------------------

I2_summary <- meta %>%
  filter(!is.na(I2)) %>%
  summarise(
    n_genes          = n(),
    mean_I2          = mean(I2),
    median_I2        = median(I2),
    I2_q05           = quantile(I2, 0.05),
    I2_q25           = quantile(I2, 0.25),
    I2_q75           = quantile(I2, 0.75),
    I2_q95           = quantile(I2, 0.95),
    n_I2_lt25        = sum(I2 < 25),
    n_I2_25_75       = sum(I2 >= 25 & I2 <= 75),
    n_I2_gt75        = sum(I2 > 75)
  )

out_I2 <- file.path(int_dir, "meta_DE_I2_summary.tsv")
readr::write_tsv(I2_summary, out_I2)
message(">> Resumen de I2 guardado en: ", out_I2)

# ------------------------------------------------------------
# 4) Distribución de k (nº datasets)
# ------------------------------------------------------------

k_dist <- meta %>%
  count(k, name = "n_genes") %>%
  arrange(k)

out_k <- file.path(int_dir, "meta_DE_k_distribution.tsv")
readr::write_tsv(k_dist, out_k)
message(">> Distribución de k guardada en: ", out_k)

# ------------------------------------------------------------
# 5) Concordancia y firmas adicionales
# ------------------------------------------------------------
out_conc <- file.path(int_dir, "meta_gene_concordance.tsv")
readr::write_tsv(concord, out_conc)
message(">> Concordancia guardada en: ", out_conc)

meta_core_highconf <- meta %>%
  filter(
    !is.na(FDR_meta), FDR_meta < 0.05,
    !is.na(k), k >= 3,
    !is.na(concordance), concordance == 1,
    !is.na(I2), I2 <= 50
  )

meta_sig_hetero <- meta %>%
  filter(!is.na(FDR_meta), FDR_meta < 0.05, !is.na(I2), I2 > 75)

out_core <- file.path(int_dir, "meta_core_highconf.tsv")
out_sighet <- file.path(int_dir, "meta_sig_heterogeneous.tsv")
readr::write_tsv(meta_core_highconf, out_core)
readr::write_tsv(meta_sig_hetero, out_sighet)
message(">> Firmas adicionales guardadas en: ", out_core, " y ", out_sighet)

message(">> 220_export_meta_results.R COMPLETADO.")
problems()
