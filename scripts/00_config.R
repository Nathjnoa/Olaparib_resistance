#!/usr/bin/env Rscript
# ============================================================
# 00_config.R — Centralized configuration for olaparib_resistance pipeline
#
# Source at the top of every script:
#   source(file.path(dirname(sys.frame(1)$ofile), "00_config.R"))
# or from within the project directory:
#   source("scripts/00_config.R")
# ============================================================

# ============================================================
# Project paths
# ============================================================
proj_dir <- Sys.getenv("OLAPARIB_RESISTANCE_DIR", unset = "")
if (!nzchar(proj_dir)) {
  proj_dir <- "~/bioinfo/projects/olaparib_resistance"
}
proj_dir <- path.expand(proj_dir)
if (!dir.exists(proj_dir)) {
  stop("No existe el directorio del proyecto: ", proj_dir,
       "\nDefine OLAPARIB_RESISTANCE_DIR o ejecuta desde el directorio correcto.")
}

raw_dir       <- file.path(proj_dir, "data", "raw")
processed_dir <- file.path(proj_dir, "data", "processed")
tables_dir    <- file.path(proj_dir, "results", "tables")
figures_dir   <- file.path(proj_dir, "results", "figures")
logs_dir      <- file.path(proj_dir, "logs")

dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(logs_dir,      showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Analysis thresholds
# ============================================================
FDR_CUTOFF     <- 0.05   # primary significance cutoff
FDR_CUTOFF_VIZ <- 0.10   # relaxed cutoff for visualizations / GSEA display
LOGFC_CUTOFF   <- 0.584  # |log2FC| >= 0.584  <=> FC >= 1.5

MIN_GS       <- 15    # minimum gene set size for GSEA/GSVA
MAX_GS       <- 500   # maximum gene set size for GSEA/GSVA
TOP_N_GENES  <- 20    # top genes per direction for heatmaps
TOP_N_PLOT   <- 10    # top pathways per direction in lollipop plots

# ============================================================
# Dataset registry (3 datasets — BRCAprof excluded)
# ============================================================
suppressPackageStartupMessages(library(tibble))

DATASETS <- tibble(
  dataset_id = c("A2780",        "PEO1",                                 "UWB_BRCAdef"),
  gse        = c("GSE153867",    "GSE117765",                            "GSE235980"),
  de_method  = c("limma",        "limma",                                "DESeq2"),
  expr_type  = c("FPKM",         "TPM",                                  "counts"),
  de_file    = c(
    "de/GSE153867_A2780_limma_DE_FPKM.tsv",
    "de/GSE117765_PEO1_limma_TPM_Res_vs_Par.tsv",
    "de/GSE235980_BRCAdef_Res_vs_Par_DESeq2.tsv"
  ),
  n_parental  = c(8L, 2L, 2L),
  n_resistant = c(8L, 4L, 2L)
)

# ============================================================
# Sample ID mappings
# ============================================================

# A2780: parental = O-1..8, resistant = C-1..8
A2780_PAR <- paste0("O-", 1:8)
A2780_RES <- paste0("C-", 1:8)

# UWB1.289 BRCA-deficient: parental first, then resistant
UWB_DEF_PAR <- c("AN14028320", "AN14028322")
UWB_DEF_RES <- c("AN14028316", "AN14028318")

# ============================================================
# Color palettes
# ============================================================

# Volcano plots
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

# Lollipop / direction (UP/DOWN)
pal_dir <- c("UP" = "#D81B60", "DOWN" = "#1E88E5")

# Venn diagram sets (A2780=blue, UWB=green, PEO1=red)
venn_palette <- c(A2780 = "#1f78b4", UWB = "#33a02c", PEO1 = "#e31a1c")

# ============================================================
# Duplicate gene handling policy:
#   - DE tables for meta-analysis: keep row with smallest padj
#   - GSEA rankings: keep entry with max |stat|
#   - GSVA expression input: keep first occurrence (distinct)
# ============================================================

# ============================================================
# Shared helpers
# ============================================================

# Clean Hallmark pathway labels (strip prefix, title-case, fix acronyms)
clean_hallmark <- function(x) {
  y <- x %>%
    stringr::str_replace_all("^HALLMARK_", "") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("\\s+", " ") %>%
    stringr::str_trim() %>%
    stringr::str_to_title()
  y <- y %>%
    stringr::str_replace_all("\\bE2f\\b",         "E2F") %>%
    stringr::str_replace_all("\\bNfkb\\b",        "NFKB") %>%
    stringr::str_replace_all("\\bTnfa\\b",        "TNFA") %>%
    stringr::str_replace_all("\\bIl6 Jak Stat3\\b","IL6 JAK STAT3") %>%
    stringr::str_replace_all("\\bMtorc1\\b",      "MTORC1") %>%
    stringr::str_replace_all("\\bUv\\b",          "UV") %>%
    stringr::str_replace_all("\\bG2m\\b",         "G2M") %>%
    stringr::str_replace_all("\\bMyc\\b",         "MYC")
  y
}

# Collapse duplicate gene symbols: keep entry with max |stat|
collapse_duplicate_symbols <- function(stats_named) {
  tibble::tibble(symbol = names(stats_named), stat = as.numeric(stats_named)) %>%
    dplyr::filter(!is.na(symbol), symbol != "", !is.na(stat), is.finite(stat)) %>%
    dplyr::group_by(symbol) %>%
    dplyr::slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    tibble::deframe()
}

# Map Ensembl IDs to SYMBOL (pass-through if already symbols)
map_to_symbol_if_needed <- function(ids) {
  ids2 <- gsub("\\..*$", "", ids)
  is_ens <- grepl("^ENSG", ids2)
  if (!any(is_ens)) return(ids2)
  sym <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys    = unique(ids2),
    keytype = "ENSEMBL",
    column  = "SYMBOL",
    multiVals = "first"
  )
  unname(sym[ids2])
}

# ============================================================
# Logging helper
# ============================================================
start_log <- function(script_name) {
  ts     <- format(Sys.time(), "%Y%m%d_%H%M%S")
  logfile <- file.path(logs_dir, paste0(script_name, "_", ts, ".log"))
  con    <- file(logfile, open = "wt")
  sink(con, append = FALSE, type = "output")
  sink(con, append = FALSE, type = "message")
  message("=== ", script_name, " started: ", Sys.time(), " ===")
  invisible(list(file = logfile, con = con))
}

stop_log <- function(log_handle) {
  message("=== Done: ", Sys.time(), " ===")
  print(sessionInfo())
  sink(type = "message")
  sink(type = "output")
  close(log_handle$con)
  invisible(log_handle$file)
}
