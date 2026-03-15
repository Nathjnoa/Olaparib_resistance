# olaparib_resistance

Random-effects meta-analysis of transcriptomic olaparib resistance across 3 GEO datasets
(A2780, PEO1, UWB1.289 BRCA-deficient). Differential expression, GSEA Hallmarks, GSVA,
and gene-level meta-analysis with forest plots and meta-GSEA.

## Environment

```bash
conda activate omics-R
cd ~/bioinfo/projects/olaparib_resistance
```

## Pipeline (10 scripts + config)

```bash
# Config: sourced by all scripts (thresholds, sample IDs, helpers)
# scripts/00_config.R

# Run in order:
Rscript scripts/01_differential_expression.R
Rscript scripts/02_1_deg_intersections.R
Rscript scripts/02_2_heatmaps_topDEGs.R
Rscript scripts/03_gsea_hallmarks.R
Rscript scripts/04_1_gsva_analysis.R
Rscript scripts/04_2_gsva_heatmaps.R
Rscript scripts/05_1_meta_analysis.R
Rscript scripts/05_2_meta_plots.R        # args optional: 0.05 10 0.5 12
Rscript scripts/05_3_meta_gsea.R         # args optional: z 0.10 10 1
```

## Script descriptions

| Script | Description |
|--------|-------------|
| `00_config.R` | Centralized config: paths, thresholds, sample IDs, color palettes, helper functions |
| `01_differential_expression.R` | limma (A2780, PEO1) + DESeq2 (UWB BRCA-def); volcanos, PCA |
| `02_1_deg_intersections.R` | Venn diagrams + UpSet plots across 3 datasets |
| `02_2_heatmaps_topDEGs.R` | ComplexHeatmap: top 20 UP + 20 DOWN DEGs per dataset |
| `03_gsea_hallmarks.R` | GSEA Hallmarks per dataset + integrated NES heatmap |
| `04_1_gsva_analysis.R` | GSVA scores + limma DE on pathway scores |
| `04_2_gsva_heatmaps.R` | Fig3A (per-sample Z-score) + Fig3B (integrated logFC) |
| `05_1_meta_analysis.R` | Build input + REML RE meta-analysis + export signatures |
| `05_2_meta_plots.R` | Meta volcano + gene lollipop + forest plots |
| `05_3_meta_gsea.R` | fGSEA Hallmarks + Reactome on meta-rankings |

## Key paths

- `data/raw/` — raw data per dataset (never modify)
- `results/figures/` — all output figures
- `results/tables/` — all output tables
- `logs/` — execution logs with timestamps
- `docs/RUNBOOK.md` — step-by-step execution guide
- `docs/METHODS.md` — methods section for manuscript
- `docs/OUTPUTS.md` — map of all output files

## Design decisions

- **3 datasets only**: A2780 (GSE153867), PEO1 (GSE117765), UWB BRCA-def (GSE235980).
  BRCA-proficient samples from GSE235980 are excluded (different resistance mechanism).
- **limma on log2(FPKM+1) and log2(TPM+1)**: correct for pre-normalized data. voom is NOT used
  (voom is for raw counts only; Ritchie et al. 2015, Law et al. 2014).
- **DESeq2 `~ resistance`** for UWB: 4 BRCA-def samples only; BRCA_status is constant.
- **SE approximation**: SE = |log2FC| / |t| (exact for limma t-statistic).
- **REML meta-analysis**: DL fallback if Fisher scoring fails to converge.
- **Duplicate gene policy**: smallest padj (DE/meta), max |stat| (GSEA rankings).
- **set.seed(1)** in 05_3_meta_gsea.R for fgsea reproducibility.

## Optional environment variable

```bash
export OLAPARIB_RESISTANCE_DIR=/home/jcarvajalv/bioinfo/projects/olaparib_resistance
```
