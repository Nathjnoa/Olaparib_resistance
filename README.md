# olaparib_resistance

R pipeline for transcriptomic meta-analysis of olaparib resistance across 3 ovarian cancer cell line datasets (A2780, PEO1, UWB1.289 BRCA-deficient). Includes differential expression, GSEA Hallmarks, GSVA, and random-effects meta-analysis.

## Datasets

| GEO | Cell line | BRCA1 status | Expression | DE method |
|-----|-----------|--------------|------------|-----------|
| GSE153867 | A2780 | WT | FPKM | limma |
| GSE117765 | PEO1 | BRCA1-mut | TPM | limma |
| GSE235980 | UWB1.289 | BRCA1-def | Counts | DESeq2 |

## Repository structure

```
data/raw/           — original data per dataset (never modified)
scripts/            — 10 numbered R scripts + 00_config.R
results/figures/    — all output figures (PNG + PDF)
results/tables/     — all output tables (TSV)
logs/               — execution logs with timestamps
env/                — conda environment specification
docs/               — RUNBOOK, METHODS, OUTPUTS documentation
```

## Pipeline

```
00_config.R                     — centralized config (thresholds, sample IDs, helpers)
01_differential_expression.R    — DE per dataset (limma/DESeq2)
02_1_deg_intersections.R        — Venn + UpSet (3 datasets)
02_2_heatmaps_topDEGs.R         — ComplexHeatmap top DEGs
03_gsea_hallmarks.R             — GSEA per dataset + integrated NES heatmap
04_1_gsva_analysis.R            — GSVA scores + limma DE
04_2_gsva_heatmaps.R            — Fig3A (per-sample) + Fig3B (logFC)
05_1_meta_analysis.R            — REML RE meta-analysis + signatures
05_2_meta_plots.R               — volcano + lollipop + forest plots
05_3_meta_gsea.R                — fGSEA Hallmarks + Reactome on meta-rankings
```

## Quickstart

```bash
conda activate omics-R
cd ~/bioinfo/projects/olaparib_resistance
export OLAPARIB_RESISTANCE_DIR=$(pwd)

Rscript scripts/01_differential_expression.R
Rscript scripts/02_1_deg_intersections.R
Rscript scripts/02_2_heatmaps_topDEGs.R
Rscript scripts/03_gsea_hallmarks.R
Rscript scripts/04_1_gsva_analysis.R
Rscript scripts/04_2_gsva_heatmaps.R
Rscript scripts/05_1_meta_analysis.R
Rscript scripts/05_2_meta_plots.R
Rscript scripts/05_3_meta_gsea.R
```

## Documentation

- `docs/RUNBOOK.md` — step-by-step execution guide
- `docs/OUTPUTS.md` — map of all output files
- `docs/METHODS.md` — methods section for manuscript
