# RUNBOOK

Step-by-step guide to reproduce the full pipeline (3-dataset analysis: A2780, PEO1, UWB1.289 BRCA-deficient).

## Requirements

- Conda with `env/environment.yml` (activate as `omics-R`)
- `Rscript` available in PATH
- Raw data in `data/raw/` with expected filenames (see `docs/OUTPUTS.md`)

## Environment setup

```bash
conda env create -f env/environment.yml   # first time only
conda activate omics-R
cd ~/bioinfo/projects/olaparib_resistance

# Optional: set project root explicitly
export OLAPARIB_RESISTANCE_DIR=$(pwd)
```

## Pipeline execution order

### 1. Differential expression per dataset

```bash
Rscript scripts/01_differential_expression.R
```

Outputs: 3 DE tables in `results/tables/de/`, 6 figures in `results/figures/de/` (volcano + PCA per dataset).

### 2. DEG intersections and heatmaps

```bash
Rscript scripts/02_1_deg_intersections.R    # Venn + UpSet (3 datasets)
Rscript scripts/02_2_heatmaps_topDEGs.R    # ComplexHeatmap top 20 UP + 20 DOWN
```

### 3. GSEA Hallmarks per dataset + integrated NES heatmap

```bash
Rscript scripts/03_gsea_hallmarks.R
```

Outputs: per-dataset TSVs + lollipop plots + integrated NES heatmap.

### 4. GSVA + differential pathway activity

```bash
Rscript scripts/04_1_gsva_analysis.R   # GSVA scores + limma DE per dataset
Rscript scripts/04_2_gsva_heatmaps.R   # Fig3A (per-sample) + Fig3B (logFC)
```

### 5. Meta-analysis, plots, and meta-GSEA

```bash
Rscript scripts/05_1_meta_analysis.R          # Build input + REML + export signatures
Rscript scripts/05_2_meta_plots.R             # Volcano + lollipop + forest plots
# Optional args: Rscript scripts/05_2_meta_plots.R 0.05 10 0.5 12
Rscript scripts/05_3_meta_gsea.R              # fGSEA Hallmarks + Reactome
# Optional args: Rscript scripts/05_3_meta_gsea.R z 0.10 10 1
```

## Full one-liner (sequential)

```bash
conda activate omics-R
cd ~/bioinfo/projects/olaparib_resistance
for script in \
    scripts/01_differential_expression.R \
    scripts/02_1_deg_intersections.R \
    scripts/02_2_heatmaps_topDEGs.R \
    scripts/03_gsea_hallmarks.R \
    scripts/04_1_gsva_analysis.R \
    scripts/04_2_gsva_heatmaps.R \
    scripts/05_1_meta_analysis.R \
    scripts/05_2_meta_plots.R \
    scripts/05_3_meta_gsea.R; do
  echo "=== Running $script ==="
  Rscript $script || { echo "FAILED: $script"; exit 1; }
done
```

## Checklist: verifying a successful run

- [ ] `results/tables/de/` contains 3 DE tables (no BRCAprof file)
- [ ] UWB DESeq2 design is `~ resistance` on 4 samples
- [ ] Venn/UpSet show 3 groups (A2780, UWB, PEO1)
- [ ] Heatmaps: 3 panels only
- [ ] GSEA lollipops and NES heatmap: 3 datasets
- [ ] GSVA: 3 datasets in all outputs
- [ ] `meta_DE_input_long.tsv` has max k=3 per gene
- [ ] Forest plots show 3 study estimates + META(RE)
- [ ] All logs go to `logs/` (no VennDiagram.log in root)
- [ ] No `BRCAprof` string in any output filename

## Notes

- BRCA-proficient samples from GSE235980 are excluded from all analyses. See `docs/METHODS.md` for biological rationale.
- All thresholds, sample IDs, and palettes are centralized in `scripts/00_config.R`.
- Logs with timestamps are written to `logs/` by each script.
- Outputs are described in `docs/OUTPUTS.md`.
