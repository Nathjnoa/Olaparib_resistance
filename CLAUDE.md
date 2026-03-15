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
- **set.seed(42)** in scripts 01–04 (after start_log); **set.seed(1)** in 05_3_meta_gsea.R for fgsea reproducibility.
- **sessionInfo()** printed automatically by stop_log() in 00_config.R — captured in every script's log.

## Plan del Manuscrito

Estado actual: análisis completo, en etapa de construcción de figuras de publicación.
Figura en progreso: **Figura 2** (próxima a construir).

### Narrativa del artículo

La historia del paper sigue esta secuencia lógica:

1. **Fig 1** — Validación individual: ¿tienen señal los 3 datasets por separado?
2. **Fig 2** — Solapamiento: ¿qué genes son compartidos entre datasets?
3. **Fig 3** — Pathways por dataset: ¿qué vías están consistentemente alteradas (GSEA Hallmarks)?
4. **Fig 4** — GSVA: actividad de pathways a nivel de muestra individual (Fig3A) e integrada (Fig3B)
5. **Fig 5** — Meta-análisis: síntesis cuantitativa de la evidencia génica entre los 3 estudios
6. **Fig 6** — Meta-GSEA: pathways enriquecidos a partir del ranking meta-analítico

### Figura 1 — COMPLETADA

**Script:** `scripts/06_fig1_study_design.R`
**Output:** `results/figures/fig1/Fig1.pdf` + `Fig1.png`
**Pregunta:** *¿Son los 3 datasets individualmente válidos y muestran evidencia transcriptómica de resistencia a olaparib?*
**Paneles:** A = diagrama de diseño | B = 3 PCAs | C = 3 volcanos | D = conteo de DEGs por dataset
**Decisiones de estilo:** double_col 200mm, leyendas fuera a la derecha del último panel de cada fila, sin columna "Method" en Panel A (va en la sección Methods del paper). Ver lógica completa en `docs/GUIDE.md` → sección "Figura 1".

### Figura 2 — PENDIENTE (próxima sesión)

**Script a crear:** `scripts/07_fig2_overlap.R`
**Pregunta:** *¿Existe una firma transcriptómica compartida entre las 3 líneas, o los cambios son mayoritariamente específicos de cada contexto?*

**Paneles planificados:**

| Panel | Fuente | Qué muestra |
|-------|--------|-------------|
| A | `results/figures/upset/DEGs_UP_UpSet_A2780_UWB_PEO1.png` (re-generar con tema unificado) | Estructura de solapamiento UP: exclusivos de 1 dataset, pares, triple |
| B | `results/figures/upset/DEGs_DOWN_UpSet_A2780_UWB_PEO1.png` (re-generar) | Ídem DOWN |
| C | **NUEVO — generar** | Heatmap ComplexHeatmap de los ~65 genes de la triple intersección, con todas las muestras de los 3 datasets, anotadas por dataset + condición |

**Datos disponibles para Panel C:**

- Genes: `results/tables/venn/DEGs_UP_triple_intersection_A2780_UWB_PEO1.tsv` + `DEGs_DOWN_triple_intersection_A2780_UWB_PEO1.tsv`
- Expresión: `data/processed/GSE153867_A2780_FPKM.tsv` (A2780); re-leer raw para PEO1 y UWB (misma lógica que `01_differential_expression.R`)

**Estilo:** mismo `theme_pub()` de Fig1, mismo preset double_col 200mm, leyendas fuera a la derecha. Los heatmaps individuales por dataset (script `02_2`) van a **suplementario**, no a figura principal.

**Nota UpSet:** el UpSet plots actualmente usan el paquete `UpSetR`. Re-generar con `ComplexHeatmap::upset()` o `ggupset` para tener control total sobre estilo y que sea consistente con el resto de las figuras.

### Figuras 3–6 — Por planificar en sesiones futuras

Fuentes disponibles (ya generadas por el pipeline):

- Fig 3 (GSEA Hallmarks): `results/figures/gsea_hallmarks/` — lollipops por dataset + heatmap NES integrado
- Fig 4 (GSVA): `results/figures/gsva_heatmaps/Fig3A_*` + `Fig3B_*`
- Fig 5 (Meta-análisis): `results/figures/meta_analysis/` — MetaVolcano + MetaLollipop + forest plots
- Fig 6 (Meta-GSEA): `results/figures/meta_gsea/` — lollipops Hallmarks + Reactome

## Optional environment variable

```bash
export OLAPARIB_RESISTANCE_DIR=/home/jcarvajalv/bioinfo/projects/olaparib_resistance
```
