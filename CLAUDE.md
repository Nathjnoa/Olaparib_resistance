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

### Figura 2 — COMPLETADA

**Script:** `scripts/07_fig2_overlap.R`
**Output:** `results/figures/fig2/Fig2.pdf` + `Fig2.png`
**Pregunta:** *¿Existe una firma transcriptómica compartida entre las 3 líneas, o los cambios son mayoritariamente específicos de cada contexto?*
**Paneles:** A = UpSet DEGs UP | B = UpSet DEGs DOWN | C = Heatmap 65 genes triple intersección (26 muestras, Z-score within-dataset)
**Decisiones de estilo:** double_col 200×260 mm, UpSet con `ComplexHeatmap::upset()` (colores por dataset en set names), heatmap con anotaciones Dataset + Condition + BRCA status. Ver lógica completa en `docs/GUIDE.md` → sección "Figura 2".

### Figura 3 — COMPLETADA

**Script:** `scripts/08_fig3_hallmarks.R`
**Output:** `results/figures/fig3/Fig3.pdf` + `Fig3.png`
**Pregunta:** *¿Convergen los cambios transcriptómicos en los mismos programas biológicos a través de los 3 modelos?*
**Paneles:** A/B/C = lollipops GSEA Hallmarks por dataset (top 5 por dirección, FDR<0.05, ejes X unificados) | D = heatmap NES integrado transpuesto (datasets=filas, pathways=columnas), 17 Hallmarks con FDR<0.1 en ≥2 datasets, asteriscos de significancia, anotación BRCA status
**Decisiones de estilo:** double_col 200×190 mm, pathways abreviados, heatmap transpuesto para aprovechar ancho horizontal. Ver lógica completa en `docs/GUIDE.md` → sección "Figura 3".

### Figura 4 — COMPLETADA

**Script:** `scripts/09_fig4_gsva.R`
**Output:** `results/figures/fig4/Fig4.pdf` + `Fig4.png`
**Pregunta:** *¿La activación diferencial de pathways observada en GSEA (Fig 3) es consistente a nivel de muestra individual, y se confirma con el marco estadístico independiente GSVA + limma?*
**Paneles:** A = heatmap per-sample Z(GSVA) (26 muestras, anotaciones Dataset/Condition/BRCA, sin column names, columnas sin clustering, filas con hclust) | B = ΔGSVA logFC integrado (3 datasets como columnas, mismos pathways, asteriscos FDR<0.05, puntos FDR<0.10)
**Decisiones de estilo:** double_col 200×150 mm, layout 65%/35% viewports, row order pre-computado con hclust para sincronizar A y B. Ver lógica en `docs/GUIDE.md` → sección "Figura 4".
**Hallazgo clave visible:** IFN-α e IFN-γ responden en direcciones opuestas (UP en A2780 *, DOWN en UWB1.289 *) — heterogeneidad de dirección directamente visible. PEO1 muestra logFC cercanos a 0 para la mayoría de pathways recurrentes (los convergentes son principalmente A2780 + UWB).

### Figura 5 — COMPLETADA

**Script:** `scripts/10_fig5_meta.R`
**Output:** `results/figures/fig5/Fig5.pdf` + `Fig5.png`
**Pregunta:** *¿Cuál es la firma transcriptómica integrada de resistencia a olaparib a través de los 3 modelos, y qué programas biológicos implica?*
**Paneles:** A = MetaVolcano (beta_meta vs -log10(FDR), etiquetados los mismos 20 genes del lollipop para coherencia visual) | B = Lollipop top genes high-confidence (k=3, I²<50%, n=361 total, mostrados top 10 UP + 10 DOWN) | C = Meta-GSEA Hallmarks full width (ranking beta_meta, 8 pathways FDR<0.05)
**Decisiones clave:**
- Meta-GSEA usa ranking **beta_meta** (NO z_meta): z solo da 1 pathway (KRAS_DN, FDR=0.01) por alta heterogeneidad (I²=96%); beta da 8 pathways (FDR<0.05) que capturan la dirección promedio del efecto
- Labels del volcano = mismos 20 genes del lollipop → el lector puede cruzar información entre paneles A y B
- Forest plots → Supplementary (no en figura principal)
**Hallazgos clave:**
- n=613 UP / 613 DOWN genes (FDR<0.05); 361 high-confidence (k=3, I²<50%)
- Panel C: Estrogen response, KRAS signaling DN, Xenobiotic metabolism, E2F/P53 UP en resistentes; EMT e IL2-STAT5 DOWN (dirección ponderada — heterogeneidad alta entre modelos)
- Panel B DOWN más destacados: LTBP1, FAM167A, USP18, DNAJB11, BAMBI; UP: URB1, JUN, NEDD9

Fuentes disponibles (ya generadas por el pipeline):

- Fig 4 (GSVA): `results/figures/gsva_heatmaps/Fig3A_*` + `Fig3B_*`
- Fig 5 (Meta-análisis): `results/figures/meta_analysis/` — MetaVolcano + MetaLollipop + forest plots
- Fig 6 (Meta-GSEA): `results/figures/meta_gsea/` — lollipops Hallmarks + Reactome

## Optional environment variable

```bash
export OLAPARIB_RESISTANCE_DIR=/home/jcarvajalv/bioinfo/projects/olaparib_resistance
```
