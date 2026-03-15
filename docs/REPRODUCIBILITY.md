# REPRODUCIBILITY

Guía para asegurar reproducibilidad y diagnosticar discrepancias.

## Entorno

- Usar el entorno conda definido en `env/environment.yml`.
- Activar con `conda activate omics-R` antes de ejecutar cualquier script.
- Registrar versión de R y paquetes: en R ejecutar `sessionInfo()` y guardar el output en `logs/`.

## Determinismo — estado actual

**⚠️ Advertencia conocida:** solo `05_3_meta_gsea.R` fija `set.seed(1)`. Los demás scripts
que generan PCA, clustering de heatmaps (ComplexHeatmap) y puntuaciones GSVA no tienen seed fijo.
Las tablas de resultados son deterministas (no dependen de aleatoriedad). Las figuras con
clustering o layout pueden variar levemente entre ejecuciones en diferentes plataformas.

Esto no afecta la validez científica de los resultados tabulares, pero sí la reproducibilidad
exacta de las figuras.

## Checklist: qué revisar si los resultados difieren

- [ ] Versiones de R/Bioconductor y paquetes (correr `sessionInfo()` y comparar).
- [ ] Nombres y contenidos de los archivos en `data/raw/`.
- [ ] Variables de entorno (`OLAPARIB_RESISTANCE_DIR` si se usa).
- [ ] Orden de ejecución de scripts (ver `docs/RUNBOOK.md`).
- [ ] Cambios en thresholds (`FDR_CUTOFF`, `LOGFC_CUTOFF`) en `scripts/00_config.R`.
- [ ] Que los archivos en `results/tables/de/` existan antes de correr scripts 02+.

## Manifest de corrida (recomendado)

Guarda un manifest por corrida con inputs, parámetros y outputs.

Ejemplo (JSON):

```json
{
  "run_id": "YYYYMMDD-HHMMSS",
  "inputs": [
    "data/raw/GSE153867/GSE153867_fpkm.txt",
    "data/raw/GSE235980/GSE235980_CountReads.txt.gz"
  ],
  "parameters": {
    "fdr": 0.05,
    "log2fc_cutoff": 0.584,
    "meta_stat": "z"
  },
  "outputs": [
    "results/tables/de/GSE153867_A2780_limma_DE_FPKM.tsv",
    "results/tables/integrated/meta_DE_random_effects.tsv"
  ],
  "environment": {
    "r_version": "{sessionInfo()$R.version$version.string}",
    "conda_env": "omics-R"
  }
}
```

## Problemas conocidos

| Problema | Impacto | Estado |
| --- | --- | --- |
| Seeds no fijados en 01–04 | Figuras no exactamente reproducibles | Abierto |
| SE limma sin validación de cero | Genes con t≈0 descartados silenciosamente | Abierto |
| `sessionInfo()` no exportado automáticamente | Versiones de paquetes no registradas | Abierto |
