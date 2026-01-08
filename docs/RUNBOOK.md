# RUNBOOK

Guia paso a paso para reproducir el flujo completo del proyecto.

## Requisitos
- Conda con el entorno de `env/environment.yml`.
- Rscript disponible.
- Datos crudos ubicados en `data/raw/` con los nombres esperados.

## Preparacion del entorno
1. Crear y activar el entorno:
   - `conda env create -f env/environment.yml`
   - `conda activate proyecto_env`
2. (Opcional) definir la ruta del proyecto:
   - `export OLAPARIB_RESISTANCE_DIR=/home/jcarvajalv/bioinfo/projects/olaparib_resistance`

## Ejecucion del flujo
Ejecuta los scripts en este orden para reproducir tablas y figuras:

1. Expresion diferencial por dataset (limma/DESeq2) y figuras basicas:
   - `Rscript scripts/01_1_run_olaparib_all.R`

2. Intersecciones de DEGs y heatmaps por dataset:
   - `Rscript scripts/02_2_venn_DEGs_up_down.R`
   - `Rscript scripts/02_2_UpSet.R`
   - `Rscript scripts/02_3_heatmaps_ComplexHeatmap_topDEGs.R`

3. GSEA Hallmarks por dataset e integracion:
   - `Rscript scripts/03_1_gsea_hallmarks.R`
   - `Rscript scripts/03_2_gsea_hallmarks_heatmap_integrated.R`

4. GSVA Hallmarks y analisis diferencial de GSVA:
   - `Rscript scripts/04_1_GSVA.R`
   - `Rscript scripts/04_2_GSVA_limma_DE.R`
   - `Rscript scripts/04_3_gsva_heatmap_per_sample_Fig3A.R`
   - `Rscript scripts/04_4_gsva_logFC_heatmap_integrated_Fig3B.R`

5. Metaanalisis y reportes:
   - `Rscript scripts/05_1_build_meta_input.R`
   - `Rscript scripts/05_2_meta_analysis_random_effects.R`
   - `Rscript scripts/05_3_export_meta_results.R`
   - `Rscript scripts/05_4_meta_plots_volcano_lollipop.R`
     - Opcional: `Rscript scripts/05_4_meta_plots_volcano_lollipop.R 0.05 10 0.5 12`
   - `Rscript scripts/05_5_metaGSEA_hallmarks_reactome.R`
     - Opcional: `Rscript scripts/05_5_metaGSEA_hallmarks_reactome.R z 0.10 10 1`
   - `Rscript scripts/05_6_meta_forest_plots.R`

## Checklist: como reproducir esta corrida
- [ ] Verificar que los archivos en `data/raw/` existan y tengan el nombre esperado.
- [ ] Activar el entorno `proyecto_env`.
- [ ] Exportar `OLAPARIB_RESISTANCE_DIR` si se ejecuta fuera del root del repo.
- [ ] Ejecutar los scripts en el orden del runbook.
- [ ] Confirmar outputs en `results/` (ver `docs/OUTPUTS.md`).

