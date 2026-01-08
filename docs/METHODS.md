# METHODS

## Diseno del estudio y datos
Se analizaron datasets publicos de expresion para resistencia a olaparib en lineas celulares, incluyendo A2780 (FPKM), PEO1 (TPM) y UWB1.289 con estados BRCA1 deficiente y proficiente (counts). Los archivos crudos se organizaron por dataset en `data/raw/` y se procesaron de manera consistente para comparaciones Resistant vs Parental. Las fechas de descarga y fuentes exactas deben registrarse como `{FECHA_DESCARGA}` y `{FUENTE_DATOS}`.

## Preprocesamiento y control de calidad
Para matrices FPKM/TPM se aplico transformacion log2(x + 1). Se verifico la consistencia de columnas de muestra y se filtraron genes sin varianza antes de los modelos lineales. Para counts se uso un filtrado minimo por conteo y transformacion de estabilidad de varianza (VST) para visualizaciones. Se generaron PCA y graficas volcán para control visual de separacion entre condiciones.

## Expresion diferencial
- A2780 y PEO1: se uso limma con diseno sin intercepto y contraste Resistant vs Parental.
- UWB1.289: se uso DESeq2 con el diseño correspondiente y contrastes especificos por estado BRCA1.
Los resultados se reportaron como log2FC, estadisticos t/z y FDR (BH). Se usaron umbrales FDR < 0.05 para listas principales y FDR < 0.10 para visualizaciones ampliadas.

## Analisis funcional
Se realizo GSEA con Hallmarks de MSigDB para cada dataset a partir de rankings por estadistico de prueba (t o z). Se genero un heatmap integrado de NES para terminos recurrentes en multiples datasets. 

## GSVA y analisis diferencial de GSVA
Se calcularon scores GSVA con colecciones Hallmarks y se integraron por dataset. Luego se aplico limma a los scores GSVA para estimar diferencias Resistant vs Parental, y se generaron heatmaps a nivel de muestra y a nivel integrado de logFC.

## Metaanalisis
Se construyo una tabla long con log2FC y error estandar por gen y dataset, y se aplico un metaanalisis de efectos aleatorios (REML), con estimaciones de heterogeneidad (I2, tau2) y FDR a nivel global. Se exportaron firmas UP/DOWN, rankings y tablas de concordancia.

## Visualizacion
Se generaron volcanos, PCA, Venn/UpSet, heatmaps (ComplexHeatmap), lollipop plots y forest plots. Las figuras se exportaron principalmente en PNG y PDF con tamaños definidos en cada script.

## Software y entorno computacional
- R `{R_VERSION}` con Bioconductor `{BIOC_VERSION}`.
- Paquetes clave: limma `{LIMMA_VERSION}`, DESeq2 `{DESEQ2_VERSION}`, GSVA `{GSVA_VERSION}`, fgsea `{FGSEA_VERSION}`, clusterProfiler `{CLUSTERPROFILER_VERSION}`, ComplexHeatmap `{COMPLEXHEATMAP_VERSION}`.
- El entorno conda se define en `env/environment.yml`.

## Reproducibilidad
El pipeline se ejecuta mediante scripts R reproducibles con rutas controladas por la variable `OLAPARIB_RESISTANCE_DIR`. Se recomienda registrar el entorno y conservar un manifest con inputs, parametros y outputs por corrida.

### Checklist de reproducibilidad
- [ ] Registrar `sessionInfo()` y versiones de paquetes.
- [ ] Guardar los paths exactos de `data/raw/`.
- [ ] Registrar parametros (FDR, log2FC, top_n).
- [ ] Ejecutar los scripts en el orden del `docs/RUNBOOK.md`.

### Template de comandos
```bash
export OLAPARIB_RESISTANCE_DIR=/path/to/olaparib_resistance
Rscript scripts/01_1_run_olaparib_all.R
Rscript scripts/02_2_venn_DEGs_up_down.R
Rscript scripts/02_2_UpSet.R
Rscript scripts/02_3_heatmaps_ComplexHeatmap_topDEGs.R
Rscript scripts/03_1_gsea_hallmarks.R
Rscript scripts/03_2_gsea_hallmarks_heatmap_integrated.R
Rscript scripts/04_1_GSVA.R
Rscript scripts/04_2_GSVA_limma_DE.R
Rscript scripts/04_3_gsva_heatmap_per_sample_Fig3A.R
Rscript scripts/04_4_gsva_logFC_heatmap_integrated_Fig3B.R
Rscript scripts/05_1_build_meta_input.R
Rscript scripts/05_2_meta_analysis_random_effects.R
Rscript scripts/05_3_export_meta_results.R
Rscript scripts/05_4_meta_plots_volcano_lollipop.R
Rscript scripts/05_5_metaGSEA_hallmarks_reactome.R
Rscript scripts/05_6_meta_forest_plots.R
```

### Parametros clave
| Parametro | Significado | Valor |
| --- | --- | --- |
| `FDR` | Umbral para significancia | 0.05 (principal), 0.10 (visualizaciones) |
| `log2FC` | Corte para DEGs en intersecciones | 0.584 (equivalente a FC 1.5) |
| `minGS` / `maxGS` | Tamano de gene set en GSEA/GSVA | 15 / 500 |
| `top_n_plot` | Terminos mostrados en lollipop | 10 (por direccion) |
| `metodo_meta` | Modelo de metaanalisis | RE (REML) |

