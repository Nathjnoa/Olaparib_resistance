# OUTPUTS

Mapa de salidas generadas por el pipeline. Los nombres pueden incluir sufijos o variar segun parametros.

## Tablas principales (`results/tables/`)
- `GSE153867_A2780_limma_DE_FPKM.tsv`: DE limma (A2780).
- `GSE117765_PEO1_limma_TPM_Res_vs_Par.tsv`: DE limma (PEO1).
- `GSE235980_BRCAdef_Res_vs_Par_DESeq2.tsv`: DESeq2 (UWB BRCA-def).
- `GSE235980_BRCAprof_Res_vs_Par_DESeq2.tsv`: DESeq2 (UWB BRCA-prof).
- `DEG_counts.tsv` / `venn/*.tsv`: conteos e intersecciones de DEGs.

## GSEA Hallmarks (`results/tables/gsea_hallmarks/`)
- `*_GSEA_Hallmarks.tsv`: resultados por dataset.

## GSVA (`results/tables/gsva/` y `results/tables/gsva_DE/`)
- `GSVA_*.tsv`: scores GSVA por dataset y merged.
- `GSVA_metadata.tsv`: metadata de muestras.
- `GSVA_DE_*.tsv`: diferenciales por dataset y tabla integrada.
- `GSVA_core_pathways_summary.tsv`: resumen de pathways recurrentes.

## Metaanalisis (`results/tables/integrated/`)
- `meta_DE_input_long.tsv` y `meta_DE_input_wide.tsv`: insumos del metaanalisis.
- `meta_DE_random_effects.tsv`: resultados RE.
- `meta_DE_meta_up.tsv` / `meta_DE_meta_down.tsv`: firmas por direccion.
- `meta_DE_meta_rank_*.tsv`: rankings.
- `meta_DE_I2_summary.tsv`, `meta_DE_k_distribution.tsv`, `meta_gene_concordance.tsv`.

## Figuras (`results/figures/`)
- `*_volcano.png`, `*_PCA.png`: DE por dataset.
- `venn/*.png`, `upset/*.png`: intersecciones de DEGs.
- `heatmaps/*.png`: heatmaps top DEGs.
- `gsea_hallmarks/*`: lollipop, emap, y heatmap integrado.
- `gsva_heatmaps/*`: Fig3A/ Fig3B.
- `meta_analysis/*`: volcano y lollipop meta.
- `meta/*`: forest plot meta.

## Logs
- Archivos `*.log` en el root (p.ej. VennDiagram) reflejan ejecuciones previas.

