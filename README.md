# Olaparib_resistance

Pipeline de bioinformática para analizar resistencia a olaparib en líneas celulares
(A2780, PEO1 y UWB1.289 BRCA-def/prof). Incluye expresión diferencial, GSEA de
Hallmarks, GSVA y metaanálisis con reportes reproducibles.

## Estructura del repositorio
- data/raw: datos originales por dataset
- data/processed: datos derivados intermedios
- scripts: scripts principales en R
- notebooks: análisis interactivo (si aplica)
- results: tablas y figuras generadas
- env: entorno reproducible (conda)
- logs: salidas de ejecución
- docs: documentación de reproducibilidad

## Requisitos
- Conda (recomendado) con `env/environment.yml`
- Rscript disponible en el PATH

## Quickstart
1. Crear el entorno:
   - `conda env create -f env/environment.yml`
   - `conda activate proyecto_env`
2. (Opcional) definir ruta del proyecto:
   - `export OLAPARIB_RESISTANCE_DIR=/home/jcarvajalv/bioinfo/projects/olaparib_resistance`
3. Ejecutar el flujo completo siguiendo `docs/RUNBOOK.md`.

## Documentación
- `docs/RUNBOOK.md`: pasos reproducibles
- `docs/OUTPUTS.md`: mapa de salidas
- `docs/REPRODUCIBILITY.md`: checklist y troubleshooting
- `docs/METHODS.md`: métodos redactados para manuscrito
