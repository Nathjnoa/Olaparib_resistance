# REPRODUCIBILITY

Guia para asegurar reproducibilidad y diagnosticar discrepancias.

## Entorno
- Usar el entorno conda definido en `env/environment.yml`.
- Registrar version de R y paquetes:
  - En R: `sessionInfo()`
- Registrar version de Python (si aplica): `python --version` y `pip freeze`/`conda list`.

## Determinismo
- Los scripts estan diseñados para ser deterministas bajo los mismos inputs.
- Si se actualizan paquetes (limma/DESeq2/GSVA/fgsea), los resultados pueden variar levemente.

## Manifest de corrida (recomendado)
Guarda un manifest por corrida con inputs, parametros y outputs.

Ejemplo (JSON):
```json
{
  "run_id": "YYYYMMDD-HHMMSS",
  "inputs": ["data/raw/GSE153867/GSE153867_fpkm.txt", "data/raw/GSE235980/GSE235980_CountReads.txt.gz"],
  "parameters": {"fdr": 0.05, "log2fc_cutoff": 0.584},
  "outputs": ["results/tables/GSE153867_A2780_limma_DE_FPKM.tsv"],
  "environment": {"r": "{R_VERSION}", "conda": "proyecto_env"}
}
```

## Checklist: que revisar si los resultados difieren
- [ ] Versiones de R/Bioconductor y paquetes.
- [ ] Nombres y contenidos de los archivos en `data/raw/`.
- [ ] Variables de entorno (p.ej. `OLAPARIB_RESISTANCE_DIR`).
- [ ] Orden de ejecucion de scripts.
- [ ] Cambios en thresholds (FDR/log2FC) usados en graficas.

