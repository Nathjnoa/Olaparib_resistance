# METHODS

## Study design and datasets

Three publicly available transcriptomic datasets of olaparib resistance in ovarian cancer cell lines were analyzed:

| Dataset | Cell line | BRCA1 status | Expression type | DE method | Parental (n) | Resistant (n) |
|---------|-----------|--------------|-----------------|-----------|-------------|---------------|
| GSE153867 | A2780 | WT | FPKM | limma | 8 | 8 |
| GSE117765 | PEO1 | BRCA1-mut | TPM | limma | 2 | 4 clones |
| GSE235980 | UWB1.289 | BRCA1-def | Counts | DESeq2 | 2 | 2 |

All comparisons contrast Resistant vs Parental (olaparib-resistant subline/clone vs parental sensitive line). Raw files are stored in `data/raw/` organized by dataset.

**Note on UWB1.289**: Only BRCA1-deficient samples (n=2 parental + 2 resistant) were included. BRCA1-proficient sublines from the same GEO entry were excluded because they represent a distinct resistance mechanism (BRCA1 reversion) that is biologically incoherent with the BRCA1-deficient context studied here. Including proficient samples would double-count UWB biology and introduce a heterogeneous mixture of resistance mechanisms into the meta-analysis.

## Preprocessing and quality control

For FPKM/TPM matrices, log2(x + 1) transformation was applied before all linear modeling. Genes with zero or near-zero variance were filtered before fitting. For count data (UWB1.289), minimum count filtering was applied and variance-stabilizing transformation (VST) was used for visualization. PCA and volcano plots were generated for visual confirmation of condition separation.

## Differential expression

### limma (A2780, PEO1)

Standard limma was applied to log2-transformed FPKM (A2780) and log2-transformed TPM (PEO1). A design matrix without intercept was constructed and a Resistant-vs-Parental contrast was estimated. Empirical Bayes shrinkage (`eBayes`) was applied; statistical output includes moderated t-statistics, log2FC, and FDR (Benjamini-Hochberg).

**Rationale for limma on log2(FPKM+1) and log2(TPM+1)**: The voom transformation is specifically designed for raw count data, where mean-variance relationships follow a Poisson/negative-binomial model. FPKM and TPM are pre-normalized abundance estimates on a continuous scale; applying voom to these values is methodologically incorrect. Standard limma on log-transformed continuous expression data is the recommended approach for pre-normalized data (Ritchie et al., 2015 *Nucleic Acids Res*; Law et al., 2014 *Genome Biol*).

### DESeq2 (UWB1.289 BRCA-deficient)

DESeq2 was run on the 4 BRCA-deficient samples (2 parental, 2 resistant). The design formula `~ resistance` was used, with `resistance` as a two-level factor (Parental/Resistant). Log2FC shrinkage was applied using `lfcShrink()` with the apeglm method (`coef = "resistance_Resistant_vs_Parental"`). The SE output from `lfcShrink()` was used directly as input SE for the meta-analysis.

## Functional enrichment — GSEA

Gene Set Enrichment Analysis (GSEA) was performed per dataset using the Hallmarks collection (MSigDB H, n=50 gene sets). Gene lists were ranked by: t-statistic for limma datasets; z = log2FC / lfcSE for DESeq2. GSEA was run with `clusterProfiler::GSEA()` (minGS=15, maxGS=500, eps=0). Lollipop plots show top significant terms (FDR < 0.10) and an integrated normalized enrichment score (NES) heatmap displays terms recurrent in ≥2 datasets.

## Functional enrichment — GSVA

Pathway activity scores were computed per sample using GSVA (`GSVA::gsva()`, Gaussian kernel for continuous expression, Hallmarks collection). Per-dataset limma was then applied to the GSVA score matrix to estimate differential pathway activity (Resistant vs Parental). Pathways recurrent in ≥2 datasets (FDR < 0.10) were selected for integrated visualizations.

## Meta-analysis

A gene-level random-effects meta-analysis was performed across the 3 datasets using `metafor::rma.uni()` with REML estimation. For each gene, the effect size (log2FC) and its SE were extracted:

- **limma datasets**: SE = |log2FC| / |t| (exact relationship for the moderated t-statistic)
- **DESeq2 dataset**: SE = lfcSE from `lfcShrink()` (direct output)

For genes present in ≥2 datasets, the REML model was fit; a DerSimonian-Laird fallback was used if Fisher scoring failed to converge. Heterogeneity was quantified via I² and τ². Fixed-effects estimates were also computed for sensitivity comparison.

**Interpretation of I²**: I² values ≥50% indicate substantial heterogeneity, which may reflect dataset-specific biology, different cell line backgrounds, or technical variation. In this 3-study meta-analysis (k=2–3 per gene), I² estimates are inherently unstable; heterogeneity statistics should be interpreted with caution.

Downstream signatures (UP/DOWN gene lists, high-confidence set: FDR < 0.05 + I² < 50%, concordance scores) were exported. Meta-GSEA on the gene-level z-ranking was performed using `fgsea::fgseaMultilevel()` with Hallmarks and Reactome gene sets (minGS=15, maxGS=500, seed=1).

**SE derivation note**: The approximation SE ≈ |log2FC| / |t| is exact for the moderated t-statistic in limma, where t = log2FC / SE by definition. This is the standard derivation used in multi-study transcriptomic meta-analyses.

## Meta-GSEA ranking — beta_meta vs z_meta

Gene-level meta-GSEA was performed using two alternative rankings: (1) z_meta = beta_meta / SE_meta (standardized effect size), and (2) beta_meta (unstandardized pooled log2FC). Given the high heterogeneity of the meta-analysis (median I² = 96%), the z_meta ranking is dominated by between-study variance and produces only 1 significant Hallmark pathway (KRAS_SIGNALING_DN, FDR = 0.0099). The beta_meta ranking captures the average direction of the effect and yields 8 significant Hallmark pathways at FDR < 0.05. The beta_meta-based meta-GSEA is reported in the main figures; both rankings are available in `results/tables/meta_gsea/`.

## Publication figures

Five publication-quality figures (double-column format, 200 mm wide, PDF + PNG at 600 DPI) were generated using ggplot2, patchwork, ComplexHeatmap, and grid (R):

| Figure | Script | Question answered |
|--------|--------|-------------------|
| Fig 1 | `06_fig1_study_design.R` | Are the 3 datasets individually valid and show transcriptomic evidence of olaparib resistance? |
| Fig 2 | `07_fig2_overlap.R` | Is there a shared transcriptomic signature across all 3 lines, or are changes mostly context-specific? |
| Fig 3 | `08_fig3_hallmarks.R` | Do transcriptomic changes converge on the same biological programs across the 3 models? |
| Fig 4 | `09_fig4_gsva.R` | Is the differential pathway activation from GSEA consistent at the individual sample level? |
| Fig 5 | `10_fig5_meta.R` | What is the integrated transcriptomic signature of olaparib resistance and what biological programs does it implicate? |

All figure scripts source `scripts/00_config.R` for centralized parameters and color palettes.

## Visualizations

Volcano plots, PCA, Venn/UpSet diagrams, ComplexHeatmap heatmaps (Z-score and logFC), lollipop plots, and forest plots were generated. Publication figures were exported as PDF (vector) and PNG (600 DPI). Thresholds used: FDR < 0.05 (primary significance), FDR < 0.10 (visualization), |log2FC| ≥ 0.584 (FC ≥ 1.5) for intersection lists.

## Software and computational environment

- R version: see `logs/sessionInfo_*.txt`
- Bioconductor packages: limma, DESeq2, GSVA, clusterProfiler, fgsea, metafor, ComplexHeatmap, msigdbr
- Full package versions: captured in `logs/sessionInfo_*.txt` at pipeline execution
- Conda environment: `env/environment.yml`
- All scripts source `scripts/00_config.R` for centralized parameters, thresholds, and sample IDs

## Reproducibility checklist

- [ ] Record `sessionInfo()` and package versions
- [ ] Verify paths in `data/raw/` match expected filenames
- [ ] Export `OLAPARIB_RESISTANCE_DIR` if running outside project root
- [ ] Run scripts in the order listed in `docs/RUNBOOK.md`

## Key parameters

| Parameter | Description | Value |
|-----------|-------------|-------|
| `FDR_CUTOFF` | Primary significance threshold | 0.05 |
| `FDR_CUTOFF_VIZ` | Visualization threshold | 0.10 |
| `LOGFC_CUTOFF` | log2FC cutoff for intersection lists | 0.584 (FC ≥ 1.5) |
| `MIN_GS` / `MAX_GS` | Gene set size limits (GSEA/GSVA) | 15 / 500 |
| `TOP_N_GENES` | Top DEGs per direction (heatmaps) | 20 |
| `TOP_N_PLOT` | Terms shown in lollipop plots | 10 |
| Meta model | Meta-analysis estimator | RE (REML), DL fallback |

