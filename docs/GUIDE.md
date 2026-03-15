# Guía Académica del Análisis

Una explicación de qué se hizo, por qué se hizo de esa manera, y qué significa cada paso.
Escrita para cualquier lector: desde alguien que apenas está aprendiendo bioinformática hasta un revisor de revista.

---

## ¿Cuál es la pregunta biológica?

El olaparib es un inhibidor de PARP, una enzima que repara roturas de cadena simple en el ADN. Las células con mutaciones en BRCA1 o BRCA2 tienen un defecto en la reparación por recombinación homóloga (HR). El olaparib aprovecha eso: si una célula no puede reparar por HR y además le bloqueas PARP, acumula daño en el ADN y muere. A esto se le llama **letalidad sintética**.

El problema clínico es que muchas pacientes con cáncer de ovario BRCA-deficiente responden inicialmente a olaparib pero desarrollan **resistencia** con el tiempo. ¿Qué cambia en esas células para que sobrevivan?

Este proyecto busca responder esa pregunta a nivel transcriptómico: **¿qué genes y vías están consistentemente alterados cuando las células adquieren resistencia a olaparib?**

---

## ¿Por qué usar 3 datasets en lugar de 1?

Cualquier resultado de un único experimento puede ser artefacto de esa línea celular, ese laboratorio, o ese protocolo. Si el mismo gen aparece consistentemente alterado en 3 experimentos independientes, con 3 líneas celulares distintas, en 3 laboratorios distintos, la probabilidad de que sea una casualidad es mucho menor.

Los 3 datasets usados son:

| Dataset | Línea | Contexto BRCA |
| --- | --- | --- |
| GSE153867 | A2780 | BRCA1 silvestre |
| GSE117765 | PEO1 | Mutación BRCA1 C61G |
| GSE235980 | UWB1.289 | Deleción BRCA1 (solo muestras def.) |

Las tres son líneas de cáncer de ovario epitelial con resistencia adquirida in vitro a olaparib, comparadas contra su línea parental sensible.

---

## Paso 1: Expresión diferencial por dataset

### ¿Qué se hace?

Se compara la expresión de genes entre células **resistentes** vs. células **parentales (sensibles)** en cada dataset por separado. Para cada gen se calcula:

- **log2FC** (fold-change en escala logarítmica): cuánto más o menos se expresa el gen en resistentes
- **p-valor ajustado (FDR)**: qué tan probable es que esa diferencia sea por azar, corregido por el hecho de que se están haciendo miles de pruebas a la vez

### ¿Por qué métodos distintos (limma vs DESeq2)?

Los datos de A2780 y PEO1 ya vienen en formato **FPKM/TPM**, que son estimados de abundancia normalizados (counts por kilobase por millón). Para este tipo de datos continuos se usa **limma** en escala log2, que asume distribución normal y es el método correcto.

Para UWB1.289 se tienen **counts crudos** (número de lecturas por gen). Los counts siguen una distribución binomial negativa, y para este caso se usa **DESeq2**, que modela correctamente esa distribución. Usar DESeq2 sobre FPKM/TPM sería incorrecto — y viceversa para limma sobre counts crudos.

**Nota técnica:** voom (de limma) transforma counts para que se puedan tratar como datos gaussianos. Pero voom también es para counts, no para FPKM/TPM. Por eso A2780 y PEO1 usan limma estándar sobre log2(x+1) y no voom.

### ¿Por qué FDR y no p-valor simple?

Se prueba simultáneamente la diferencia en ~20,000 genes. Si se usa p < 0.05 sin corrección, el 5% de los genes (≈1,000) darían "significativo" por azar. La corrección de Benjamini-Hochberg controla la **tasa de falsos descubrimientos (FDR)**: si FDR < 0.05, se acepta que hasta el 5% de los genes declarados significativos podrían ser falsos positivos.

### Figura 1: ¿Por qué es la primera figura y qué pregunta responde?

**Pregunta:** *¿Son los tres datasets individualmente válidos, y cada uno muestra evidencia transcriptómica de resistencia a olaparib?*

Esta es la pregunta obligatoria de cualquier meta-análisis antes de integrar nada. Si el lector no confía en los datos individuales, no confiará en el meta-resultado. El paper tiene un problema de credibilidad inherente: tres datasets con métodos distintos (limma-FPKM, limma-TPM, DESeq2-counts), tamaños de muestra muy diferentes (n=16 vs n=6 vs n=4) y fondos genéticos distintos (BRCA-WT, BRCA1-mut, BRCA1-def). La Figura 1 resuelve ese problema antes de avanzar al análisis integrado.

**Por qué cada panel:**

- **Panel A (diseño del estudio):** orienta al lector — tres líneas celulares, tres contextos BRCA distintos, mismo paradigma experimental (parental → resistente). Sin este mapa, los paneles siguientes se leen sin contexto.

- **Panel B (PCA por dataset):** el PCA es la prueba visual de que la resistencia domina la varianza en cada dataset. Si parental y resistente no se separan en PC1 o PC2, el experimento no tiene señal suficiente para ningún análisis posterior. Que los tres datasets muestren separación — incluso UWB con solo n=4 — valida que la señal biológica existe.

- **Panel C (volcano plots):** muestra la magnitud y la escala del cambio diferencial en cada dataset. La diferencia de escalas entre A2780 (>7,500 DEGs), PEO1 (~1,800) y UWB (~4,700) no es falla — es información. Le dice al lector que PEO1 tiene menos poder estadístico (n=2 vs n=4), anticipando por qué la intersección simple da pocos genes y por qué el meta-análisis es necesario.

- **Panel D (conteo de DEGs):** cuantifica explícitamente esa diferencia de escala. Hace evidente que una intersección estricta (gen significativo en los 3) descartaría la mayoría de la señal real de PEO1. Justifica metodológicamente el meta-análisis de efectos aleatorios que viene después.

**Lógica de secuencia:** solo después de que el lector acepta que los tres datasets individualmente tienen señal coherente tiene sentido mostrarles la intersección (Fig. 2), el análisis de pathways (Fig. 3/4) y el meta-análisis (Fig. 5). Si la Figura 1 no estuviera, el lector llegaría a los resultados integrados sin base para confiar en ellos.

**Script:** `scripts/06_fig1_study_design.R`
**Output:** `results/figures/fig1/Fig1.pdf` + `Fig1.png`

---

## Paso 2: Intersección de DEGs

### ¿Qué se hace?

Se buscan los genes que son diferencialmente expresados (FDR < 0.05, |log2FC| ≥ 0.584 ≈ FC ≥ 1.5) en **los 3 datasets simultáneamente**, separando los que suben y los que bajan.

### Visualización: Venn y UpSet

- **Diagrama de Venn:** muestra visualmente cuántos genes son únicos de cada dataset y cuántos solapan en 2 o 3 datasets.
- **UpSet plot:** alternativa más legible para conjuntos grandes. Cada barra vertical es una combinación de datasets, y la altura indica cuántos genes pertenecen a esa combinación exacta.

### ¿Por qué el solapamiento es pequeño (65 genes de >7,500)?

PEO1 tiene muy pocos DEGs (1,823) comparado con A2780 (7,532). Con n=4 muestras y una variabilidad natural entre clones de resistencia, el poder estadístico para detectar genes diferencialmente expresados es menor. Los 65 genes que solapan son los que superaron el umbral en **los 3** datasets — un criterio muy estricto. El meta-análisis (paso 5) aborda este problema de manera más sofisticada.

---

## Paso 3: GSEA Hallmarks por dataset

### ¿Qué son los gene sets de Hallmarks?

MSigDB (Molecular Signatures Database) compiló 50 "Hallmarks" del cáncer: colecciones curadas de genes que definen procesos biológicos bien conocidos como proliferación (E2F_TARGETS), respuesta a hipoxia (HYPOXIA), transición epitelio-mesénquima (EMT), etc.

### ¿Qué hace GSEA?

Gene Set Enrichment Analysis evalúa si los genes de un pathway determinado tienden a estar **colectivamente** en los extremos del ranking de expresión diferencial (muy up o muy down), en lugar de distribuirse aleatoriamente. No requiere un corte arbitrario de significancia — usa toda la lista de genes rankeados.

El resultado principal es el **NES (Normalized Enrichment Score):** positivo = el pathway tiende hacia arriba en resistentes, negativo = tiende hacia abajo.

### ¿Por qué Hallmarks y no KEGG o GO?

- **KEGG** tiene pathways definidos por reacciones bioquímicas, no siempre corresponden a módulos de expresión coherentes.
- **GO** tiene miles de términos, muchos muy generales (p. ej. "proceso biológico") o muy redundantes.
- **Hallmarks** son 50 procesos curados con alta coherencia interna — una primera lectura biológica limpia. En el meta-GSEA (paso 5) se usan también Reactome para mayor resolución.

---

## Paso 4: GSVA

### ¿En qué se diferencia GSVA de GSEA?

GSEA da una puntuación por **comparación** (resistentes vs. parentales como grupo). GSVA da una puntuación por **muestra individual**: cuán activo está ese pathway en cada muestra.

Esto permite:
1. Ver la heterogeneidad dentro de un grupo (¿todas las células resistentes tienen el mismo perfil?)
2. Aplicar después un modelo estadístico (limma) a las puntuaciones GSVA para comparar grupos

Es como pasar de "¿cuántos géneros musicales están representados en la playlist?" (GSEA) a "¿qué géneros escucha cada persona?" (GSVA).

### Figuras generadas

- **Fig3A:** heatmap de puntuaciones GSVA por muestra — permite ver si cada muestra tiene el patrón esperado
- **Fig3B:** heatmap de logFC de GSVA entre resistentes y parentales — integra los 3 datasets en una sola visualización

---

## Paso 5: Meta-análisis de genes

### ¿Por qué meta-análisis y no simplemente intersección?

La intersección requiere que el gen sea significativo en **todos** los datasets. Con n=4 en PEO1, muchos genes reales no alcanzan significancia estadística en ese dataset solo. El meta-análisis combina la **evidencia estadística** de los 3 datasets en lugar de requerir significancia individual en cada uno.

### Modelo de efectos aleatorios (REML)

Se usa un modelo de **efectos aleatorios** (Random Effects, REML = Restricted Maximum Likelihood):

- **¿Por qué efectos aleatorios y no fijos?** Los efectos fijos asumen que hay un único efecto verdadero que todos los estudios están midiendo con ruido. Los efectos aleatorios permiten que el efecto verdadero sea distinto entre estudios (cada línea celular puede tener un efecto real distinto), y estiman la **distribución** de esos efectos. Como las 3 líneas son biológicamente distintas (BRCA1 WT, mut, def), el modelo de efectos aleatorios es más honesto.

- **REML vs. DerSimonian-Laird (DL):** REML es más preciso para estimar la varianza entre estudios (τ²) cuando hay pocos estudios (k=3). DL tiende a subestimar τ² con k pequeño. Se usa DL como fallback si REML no converge.

### ¿Qué es I²?

I² mide qué proporción de la variación total entre los efectos observados se debe a **heterogeneidad real** (diferencias entre estudios) en lugar de azar (ruido de muestreo):

- I² < 25%: heterogeneidad baja → el gen responde de manera similar en los 3 estudios
- I² 25–50%: moderada
- I² > 75%: alta → el gen responde de manera muy distinta en cada estudio

En este meta-análisis, la **mediana de I² es 96%**. Esto significa que la mayoría de los genes no tienen un efecto homogéneo entre líneas — la biología de resistencia es específica de contexto.

### ¿Cómo se calcula el SE para limma?

DESeq2 reporta directamente el error estándar (lfcSE). Para limma, el SE se aproxima como:

```
SE = |log2FC| / |t_statistic|
```

Esta es la relación exacta del estadístico t de limma (t = log2FC / SE por definición). No es una aproximación — es algebraicamente exacta dado cómo limma define su estadístico.

### Firma high-confidence

Se definen como **high-confidence** los genes con:
- k = 3 (en todos los datasets)
- FDR meta < 0.05
- I² < 50% (heterogeneidad baja-moderada)
- Concordancia de dirección = 1 (misma dirección en todos)

Resultado: **361 genes** que son los candidatos más robustos para validación o uso como biomarcadores.

---

## Paso 5b: Meta-GSEA

### ¿Qué es y cómo se diferencia del GSEA por dataset?

En lugar de rankear genes por el estadístico t de un solo dataset, se rankea por el **z-score del meta-análisis** (z = beta/SE_meta). Este ranking integra la evidencia de los 3 estudios y se usa como input para fGSEA (fast GSEA, implementación eficiente de GSEA con permutaciones).

Se usa también el ranking por **beta_meta** (log2FC del meta) como alternativa — los dos deberían dar resultados similares si no hay heterogeneidad fuerte.

### ¿Por qué Reactome además de Hallmarks?

Reactome tiene >2,000 pathways con mayor resolución mecanística que los 50 Hallmarks. Sin embargo, más pathways = más corrección múltiple = más difícil alcanzar significancia. En este caso, ningún pathway Reactome alcanza FDR < 0.05 en el meta-GSEA, lo cual es coherente con la alta heterogeneidad del meta-análisis.

---

## ¿Por qué los resultados son tan heterogéneos?

I² mediana = 96% significa que **las 3 líneas celulares usan rutas distintas para resistir olaparib**. Esto tiene sentido biológico: hay múltiples mecanismos de resistencia a PARPi documentados en la literatura:

1. Restauración de HR (reversión de mutación BRCA, fusiones alternativas, pérdida de 53BP1)
2. Protección de horquillas de replicación
3. Eflujo de fármaco (PGP/MDR1)
4. Señalización de supervivencia alternativa (PI3K, NF-κB, STAT3)
5. Silenciamiento epigenético de genes que sensibilizan

Cada línea puede haber escogido una o varias de estas rutas. La intersección (65 genes, 361 high-confidence) captura los elementos **transversales** a esas rutas distintas.

---

## Interpretación global: tres ejes de resistencia

Los resultados sugieren que la resistencia a olaparib en estas células implica tres procesos convergentes pero con expresión variable entre líneas:

### Eje 1: EMT parcial / programa de invasión

SNAI2, JUN, NEDD9, FSCN1 suben consistentemente. Las células resistentes adquieren características mesenquimales selectivas (migración, supervivencia) sin completar la transición EMT completa. BAMBI y LTBP1 bajan, desreprimiendo TGF-β, pero sin activación canónica de Smads (señalización no canónica de TGF-β).

### Eje 2: Supresión de inmunidad innata

USP18 e IFITM1 bajan en los 3 datasets. UWB suprime masivamente los programas de IFN-gamma e IFN-alfa. Esta supresión podría facilitar la evasión inmune in vivo y se conecta con el eje cGAS-STING-TBK1, que está emergiendo como un nodo central en la biología de resistencia a PARPi.

### Eje 3: Reconfiguración del aparato secretor y síntesis proteica

URB1 (biogénesis de ribosomas) sube. DNAJB11, TMED4 (tráfico ER-Golgi) bajan. El hallmark UPR (estrés del RE) tiende a bajar en el meta-GSEA. Las células resistentes parecen reconfigurar su ruta secretora y reducir la respuesta al estrés del RE, posiblemente para aumentar capacidad de síntesis proteica.

---

## Limitaciones del estudio

1. **Solo líneas celulares in vitro:** no capturan microambiente tumoral, heterogeneidad clonal, ni efectos inmunes del sistema in vivo.
2. **n pequeño:** especialmente PEO1 (n=4). Los estimados de varianza son menos estables.
3. **Datos pre-normalizados:** FPKM/TPM en A2780 y PEO1 limitan algunas opciones metodológicas.
4. **Causalidad:** todos los hallazgos son asociativos. Un gen alterado en células resistentes no necesariamente *causa* la resistencia.
5. **Heterogeneidad real:** I² = 96% implica que la mayoría de los efectos no se replican de forma homogénea — los resultados son válidos como hipótesis pero requieren validación en sistemas adicionales.

---

## Glosario rápido

| Término | Definición |
| --- | --- |
| log2FC | Logaritmo en base 2 del fold-change. log2FC = 1 → gen 2x más alto. log2FC = -1 → gen 2x más bajo. |
| FDR | False Discovery Rate. Proporción esperada de falsos positivos entre los genes declarados significativos. |
| NES | Normalized Enrichment Score. Puntuación de enriquecimiento de GSEA normalizada por el tamaño del gene set. |
| I² | Porcentaje de variación entre estudios que no se explica por azar. I² alto = heterogeneidad real. |
| τ² | Varianza entre estudios en el modelo de efectos aleatorios. |
| beta_meta | Estimado del efecto (log2FC) combinado del meta-análisis. |
| SE | Error estándar del estimado de efecto. |
| REML | Restricted Maximum Likelihood. Método de estimación de parámetros en modelos mixtos. |
| GSVA | Gene Set Variation Analysis. Calcula actividad de pathways por muestra. |
| EMT | Epithelial-Mesenchymal Transition. Proceso por el que células epiteliales adquieren características mesenquimales (invasión, migración). |
| ISG | Interferon-Stimulated Gene. Gen cuya expresión es inducida por interferones. |
| PARPi | Inhibidor de PARP (Poly ADP-Ribose Polymerase). |
| HR | Homologous Recombination. Reparación del ADN de alta fidelidad, deficiente en células BRCA-mutadas. |
| cGAS-STING | Vía de sensado de ADN citosólico que activa la respuesta inmune innata. |
