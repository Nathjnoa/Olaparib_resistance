# Análisis Preliminar de Resultados

**Fecha:** 2026-03-15
**Sesión:** Revisión post-pipeline completo (scripts 01–05_3)
**Estado:** Preliminar — pendiente validación funcional y publicación

---

## Contexto

Meta-análisis transcriptómico de resistencia a olaparib en cáncer de ovario usando 3 datasets públicos de GEO:

| Dataset | Línea | BRCA1 | Método DE | Muestras |
| --- | --- | --- | --- | --- |
| GSE153867 | A2780 | WT | limma (FPKM) | 8 par / 8 res |
| GSE117765 | PEO1 | BRCA1-mut | limma (TPM) | 2 par / 4 res |
| GSE235980 | UWB1.289 | BRCA1-def | DESeq2 (counts) | 2 par / 2 res |

---

## 1. DEGs por dataset y reproducibilidad

| Dataset | UP | DOWN | Total |
| --- | --- | --- | --- |
| A2780 | 3,292 | 4,240 | 7,532 |
| UWB BRCA-def | 2,684 | 2,077 | 4,761 |
| PEO1 | 814 | 1,009 | 1,823 |

**Triple intersección (FDR < 0.05, |log2FC| ≥ 0.584 en los 3):**
- UP en los 3: **30 genes**
- DOWN en los 3: **35 genes**
- Core total: **65 genes**

PEO1 es el cuello de botella (n=4 muestras, TPM+limma = menor poder estadístico). Los 65 genes de triple intersección representan un **piso conservador** de la biología compartida — el meta-análisis recupera mucha más señal al combinar estadísticamente los tamaños de efecto en lugar de requerir intersección estricta.

**Caveat central: heterogeneidad extrema.**
La mediana de I² entre los genes del meta-análisis es **96%** (media: 82.6%). Solo el 10.1% de los genes tiene I² < 25%. Esto no es ruido: las 3 líneas celulares adquieren resistencia a olaparib por rutas transcriptómicas parcialmente distintas. La heterogeneidad es un resultado biológico.

---

## 2. Core de alta confianza (meta-análisis)

**361 genes high-confidence:** k=3 datasets, FDR < 0.05, I² < 50%, concordancia de dirección = 1.

Estos son los candidatos más robustos a biomarcadores generalizables de resistencia.

### Top genes DOWN (más reproducibles)

| Gen | beta_meta | z_meta | FDR | I² | Función |
| --- | --- | --- | --- | --- | --- |
| LTBP1 | −0.84 | −29.1 | 4.9e-183 | 0% | Secuestra TGF-β en matriz extracelular |
| FAM167A | −1.09 | −27.5 | 2.7e-163 | 0% | Función desconocida; candidato novel |
| AK5 | −1.02 | −26.8 | 1.5e-155 | 0% | Metabolismo de adenosina |
| USP18 | −0.94 | −22.6 | 3.9e-110 | 0% | Regulador negativo de señalización IFN tipo I |
| DNAJB11 | −0.77 | −18.8 | 1.9e-76 | 19.5% | Chaperón del RE; sensor de estrés ER |
| BAMBI | −0.79 | −16.2 | 1.7e-56 | 0% | Pseudo-receptor inhibidor de TGF-β/BMP |

### Top genes UP (más reproducibles)

| Gen | beta_meta | z_meta | FDR | I² | Función |
| --- | --- | --- | --- | --- | --- |
| URB1 | 0.67 | +25.3 | 9.1e-138 | 0% | Biogénesis de ribosomas |
| NEDD9 | 0.73 | +19.3 | 2.1e-80 | 0% | Dinámica de adhesiones focales, invasión |
| TUBB6 | 0.68 | +17.4 | 6.5e-65 | 0% | Tubulina β6; reorganización del citoesqueleto |
| **JUN** | **0.91** | **+16.1** | **1.5e-55** | **0.04%** | **Factor de transcripción AP-1** |
| UCK1 | 0.63 | +15.7 | 5.2e-53 | 0% | Metabolismo de pirimidinas |
| TOP3A | 0.38 | +12.6 | 2.9e-34 | 0% | Topoisomerasa III alfa |

### Core triple intersección (30 UP + 35 DOWN) — destacados

**UP coherentes con resistencia a PARPi:**

- **SNAI2** (Slug) — factor de transcripción maestro de la EMT
- **JUN** (c-Jun/AP-1) — señalización de supervivencia y resistencia a fármacos
- **NOTCH1** — promueve stemness y quimiorresistencia en ovario
- **IER3** — gen antiapoptótico inducido por NF-κB
- **FSCN1** (Fascin-1) — reorganización actínica en migración invasiva
- **NEDD9** — dinámica de adhesiones focales, invasión
- **MDK** (Midkine) — factor de crecimiento pro-supervivencia

**DOWN coherentes con resistencia a PARPi:**

- **IFITM1** — proteína inducida por interferón; su bajada sugiere supresión de inmunidad innata
- **DNAJB11** — chaperón ER; su bajada sugiere reconfiguración del UPR
- **LTBP1** — secuestra TGF-β en la matriz; su bajada libera TGF-β disponible
- **BAMBI** — inhibe TGF-β/BMP; su bajada desinhibiría la señalización TGF-β
- **NLRP2** — componente del inflamasoma; supresión de muerte celular inflamatoria

---

## 3. Hallmarks y pathways

### GSEA Hallmarks por dataset — discordancia significativa

| Hallmark | A2780 NES | PEO1 NES | UWB NES |
| --- | --- | --- | --- |
| E2F_TARGETS | +2.35* | NS | −1.54* |
| G2M_CHECKPOINT | +1.78* | NS | −1.80* |
| MYC_TARGETS_V1 | +1.79* | NS | −1.63* |
| INTERFERON_GAMMA | NS | NS | −2.37* |
| INTERFERON_ALPHA | NS | NS | −2.29* |
| EMT | −1.69* | NS | +1.54* |
| TNFA_VIA_NFKB | +1.58* | +1.61* | −1.45* |

*padj < 0.05

La discordancia de dirección en E2F, EMT, MYC, y TNFA entre líneas confirma que la heterogeneidad de I² es biológicamente real. No hay un pathway único "de resistencia" — cada línea usa rutas diferentes.

### Meta-GSEA (sobre ranking z del meta-análisis)

**Hallmarks — único significativo a padj < 0.05:**
- **KRAS_SIGNALING_DN** (NES = +2.14, padj = 0.0099): los genes reprimidos por KRAS activo están upregulados en células resistentes.

**Reactome — ninguno a padj < 0.05.** Los más enriquecidos nominalmente:
- ROS_AND_RNS_PRODUCTION_IN_PHAGOCYTES (NES = +2.30, p = 0.002)
- ION_HOMEOSTASIS (NES = +2.25, p = 0.006)
- TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS (NES = −1.95, p = 0.003)
- COLLAGEN_BIOSYNTHESIS (NES = −1.77, p = 0.008)

La pobreza de señal en meta-GSEA Reactome refleja directamente la heterogeneidad genómica entre líneas — el ranking z combinado diluye las señales discordantes.

### GSVA — pathways recurrentes

De 10 pathways significativos en ≥2 datasets (FDR < 0.10), solo **APICAL_JUNCTION** fue consistente en dirección. Los demás (IFN-alpha, IFN-gamma, E2F, G2M, mTORC1) fueron significativos en 2 datasets pero con **direcciones opuestas**.

---

## 4. Hallazgos más interesantes

### A. Supresión del eje IFN en UWB — paradoja USP18

UWB muestra supresión masiva de INTERFERON_GAMMA (NES = −2.37) e INTERFERON_ALPHA (NES = −2.29). Pero USP18 — un regulador **negativo** de la vía IFN tipo I — también baja consistentemente en los 3 datasets (I² = 0%). Si USP18 baja, la señalización IFN debería *aumentar* (menos inhibición), pero los genes blanco de IFN están suprimidos. Esto es paradójico e implica que la vía IFN está apagada en un punto **upstream** de USP18 (posiblemente epigenéticamente), y que USP18 baja como efecto secundario.

Esto conecta con literatura reciente:
- Guo et al. 2024 (PMID 39660334): cGAS-STING-TBK1-IRF3 promueve resistencia a olaparib vía NF-κB/STAT3
- Ding et al. 2023 (PMID 36609487): STAT3 media inmunosupresión por macrófagos en resistencia adaptativa a PARPi en tumores BRCA1-def
- Liu et al. 2024 (PMID 39412086): agonistas de STING combinados con olaparib revierten resistencia

**Hipótesis:** el silenciamiento de ISGs (IFITM1, USP18, etc.) en células resistentes es epigenético y podría revertirse con inhibidores de DNMT (5-azacitidina) combinados con olaparib.

### B. JUN/AP-1 — señal más consistente de todo el meta-análisis

JUN: beta = 0.91, z = +16.1, I² = **0.04%**. Está en la triple intersección y es high-confidence. AP-1 promueve supervivencia, evasión apoptótica y resistencia a múltiples fármacos. Sin embargo, no hay literatura que lo vincule específicamente a resistencia a PARPi en ovario — esto podría ser un hallazgo novel. Candidato prioritario para validación funcional con el inhibidor T-5224 (AP-1 inhibitor) en combinación con olaparib.

### C. Remodelación TGF-β sin activación Smad canónica

BAMBI (I² = 0%) y LTBP1 (I² = 0%) bajan consistentemente → más TGF-β disponible. Pero en meta-GSEA Reactome, TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS está negativamente enriquecido (NES = −1.95). Señalización TGF-β desinhibida pero sin activación canónica Smad → señalización no canónica (MAPK, PI3K, Rho-GTPasas) o contexto-específica.

### D. EMT selectiva: invasión sí, programa completo no

SNAI2, FSCN1, NEDD9 y HSPG2 están en la triple intersección UP — componentes de invasión y migración. Sin embargo, el hallmark EMT completo es contradictorio entre líneas (DOWN en A2780, UP en UWB). Las células resistentes activan selectivamente genes de invasión y supervivencia del programa EMT sin completar la transición completa — EMT parcial o "mesenchymal shift" selectivo.

---

## 5. Caveats

- **Heterogeneidad extrema (mediana I² = 96%):** la mayoría de los genes no replican entre las 3 líneas. Solo los 361 high-confidence genes son candidatos a marcadores generalizables.
- **Líneas celulares en cultivo:** no capturan heterogeneidad tumoral, microambiente o efectos inmunes in vivo.
- **Datos pre-normalizados:** FPKM/TPM en A2780 y PEO1 limitan el rigor estadístico.
- **Causalidad no establecida:** todos los hallazgos son asociativos. Se requiere validación funcional.
- **n pequeño:** PEO1 tiene n=4 muestras — estimados de varianza inestables.

---

## 6. Preguntas abiertas prioritarias

1. **¿El silenciamiento IFN es epigenético?** Buscar datos de metilación de promotores de ISGs en estas líneas (GEO/ENCODE). Si sí: hipótesis terapéutica con 5-azacitidina + olaparib.

2. **¿JUN/AP-1 conduce funcionalmente la resistencia?** Validar con knockdown o inhibidor T-5224 en A2780, PEO1 y UWB resistentes.

3. **¿Los 361 genes high-confidence predicen resistencia en tumores de pacientes?** Validar en cohortes clínicas con datos transcriptómicos y respuesta a PARPi (TCGA ovario, ensayos ARIEL/NOVA).

4. **¿Cuál es la naturaleza del señalamiento TGF-β no canónico?** Verificar fosforilación de ERK/p38/AKT en células resistentes tras estimulación con TGF-β.

---

## 7. Próximos pasos analíticos sugeridos

- [ ] Enriquecimiento de los 361 genes high-confidence en bases de datos de fármacos (DGIdb, ChEMBL) para drug repurposing
- [ ] Validación de la firma en dataset externo (si está disponible o se puede obtener)
- [ ] Análisis de conectividad (LINCS/CMap) con la firma UP/DOWN para identificar compuestos que reviertan el perfil
- [ ] Búsqueda de datos de metilación en GEO para las mismas líneas (correlación con silenciamiento IFN)
- [ ] Figura integradora: heatmap de los 361 genes high-confidence con anotaciones de vía y I²
