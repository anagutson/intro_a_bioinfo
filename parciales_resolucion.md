# Resolución de Parciales — Bioinformática

> Todos los parciales (2018, 2022, 2023, 2024) tienen la misma estructura de 4 problemas.
> Este documento está dividido en dos partes: **teoría con respuestas modelo** y **ejercicios resueltos por año**.

---

## Mapa del parcial

| Problema | Tema central | Subtemas que siempre aparecen |
|----------|-------------|-------------------------------|
| **P1** | Bases de datos primarias | Clasificar BDs (curada/no curada/primaria/secundaria/redundante), BLAST paso a paso, PSI-BLAST, registro/campo/índice, PubMed |
| **P2** | Alineamiento + MSA | NW o SW (llenar la matriz), árbol guía UPGMA, ClustalW con perfiles, SP score |
| **P3** | NGS / Ensamblado / Anotación | Profundidad, cobertura, Lander-Waterman, De Bruijn, contigs, scaffolds, missing genes, proteínas hipotéticas |
| **P4** | BD secundarias (HMM) | 3 problemas del HMM, gathering cutoff, VP/VN/FP/FN, leer WebLogo, expresión regular vs HMM |

---

# PARTE 1 — Teoría y respuestas modelo

---

## P1 — Bases de Datos y BLAST

### Criterios de clasificación

| Eje | Definición | Ejemplo curada | Ejemplo no curada |
|-----|-----------|---------------|-------------------|
| **Curada / No curada** | Curada = revisión manual por expertos antes de aceptar el registro. No curada = el autor deposita directamente sin verificación sistemática. | Swiss-Prot, RefSeq, PDB | GenBank, TrEMBL, GEO |
| **Primaria / Secundaria** | Primaria = dato crudo de experimentos (secuencias, estructuras, expresión). Secundaria = dato derivado de análisis bioinformático de fuentes primarias. | GenBank, Swiss-Prot, GEO | SGD, FlyBase, Pfam, KEGG |
| **Redundante / No redundante** | Redundante = múltiples registros para la misma entidad. No redundante = una entrada canónica por entidad. | GenBank, TrEMBL | RefSeq, Swiss-Prot |

> **Caso especial UniProt:** las secuencias vienen de GenBank/PDB (→ primaria), pero las anotaciones funcionales las agregan curadores revisando literatura (→ secundaria). En el parcial se clasifica como **primaria**.

---

### Tabla de BDs por tipo

| Tipo | Base de datos | Prim/Sec | Curada | Redundante | Nota clave |
|------|--------------|----------|--------|------------|------------|
| **Nucleótidos** | GenBank (NCBI) | Prim | No | Sí | Repositorio NCBI; INSDC |
| | EMBL / ENA | Prim | No | Sí | Equivalente europeo; INSDC |
| | RefSeq | Prim | Parcial | No | Una secuencia canónica por molécula |
| **Proteínas** | Swiss-Prot | Prim | **Sí** | No | Revisión manual; 1 entrada por proteína |
| | TrEMBL | Prim | No | Sí | Traducción automática de GenBank |
| | PDB | Prim | Sí | Sí | Estructuras 3D; redundante porque la misma proteína tiene cientos de estructuras |
| **Expresión** | GEO (NCBI) | Prim | No | Sí | Microarrays y RNA-seq |
| | ArrayExpress | Prim | No | Sí | Equivalente europeo de GEO |
| **Organismo** | FlyBase | Sec | Sí | No | Drosophila melanogaster |
| | WormBase | Sec | Sí | No | C. elegans |
| | SGD | Sec | Sí | No | S. cerevisiae |
| | EcoCyc | Sec | Sí | No | E. coli |
| **Compuestos** | PubChem | Prim | No | Sí | La más grande de libre acceso |
| | DrugBank | Sec | Sí | No | Fármacos con targets y mecanismo |
| | ChEMBL | Sec | Sí | No | Compuestos bioactivos (EBI) |
| **Literatura** | Europe PMC | Prim | No | No | Texto completo open access + preprints; tiene API |
| | Web of Science | Prim | No | No | Comercial; amplia cobertura |
| | Scopus | Prim | No | No | Comercial; alternativa a WoS |
| **Dominios** | Pfam | Sec | Sí | No | Familias de dominios con HMMs |
| **Enfermedades** | OMIM | Sec | Sí | No | Genes y enfermedades mendelianas |
| **Vías metab.** | KEGG | Sec | Sí | No | Pathways metabólicos y de señalización |

---

### Respuestas por categoría — formato parcial

**a) Nucleótidos → GenBank**
Primaria, no curada, redundante. Cada laboratorio deposita directamente. Alternativa si quiero calidad: RefSeq (no redundante, parcialmente curada).

**b) Proteínas → Swiss-Prot**
Primaria, curada, no redundante. Revisión manual; una entrada por proteína con función verificada. Si necesito cobertura máxima: TrEMBL (no curada, redundante).

**c) Expresión génica → GEO**
Primaria, no curada, redundante. Repositorio NCBI de experimentos de microarrays y RNA-seq. Permite reanálisis de experimentos publicados.

**d) Organismo → SGD / FlyBase / WormBase**
Secundarias, curadas, no redundantes. Integran GenBank + UniProt + literatura para un organismo modelo. Son secundarias porque el contenido se genera por análisis bioinformático de fuentes primarias, no por nuevos experimentos.

**e) Compuestos químicos → PubChem o DrugBank**
- PubChem: primaria, no curada, redundante. La más grande de libre acceso.
- DrugBank: secundaria, curada, no redundante. Mejor para información farmacológica confiable.

**f) Literatura (excluyendo PubMed) → Europe PMC**
Primaria, no curada. Incluye texto completo de artículos open access y preprints; tiene API para minería de texto. Alternativas comerciales: Web of Science, Scopus.

**g) Libre → OMIM / Pfam / KEGG**
Las tres son secundarias y curadas:
- OMIM: genes y enfermedades genéticas humanas.
- Pfam: familias de dominios proteicos con modelos HMM.
- KEGG: vías metabólicas; útil para análisis de enriquecimiento.

---

### BLAST — pasos del algoritmo

**Pregunta: "¿Cómo buscaría homólogos de la proteína X?"**

1. Obtener la secuencia en formato FASTA desde UniProt (Swiss-Prot si quiero curada).
2. Usar **BLASTp** contra la BD adecuada (nr, UniProtKB/Swiss-Prot).
3. Aplicar filtro taxonómico si corresponde.

**Los 4 pasos de BLASTp:**

| Paso | Qué hace | Output |
|------|---------|--------|
| 1. Word list | Fragmenta la query en words de largo w=3. Genera vecinas con score BLOSUM62 ≥ T | Lista de words + vecinas |
| 2. Búsqueda en índice | Busca cada word en el índice pre-construido de la BD | Lista de hits (registro, posición) |
| 3. Extensión sin gaps (HSP) | Extiende cada hit izquierda/derecha acumulando score BLOSUM62; para cuando cae bajo umbral | HSPs con score ≥ S |
| 4. Extensión con gaps + E-value | Aplica SW con gaps a los HSPs; calcula E = K·m·N·e^(−λS); reporta los con E < umbral | Hits finales con E-value |

**Tamaño de palabra (word size):**

| W | Sensibilidad | Especificidad | Velocidad |
|---|-------------|---------------|-----------|
| Pequeño (W=3) | Alta (detecta homólogos lejanos) | Baja (más FP) | Lenta |
| Grande (W=6) | Baja (pierde divergentes) | Alta (menos FP) | Rápida |

Cantidad de words para secuencia de largo L: **L − W + 1**

**¿Cómo encontrar homólogos más divergentes? → PSI-BLAST**
Ronda 1 con BLOSUM62 → construye PSSM con los hits confiables → rondas siguientes buscan con esa PSSM posición-específica. Detecta homólogos que BLAST estándar no ve porque las posiciones conservadas tienen mayor puntaje y las variables menor penalidad.

**Registro, campo e índice:**
- **Registro:** unidad completa de información (ej: toda la entrada de una proteína en UniProt).
- **Campo:** categoría dentro del registro (ej: "organismo", "secuencia", "función").
- **Índice:** estructura pre-calculada que permite encontrar registros rápido a partir de un campo. Sin índice, BLAST sería O(n·m) por registro → inviable.

**Identities vs Positives en el output de BLAST:**
- **Identities:** posiciones donde los dos aminoácidos son exactamente iguales.
- **Positives:** posiciones donde el score BLOSUM62 del par es > 0. Incluye cambios conservativos (Ser→Ala = +1 → Positive pero no Identity).

**Búsqueda en PubMed:**
```
("Wilson Disease"[Title/Abstract])
AND ("ATP7B"[Title/Abstract] OR "CopA"[Title/Abstract])
AND ("2020/01/01"[Date - Publication] : "2025/12/31"[Date - Publication])
```
- `[Title/Abstract]` → más específico que buscar en todo el registro.
- Incluir sinónimos del gen/enfermedad.

---

## P2 — Alineamiento y MSA

### Needleman-Wunsch — procedimiento completo

**Parámetros:** Match = +X, Mismatch = −Y, Gap = −Z (los da el enunciado).

**Pasos:**
1. Crear matriz de (|seqA|+1) × (|seqB|+1).
2. **Inicializar:** F(i,0) = i × gap; F(0,j) = j × gap.
3. **Llenar** cada celda con el máximo de:
   - **Diagonal** ← F(i−1, j−1) + score(aᵢ, bⱼ)  → match o mismatch
   - **Arriba** ← F(i−1, j) + gap  → gap en la secuencia de la fila
   - **Izquierda** ← F(i, j−1) + gap  → gap en la secuencia de la columna
4. **Traceback:** desde la celda (m,n) seguir flechas hasta (0,0).

> Si dos o tres movimientos dan el mismo máximo → empate → múltiples traceback posibles → múltiples alineamientos óptimos, todos con el mismo score. Todos son correctos.

---

### UPGMA + ClustalW — procedimiento completo

**Paso 1 — Scores de pares:** alinear cada par de secuencias (NW o SW) y anotar el score.

**Paso 2 — Árbol UPGMA:**
- Encontrar el par con **mayor score** (= más similares).
- Unirlos en un nodo compuesto (A+B).
- Actualizar la distancia del nodo a las demás: d(A+B, C) = [d(A,C) + d(B,C)] / 2
- Repetir hasta unir todo.

**Paso 3 — MSA con perfiles:**
- Primero alinear las dos más similares → perfil 1.
- Agregar la siguiente secuencia al perfil (no re-alinear, solo insertar gaps donde el perfil los tiene).
- Puntuar con **SP score** = suma de scores de todos los pares en cada columna.

---

## P3 — NGS, Ensamblado y Anotación

### Profundidad y cobertura

**Profundidad (depth):** número de veces que una base fue leída, en promedio.
```
Profundidad target = (N_reads × L_read) / L_genoma
```
- Se calcula **antes** del experimento para saber cuántas lecturas pedir.
- **Lander-Waterman:** P(base sin leer) = e^(−D)

| Profundidad | % bases sin leer |
|-------------|-----------------|
| 1× | 37% |
| 5× | 0.67% |
| 10× | ~0% |
| 30× | estándar clínico para variant calling |

**Cobertura (coverage):** % del genoma cubierto con profundidad ≥ umbral. Distinto a profundidad.

---

### De Bruijn → contigs → scaffolds

**Contigs** (solapamiento de lecturas):
1. Fragmentar lecturas en k-mers.
2. Construir el grafo de De Bruijn: nodos = k-mers, aristas = lecturas.
3. Resolver el grafo: caminos de bajo peso → errores de secuenciación (eliminar); ciclos → repeticiones.
4. Caminos lineales resueltos = contigs.

**Scaffolds** (ordenar contigs):
- Usar lecturas **paired-end**: dos lecturas del mismo fragmento de tamaño conocido.
- Si la lectura 1 mapea en el contig A y la lectura 2 en el contig B → A y B están orientados y separados por ≈ (tamaño_inserto − 2×largo_lectura).
- Con genoma de referencia: mapear contigs contra la referencia → inferir posición y orden directamente.

**¿Por qué las repeticiones complican el ensamblado?**
Los k-mers de la región repetida se comparten entre las distintas copias → generan ciclos/bifurcaciones en el grafo → el ensamblador no puede resolver el camino → corta el contig en los bordes de la repetición → más contigs, menor N50.
**Solución:** lecturas largas (PacBio, Nanopore) que abarcan la repetición entera.

---

### Missing genes vs proteínas hipotéticas

| | Proteína hipotética | Missing gene |
|--|---------------------|--------------|
| **Qué es** | ORF predicho (tiene start/stop y largo adecuado) pero sin homólogos en BD → función desconocida | Gen que debería estar según la BD de referencia (RAST/FigFams) pero no se encontró en el ensamblado |
| **¿Se anotó?** | Sí, como ORF sin función | No se encontró ningún ORF |
| **Causa probable** | Proteína nueva o muy divergente | Cobertura baja → hueco en el ensamblado → no se formó contig → no hay ORFs |

**¿Por qué puede haber baja cobertura en algunas regiones?**
- Distribución de lecturas es una Poisson → por azar algunas zonas se leen menos.
- GC bias: PCR amplifica peor regiones con >65% GC (problema en Illumina).
- Estructuras secundarias del ADN interfieren con la síntesis.
- Fragmentación sesgada en la preparación de la librería.

**¿Cómo mejorar el ensamblado?**
1. Lecturas largas (PacBio o Nanopore) → resuelven repeticiones.
2. Mate-pair con inserto de 10–50 kb → conecta contigs más separados.
3. Mayor profundidad → reduce bases sin leer (Lander-Waterman).
4. Ensamblado híbrido: cortas (alta precisión) + largas (mayor continuidad).

---

## P4 — HMM y Bases de Datos Secundarias

### Los 3 problemas del HMM

```
┌─────────────────────────────────────────────────────────────┐
│  1. ENTRENAMIENTO                                           │
│     Input:  MSA semilla (secuencias que SÉ que pertenecen) │
│     Output: parámetros del modelo                          │
│                                                             │
│     Por cada columna del MSA:                              │
│     - contar frecuencias de aa → P(emisión) del estado MATCH│
│     - contar transiciones M→M, M→I, M→D → P(transición)   │
├─────────────────────────────────────────────────────────────┤
│  2. SCORING (puntaje)                                       │
│     Input:  secuencia nueva + modelo entrenado             │
│     Output: score / bit-score                              │
│                                                             │
│     Calcular P(modelo generó esa secuencia).               │
│     Comparar con el umbral → pertenece o no a la familia.  │
├─────────────────────────────────────────────────────────────┤
│  3. ALINEAMIENTO                                            │
│     Input:  secuencia nueva + modelo entrenado             │
│     Output: alineamiento de la secuencia contra el perfil  │
│                                                             │
│     Algoritmo de Viterbi: encuentra el camino más probable │
│     a través de los estados que genera la secuencia.       │
│     (Programación dinámica, análogo a NW pero contra HMM)  │
└─────────────────────────────────────────────────────────────┘
```

---

### Flujo completo — desde cero hasta buscar en un genoma

```
Secuencias semilla (verdaderos positivos confirmados)
         ↓
Alineamiento múltiple (MSA)  ←  define la estructura del HMM
         ↓                       (columnas del MSA = estados MATCH)
Calcular parámetros:
  - P(aa | posición) por columna → probabilidades de emisión
  - P(M→M, M→I, M→D) → probabilidades de transición
         ↓
Modelo HMM entrenado
         ↓ ─── SCORING ───
Evaluar secuencias semilla con el modelo → distribución de scores de VP
Evaluar secuencias negativas o aleatorias → distribución de scores de VN
         ↓
Graficar histograma: determinar umbral (gathering cutoff)
en el punto de separación entre las dos distribuciones
         ↓ ─── ALINEAMIENTO ───
Evaluar secuencias nuevas (ej: proteoma completo)
Score ≥ cutoff → pertenece a la familia → candidato
Score < cutoff → no pertenece
```

---

### Cutoffs en Pfam

| Cutoff | Definición |
|--------|-----------|
| **GA (Gathering)** | Score mínimo para incluir en la familia completa |
| **TC (Trusted)** | Score mínimo de la secuencia conocida con score más bajo que SÍ pertenece |
| **NC (Noise)** | Score máximo de una secuencia que NO pertenece |

> Idealmente TC > GA > NC. Si se solapan → FP y FN inevitables.

Para hacer la predicción **más específica:** subir el umbral → los FP quedan por debajo → menos FP, pero más FN.

---

### VP / VN / FP / FN

```
                     Predicción del modelo
                    Positivo    Negativo
               ┌───────────────────────────┐
Realidad   +   │    VP     │     FN        │
           −   │    FP     │     VN        │
               └───────────────────────────┘
```

| | Efecto sobre métricas |
|--|----------------------|
| **VP** (verdadero positivo) | Bien clasificado; contribuye a sensibilidad |
| **VN** (verdadero negativo) | Bien clasificado; contribuye a especificidad |
| **FP** (falso positivo) | Predijo pertenece, pero no → reduce especificidad |
| **FN** (falso negativo) | Predijo no pertenece, pero sí → reduce sensibilidad |

**Sensibilidad** = VP / (VP + FN) → qué fracción de los positivos reales detecté.
**Especificidad** = VN / (VN + FP) → qué fracción de los negativos reales rechacé.

Subir el umbral → ↑ especificidad, ↓ sensibilidad.
Bajar el umbral → ↑ sensibilidad, ↓ especificidad.

---

### Leer un WebLogo

- **Columna alta** = conservada → esa posición tiene función crítica (mutarla rompe la función).
- **Letra grande** = frecuente en esa posición = residuo más probable en ese estado MATCH.
- **Columna baja** = variable → el modelo tolera sustituciones ahí.

**Interpretación funcional por residuo conservado:**

| Residuo | Interpretación probable |
|---------|------------------------|
| D / E | Sitio catalítico ácido o coordinación de iones |
| K / R | Interacción con ADN (carga positiva) |
| C | Puente disulfuro o sitio activo |
| W / F | Empaquetamiento aromático o unión a ligando |
| G | Loop de alta flexibilidad (Gly sin cadena lateral) |

**Contar estados MATCH:** 1 columna del logo = 1 estado MATCH. Logo de 15 columnas = HMM con 15 estados MATCH.

**Inserción vs deleción respecto del modelo:**
- **Inserción:** la secuencia tiene residuos EXTRA que el modelo no espera.
- **Deleción:** a la secuencia le FALTA algo que el modelo espera (gap en la secuencia query al alinearla contra el perfil).

---

### Expresión regular (PROSITE) vs HMM

| | Expresión regular | HMM |
|--|-------------------|-----|
| **Tipo** | Determinístico | Probabilístico |
| **Gaps/indels** | No modela bien | Sí (estados I y D) |
| **Variabilidad** | Fija (lista corta) | Continua (prob. de emisión) |
| **Sensibilidad** | Baja | Alta |
| **Especificidad** | Alta | Menor |
| **Velocidad** | Muy rápida | Más lenta |

La expresión regular es más específica porque es binaria (cumple o no cumple). No hay "score intermedio" → menos FP pero más FN.
El HMM es más sensible porque tolera variaciones → detecta homólogos divergentes.

**Cómo escribir una expresión regular desde un logo:**
- Posición con una letra dominante → escribir esa letra: `G`
- Posición con dos opciones → lista: `[AG]`
- Posición variable → comodín: `x`
- N posiciones variables → `x(N)` o `x(min,max)`

---

## Tips y errores comunes

### Temas garantizados

1. **Pasos de BLAST** — los 4: word list → búsqueda en índice → extensión sin gaps → extensión con gaps + E-value.
2. **Identities vs Positives** — Identities = exactamente iguales; Positives = BLOSUM62 > 0 (incluye cambios conservativos).
3. **3 problemas del HMM** — Entrenamiento / Scoring / Alineamiento. Con input y output de cada uno.
4. **MSA + UPGMA** — construir la matriz, hacer UPGMA, agregar secuencias al perfil en orden del árbol.
5. **VP/VN/FP/FN + efecto del umbral** — tabla de 2×2, qué métrica afecta cada tipo.
6. **Contigs vs scaffolds vs missing genes vs proteínas hipotéticas** — diferencias claras.

### Confusiones frecuentes

| Error | Corrección |
|-------|-----------|
| Confundir **profundidad** con **cobertura** | Profundidad = lecturas por base; cobertura = % del genoma sobre un umbral |
| Confundir **PSI-BLAST** con BLAST estándar | PSI-BLAST itera con PSSM; detecta homólogos lejanos. BLAST usa BLOSUM62 fija |
| Confundir **expresión regular** con HMM | Expr. reg. = determinístico, más específico. HMM = probabilístico, más sensible |
| Confundir **inserción** con **deleción** | Inserción = la secuencia tiene de más. Deleción = le falta algo comparado con el modelo |
| Confundir **UPGMA** con otros métodos filogenéticos | UPGMA se usa para el árbol guía del MSA en ClustalW. No es lo mismo que NJ o máxima verosimilitud |
| Confundir **missing gene** con **proteína hipotética** | Missing = no se encontró el ORF. Hipotética = se encontró el ORF pero sin función conocida |

---

# PARTE 2 — Ejercicios resueltos por año

---

## Parcial 2018 — Resolución completa

### Problema 1 — NW (alineamiento global)

**Secuencias:** CTGGCT (filas) y ATGCTG (columnas)
**Scoring:** Match = +2, Mismatch = 0, Gap = −1

**Matriz completa:**

|   | − | A | T | G | C | T | G |
|---|---|---|---|---|---|---|---|
| **−** | 0 | −1 | −2 | −3 | −4 | −5 | −6 |
| **C** | −1 | 0 | −1 | −2 | −1 | −2 | −3 |
| **T** | −2 | −1 | 2 | 1 | 0 | 1 | 0 |
| **G** | −3 | −2 | 1 | 4 | 3 | 2 | 3 |
| **G** | −4 | −3 | 0 | 3 | 3 | 2 | 4 |
| **C** | −5 | −4 | −1 | 2 | 5 | 4 | 3 |
| **T** | −6 | −5 | −2 | 1 | 4 | 7 | 6 |

**Score final = 6** (esquina inferior derecha → celda T/G)

Justificación de celdas del traceback:

- F(T,T) fila 6, col 5: T=T match → diagonal F(5,4)+2 = 5+2 = **7** ← elegido
- F(T,G) fila 6, col 6: mismatch → arriba F(5,6)−1 = 7−1 = **6** ← elegido
- F(G,G) fila 4, col 3: empate diagonal (F(3,2)+2=1+2=3) e izquierda (F(4,2)−1=4−1=3) → **dos caminos**

**¿Existe más de un alineamiento óptimo? Sí — exactamente dos:**

```
Alineamiento A (diagonal en el empate):
  CTGGCT−
  AT−GCTG

Alineamiento B (izquierda en el empate):
  CTG−GCT−
  A T G C T G

Score A: 0+2−1+2+2+2−1 = 6 ✓
Score B: 0+2+2−1+2+2−1 = 6 ✓
```

---

### Problema 2 — Bases de Datos Primarias

Clasificar una BD de cada tipo según curación, redundancia y tipo (primaria/secundaria):

| Categoría | BD elegida | Curada | Redundante | Prim/Sec | Justificación |
|-----------|-----------|--------|------------|----------|---------------|
| **a) Nucleótidos** | RefSeq | Sí (parcial) | No | Primaria | Un registro canónico por molécula/organismo; NCBI revisa manualmente |
| **b) Proteínas** | Swiss-Prot | Sí | No | Primaria | Revisión manual experta; 1 entrada por proteína; función verificada |
| **c) Expresión génica** | GEO | No | Sí | Primaria | Autores depositan directamente; mismo gen aparece en miles de experimentos |
| **d) Organismo** | SGD | Sí | No | **Secundaria** | Integra GenBank + UniProt + literatura para S. cerevisiae; no genera experimentos propios |
| **e) Compuestos** | PubChem | No | Sí | Primaria | La más grande de libre acceso; sin revisión manual sistemática |
| **f) Literatura** | Europe PMC | No | No | Primaria | Texto completo open access + preprints; tiene API para minería de texto |
| **g) Libre** | PDB | Sí | Sí | Primaria | Estructuras 3D; curada (validación R-factor/Ramachandran obligatoria); redundante (cientos de estructuras por proteína) |

---

### Problema 3 — MSA, SP score, E-values, Prolina

**MSA dado** — motivo N-glicosilación N−X(≠P)−S/T, gap affine: apertura=−8, extensión=−3 (gap de largo 1 = −8+(−3)= −11):

```
Prot1: L  E  N  K  T  V  A
Prot2: V  −  N  E  S  Y  A
Prot3: I  D  N  Q  T  I  A
```

**Scores por par (BLOSUM62):**

| Par | L/V I/V | E/− E/D | N/N | K/E K/Q | T/S T/T | V/Y V/I | A/A | **Total** |
|-----|---------|---------|-----|---------|---------|---------|-----|-----------|
| P1–P2 | L/V = +1 | E/− = −11 | N/N = +6 | K/E = +1 | T/S = +1 | V/Y = −1 | A/A = +4 | **1** |
| P1–P3 | L/I = +2 | E/D = +2 | N/N = +6 | K/Q = +1 | T/T = +5 | V/I = +3 | A/A = +4 | **23** |
| P2–P3 | V/I = +3 | −/D = −11 | N/N = +6 | E/Q = +2 | S/T = +1 | Y/I = −1 | A/A = +4 | **4** |

**SP score total = 1 + 23 + 4 = 28**

**Árbol guía UPGMA:** mayor score = P1–P3 (23) → se unen primero → árbol: `((Prot1, Prot3), Prot2)`.

**Orden de E-values:** S₁₃=23 > S₂₃=4 > S₁₂=1  →  **E₁₃ < E₂₃ < E₁₂**
(mayor score = menor probabilidad de verlo por azar = E-value más bajo = más significativo)

**¿Por qué la Prolina no se acepta en posición X del motivo N−X(≠P)−S/T?**

La Prolina tiene un anillo pirrolidínico que une la cadena lateral al nitrógeno del backbone:
- Elimina el H del N amídico → no puede formar puentes de hidrógeno con el backbone.
- Restringe severamente el ángulo φ → fuerza un kink en la cadena.
- Impide que el tripéptido adopte la conformación de β-turn que necesita la enzima OST (oligosacariltransferasa) para reconocer la Asn y añadir el oligosacárido.

BLOSUM62 confirma: casi todos los scores de la fila P son negativos → la Pro es evolutivamente muy poco intercambiable, lo que refleja su excepcionalidad estructural.

---

### Problema 4 — HMM y PFAM (quinasa KC)

**Pipeline completo (KC quinasa, búsqueda en proteoma de Av):**

```
Sustratos confirmados de KC (verdaderos positivos)
         ↓
MSA de las regiones de reconocimiento (CLUSTAL / MAFFT)
         ↓
hmmbuild → KC_motivo.hmm          ← ENTRENAMIENTO
  estima P(aa|posición), P(M→I), P(M→D)
         ↓
Calibración (automática) → parámetros para calcular E-values
         ↓
hmmsearch en proteoma de Av → lista de candidatos con E-value   ← SCORING
         ↓
Validar experimentalmente 100 candidatos
         ↓
Evaluación: sensibilidad, especificidad, curva ROC
```

**Interpretación del logo:**
- Columna alta → determinante de especificidad de KC; mutar ese residuo elimina el reconocimiento.
- Letra dominante revela tipo fisicoquímico: K/R → basofílica; D/E → acidofílica.
- Columna baja → KC no impone restricción → equivalente a "X" en la expresión regular.
- Número de columnas = estados MATCH del HMM = largo del motivo.

**Resultado experimental:** 5 de 100 candidatos no son sustratos → **Falsos Positivos**.
- Precisión = 95/100 = 95%.
- Causas: motivo enterrado en estructura 3D, sobreajuste del modelo con set de entrenamiento pequeño, dependencia de contexto celular no capturado por secuencia lineal.

Para aumentar especificidad: subir el umbral de E-value → los FP (scores más bajos) quedan por debajo del umbral → menos FP, pero se pierde algún VP marginal (más FN).

---

## Parcial 2022 — Problemas 2 y 4

### Problema 2 — MSA (SeqA=TCTCA, SeqB=TGTGATTGT, SeqC=ACTA)

Scoring: Match=+1, Mismatch=−1, Gap=−2

**Scores por par:**

| Par | Score |
|-----|-------|
| A–B | −1 |
| A–C | 1 |
| B–C | −4 |

**UPGMA:** mayor score = A–C (1) → unir A+C primero.
- d(A+C, B) = [d(A,B) + d(C,B)] / 2 = [−1 + (−4)] / 2 = −2.5
- Árbol guía: `((A, C), B)`

**MSA resultante:**
```
Primero alinear A con C:
A: T C T − C A G T
C: A C T A C A − −

Luego agregar B al perfil (A+C):
A: T C T − C A − − G T
C: A C T A C A − − − −
B: T G T − G A T T G T

SP score = suma de scores de todos los pares en cada columna
```

---

### Problema 4 — HMM e inserción/deleción (hNeu1)

**¿Y192 de hNeu1 es inserción o deleción respecto del modelo HMM?**

Es una **inserción**: la Y192 no aparece en ninguna columna MATCH del logo ni en las otras secuencias del alineamiento. hNeu1 tiene un residuo "extra" que el modelo no espera → estado Inserción del HMM.

> Regla: **inserción** = la secuencia tiene algo que el modelo NO tiene. **Deleción** = le falta algo que el modelo espera.

**Clasificación de cNeu1 y gNeu1:**
- cNeu1 tiene el sitio G-Y-x(1-2)-[VIL] pero no se fosforila → **Falso Positivo (FP)**.
- gNeu1 no tiene el sitio y no se fosforila → **Verdadero Negativo (VN)**.

**¿Qué 3 problemas resolvió PFAM para encontrar dominios en hNeu1?**
1. **Entrenamiento:** construyó un HMM por familia (Activin_recp, TGF_beta_GS, PK_Tyr_Ser-Thr) desde alineamientos semilla.
2. **Scoring:** evaluó la secuencia de hNeu1 contra cada HMM y calculó scores.
3. **Alineamiento:** alineó hNeu1 contra los modelos para reportar qué regiones corresponden a cada dominio y en qué posiciones.

---

## Parcial 2023 — Problema 3 (NGS/Ensamblado)

**a) Pasos para formar contigs y scaffolds:**

Contigs:
1. Fragmentar las lecturas en k-mers.
2. Construir el grafo de De Bruijn (nodos = k-mers, aristas = lecturas).
3. Resolver el grafo: caminos de bajo peso = errores (eliminar); ciclos = repeticiones.
4. Los caminos lineales resueltos son los contigs.

Scaffolds:
- Lecturas paired-end: dos lecturas del mismo fragmento de tamaño conocido.
- Si lectura₁ mapea en contig A y lectura₂ en contig B → A y B están orientados y separados por ≈ inserto.
- Con genoma de referencia: mapear contigs contra la referencia → orden y posición directamente.

**b) Regiones repetitivas:**
Los k-mers de la repetición son compartidos entre las copias → bifurcaciones/ciclos en el grafo → el ensamblador no puede resolverlos → corta los contigs en los bordes → más contigs, menor N50.
Solución: lecturas largas (PacBio, Nanopore) que abarcan la repetición entera; mate-pair con inserto grande.

**c) ¿Qué significa profundidad de 50×?**
Cada base fue leída 50 veces en promedio (distribución Poisson).
Baja cobertura en algunas regiones: GC bias, estructuras secundarias del ADN, fragmentación sesgada, varianza estadística de la Poisson.
Missing genes con RAST: genes esperados que no se detectaron → probablemente en zonas de baja cobertura → no se formaron contigs → no hay ORFs en esa región.

**d) Proteínas hipotéticas vs missing genes:**
- **Hipotéticas:** hay un ORF predicho (start/stop, largo adecuado) pero BLAST no encuentra ningún homólogo → función desconocida.
- **Missing genes:** genes esperados (por comparación con la BD de RAST/FigFams) que no se encontraron en el ensamblado → ausentes del genoma ensamblado.

---

## Parcial 2024 — Problema 2 (NW + MSA)

SeqA = AAGTCA, SeqB = AAATC, SeqC = AAATA
Match=+1, Mismatch=−1, Gap=−1

**Scores por par:**
```
A vs B:  A A G T C A           Score = 1+1−1+1+1−1 = 2
         A A A T C −

A vs C:  A A G T C A           Score = 1+1−1+1−1+1 = 2
         A A A T − A

B vs C:  A A A T C             Score = 1+1+1+1−1 = 3  ← mayor (sin gaps)
         A A A T A
```

**UPGMA:** mayor score = B–C (3) → unir B+C primero.
- d(B+C, A) = [d(B,A) + d(C,A)] / 2 = [2 + 2] / 2 = 2
- Árbol: `((B, C), A)`

**MSA:** primero alinear B con C, luego agregar A al perfil (B+C).
