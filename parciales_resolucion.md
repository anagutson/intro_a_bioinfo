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

**Secuencias:** CTGGCT (columnas, horizontal) y ATGCTG (filas, vertical)
**Scoring:** Match = +2, Mismatch = 0, Gap = −1

**Matriz completa:**

|   | − | C | T | G | G | C | T |
|---|---|---|---|---|---|---|---|
| **−** | 0 | −1 | −2 | −3 | −4 | −5 | −6 |
| **A** | −1 | 0 | −1 | −2 | −3 | −4 | −5 |
| **T** | −2 | −1 | 2 | 1 | 0 | −1 | −2 |
| **G** | −3 | −2 | 1 | 4 | 3 | 2 | 1 |
| **C** | −4 | −1 | 0 | 3 | 4 | 5 | 4 |
| **T** | −5 | −2 | 1 | 2 | 3 | **4** | **7** |
| **G** | −6 | −3 | 0 | 3 | 4 | **3** | **6** |

**Score final = 6** (celda esquina inferior derecha: fila G, columna T)

**Justificación de celdas clave del traceback:**

F(T,T) — fila T (pos 5), col T (pos 6): T = T → **match**
- Diagonal: F(4,5) + 2 = 5 + 2 = **7** ← elegido
- Arriba:    F(4,6) − 1 = 4 − 1 = 3
- Izquierda: F(5,5) − 1 = 4 − 1 = 3

F(G,T) — fila G (pos 6), col T (pos 6): G ≠ T → **mismatch**
- Diagonal:  F(5,5) + 0 = 4
- Arriba:    F(5,6) − 1 = 7 − 1 = **6** ← elegido
- Izquierda: F(6,5) − 1 = 3 − 1 = 2

F(G,G) — fila G (pos 3), col G (pos 3): G = G → **empate**
- Diagonal:  F(2,2) + 2 = 1 + 2 = **3**
- Izquierda: F(3,2) − 1 = 4 − 1 = **3**
→ dos caminos posibles → dos alineamientos óptimos

**¿Existe más de un alineamiento óptimo? Sí — exactamente dos:**

**Alineamiento A** (traceback diagonal en el empate de F(G,G)):
```
seq1 (CTGGCT): C  T  G  G  C  T  −
seq2 (ATGCTG): A  T  −  G  C  T  G
```

Verificación columna por columna:

| Col | seq1 | seq2 | Score |
|-----|------|------|-------|
| 1 | C | A | mismatch = 0 |
| 2 | T | T | match = +2 |
| 3 | G | − | gap = −1 |
| 4 | G | G | match = +2 |
| 5 | C | C | match = +2 |
| 6 | T | T | match = +2 |
| 7 | − | G | gap = −1 |
| **Total** | | | **0+2−1+2+2+2−1 = 6 ✓** |

**Alineamiento B** (traceback izquierda en el empate de F(G,G)):
```
seq1 (CTGGCT): C  T  G  G  C  T  −
seq2 (ATGCTG): A  T  G  −  C  T  G
```

Verificación columna por columna:

| Col | seq1 | seq2 | Score |
|-----|------|------|-------|
| 1 | C | A | mismatch = 0 |
| 2 | T | T | match = +2 |
| 3 | G | G | match = +2 |
| 4 | G | − | gap = −1 |
| 5 | C | C | match = +2 |
| 6 | T | T | match = +2 |
| 7 | − | G | gap = −1 |
| **Total** | | | **0+2+2−1+2+2−1 = 6 ✓** |

Ambos tienen score = 6. Los dos son alineamientos globales óptimos válidos.

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

**MSA dado** — motivo N-glicosilación N−X(≠P)−S/T
Gap affine: apertura = −8, extensión = −3 → gap de largo 1 = −8 + (−3×1) = **−11**

```
Prot1: L  E  N  K  T  V  A
Prot2: V  −  N  E  S  Y  A
Prot3: I  D  N  Q  T  I  A
```

**Scores por par (BLOSUM62):**

| Pos | P1–P2 | P1–P3 | P2–P3 |
|-----|-------|-------|-------|
| 1 | L/V = +1 | L/I = +2 | V/I = +3 |
| 2 | E/− = −11 | E/D = +2 | −/D = −11 |
| 3 | N/N = +6 | N/N = +6 | N/N = +6 |
| 4 | K/E = +1 | K/Q = +1 | E/Q = +2 |
| 5 | T/S = +1 | T/T = +5 | S/T = +1 |
| 6 | V/Y = −1 | V/I = +3 | Y/I = −1 |
| 7 | A/A = +4 | A/A = +4 | A/A = +4 |
| **Total** | **1** | **23** | **4** |

**SP score total = S(1,2) + S(1,3) + S(2,3) = 1 + 23 + 4 = 28**

Verificación columna a columna:

| Pos | S(1,2) | S(1,3) | S(2,3) | Subtotal |
|-----|--------|--------|--------|----------|
| 1 | 1 | 2 | 3 | 6 |
| 2 | −11 | 2 | −11 | −20 |
| 3 | 6 | 6 | 6 | 18 |
| 4 | 1 | 1 | 2 | 4 |
| 5 | 1 | 5 | 1 | 7 |
| 6 | −1 | 3 | −1 | 1 |
| 7 | 4 | 4 | 4 | 12 |
| **Total** | | | | **28 ✓** |

**Árbol guía UPGMA:** S(1,3)=23 es el mayor → Prot1 y Prot3 se unen primero.
```
        ┌── Prot1
   ┌────┤
   │    └── Prot3
───┤
   └──────── Prot2
```
Prot1 y Prot3 comparten sustituciones conservativas sin gaps: L↔I (alifáticos), E↔D (ácidos), T↔T, V↔I. Prot2 se aleja por el gap en posición 2 (−11) y Y↔V (−1).

**Orden de E-values:** S₁₃=23 > S₂₃=4 > S₁₂=1  →  **E₁₃ < E₂₃ < E₁₂**

E-value = K·m·n·e^(−λS). A mayor score → menor exponente → menor E-value → más significativo. El par 1–3 es el menos probable de verse por azar.

**¿Por qué la Prolina no se acepta en posición X del motivo N−X(≠P)−S/T?**

La Prolina tiene un anillo pirrolidínico que forma un enlace covalente entre la cadena lateral y el nitrógeno del backbone. Esto produce tres efectos:

1. **Elimina el H del N amídico** → no puede formar puentes de hidrógeno con el backbone.
2. **Restringe severamente el ángulo φ** → fuerza un kink (quiebre) en la cadena polipeptídica.
3. **Impide la conformación de β-turn** que necesita la enzima oligosacariltransferasa (OST) para reconocer el motivo N-X-S/T y transferir el oligosacárido a la Asn.

BLOSUM62 confirma la excepcionalidad de la Pro: casi todos los scores de su fila son negativos (P–L = −3, P–I = −3, P–F = −4, solo P–P = +7). Esto refleja que la Pro es evolutivamente muy poco intercambiable, exactamente porque sus sustituciones destruyen la estructura local.

---

### Problema 4 — HMM y PFAM (quinasa KC)

**Pipeline completo para buscar sustratos en el proteoma de Av:**

```
Sustratos confirmados de KC (verdaderos positivos experimentales)
         ↓
MSA de las regiones de reconocimiento (CLUSTAL / MAFFT)
         ↓                      ← PROBLEMA DE ENTRENAMIENTO
hmmbuild → KC_motivo.hmm
  Estima: P(aa | posición) por columna   → prob. de emisión
          P(M→M), P(M→I), P(M→D)        → prob. de transición
         ↓
Calibración automática → parámetros λ y K para calcular E-values
         ↓                      ← PROBLEMA DE SCORING
hmmsearch -E 0.001 KC_motivo.hmm proteoma_Av.fasta
  Produce lista de candidatos ordenados por score y E-value
         ↓
Validación experimental: 100 candidatos → ensayos de unión y fosforilación
         ↓
Evaluación: curva ROC (sensibilidad vs 1−especificidad), ajuste del umbral
```

**Análisis del logo HMM (KC quinasa) — respuesta completa de parcial:**

**1. Largo del motivo y estados MATCH**

El logo tiene ~65 columnas → el HMM tiene 65 estados MATCH → el motivo reconocido por KC mide ~65 aa.
Sin embargo, los extremos del logo tienen columnas muy bajas (poca información) → son posiciones poco conservadas que se pueden recortar. El núcleo funcional es de ~30 aa (las columnas con información apreciable en el centro del logo).

**2. Triptófanos (W) fuertemente conservados**

Hay varias posiciones con W (triptófano) como residuo dominante y columna alta. Esto es muy llamativo porque W es el aminoácido más grande, aromático, con anillo indol.

*¿Qué significa?*
- Las posiciones con W conservado probablemente forman una **superficie aromática** que encaja en el dominio beta-glove de KC (interacciones de apilamiento π-π o CH-π entre el W del sustrato y residuos aromáticos de KC).
- W en el núcleo de la proteína estabiliza el plegamiento por empaquetamiento hidrofóbico aromático.
- La conservación extrema de W indica que **no tolera sustitución**: si se reemplaza W por otro aminoácido (incluso F o Y, también aromáticos pero más pequeños), KC probablemente pierde la unión. Esto se puede validar por mutagénesis dirigida.

*Experimento de validación:* mutar W → F (sustitución conservativa, mismo anillo aromático pero sin el N del indol) y luego W → A (sustitución disruptiva, elimina el anillo) → medir fosforilación in vitro. Si W→F mantiene actividad y W→A la elimina, confirma el rol del anillo; si ambas eliminan actividad, confirma que el indol específico es necesario.

**3. Dos argininas (R) conservadas**

Dos columnas con R (arginina) bien definidas y de altura considerable.

*¿Qué significa?*
- R tiene carga positiva (pKa ≈ 12.5, siempre cargada a pH fisiológico) y puede:
  - Interactuar con residuos ácidos (D/E) del sitio activo o de superficie de KC (puentes salinos).
  - Interactuar con el grupo fosfato del ATP durante la transferencia (mecanismo catalítico).
  - Formar la "señal de reconocimiento básico" si KC es una quinasa basofílica (como PKA, que reconoce R-R-X-S/T).
- La posición relativa de las dos R en el motivo puede ser clave: si están a distancia fija del sitio de fosforilación (S/T), KC las usa como ancla de posicionamiento.

**4. Distancias fijas entre residuos conservados**

La disposición regular de columnas altas a lo largo del logo indica que el espaciado entre W y R es constante en todos los sustratos.

*¿Qué significa?*
- El reconocimiento no es solo por identidad de aminoácido sino por **patrón espacial**: KC necesita que W esté a X posiciones de R.
- Esto implica que la estructura 3D del sustrato en la región del motivo debe adoptar una conformación particular (ej: una hélice α posiciona los residuos en la misma cara cada ~3.5 posiciones; una lámina β los posiciona alternos).
- Una expresión regular podría reflejar esto: `W-x(n)-R-x(m)-R-x(p)-[ST]` donde n, m, p son distancias fijas.

**5. Región central variable con probabilidad de inserción/deleción**

En el centro del logo hay una región con columnas de poca altura y, si se hace zoom, se ve la señal de inserción/deleción del HMM.

*¿Qué significa?*
- Esta región **no es reconocida directamente por KC**: es un loop o región flexible entre los elementos estructurales que sí importan.
- El HMM la modela con estados I (inserción) y D (deleción): algunos sustratos tienen residuos extra aquí, otros tienen menos → el largo de esa región varía entre sustratos.
- Esto es una ventaja del HMM sobre la expresión regular: la expresión regular no puede capturar bien esta variabilidad de largo; el HMM sí (estados I/D con sus probabilidades).

**Síntesis — ¿qué dice el logo sobre el mecanismo de reconocimiento de KC?**

KC reconoce un motivo estructurado de ~30 aa (núcleo) con:
- Una o varias posiciones W que forman la interfaz hidrofóbica aromática con el dominio beta-glove de KC.
- Dos posiciones R que anclan el sustrato por interacciones electrostáticas y posicionan el sitio S/T para la transferencia del fosfato.
- Un espaciado fijo entre W y R que refleja la estructura 3D requerida.
- Un loop central de largo variable que no contacta KC directamente.

Esto caracteriza a KC como una quinasa con **especificidad de estructura terciaria**, no solo de secuencia primaria — lo que explica los falsos positivos: un sustrato puede tener la secuencia correcta pero si el loop o una estructura secundaria entierra el motivo, KC no puede acceder.

**Resultado experimental:** 5 de 100 candidatos no son sustratos → **Falsos Positivos (FP)**
- Precisión = 95/100 = 95%.
- Sensibilidad = VP / (VP + FN); Especificidad = VN / (VN + FP).

Causas de los FP:
1. El motivo existe en la secuencia pero está **enterrado en la estructura 3D** → KC no puede acceder in vivo.
2. Set de entrenamiento pequeño → el modelo sobreajustó y aprendió señales espurias.
3. KC requiere co-localización, andamiajes o fosforilaciones previas no capturables por secuencia lineal.
4. Umbral de E-value permisivo → hits marginalmente similares al modelo.

Para aumentar especificidad: subir el umbral de E-value → los FP (scores más bajos) quedan por debajo del corte → menos FP, pero algunos VP marginales pasan a ser FN.

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

---

## Parcial (año desconocido) — Ensamblado, Mapeo y Anotación

### Ejercicio 1 — Grafo de De Bruijn con k=3

**Lecturas:**

| Lectura | Veces |
|---------|-------|
| ACTTTC | 15 |
| ACTGTC | 1 |
| TTCTGG | 14 |
| TCTGGT | 1 |
| TCTGGA | 8 |
| ACTGGA | 7 |

---

#### a) Construcción del grafo de De Bruijn (k=3)

Con k=3: nodos = 2-mers (k−1), aristas = 3-mers (k).
Cada lectura de largo 6 genera 4 k-mers. La frecuencia de cada arista = suma de frecuencias de las lecturas que la generan.

**Todos los k-mers y sus frecuencias:**

| k-mer | Arista | Frecuencia | De qué lecturas viene |
|-------|--------|------------|----------------------|
| CTG | CT → TG | 31 | ACTGTC(1) + TTCTGG(14) + TCTGGT(1) + TCTGGA(8) + ACTGGA(7) |
| TGG | TG → GG | 30 | TTCTGG(14) + TCTGGT(1) + TCTGGA(8) + ACTGGA(7) |
| TTC | TT → TC | 29 | ACTTTC(15) + TTCTGG(14) |
| ACT | AC → CT | 23 | ACTTTC(15) + ACTGTC(1) + ACTGGA(7) |
| TCT | TC → CT | 23 | TTCTGG(14) + TCTGGT(1) + TCTGGA(8) |
| CTT | CT → TT | 15 | ACTTTC(15) |
| TTT | TT → TT | 15 | ACTTTC(15) |
| GGA | GG → GA | 15 | TCTGGA(8) + ACTGGA(7) |
| TGT | TG → GT | **1** | ACTGTC(1) ← posible error |
| GTC | GT → TC | **1** | ACTGTC(1) ← posible error |
| GGT | GG → GT | **1** | TCTGGT(1) ← posible error |

**Grafo completo:**

```
                 CTT(15)      TTT(15, self)
         ┌──────────────► TT ◄────────┐
         │                │           │
AC ─ACT(23)─► CT           └─TTC(29)─► TC ─TCT(23)─► CT
                │                                      ▲
                └─CTG(31)─► TG ─TGG(30)─► GG ─GGA(15)─► GA
                                │              │
                               TGT(1)         GGT(1)
                                │              │
                                ▼              ▼
                               GT ─GTC(1)─► TC
```

---

#### b) Camino óptimo y secuencia consenso

**Criterio:** eliminar aristas con frecuencia = 1 (TGT, GTC, GGT) → son errores de secuenciación (ver inciso d).

**Grafo simplificado** (sin aristas de frecuencia 1):

```
AC ─ACT─► CT ─CTT─► TT ─TTT(self)─► TT ─TTC─► TC ─TCT─► CT ─CTG─► TG ─TGG─► GG ─GGA─► GA
```

Grados en el grafo simplificado:
- AC: entrada=0, salida=1 → **nodo inicio**
- CT: entrada=2 (ACT, TCT), salida=2 (CTT, CTG) → balanceado
- TT: entrada=2 (CTT, TTT-self), salida=2 (TTT-self, TTC) → balanceado
- TC: entrada=1 (TTC), salida=1 (TCT) → balanceado
- TG: entrada=1 (CTG), salida=1 (TGG) → balanceado
- GG: entrada=1 (TGG), salida=1 (GGA) → balanceado
- GA: entrada=1 (GGA), salida=0 → **nodo fin**

→ Existe exactamente un **camino euleriano** (inicio en AC, fin en GA):

```
AC →ACT→ CT →CTT→ TT →TTT→ TT →TTC→ TC →TCT→ CT →CTG→ TG →TGG→ GG →GGA→ GA
```

**Secuencia consenso:** cada nodo aporta su primer carácter + el último carácter de la última arista:

```
AC + T + T + T + C + T + G + G + A = ACTTTCTGGA
```

**Verificación:** los reads de alta frecuencia son subcadenas de ACTTTCTGGA:
- ACTTTC(15): posiciones 1–6 ✓
- TTCTGG(14): posiciones 4–9 ✓
- TCTGGA(8): posiciones 5–10 ✓

**Justificación de la elección:** este es el único camino euleriano posible en el grafo simplificado. Traversa todas las aristas de alta frecuencia exactamente una vez. Maximiza la cobertura de las lecturas confiables.

---

#### c) ¿Cuántos nodos tiene el grafo vs el total posible?

Con k=3, los nodos son 2-mers. Total posible: 4² = **16**.
Nodos presentes en el grafo: AC, CT, TT, TC, TG, GT, GG, GA = **8 nodos**.

El grafo usa 8/16 = 50% de los nodos posibles.

**Ventaja bioinformática:** el grafo De Bruijn solo almacena los k-mers que realmente aparecen en las lecturas. Para genomas grandes con k=31 (estándar en ensambladores como SPAdes), existen 4³¹ ≈ 4.6×10¹⁸ k-mers posibles pero solo una fracción ínfima aparece en los datos. Esto hace que el algoritmo sea **eficiente en memoria**: la tabla hash solo guarda los k-mers observados, no el espacio completo. Cuanto mayor el k, menor la proporción de k-mers posibles que están presentes → mayor ahorro relativo.

---

#### d) Histograma de frecuencias y eliminación de errores

**Histograma:**

```
Frecuencia
   31 |  █
   30 |  █
   29 |  █
   23 |  █  █
   15 |  █  █  █
    1 |  █  █  █
       CTG TGG TTC ACT TCT CTT TTT GGA TGT GTC GGT
```

| Frecuencia | k-mers |
|-----------|--------|
| 31 | CTG |
| 30 | TGG |
| 29 | TTC |
| 23 | ACT, TCT |
| 15 | CTT, TTT, GGA |
| **1** | **TGT, GTC, GGT** ← eliminar |

**K-mers a eliminar:** TGT, GTC, GGT (frecuencia = 1).

**¿Por qué?** Con una profundidad promedio de 15–31×, un k-mer verdadero debería aparecer muchas veces. Un k-mer con frecuencia = 1 es casi con certeza un **error de secuenciación**: si una base fue leída incorrectamente en una sola lectura, genera k-mers únicos que no aparecen en ninguna otra lectura. La probabilidad de que exactamente el mismo error ocurra dos veces es mínima.

**Lecturas erróneas correspondientes:** ACTGTC(1) y TCTGGT(1) — ambas aparecen exactamente una vez. Comparten sus k-mers correctos (ACT, CTG, etc.) con otras lecturas, pero los k-mers únicos que aportan (TGT, GTC, GGT) delatan que son errores.

**¿Se simplifica el grafo?** Sí:
- Al eliminar TGT, GTC, GGT, el nodo GT queda sin aristas → se elimina.
- El grafo pasa de 11 aristas y 8 nodos a **8 aristas y 7 nodos**.
- Más importante: los nodos TG, GG y TC pasan a ser balanceados → el grafo tiene un único camino euleriano (AC → GA), lo que permite reconstruir la secuencia consenso sin ambigüedad.

---

### Ejercicio 2 — Definiciones

**a) Genoma Humano de Referencia**

Es una secuencia consenso del genoma humano diploide construida a partir de un número pequeño de individuos anónimos. No representa el genoma de ninguna persona en particular, sino una referencia coordinada y estable. Se usa como sistema de coordenadas universal para: mapear lecturas de secuenciación, reportar la posición de variantes, anotar genes, y comparar estudios entre laboratorios. La versión actual es GRCh38/hg38. Fue construido por el IHGSC (Human Genome Project) y es mantenido por el Genome Reference Consortium (GRC). Su limitación principal: no captura la diversidad de haplotipos humanos (es una secuencia mosaico de muy pocos individuos).

**b) Variante**

Diferencia entre la secuencia de un individuo y el genoma de referencia (o entre individuos de la misma especie). Tipos principales:

| Tipo | Descripción | Ejemplo |
|------|-------------|---------|
| SNV / SNP | Cambio de una sola base | A→G en posición 12.345 |
| Indel | Inserción o deleción pequeña (1–50 bp) | +AT en posición 6.789 |
| CNV | Variación en el número de copias de una región | Duplicación del cromosoma 17q |
| Variante estructural | Inversión, translocación, fusión | Fusión BCR-ABL en leucemia |

Las variantes se clasifican clínicamente como: Patogénica, Probablemente patogénica, VUS (Variant of Uncertain Significance), Probablemente benigna, Benigna. Esta clasificación usa bases de datos como ClinVar, gnomAD y OMIM.
