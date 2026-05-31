# Apunte Bioinformática — Examen a Libro Abierto

---

# ÍNDICE RÁPIDO

1. [Marco conceptual y estructura del parcial](#1-marco-conceptual)
2. [Alineamiento por pares: NW / SW](#2-alineamiento-nw-sw)
3. [BLAST y PSI-BLAST](#3-blast)
4. [MSA — CLUSTALW y SP score](#4-msa)
5. [Bases de datos primarias](#5-bases-de-datos-primarias)
6. [Bases de datos secundarias — HMM, PROSITE, InterPro](#6-bases-de-datos-secundarias)
7. [Matrices de sustitución — PAM vs BLOSUM](#7-matrices-pam-blosum)
8. [Aminoácidos](#8-aminoacidos)
9. [NGS, Ensamblado, Genómica Humana](#9-ngs)
10. [Resolución del Parcial 2018](#10-parcial-2018)

---

## 1. Marco Conceptual

La bioinformática existe en la intersección entre biología y ciencias de la computación. Antes de entrar a los algoritmos, conviene entender cómo se organizan los niveles de abstracción de la disciplina y qué tipo de preguntas se hacen en el parcial.

### Jerarquía: Teoría → Algoritmo → Programa

| Nivel | Definición | Ejemplo |
|-------|-----------|---------|
| **Teoría** | Modelo que describe el fenómeno biológico | Evolución por sustitución |
| **Algoritmo** | Serie de pasos lógico-matemáticos, independiente del lenguaje | Needleman-Wunsch |
| **Programa** | Implementación en un lenguaje específico | BLAST, CLUSTALW, hmmbuild |
| **Dato** | Lo medido experimentalmente; input del algoritmo | Secuencia FASTA |

> En el parcial se puede responder con el nombre del **algoritmo** O con el nombre del **programa** — ambas respuestas son equivalentes. "Haría un alineamiento de pares" = "usaría NW/SW" = "usaría BLAST".

### Estructura típica del parcial

| Problema | Tema |
|----------|------|
| **P1** | Alineamiento / Programación Dinámica (NW o SW) |
| **P2** | Bases de datos primarias (criterios: curación, redundancia, primaria/secundaria) |
| **P3** | MSA y árbol guía (CLUSTALW, SP score, E-value, PROSITE) |
| **P4** | Bases de datos secundarias (HMM, logos, sensibilidad/especificidad) |

### Tips del profe (de la clase de repaso)

- Leer el parcial **completo** antes de empezar.
- Responder de lo más **fácil a lo más difícil**.
- La respuesta **simple es la correcta** — el parcial no tiene trampa.
- Cuando veas un logo, fijarse si las letras son aminoácidos (A/R/N...) o nucleótidos (A/T/G/C). ¡Son distintos!
- Simple answer = correct answer. No complicar de más.

### Principio de "culpa por asociación"

Si una secuencia desconocida se parece mucho a un grupo de genes conocidos, se predice que tiene la misma función **sin necesidad de experimentar directamente**. Es el principio fundamental detrás de las búsquedas en BD.

---

## 2. Alineamiento NW / SW

### Introducción: el problema de comparar dos secuencias

Dos secuencias homólogas (con un ancestro común) suelen diferir porque a lo largo de la evolución se acumulan sustituciones, inserciones y deleciones. Para detectar esa homología hay que **alinear** las secuencias: colocarlas una sobre la otra de manera que las posiciones equivalentes queden en la misma columna, introduciendo guiones (gaps) donde hubo inserciones/deleciones.

El problema de encontrar el **mejor alineamiento posible** se resuelve con programación dinámica. La idea clave es que el alineamiento óptimo de dos prefijos de longitud i y j puede calcularse a partir del alineamiento óptimo de prefijos más cortos → se construye la solución celda a celda en una matriz.

Existen dos variantes del mismo principio según qué pregunta se quiere responder:
- **NW** (Needleman-Wunsch, 1970): alineamiento **global** — fuerza a alinear las secuencias de punta a punta. Útil cuando las dos secuencias son similares en largo y función.
- **SW** (Smith-Waterman, 1981): alineamiento **local** — encuentra la región de mayor similitud. Útil cuando las secuencias son divergentes o tienen dominios distintos.

### Comparación rápida

| | **Needleman-Wunsch** | **Smith-Waterman** |
|--|--|--|
| Tipo | **Global** (toda la secuencia) | **Local** (mejor región) |
| Cuándo usar | Secuencias similares en longitud | Dominios en secuencias divergentes |
| Diferencia en recurrencia | — | **+ opción 0** (evita scores negativos) |
| Inicio traceback | Celda (n, m) — esquina inferior derecha | Celda con **valor máximo** de toda la matriz |
| Fin traceback | Celda (0, 0) | Celda con **valor 0** |
| Complejidad | O(nm) tiempo y espacio | O(nm) tiempo y espacio |

### Inicialización

```
NW:  F(i, 0) = i · gap_penalty   (para todo i)
     F(0, j) = j · gap_penalty   (para todo j)
     F(0, 0) = 0

SW:  F(i, 0) = 0   y   F(0, j) = 0   (siempre)
     F(0, 0) = 0
```

### Recurrencia NW

```
F(i,j) = max {
    F(i-1, j-1) + s(xᵢ, yⱼ)    ← diagonal (match/mismatch)
    F(i-1, j)   + gap            ← arriba    (gap en Y)
    F(i, j-1)   + gap            ← izquierda (gap en X)
}
```

### Recurrencia SW (idéntica + opción 0)

```
F(i,j) = max {
    0                             ← NUNCA baja de cero
    F(i-1, j-1) + s(xᵢ, yⱼ)
    F(i-1, j)   + gap
    F(i, j-1)   + gap
}
```

### Gap lineal vs gap affine

```
Gap lineal:   W_k = gap · k
Gap affine:   W_k = a + b·k      (a = costo apertura, b = costo extensión)
Biológico: a >> b — abrir es caro, extender es barato
```

Con gap affine: venir desde la izquierda (misma fila) = extender gap (costo b). Venir desde la diagonal = abrir gap nuevo (costo a).

### Múltiples alineamientos óptimos

Ocurren cuando **dos opciones de la recurrencia dan el mismo valor** en alguna celda del traceback. Hay que enumerar todos los caminos posibles — todos tienen el mismo score final.

Para encontrarlos:
1. Llenar la matriz de scores normalmente.
2. Al hacer traceback, en cada celda empate → bifurcar el camino.
3. Cada camino completo desde esquina a esquina = un alineamiento óptimo distinto.

### Cómo leer el alineamiento desde el traceback

| Movimiento | Significado |
|-----------|-------------|
| ↖ diagonal | match o mismatch entre xᵢ e yⱼ |
| ↑ arriba | gap en la secuencia horizontal (Y) |
| ← izquierda | gap en la secuencia vertical (X) |

---

## 3. BLAST

### Introducción: el problema de escala

NW/SW encuentran el alineamiento óptimo entre dos secuencias. Pero si quiero saber si mi proteína nueva tiene homólogos en una base de datos con millones de secuencias, aplicar SW a cada par es computacionalmente inviable. Se necesita un método que sea rápido a escala de base de datos completa, aunque sacrifique el óptimo garantizado. Esa es la razón de ser de BLAST: una **heurística** que usa una estrategia de "sembrar y extender" para lograr velocidad O(m) en la práctica, donde m es el largo de la query.

### Por qué existe BLAST

NW/SW son O(nm). Buscar una proteína de 200 aa contra 1M de proteínas = 2×10¹¹ operaciones → **inviable**. BLAST es una **heurística** que sacrifica el óptimo garantizado por velocidad, logrando O(m) en la práctica.

**Idea:** primero buscar "semillas" exactas cortas (words) que probablemente estén en el alineamiento óptimo, y solo desde ahí extender.

### Los 4 pasos del algoritmo BLAST

#### Paso 1 — Word list

Para una query de largo m y word size w=3, se extraen m−w+1 words solapadas:
```
Query: A L V G T T Y H H V D R R
Words: ALV, LVG, VGT, GTT, ...   (11 words para m=13)
```

Para cada word, se buscan **todas las palabras vecinas** de largo w con score BLOSUM62 ≥ T (default T=11):
```
Vecinas de "ALV": ALV, ALL, AVL, ... (unas ~50 para w=3, T=11)
```
→ Resultado: **word list** con ~50×(m−2) entradas.

> Efecto del umbral T: T↑ → lista más corta → más rápido, menos sensible.  
> Efecto de w: w↑ → más rápido, menos sensible.

#### Paso 2 — Búsqueda de hits exactos

La BD está **pre-indexada**: tabla que, para cada word, lista todas las posiciones de todas las secuencias donde aparece **exactamente**.

BLAST busca cada word de la word list en ese índice → **hit exacto** = par (posición en query, posición en secuencia BD).

> Innovación clave: la estructura de datos permite búsqueda en **O(1)** → total O(m) en lugar de O(nm).

#### Paso 3 — Extensión sin gaps (ungapped extension)

Desde cada hit, se extiende en ambas direcciones sumando scores BLOSUM62, **sin gaps**.

```
         ← izq  |  hit  |  der →
query:   L  V  [A  L  V]  G  T
BD seq:  L  V  [A  L  L]  G  T
                ↑ extiende acumulando score
```

Se detiene cuando el score cae demasiado. Resultado: **HSP** (High-Scoring Segment Pair). Solo se guardan HSPs con score ≥ umbral S.

#### Paso 4 — Extensión con gaps (gapped extension)

Solo para los HSPs que superaron el umbral, se hace alineamiento local con gaps (tipo SW) centrado en el HSP. Se calcula el **E-value** y se filtra.

### E-value — qué es y cómo interpretarlo

```
E = K · m · N · e^(−λS)

  m = largo de la query
  N = tamaño total de la BD (número de residuos)
  S = score del HSP
  K, λ = parámetros estadísticos (dependen de la matriz y gaps)
```

**Significado:** número esperado de hits con score ≥ S que aparecerían **por azar** en una BD de ese tamaño. No es una probabilidad, es un conteo esperado.

```
E = 1e-10  →  casi imposible por azar  → hit muy significativo
E = 0.01   →  1 hit cada 100 búsquedas al azar → aceptable
E = 1      →  1 hit por azar en esta búsqueda  → dudoso
E = 100    →  100 hits por azar                → ruido
```

> E-value depende de N (tamaño de BD) → el mismo score da distinto E-value en distintas búsquedas. Por eso se reporta E-value y no solo el score raw.
>
> A mayor score → **menor** E-value → más significativo.

### Tipos de BLAST

| Programa | Query | BD | Cuándo usarlo |
|----------|-------|----|---------------|
| **blastp** | proteína | proteínas | buscar homólogos de una proteína |
| **blastn** | nucleótidos | nucleótidos | buscar homólogos de un gen |
| **blastx** | nucl → 6 marcos | proteínas | gen nuevo, buscar proteínas homólogas |
| **tblastn** | proteína | nucl → 6 marcos | proteína conocida, buscar en genoma no anotado |
| **tblastx** | nucl 6 marcos | nucl 6 marcos | comparar genomas no anotados (muy costoso) |

> Sensibilidad para detectar homólogos remotos: **blastp > blastx > tblastn**.

### PSI-BLAST — BLAST iterativo

**Problema que resuelve:** BLAST usa BLOSUM62 genérica. PSI-BLAST aprende que ciertas posiciones de la familia son muy conservadas y las trata distinto.

```
Ronda 1:
  query → BLAST normal (BLOSUM62) → hits con E < umbral

Construcción de PSSM:
  MSA de los hits → para cada posición i, para cada aa a:
    PSSM[i][a] = log( frecuencia observada / frecuencia esperada al azar )
  → Posiciones muy conservadas: score alto para el aa correcto, muy negativo para otros
  → Posiciones variables: scores similares para muchos aa

Ronda 2:
  query → buscar en BD con PSSM[i][a] (en lugar de BLOSUM62)
  → Detecta homólogos remotos que ronda 1 perdió

Repetir hasta convergencia
```

> **Peligro:** si un falso positivo entra al MSA, contamina la PSSM y en la siguiente ronda aparecen más FP → "PSI-BLAST drift". Umbral típico para incluir en MSA: E < 0.001.

### Complejidad comparada

| Método | Complejidad | Tipo | Sensibilidad |
|--------|-------------|------|--------------|
| NW / SW | **O(nm)** | Exacto | Máxima |
| BLAST | **O(m)** | Heurístico | Alta |
| PSI-BLAST | **O(m) × iter** | Heurístico iterativo | Muy alta |
| HMMer | O(m × L) | Probabilístico | Máxima para familia |

---

## 4. MSA — CLUSTALW y SP score

### Introducción: ¿por qué alinear múltiples secuencias?

Un alineamiento de pares (NW/SW) compara dos secuencias. Pero para entender una **familia de proteínas** — detectar qué posiciones están conservadas evolutivamente, qué residuos son funcionalmente críticos, o construir un árbol filogenético — necesito alinear simultáneamente decenas o cientos de secuencias. Eso es un **MSA** (Multiple Sequence Alignment).

El problema: el MSA óptimo por programación dinámica tiene complejidad O(Lᴺ) siendo N el número de secuencias → intratable para N > 3. Se necesita una heurística. La más usada es el enfoque **progresivo**: en lugar de alinear todo junto, primero alinear los pares más parecidos y luego ir agregando los más distantes, guiándose por un árbol de similitud.

**CLUSTALW** es el algoritmo estándar de MSA progresivo. Su salida (el MSA) es la entrada para construir perfiles HMM (§6), calcular distancias evolutivas y definir árboles filogenéticos.

### CLUSTALW: algoritmo de MSA progresivo

```
Paso 1: Alinear todos los pares (pairwise) → scores → matriz de distancias N×N
Paso 2: Construir árbol guía (Neighbor-Joining) desde la matriz de distancias
Paso 3: Alinear progresivamente siguiendo el árbol
        → primero los más similares (rama más corta del árbol)
        → luego agregar los más distantes
        → los gaps introducidos en pasos previos se preservan
```

**Limitación clave:** errores de pasos tempranos se propagan y **NO se corrigen** en pasos posteriores.

### Árbol guía

- El par con **mayor score de alineamiento** (más similares → menor distancia) se une **primero**.
- Si hay 3 secuencias y S(1,3) > S(2,3) > S(1,2): el árbol une primero 1 y 3, luego agrega 2.

```
Ejemplo árbol con Prot1, Prot2, Prot3 (S13=23, S23=4, S12=1):

        ┌── Prot1
   ┌────┤           ← primeros en unirse (mayor score)
   │    └── Prot3
───┤
   └──────── Prot2  ← se agrega después (score menor)
```

### SP score (Sum of Pairs)

El SP score es la suma de los scores de **todos los pares posibles** del MSA.

```
SP = Σ S(seqᵢ, seqⱼ)   para todo i < j
```

Para 3 secuencias: SP = S(1,2) + S(1,3) + S(2,3)

**Cálculo por columna:**

Para cada columna del MSA, sumar los scores entre todos los pares de caracteres en esa columna. Luego sumar todas las columnas.

**Regla gap-gap = 0:** cuando dos gaps quedan alineados en la misma columna, su contribución al SP es 0 (no penaliza). Un gap alineado con un residuo sí penaliza.

**Ejemplo con 3 secuencias, scoring match=1, mismatch=0, gap-char=−1, gap-gap=0:**
```
Seq1: A  C  -  T
Seq2: A  C  G  T
Seq3: A  T  G  T

col 1 (A,A,A): 1+1+1 = 3
col 2 (C,C,T): 1+0+0 = 1
col 3 (-,G,G): −1+−1+1 = −1   ← gap con residuo penaliza; dos residuos iguales suman
col 4 (T,T,T): 1+1+1 = 3
SP = 6
```

### Relación SP score ↔ árbol guía

Con 3 secuencias, SP_total = S(1,2) + S(1,3) + S(2,3), independientemente del orden en que se construyó el MSA.

Con 4 secuencias agrupadas como (1&2) vs (3&4):
```
SP_total = S(1,2) + S(3,4) + S(1,3) + S(1,4) + S(2,3) + S(2,4)
```

---

## 5. Bases de Datos Primarias

### Introducción: ¿por qué existen bases de datos biológicas?

Cada vez que se secuencia un gen, se resuelve una estructura proteica o se hace un experimento de expresión génica, se genera información que tiene valor más allá del laboratorio que la produjo. Las **bases de datos biológicas** centralizan esos datos para que cualquier investigador pueda aplicar el principio de "culpa por asociación": si mi secuencia nueva se parece a una secuencia con función conocida, puedo inferir su función sin necesidad de experimentar desde cero.

No todas las bases de datos son iguales. Difieren en tres dimensiones clave que determinan cuándo conviene usarlas: **curación** (¿quién controla la calidad?), **redundancia** (¿aparece la misma secuencia varias veces?) y **origen de los datos** (¿datos crudos experimentales o resultado de análisis bioinformático?). Entender esas dimensiones es exactamente lo que se pregunta en el Problema 2 del parcial.

### Criterios de clasificación

| Dimensión | Opciones | Definición |
|-----------|----------|-----------|
| **Curación** | Curada vs No curada | Curada = revisión manual por expertos. No curada = el autor deposita sin revisión. |
| **Redundancia** | Redundante vs No redundante | Redundante = el mismo dato puede aparecer en múltiples registros. |
| **Origen** | Primaria vs Secundaria | Primaria = datos experimentales crudos. Secundaria = resultado del análisis bioinformático de datos primarios. |

> **Regla:** para clasificar una BD, responder las 3 preguntas. Se puede ser curada y primaria (RefSeq), o no curada y primaria (GenBank), o curada y secundaria (SGD).

### Bases de nucleótidos

| BD | Curación | Redundancia | Tipo | Notas |
|----|----------|-------------|------|-------|
| **GenBank** | No curada | Redundante | Primaria | Repositorio NCBI; el autor deposita y es responsable; INSDC |
| **EMBL/ENA** | No curada | Redundante | Primaria | Equivalente europeo (EBI); mismo formato INSDC |
| **DDBJ** | No curada | Redundante | Primaria | Equivalente japonés; sincronizado con INSDC |
| **RefSeq** | Curada | No redundante | Primaria | 1 registro por molécula por organismo; NCBI lo mantiene |
| **Ensembl** | Curada | No redundante | Secundaria | Browser genómico; anotación de vertebrados |

> GenBank es análogo = **paper**; RefSeq es análogo = **review**.

**Prefijos RefSeq:**
| Prefijo | Contenido |
|---------|-----------|
| NM | mRNA verificado experimentalmente |
| NP | Proteína verificada experimentalmente |
| XM | mRNA predicho in silico |
| XP | Proteína predicha in silico |
| NT | DNA genómico de contigs |

### Bases de proteínas

| BD | Curación | Redundancia | Tipo | Notas |
|----|----------|-------------|------|-------|
| **Swiss-Prot** | Curada | No redundante | Primaria | 1 registro = 1 gen; todas las isoformas juntas |
| **TrEMBL** | No curada | Redundante | Primaria | Traducción automática de GenBank/EMBL |
| **UniProtKB** | Mixta | Baja | Primaria | Swiss-Prot + TrEMBL integrados |
| **GenPept** | No curada | Redundante | Primaria | Traducción automática del INSDC (NCBI) |
| **PDB** | Curada | Redundante | Primaria | Estructuras 3D; misma proteína puede tener cientos de estructuras |

> Swiss-Prot ≠ TrEMBL: Swiss-Prot curada manual → confiable. TrEMBL automática → enorme pero con más errores.
>
> UniProt organiza por **unidad funcional** (gen) → todas las isoformas en 1 registro. RefSeq organiza por **molécula** → cada isoforma de mRNA = 1 registro separado.

### Bases de expresión génica

| BD | Curación | Redundancia | Tipo | Notas |
|----|----------|-------------|------|-------|
| **GEO** | No curada | Redundante | Primaria | Repositorio NCBI de microarrays y RNA-seq |
| **ArrayExpress** | No curada | Redundante | Primaria | Equivalente europeo (EBI) |
| **GTEx** | Curada | No redundante | Secundaria | Expresión por tejido en humanos sanos |

### Bases de variantes genéticas

| BD | Curación | Redundancia | Notas |
|----|----------|-------------|-------|
| **dbSNP** | No curada | Redundante | Repositorio de todos los SNPs; IDs `rs` |
| **ClinVar** | Curada | No redundante | Variantes con significado clínico (Patogénica/VUS/Benigna) |
| **OMIM** | Curada | No redundante | Genes y enfermedades mendelianas |
| **gnomAD** | Curada | No redundante | Frecuencias alélicas en ~125.000 genomas de referencia |
| **COSMIC** | Curada | No redundante | Mutaciones somáticas en cáncer |

### Otras bases importantes

| BD | Tipo | Curación | Notas |
|----|------|----------|-------|
| **SGD** | Secundaria | Curada | S. cerevisiae; integra GenBank + UniProt + literatura |
| **FlyBase** | Secundaria | Curada | D. melanogaster |
| **KEGG** | Secundaria | Curada | Pathways metabólicos y de señalización |
| **PubMed** | Secundaria | Curada | Literatura biomédica; indexada por MeSH |
| **PubChem** | Primaria | No curada | Estructuras químicas; redundante |
| **DrugBank** | Secundaria | Curada | Fármacos con targets, mecanismo, farmacocinética |
| **IntAct** | Secundaria | Curada | Interacciones proteína-proteína |

---

## 6. Bases de Datos Secundarias — HMM, PROSITE, InterPro

### El problema: una BD primaria no es un clasificador

Las BDs primarias (§5) guardan registros individuales: una secuencia, un experimento, una estructura. Cuando yo tengo una proteína nueva y quiero saber a qué familia pertenece, puedo hacer BLAST contra Swiss-Prot y ver qué proteínas conocidas se parecen. Pero hay un límite: dos proteínas de la misma familia pueden tener menos del 30% de identidad de secuencia (zona twilight) y BLAST las pierde.

La diferencia conceptual fundamental es esta:

- **BD primaria**: el registro = una secuencia (resultado de un experimento). La pregunta que responde: "¿hay una secuencia parecida a la mía?"
- **BD secundaria**: el registro = un **modelo matemático** que describe una familia entera. La pregunta que responde: "¿pertenece mi secuencia a esta familia?"

En una BD secundaria el dato no viene de un experimento, viene del **análisis de datos de otras BDs**. El registro es una abstracción: captura la esencia matemática de qué tienen en común todos los miembros de la familia, incluyendo los que tienen solo 20% de identidad de secuencia entre sí.

### ¿Qué es un dominio? ¿Qué es una familia?

**Familia**: grupo de proteínas evolutivamente relacionadas (mismo ancestro común) que comparten función. El término viene del árbol filogenético.

**Dominio**: segmento de secuencia proteica con tres propiedades simultáneas:
1. **Función propia** — el dominio tiene actividad biológica característica por sí solo. Ej: dominio kinasa, dominio de unión ATP, dominio de unión al ADN.
2. **Estructura independiente** — el dominio se pliega correctamente incluso si se lo separa del resto de la proteína. Esto permite estudiarlo en forma aislada.
3. **Evolución independiente** — los dominios se "mezclan y combinan" evolutivamente (domain shuffling). Las proteínas de organismos superiores son más complejas precisamente porque tienen más dominios.

Ejemplo: un anticuerpo tiene dominios variables (reconocen el antígeno) y dominio Fc (interacciona con células inmunes). Cambiando el dominio Fc se obtienen las diferentes clases de inmunoglobulinas (IgG, IgM, etc.) que tienen la misma especificidad pero diferente función efectora. Los exones suelen corresponderse con dominios — el splicing alternativo es en parte un mecanismo de ensamblar combinaciones de dominios.

### El desafío: construir el modelo matemático

Para que una BD secundaria sea útil, el modelo matemático de cada registro tiene que ser capaz de **clasificar**: dada una secuencia nueva, decir si pertenece o no. Para hacer eso hay que resolver tres problemas:

| Problema | Pregunta | Descripción |
|----------|----------|-------------|
| **Entrenamiento** | ¿Cómo estimo los parámetros del modelo? | A partir de las secuencias conocidas de la familia, determinar los parámetros que definen el modelo de ese registro |
| **Puntaje (scoring)** | ¿Qué tan bien encaja una secuencia nueva en este modelo? | Dado el modelo ya entrenado, asignar un número a la secuencia nueva |
| **Alineamiento** | ¿Cómo se posiciona la secuencia nueva dentro del modelo? | Encontrar cuál es la correspondencia óptima entre posiciones de la secuencia y posiciones del modelo |

Estos tres problemas están **íntimamente ligados**: para calcular el puntaje hay que hacer el alineamiento; para hacer el alineamiento hay que tener el modelo entrenado.

**La lógica del umbral:**
Una vez construido el modelo, se evalúan dos grupos de secuencias:
- Secuencias que **sí** deben pertenecer a la familia → distribución de scores "positivos"
- Secuencias que **no** pertenecen → distribución de scores "negativos"

Lo ideal: dos histogramas bien separados, sin solapamiento. En la práctica hay solapamiento → se elige un umbral de corte. Las secuencias que caen del lado equivocado son los falsos positivos (FP) y falsos negativos (FN).

---

### Solución 1: PROSITE — expresiones regulares (alta especificidad, baja sensibilidad)

La forma más simple de modelar un motivo es con un lenguaje formal que describa exactamente qué aminoácidos se toleran en cada posición. PROSITE usa una sintaxis inspirada en expresiones regulares:

| Constructor | Significado | Ejemplo |
|-------------|-------------|---------|
| Letra | Aminoácido fijo obligatorio | `N` = solo Asn |
| `x` | Cualquier aminoácido | `N-x-T` |
| `[AB]` | A **o** B en esa posición | `[ST]` = Ser o Thr |
| `{AB}` | Cualquiera **excepto** A o B | `{P}` = cualquiera menos Pro |
| `A(n)` | Repetición exacta n veces | `x(3)` = tres cualquieras |
| `x(n,m)` | Entre n y m repeticiones | `x(2,4)` |
| `x(0,n)` | Elemento opcional (0 a n veces) | `x(0,1)` = puede o no estar |

> Los guiones `-` son separadores opcionales para legibilidad.

**Ejemplos:**

| Patrón | Cómo leerlo |
|--------|-------------|
| `[AC]-x-V-x(4)-{ED}` | (A o C) — cualquiera — V — 4×cualquiera — (no E ni D) |
| `A-x(0,1)-{V}` | A — (un aa opcional) — (no V) |
| `[ST]-C-x(1,2)-M` | (S o T) — C — (1 o 2 aa) — M |

**Motivo de N-glicosilación:** `N-{P}-[ST]`
- Asn − (cualquier aa **excepto** Pro) − (Ser **o** Thr)
- **¿Por qué no Prolina?** La Pro tiene su nitrógeno backbone atrapado en un anillo pirrolidínico → no puede rotar → fuerza un quiebre (kink) rígido en la cadena peptídica → el tripéptido N-P-S/T no puede adoptar la conformación de β-turn que necesita la enzima OST para realizar la glicosilación. Por eso Pro está **excluida** de esa posición.

**Limitación de PROSITE:** solo captura variabilidad posición a posición mediante un lenguaje discreto. No aprende qué tan frecuente es cada aminoácido en cada posición, ni modela inserciones/deleciones con probabilidades. Resultado: alta especificidad (pocas falsas alarmas) pero baja sensibilidad (se pierde muchos miembros divergentes de la familia).

---

### Solución 2: PSSM (Position-Specific Scoring Matrix)

Una PSSM es una tabla de 20 columnas (una por aminoácido) × L filas (una por posición del motivo). Cada celda contiene `score(posición i, aa a) = log[ frecuencia observada / frecuencia de fondo ]`.

Se construye a partir del MSA de las secuencias conocidas de la familia:
- Posición muy conservada (ej: siempre Asp) → score muy alto para Asp, muy negativo para el resto
- Posición variable → scores similares para muchos aminoácidos

**Ventaja sobre PROSITE:** cuantifica la variabilidad en lugar de usar "sí/no". Es lo que usa PSI-BLAST para refinar su búsqueda iterativamente.

**Limitación de la PSSM:** trata todas las posiciones independientemente. No modela explícitamente inserciones ni deleciones — solo funciona bien si las secuencias que comparo tienen exactamente el mismo largo que el motivo.

---

### Solución 3: HMM de perfil — el modelo más potente

En 1996 se demostró que un perfil de secuencias (como la PSSM) puede describirse como un **Modelo Oculto de Markov (Hidden Markov Model)**. Esto importa porque los HMMs tienen toda una teoría matemática detrás que permite resolver los 3 problemas de forma rigurosa y eficiente.

#### ¿Qué es un HMM?

Un HMM es una caja negra que **genera secuencias**. Por dentro tiene:
- **Estados** (ocultos — no los observamos directamente)
- **Probabilidades de transición** entre estados: P(ir del estado A al estado B)
- **Probabilidades de emisión**: para cada estado, la probabilidad de que "escupa" cada símbolo (aminoácido o nucleótido)

La clave: los estados están **ocultos**. Solo vemos la secuencia de símbolos emitidos, nunca la secuencia de estados que la generó. Eso explica el nombre "oculto".

**Ejemplo simplificado** (para intuir el concepto):
Supongamos un HMM con 2 estados: Estado 1 emite A y T con frecuencia (zona AT-rica), Estado 2 emite G y C con frecuencia (zona GC-rica). Si el modelo está el 98% del tiempo en Estado 1 y 2% en Estado 2, las secuencias generadas tendrán mayoría de A y T. Si veo una región de la secuencia con muchas G y C, puedo inferir que el modelo estaba probablemente en Estado 2 ahí — aunque no puedo saberlo con certeza. Los estados permanecen ocultos.

#### Los 3 tipos de estados en un HMM de perfil de proteínas

Para modelar una familia de proteínas, el HMM usa tres tipos de estados:

```
          ┌─────────────────── Match ──────────────────┐
          │  Estado M(i): corresponde a la columna i   │
          │  del MSA. Emite los 20 aa con probabilidad │
          │  = frecuencia observada en esa columna.    │
          └────────────────────────────────────────────┘

          ┌─────────────────── Insert ─────────────────┐
          │  Estado I(i): puede aparecer entre dos      │
          │  estados Match. Modela una inserción        │
          │  (aminoácido extra que no corresponde a     │
          │  ninguna columna del MSA). Emite con        │
          │  probabilidad de fondo (uniforme).          │
          └────────────────────────────────────────────┘

          ┌─────────────────── Delete ─────────────────┐
          │  Estado D(i): estado silente — no emite     │
          │  ningún símbolo. Modela una deleción        │
          │  (un gap en esa posición). Se "salta"       │
          │  el estado Match correspondiente.           │
          └────────────────────────────────────────────┘
```

La arquitectura del modelo para una familia con L columnas en el MSA es:
```
Inicio → M(1) → M(2) → M(3) → … → M(L) → Fin
           ↕      ↕      ↕
          I(1)  I(2)  I(3)      ← estados de inserción
           ↕      ↕      ↕
          D(1)  D(2)  D(3)      ← estados de deleción (silentes)
```

Las probabilidades de transición (P(M→M), P(M→I), P(M→D)) se estiman a partir de la frecuencia de gaps observados en el MSA de entrenamiento.

#### Problema 1: Entrenamiento — construir el modelo desde el MSA

```
Paso 1: Tomar el MSA de las secuencias conocidas de la familia
Paso 2: Para cada columna i del MSA:
        → Contar la frecuencia de cada aminoácido → probabilidades de emisión de M(i)
        → Contar la frecuencia de gaps → probabilidades de ir a D(i) o I(i)
Paso 3: Suavizar las probabilidades con frecuencias de fondo
        (porque nunca tengo infinitas secuencias para cada columna)
Resultado: el perfil HMM = archivo .hmm con todos los parámetros
Herramienta: hmmbuild
```

**¿Por qué suavizar con frecuencias de fondo?** Si en 10 secuencias de entrenamiento la posición 5 siempre tiene Phe, la probabilidad estimada de Ala en esa posición sería 0. Pero probablemente existan miembros de la familia donde esa posición tiene Tyr (similar a Phe). El suavizado mezcla la observación con la frecuencia general de cada aminoácido, evitando asignar probabilidad 0 a algo solo porque no apareció en el set de entrenamiento.

#### Problema 2: Scoring — ¿qué puntaje le doy a una secuencia nueva?

El puntaje de una secuencia S contra el HMM es la **probabilidad de que el HMM haya generado esa secuencia S**. Matemáticamente: P(S | modelo).

Intuitivamente: si corriera el HMM millones de veces y contara cuántas veces genera exactamente la secuencia S, esa fracción es el puntaje. En la práctica se calcula analíticamente usando programación dinámica.

Después del scoring se compara contra el umbral:
- Score > umbral → la secuencia probablemente pertenece a la familia
- Score < umbral → probablemente no pertenece

#### Problema 3: Alineamiento — ¿cómo se "encaja" la secuencia en el modelo?

Alinear la secuencia contra el HMM significa determinar **qué estado corresponde a cada posición de la secuencia**: ¿este aminoácido corresponde al estado M(3), o al estado I(2), o al D(4)? Esa asignación es el alineamiento.

La solución óptima es la asignación de estados que maximiza la probabilidad P(S, secuencia_de_estados | modelo). Se calcula con programación dinámica igual que NW/SW, pero en lugar de llenar una matriz de scores se llena una matriz de probabilidades.

**Conexión con inserciones y deleciones:**
- Si la secuencia tiene un aminoácido "de más" en la posición 5 → ese aminoácido se alinea contra el estado I(4) o I(5) → es una inserción
- Si la secuencia tiene un gap en la posición 3 → el estado M(3) se mapea a un D(3) → es una deleción

Esta es la gran ventaja del HMM sobre la PSSM: **modela inserciones y deleciones de forma explícita y probabilística**, por eso es mucho más potente para detectar homólogos divergentes.

---

### Pipeline HMMer completo

```
1. Secuencias conocidas de la familia (experimentalmente validadas)
        ↓
2. MSA de esas secuencias (con CLUSTALW, MAFFT, etc.)
        ↓
3. hmmbuild   →   construye el perfil HMM (.hmm)
   = PROBLEMA DE ENTRENAMIENTO
        ↓
4. Calibración del umbral:
   - Evaluar las secuencias de entrenamiento con el modelo → distribución de scores "positivos"
   - Evaluar secuencias al azar → distribución de scores "negativos"
   - Elegir el umbral que separa las dos distribuciones (minimiza FP + FN)
        ↓
5. hmmsearch  →   buscar en el proteoma o en una BD con ese umbral
   = PROBLEMA DE PUNTAJE (+ ALINEAMIENTO para cada hit)
        ↓
6. Candidatos con score > umbral → validación experimental
        ↓
7. Evaluar calidad: sensibilidad / especificidad / curva ROC
```

---

### Sequence Logo — leer el modelo visualmente

El logo es la visualización del HMM. Cada columna del logo corresponde a un estado Match del HMM.

- **Altura total de la columna** = información (bits) = qué tan conservada está esa posición
  - Máximo teórico ≈ 4.32 bits (= log₂20 — posición completamente conservada en 1 aminoácido)
  - Mínimo = 0 bits — posición completamente variable (todos los aa igualmente probables)
- **Tamaño relativo de cada letra** = probabilidad de ese aminoácido en esa posición
- **Letra más grande** = aminoácido más frecuente (más probable que el HMM emita en ese estado Match)

**Qué se puede leer de un logo:**
1. **Posiciones críticas** — columnas altas: ese residuo es casi imposible de cambiar sin perder función
2. **Identidad fisicoquímica del residuo dominante** — ¿es cargado, hidrófobo, polar? Dice qué tipo de interacción necesita la proteína ahí
3. **Tolerancia a sustituciones conservativas** — si la columna tiene altura media y 2-3 letras similares (ej: L, I, V), acepta variación dentro del grupo hidrófobo
4. **Posiciones variables** — columnas bajas: no involucradas en reconocimiento específico
5. **Largo del motivo** — número total de columnas

> **¡OJO en el parcial!** Verificar si el logo es de proteínas (letras de aminoácidos: A, R, N, D...) o de ADN (A, T, G, C). Si tiene T y U con frecuencias similares es ARN. Si tiene principalmente H, I, L, V, F... es proteína.

---

### Métricas de evaluación del modelo

Una vez que se busca en el proteoma, los resultados se clasifican en:

```
                      Predicción del modelo
                    ┌──────────┬──────────┐
                    │ POSITIVO │ NEGATIVO │
       ┌────────────┼──────────┼──────────┤
Verdad │ POSITIVO   │    VP    │    FN    │
real   │ NEGATIVO   │    FP    │    VN    │
       └────────────┴──────────┴──────────┘

Sensibilidad  = VP / (VP + FN)   ← fracción de miembros reales detectados
Especificidad = VN / (VN + FP)   ← fracción de no-miembros correctamente rechazados
```

**Falso positivo:** secuencia que supera el umbral pero no es miembro real de la familia. Causas posibles:
- El motivo existe en la secuencia pero está enterrado/inaccesible estructuralmente
- El umbral está demasiado bajo (muy permisivo)
- Set de entrenamiento pequeño → el modelo aprendió ruido estadístico

**Falso negativo:** miembro real de la familia que no supera el umbral. Causas posibles:
- Secuencia muy divergente (zona twilight)
- El umbral está demasiado alto (muy estricto)
- La región de la familia está poco representada en el set de entrenamiento

**Curva ROC:** se barre el umbral de score y en cada valor se calcula sensibilidad y especificidad → eje Y = sensibilidad, eje X = 1 − especificidad. Un modelo perfecto pasa por la esquina superior izquierda (sensibilidad=1, especificidad=1). Diagonal = clasificación aleatoria.

**Tensión entre sensibilidad y especificidad:** bajar el umbral aumenta sensibilidad (se detectan más miembros) pero baja especificidad (entran más FP). Subir el umbral hace lo contrario. El umbral óptimo depende del objetivo biológico.

---

### ¿Por qué HMM > BLAST para familias?

| Criterio | BLAST | HMM de perfil |
|----------|-------|---------------|
| Modelo | Secuencia individual | Toda la familia a la vez |
| Inserciones/deleciones | Solo penaliza con gap penalty uniforme | Modela con probabilidades posición-específicas |
| Posiciones conservadas | No las distingue de las variables | Les asigna mayor peso automáticamente |
| Detección de homólogos remotos | Pierde en zona twilight (<30% id) | Detecta hasta 10-20% de identidad |
| Construcción del modelo | No necesita entrenamiento | Requiere MSA de entrenamiento curado |

---

### Bases de datos secundarias principales

| BD | Método | Contenido | Notas |
|----|--------|-----------|-------|
| **Pfam** | HMM (hmmbuild/hmmsearch) | Familias y dominios proteicos | Seeds curados manualmente; la más usada |
| **PROSITE** | Regex + PSSM | Patrones cortos y perfiles | Alta especificidad; motivos funcionales |
| **SMART** | HMM | Dominios de señalización | Incluye contexto arquitectónico (qué dominios aparecen juntos) |
| **TIGRFAMs** | HMM | Familias de función específica | Útil en genómica microbiana |
| **CDD** (NCBI) | RPS-BLAST | Dominios conservados | Integrado en BLAST de NCBI |

### InterPro — la BD de BDs secundarias

InterPro (EBI) consolida todas las BDs de dominios/familias en una sola entrada unificada. Un registro InterPro agrupa entradas equivalentes de Pfam, PROSITE, SMART, etc. que describen el mismo dominio o familia.

- **Herramienta:** **InterProScan** — busca una secuencia contra todas las bases miembro en una sola corrida
- **Output:** lista de dominios encontrados con posición, E-value, BD de origen y términos GO (Gene Ontology)

**Tipos de entradas InterPro:**

| Tipo | Descripción |
|------|-------------|
| **Family** | Proteínas evolutivamente relacionadas que comparten función |
| **Domain** | Unidad estructural/funcional que aparece en distintos contextos proteicos |
| **Homologous superfamily** | Comparten estructura 3D sin similitud de secuencia detectable |
| **Repeat** | Secuencia corta repetida en tandem dentro de la proteína |
| **Site** | Sitio activo, de unión a ligando, o de modificación post-traduccional |

---

## 7. Matrices PAM vs BLOSUM

### Introducción: ¿por qué no usar match=1 / mismatch=−1?

En los algoritmos NW/SW necesitamos asignar un score a cada par de aminoácidos alineados. La opción más simple sería match=1 / mismatch=−1, pero eso ignora la biología: no es igual que una Leucina mute a Isoleucina (ambas hidrófobas, muy tolerado evolutivamente) que mute a Arginina (cargada positiva, casi nunca ocurre). Las **matrices de sustitución** codifican esta información: para cada par de aminoácidos dicen cuánto "vale" ese alineamiento basándose en qué tan frecuentemente ocurre esa sustitución en la evolución real.

La fórmula general es: `score(a,b) = log[ P(a y b alineados) / P(a y b por azar) ]`. Si el par ocurre más seguido de lo esperado al azar → score positivo. Si ocurre menos → score negativo.

Hay dos familias de matrices construidas con distintas filosofías:

| | PAM | BLOSUM |
|--|-----|--------|
| **Construida con** | Alineamientos globales de secuencias ≥85% id | Bloques conservados locales sin gaps |
| **Número** | ↑ = mayor distancia evolutiva | ↑ = secuencias más similares usadas |
| **Default** | — | **BLOSUM62** |
| Para secuencias similares | PAM120 | BLOSUM80 |
| Para secuencias divergentes | PAM250 | BLOSUM45 |

> **Regla de oro:** BLOSUM62 para la mayoría de búsquedas. Para comparaciones muy divergentes → BLOSUM45 o PAM250.
>
> **Intuición BLOSUM:** el número = umbral de identidad de los bloques. BLOSUM80 fue construida con bloques de secuencias ≥80% idénticas → puntúa alto las sustituciones muy conservadas, penaliza fuerte cualquier cambio. BLOSUM45 fue construida con secuencias más divergentes → más permisiva.

### ¿Qué hace positivo un score BLOSUM62?

Si dos aminoácidos comparten el **mismo grupo fisicoquímico** (ver diagrama de Venn en sección 8), su score en BLOSUM62 es positivo. Cuantos más grupos comparten → score más alto.

| Par | Score BLOSUM62 | Razón |
|-----|----------------|-------|
| L ↔ I | +2 | Aliphatic + Hydrophobic |
| D ↔ E | +2 | Ambos ácidos cargados negativamente |
| K ↔ R | +2 | Ambos básicos cargados positivamente |
| F ↔ Y | +3 | Ambos aromáticos |
| V ↔ I | +3 | Aliphatic + Hydrophobic |
| S ↔ T | +1 | Ambos polares pequeños con OH |
| P ↔ L | −3 | Pro rompe estructura; grupos distintos |
| P ↔ W | −4 | Pro rompe estructura; grupos muy distintos |
| P ↔ P | +7 | Auto-match |

---

## 8. Aminoácidos

### Introducción: por qué importan las propiedades fisicoquímicas

En bioinformática, los aminoácidos no son solo letras del alfabeto. Sus propiedades fisicoquímicas determinan qué sustituciones son toleradas evolutivamente (y por lo tanto qué scores tienen en BLOSUM), qué residuos aparecen en sitios activos, y qué motivos de secuencia son funcionalmente relevantes. Esta sección es referencia rápida para responder preguntas como "¿por qué Pro no puede estar en el motivo N-{P}-[ST]?" o "¿por qué Cys aparece en puentes disulfuro?".

### Los 20 aminoácidos: códigos y propiedades

| 1 letra | 3 letras | Nombre | Grupo | Nota especial |
|:-------:|:--------:|--------|-------|---------------|
| **A** | Ala | Alanina | Hidrofóbico pequeño | El más frecuente en proteínas |
| **R** | Arg | Arginina | Básico (+) | pKa ~12.5; carga + más estable a pH 7 |
| **N** | Asn | Asparagina | Polar sin carga | Sitio de N-glicosilación (motivo N-X-S/T) |
| **D** | Asp | Ácido aspártico | Ácido (−) | pKa ~3.9; cadena corta |
| **C** | Cys | Cisteína | Polar / hidrofóbico | Puentes disulfuro; sitios activos nucleofílicos; 2do más raro |
| **Q** | Gln | Glutamina | Polar sin carga | N↔Q y D↔E son sustituciones conservativas |
| **E** | Glu | Ácido glutámico | Ácido (−) | pKa ~4.1; un CH₂ más que Asp |
| **G** | Gly | Glicina | Tiny / flexible | Sin cadena lateral; máxima flexibilidad; aparece en loops y giros |
| **H** | His | Histidina | Básico / polar | **Único con pKa ~6** → catálisis ácido-base a pH fisiológico (tríada catalítica) |
| **I** | Ile | Isoleucina | Aliphatic | Isómero de Leu; L↔I muy conservativa |
| **L** | Leu | Leucina | Aliphatic | El aminoácido más frecuente |
| **K** | Lys | Lisina | Básico (+) | pKa ~10.5; ubiquitinación; interacciona con DNA |
| **M** | Met | Metionina | Hidrofóbico | Codón de inicio (AUG) |
| **F** | Phe | Fenilalanina | Aromático | Núcleo hidrofóbico; absorbe UV a 280nm levemente |
| **P** | Pro | Prolina | Especial (rígido) | **Rompe hélices y hojas β. Anillo unido al N backbone → no dona H. No va en X de N-glicosilación.** |
| **S** | Ser | Serina | Polar pequeño | Fosforilación; O-glicosilación |
| **T** | Thr | Treonina | Polar | Fosforilación; posición S/T en motivo N-glicosilación |
| **W** | Trp | Triptófano | Aromático | **El más raro (~1%)**; casi siempre crítico funcionalmente; absorbe fuerte a 280nm |
| **Y** | Tyr | Tirosina | Aromático / polar | Fosforilación; absorbe a 280nm |
| **V** | Val | Valina | Aliphatic | Voluminoso; volúmenes: V < L = I |

### Clasificación por grupos fisicoquímicos

```
TINY:      G  A  C  S
SMALL:     G  A  C  S  N  D  T  P
ALIPHATIC: I  V  L
AROMATIC:  F  Y  W
HYDROPHOBIC (teal): I  V  L  M  A  G  C  T  F  Y  W
POLAR (blue):       S  T  C  N  Q  D  E  K  R  H
CHARGED (yellow):   D  E  K  R  H
POSITIVE (red):     K  R  H
NEGATIVE / ÁCIDOS:  D  E
```

> **Regla:** si dos aa comparten grupo → sustitución conservativa → BLOSUM62 positivo.
> L e I comparten ALIPHATIC + HYDROPHOBIC → score alto.
> K y E tienen cargas opuestas → score negativo.

### Aminoácidos especiales (importantes para el parcial)

| aa | ¿Por qué es especial? |
|----|-----------------------|
| **G** | Sin cadena → máxima flexibilidad; aparece en turns y loops |
| **P** | Anillo con N backbone → rompe estructuras secundarias; fuerza kink; **no va en posición X de N-glicosilación** |
| **C** | Tiol → puentes disulfuro; sitios activos; el 2do más raro |
| **W** | El más raro (~1%); casi siempre funcionalmente crítico; mayor absorbancia a 280nm |
| **H** | Único pKa ~6 → catálisis ácido-base a pH 7 (tríada catalítica: Ser-His-Asp) |

### Modificaciones post-traduccionales

| Modificación | aa | Función |
|-------------|-----|---------|
| **Fosforilación** | S, T, Y | Señalización; activación/inactivación enzimática |
| **N-glicosilación** | N (motivo N-{P}-[ST]) | Plegamiento, estabilidad, reconocimiento celular |
| **Ubiquitinación** | K | Degradación por proteasoma |
| **Puente disulfuro** | C-C | Estabilidad en proteínas extracelulares |
| **Acetilación** | K (histonas), Met N-terminal | Regulación génica |

---

## 9. NGS, Ensamblado, Genómica Humana

### Introducción: el problema de leer el ADN

Desde que Watson y Crick describieron la doble hélice, los biólogos entendieron que la información biológica está en la secuencia de bases del ADN. El desafío fue desarrollar métodos para leer esa secuencia. La historia de la secuenciación es la historia de dos problemas que se fueron resolviendo progresivamente:

1. **¿Cómo identificar cada base?** → necesito un mecanismo que me diga si es A, T, C o G
2. **¿Cómo escalar?** → un genoma humano tiene 3.000 millones de bases; no puedo leerlas de a una

Las soluciones a estos dos problemas evolucionaron en generaciones tecnológicas. Cada generación resolvió las limitaciones de la anterior y abrió nuevas preguntas.

---

### Solución 1: Sanger (1ª generación) — leer un fragmento a la vez

Sanger (Premio Nobel ×2) inventó el primer método práctico basado en **didesoxinucleótidos terminadores (ddNTPs)**. La idea: si la polimerasa copia una hebra molde pero de vez en cuando incorpora un nucleótido "roto" que le impide seguir, va a generar fragmentos de todos los largos posibles (1, 2, 3… N bases). Si además marco **en qué base termina cada fragmento**, puedo ordenarlos por tamaño y leer la secuencia.

**¿Qué es un ddNTP?** Un nucleótido normal tiene un -OH en el carbono 3' que permite que la polimerasa enganche el siguiente nucleótido. Un ddNTP (didesoxinucleótido) no tiene ese -OH → la polimerasa lo incorpora, pero ahí se detiene: no puede agregar nada más. Es un terminador.

**Paso a paso:**

```
Hebra molde (lo que quiero secuenciar): 5'- G A T C A G -3'

Hago 4 tubos en paralelo. En cada uno pongo:
  • polimerasa + todos los dNTPs normales (para copiar)
  • 10% de un ddNTP distinto (para terminar al azar en esa base)

  Tubo G: 10% ddGTP  → la polimerasa copia hasta que incorpora una ddG → para
  Tubo A: 10% ddATP  → para cuando incorpora una ddA
  Tubo T: 10% ddTTP  → para cuando incorpora una ddT
  Tubo C: 10% ddCTP  → para cuando incorpora una ddC

Al final de cada tubo tengo una mezcla de fragmentos de todos los largos
posibles que terminan en esa base. Por ejemplo, el Tubo G contiene:
  • el fragmento de largo 1 (termina en la primera G de la copia)
  • el fragmento de largo 4 (termina en la segunda G)
  • etc.
```

**El gel — leer la secuencia:**

Un gel de electroforesis es como una carrera de fragmentos: los más cortos migran más rápido y llegan más abajo. El gel tiene **4 calles** (carriles verticales), una por tubo:

```
               Calle G   Calle A   Calle T   Calle C
  ARRIBA        |         |         |         |
  (fragmentos   |         |         |    ●    |   ← largo 5 → termina en T → 5ª letra = T
  largos)       |         |    ●    |         |   ← largo 4 → termina en A → 4ª letra = A
                |    ●    |         |         |   ← largo 3 → termina en G → 3ª letra = G
                |         |    ●    |         |   ← largo 2 → termina en A → 2ª letra = A
  ABAJO         |    ●    |         |         |   ← largo 1 → termina en G → 1ª letra = G
  (fragmentos
  cortos)

  Leyendo de abajo hacia arriba la calle donde aparece cada banda:
  G → A → G → A → T → ...   = la secuencia de la copia (complemento de la molde)
```

**Versión moderna — Sanger fluorescente + capilar:**

El gel con 4 calles es incómodo. La solución fue darle a cada ddNTP un **color fluorescente distinto** → ahora todo cabe en 1 sola reacción (no 4 tubos). Los fragmentos se separan por tamaño en un capilar (un tubo finísimo) y un detector lee el color de cada fragmento al pasar:

```
Fragmentos pasando por el capilar (del más corto al más largo):

  tiempo →  [●rojo] [●verde] [●rojo] [●azul] [●amarillo] ...
               G        A       G       A        T        ...
              
  El detector registra: rojo=G, verde=A, rojo=G, azul=A, amarillo=T → GAGAT...
  Eso es el electroferograma: un gráfico de picos de color = la secuencia
```

**Limitación que motiva la siguiente generación:** todo esto secuencia **un solo fragmento** por reacción. Para un genoma humano (3.000 millones de bases), necesitaría millones de reacciones → décadas de trabajo. Solución: paralelizar masivamente.

---

### Herramienta clave: PCR (Mullis, Premio Nobel)

La PCR permite:
1. **Amplificación exponencial** de un fragmento específico a partir de una sola molécula
2. **Selectividad**: en una mezcla compleja (ej: genoma entero), amplifica solo el fragmento flanqueado por los primers
3. Amplificación **clonal**: todas las copias resultantes son idénticas

La PCR es la base del paso de amplificación clonal en las tecnologías de 2ª generación. Sin PCR (o algún equivalente), no hay señal detectable.

---

### Preparación de librería — paso previo a todos los NGS

Antes de entrar al secuenciador, el ADN debe acondicionarse:

```
Muestra de ADN
    ↓ fragmentación (mecánica o enzimática → fragmentos de ~150–1000 pb)
    ↓ agregar adaptadores a los extremos (secuencias cortas conocidas, iguales para todos los fragmentos)
    ↓ [opcional] enriquecer en regiones de interés:
        • Captura: sondas magnéticas complementarias a los exones (exoma) → imán retiene solo esos fragmentos
        • Amplificación selectiva por PCR: primers para la región deseada (ej: genes virales en muestra nasal)
    ↓ → librería lista para secuenciar
```

Los adaptadores son esenciales: el secuenciador los usa como sitio de anclaje y como primer para empezar a leer.

---

### Pair-End Reads — concepto fundamental

Un fragmento de ADN doble cadena tiene **dos extremos** (ends). Los secuenciadores de 2ª generación leen ambos → generan **dos lecturas apareadas**:

```
5'───[Lectura 1]──────────────────[Lectura 2]───3'
      ←~100 pb→   (interior no leído)  ←~100 pb→
      │←──────────── inserto ~1–5 kb ───────────→│
```

- Las dos lecturas vienen del **mismo fragmento** → tamaño del inserto conocido de antemano
- Si sé dónde mapea una lectura, sé a qué distancia debe estar la otra
- **Fundamental para el scaffolding** (Paso 2 del ensamblado): permiten ordenar contigs y estimar el tamaño del gap entre ellos

---

### Profundidad (Depth) vs Cobertura (Coverage)

Son dos métricas distintas que se confunden fácilmente:

**Profundidad (depth):** número de veces que una base concreta aparece en las lecturas.
- Es un valor por posición, varía a lo largo del genoma
- Se resume como **profundidad media** = promedio sobre todas las posiciones

**Cobertura (coverage):** porcentaje del genoma cubierto por al menos X lecturas (X = umbral elegido antes de calcularlo).
- Ejemplo: "cobertura del 90% a 10X" → el 90% de las bases se leyó al menos 10 veces

```
Profundidad target (teórica) = (N_reads × L_read) / L_genoma  ← se fija ANTES del experimento

Profundidad media real  ← se calcula DESPUÉS del ensamblado
Cobertura real          ← se calcula DESPUÉS (necesita definir el umbral primero)
```

Calidad del experimento: si profundidad media real ≈ profundidad target → pocas lecturas se perdieron en ruido.

**Estándar clínico:** 30X para variant calling confiable. "Billetera mata galán": más cobertura → mejor resultado.

### Phred Quality Score

```
Q = −10 · log₁₀(P_error)

Q20 → P = 0.01  (1 error cada 100 bases)   → 99%
Q30 → P = 0.001 (1 error cada 1000 bases)  → 99.9%  ← estándar mínimo aceptable
```

### Formato FASTQ

4 líneas por read:
```
@ID            ← identificador
ACGTACGT...    ← secuencia
+              ← separador
IIHHGF...      ← calidades Phred en ASCII (Phred+33)
```

---

### Solución 2: 2ª Generación — paralelizar millones de reacciones

**Idea central:** en lugar de hacer una reacción de Sanger por fragmento, hacer millones de reacciones en paralelo en el mismo chip. El costo por letra baja porque los reactivos se reparten entre millones de lecturas simultáneas.

**Restricción compartida de todas las plataformas 2ª gen:** la señal de una sola molécula de ADN es indetectable. Para amplificar la señal hay que amplificar la cantidad de moléculas → **amplificación clonal** previa a la secuenciación. Cada fragmento se transforma en un cluster/colonia de miles de copias idénticas. Esto impone lecturas cortas (~100–600 pb) porque la PCR de amplificación clonal no funciona bien con fragmentos muy largos.

#### 2ª Gen — 454 Pyrosequencing (Roche, históricamente primera)

**Analogía general:** imaginate miles de tubitos de ensayo en miniatura (pocillos), cada uno haciendo su propia reacción de Sanger en paralelo. En lugar de leer el color de los fragmentos en un gel, detectás luz cuando se incorpora un nucleótido.

**Paso 1 — Amplificación clonal en emulsión (PCR en gotitas):**

```
Tengo mis fragmentos de ADN con adaptadores en los extremos.
Los mezclo con bolitas (beads) cubiertas de primers que se pegan a esos adaptadores.
Proporciones calibradas → en promedio 1 fragmento se pega a 1 bolita.

Ahora agito esa mezcla con aceite → se forma una EMULSIÓN:
el agua queda en gotitas microscópicas suspendidas en aceite.

Aceite  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        (  gotita 1  )   (  gotita 2  )   (  gotita 3  )
        [ bolita+frag1]  [ bolita+frag2]  [ bolita+frag3]
        (  PCR aquí  )   (  PCR aquí  )   (  PCR aquí  )
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Como cada gotita es un reactor aislado (el aceite la separa de las demás),
la PCR dentro de cada gotita amplifica SOLO el fragmento que tiene esa bolita.
Al terminar: cada bolita está rodeada de miles de copias del MISMO fragmento.
```

**Paso 2 — Cargar en la PicoTiterPlate:**

```
La PicoTiterPlate es una placa con miles de pocillos, cada uno del tamaño
exacto de 1 bolita. Como 1 bolita = 1 secuencia, cada pocillo tiene 1 secuencia distinta.

 ┌──┬──┬──┬──┬──┐
 │● │● │● │● │● │  ← bolitas con fragmentos amplificados
 ├──┼──┼──┼──┼──┤
 │● │● │● │● │● │  (en realidad hay cientos de miles de pocillos)
 └──┴──┴──┴──┴──┘
```

**Paso 3 — Secuenciación por pirofosfato (detección de luz):**

```
Agrego a todos los pocillos: polimerasa + sulfurilasa + luciferasa + luciferina.

Luego ciclo nucleótidos de a UNO:

Ciclo 1: agrego T a toda la placa
  → en los pocillos cuya siguiente base a copiar es A (complemento de T),
    la polimerasa incorpora la T y libera PIROFOSFATO (PPi)
  PPi → sulfurilasa → ATP → luciferasa → LUZ  💡
  
  Foto de la placa:  pocillo 1=luz💡  pocillo 2=oscuro  pocillo 3=luz💡 ...
  → registro: pocillo 1 incorporó T, pocillo 3 incorporó T

Lavo todo el T no incorporado.

Ciclo 2: agrego G
  → nueva foto → sé cuáles incorporaron G

Ciclo 3: agrego A ... y así sucesivamente.

Leyendo la secuencia de ciclos en los que cada pocillo emitió luz → secuencia de ese fragmento.
```

**⚠ Problema de los homopolímeros:**

```
¿Qué pasa si la secuencia tiene TTTT (4 T's seguidas)?

Cuando agrego T en el ciclo correspondiente, la polimerasa pone las 4T de GOLPE
(no para entre medio) → se liberan 4× más pirofosfato → 4× más luz.

El equipo calibra: 1 luz = 1T, doble luz = 2T, triple = 3T...
Pero con 5 o más T iguales seguidas, la escala se vuelve imprecisa → errores de lectura.

AAAAA → ¿5 A's? ¿6 A's? El equipo no lo sabe bien. Eso es el problema del homopolímero.
```

#### 2ª Gen — Ion Torrent (misma arquitectura, distinta detección)

La reacción de la luciferasa (454) requiere agregar enzimas caras que se degradan. Ion Torrent dice: "¿para qué luz si puedo medir directamente la química?" Cuando la polimerasa incorpora un nucleótido, libera un H⁺. Ion Torrent pone un **sensor de pH debajo de cada pocillo** — básicamente un pHímetro microscópico.

```
Arquitectura: IDÉNTICA a 454 hasta el paso 2
  (emulsión PCR → bolitas → pocillos de la PicoTiterPlate)

Paso 3: DETECCIÓN por cambio de pH

 ┌─────────────────┐
 │  pocillo con    │
 │  bolita+poli-   │  ← agrego nucleótido T
 │  merasa         │
 ├─────────────────┤
 │  sensor de pH   │  ← semiconductor de silicio
 └────────┬────────┘
          │
     señal eléctrica

Si la polimerasa incorpora la T → libera H⁺ → el pH cae → el sensor detecta la caída
Si no había que incorporar T (no es complementario) → no hay H⁺ → pH no cambia → silencio

Ciclo de nucleótidos: igual que 454 (T → lavar → G → lavar → A → lavar → C → lavar → ...)
Mido el cambio de pH en cada ciclo → sé qué base se incorporó en cada pocillo.
```

Ventaja: más barato y más rápido (no necesita luciferasa ni luciferina).
⚠ **Mismo problema de homopolímeros**: TTTT libera 4 H⁺ → mayor caída de pH, pero difícil distinguir 4T de 5T.

#### 2ª Gen — Illumina (dominante, ~70% del mercado)

Illumina resolvió la amplificación clonal de forma más elegante: en lugar de emulsionar en gotitas de aceite, **ancla los fragmentos directamente a una superficie sólida** (la flow cell) y los amplifica ahí mismo formando colonias locales llamadas clusters.

**Paso 1 — Bridge Amplification (amplificación en puente):**

```
La flow cell es una placa de vidrio con millones de primers cortos pegados a su superficie.

1. El fragmento (con adaptadores) se pega a uno de esos primers y se copia:

   superficie  ──────────────────────────────────────
                ↑primer  ↑primer  ↑primer  ↑primer
               [frag1]

2. El fragmento copiado se dobla y pesca el primer de al lado → forma un "puente":

   superficie  ──────────────────────────────────────
               (──────)   ← puente entre los dos primers
              original    copia

3. La polimerasa copia el puente → 2 copias. Se repite → 4 → 8 → 16 → ...
   Todas las copias quedan en el MISMO punto de la superficie = CLUSTER.

   superficie  ──────────────────────────────────────
               [cluster 1]   [cluster 2]   [cluster 3]
               (frag A ×1000) (frag B ×1000) (frag C ×1000)

1 flow cell tiene ~1.600 millones de clusters = 1.600 millones de reacciones en paralelo.
```

**Paso 2 — Secuenciación por síntesis con terminadores reversibles:**

Los 4 nucleótidos tienen DOS modificaciones: un fluoroforo de color distinto por base, y un bloqueo en el -OH 3' que impide incorporar más de uno por ciclo.

```
Ciclo 1: agrego los 4 dNTPs modificados a la flow cell.
  Cada cluster incorpora exactamente 1 nucleótido (el bloqueo frena ahí).
  Saco UNA FOTO de toda la flow cell:

  [cluster1: 🔵]  [cluster2: 🔴]  [cluster3: 🟢]  [cluster4: 🟡] ...
      A                T                C                G

  Paso de química: corto el bloqueo y el fluoroforo → el nucleótido queda libre.

Ciclo 2: agrego de nuevo los 4 dNTPs → cada cluster incorpora su siguiente base → nueva foto.
  [cluster1: 🔴]  [cluster2: 🔵]  [cluster3: 🔴]  [cluster4: 🟢] ...
      T                A                T                C

Ciclo 3: foto → ciclo 4: foto → ... → N ciclos = N bases leídas por cluster.

Leyendo el color del cluster 1 en cada foto: 🔵🔴🟢... → A, T, C... = su secuencia.
```

Ventaja: altísima precisión (~0.1% de error), sin problema de homopolímeros (se incorpora de a 1 aunque haya TTTT).
Desventaja: lecturas cortas (75–300 pb) — la señal se degrada después de muchos ciclos.

---

### Solución 3: 3ª Generación — eliminar la amplificación clonal

Las plataformas de 2ª generación tienen dos costos inherentes: (1) reactivos para la amplificación clonal y (2) lecturas cortas forzadas por ella. La 3ª generación elimina ambos detectando **molécula única** — para lo cual se necesita amplificar la señal por medios físicos en lugar de multiplicar las moléculas.

#### 3ª Gen — PacBio SMRT (secuenciación por síntesis, molécula única)

**El problema a resolver:** la señal de fluorescencia de UN SOLO fluoroforo es demasiado débil para detectarla — por eso la 2ª gen amplificaba la cantidad de moléculas. PacBio resuelve esto de otra manera: en lugar de multiplicar las moléculas, **concentra la luz** en un volumen nanoscópico usando los Zero-Mode Waveguides (ZMW).

```
¿Qué es un ZMW? Un pozo extremadamente pequeño (diámetro de ~70 nm,
mucho menor que la longitud de onda de la luz visible ~500 nm).
Por las leyes de la óptica, la luz no puede propagarse dentro de algo más
pequeño que su longitud de onda → queda atrapada y concentrada en el fondo.

Dentro de cada ZMW:
  • 1 polimerasa fija en el fondo (no hay amplificación clonal)
  • 1 molécula de ADN que se introduce → la polimerasa empieza a copiarla
  • Los 4 dNTPs con fluoroforos distintos circulan en la solución encima

Cuando la polimerasa incorpora, por ejemplo, una T:
  • La T con su fluoroforo queda retenida en el fondo del ZMW durante ~ms
  • La luz concentrada en ese volumen hace que el fluoroforo emita con señal
    ~1000× más intensa que en solución → DETECTABLE a molécula única
  • Después del ms, el fluoroforo se libera y la polimerasa sigue con la siguiente base

  tiempo →  T💡  A💡  G💡  T💡  C💡 ...
              ↑    ↑    ↑    ↑    ↑
            color distinto por base → leés la secuencia en tiempo real
```

Sin necesidad de amplificar → el fragmento puede ser TAN LARGO COMO QUIERA la polimerasa copiarlo → lecturas de 1–50 kb o más.

Desventaja: la polimerasa comete más errores (~1% vs 0.1% de Illumina). Para compensar se secuencia con mayor profundidad o se hace circular el mismo fragmento varias veces.

---

#### 3ª Gen — Oxford Nanopore (sin síntesis, lectura directa)

PacBio todavía necesita polimerasa y nucleótidos. Nanopore hace algo conceptualmente más radical: **leer el ADN directamente**, sin sintetizar nada, igual que leer un libro sin tener que copiarlo.

```
Analogía: imaginá una manguera de agua. Si la tapás con el dedo, el flujo cae.
Si la tapás con algo más grande, cae más. Nanopore hace lo mismo pero con iones
y moléculas de ADN.

Configuración:
  • Una membrana con un PORO PROTEICO (derivado de una proteína de bacteria,
    modificado por ingeniería genética para que le encajen bien las 4 bases)
  • A un lado: voltaje positivo → los iones fluyen a través del poro → corriente medible
  • El ADN simple cadena, cargado negativamente, es "jalado" hacia el voltaje positivo
    → pasa a través del poro de a una base por vez

Lo que pasa mientras el ADN pasa:

  corriente
  sin ADN:  ──────────────────────────────────  (corriente constante = poro abierto)

  con ADN:  ───────┐    ┌──┐  ┌──────┐  ┌───   (cada base tapa el poro diferente)
            base A ↓T   ↓G  ↓G       ↓C  ↓A
                   (cada caída tiene amplitud característica de esa base)

  El cuello del poro tiene aminoácidos que forman puentes de hidrógeno con cada base
  de forma diferencial → A tapa diferente que T, que G, que C → caída distinta.
  El detector registra la corriente ms a ms → identifica cada base directamente.
```

Sin síntesis, sin amplificación → sin límite de largo de lectura. El mismo poro puede leer fragmentos de megabases.

Ventajas: lecturas larguísimas, portátil (tamaño USB), no gasta reactivos de síntesis.
Desventaja: los poros proteicos son descartables (el cartucho con ~5000 poros se usa una vez). Tasa de error más alta (~5–10%). Futuro: reemplazar los poros proteicos por poros sintéticos de grafeno o semiconductores → más estables, más baratos, reutilizables.

---

### Resumen comparativo de plataformas

| Plataforma | Gen | Ampl. clonal | Long. lectura | Error | Detección |
|---|---|---|---|---|---|
| 454 (Roche) | 2ª (discontinuado) | Emulsión PCR | ~400 pb | ~1% | Bioluminiscencia (PPi → luz) |
| Ion Torrent | 2ª | Emulsión PCR | 200–600 pb | ~1% | pH (H⁺) |
| Illumina | 2ª | Bridge PCR | 75–300 pb | ~0.1% | Fluorescencia + terminador reversible |
| PacBio SMRT | 3ª | No (mol. única) | 1–50 kb | ~1% | Fluorescencia en ZMW |
| Oxford Nanopore | 3ª | No (mol. única) | >100 kb | ~5–10% | Corriente iónica |

---

### Ensamblado de genomas — 3 pasos jerárquicos

```
LECTURAS (millones de fragmentos cortos de 100–300 pb)
    │
    │ Paso 1: Superponer lecturas → Contigs
    │   - Alinear todas las lecturas contra todas (alineamiento LOCAL de a pares)
    │   - Lecturas que solapan con alta identidad → unirlas
    │   - Contigs = secuencias contiguas sin huecos (se conoce la secuencia de punta a punta)
    ↓
CONTIGS (secuencias de cientos de kb, sin gaps, secuencia conocida completa)
    │
    │ Paso 2: Usar pair-ends para ordenar contigs → Scaffolds
    │   - Las lecturas apareadas vienen del mismo fragmento de ~1–5 kb
    │   - Si un pair-end tiene una lectura en Contig A y la otra en Contig B:
    │     → A y B están separados por (tamaño inserto - lecturas leídas) bases
    │   - Scaffolds = contigs ordenados con HUECOS DE TAMAÑO CONOCIDO entre ellos
    │   - (no sé qué bases van en el hueco, pero sé cuántas son)
    ↓
SCAFFOLDS (secuencias de Mb, con gaps de tamaño conocido)
    │
    │ Paso 3: Usar mapa físico (STS markers) → ordenar scaffolds → Cromosomas
    │   - STS (Sequence Tagged Sites): marcadores moleculares de posición conocida
    │   - Se mide distancia en unidades de Morgan (cM) por genética clásica de ligamiento
    │   - Si un scaffold contiene un STS que por genética sé que está en tal posición
    │     del cromosoma 7 → ubico ese scaffold en el cromosoma 7
    │   - Solo necesario para genomas grandes (humanos, plantas)
    │   - Bacterias/virus: suele alcanzar con los scaffolds o incluso con un solo contig
    ↓
GENOMA COMPLETO (secuencias por cromosoma)
```

**¿Dónde se guarda el ensamblado?**
Las bases de datos almacenan las secuencias de los **contigs** (no del genoma completo como un solo archivo). En RefSeq los contigs tienen prefijo `NT_`. Al buscar un gen, siempre se obtiene la secuencia del contig correspondiente.

#### De Bruijn Graph (estándar para NGS)

El problema del ensamblado con millones de lecturas: no se puede hacer MSA ni comparar todas contra todas (O(n²)).

**Solución Overlap-Layout-Consensus (OLC):** generar un grafo donde los nodos son las lecturas y las aristas son solapamientos. Para genomas con millones de reads → grafo gigante e inmanejable.

**Solución De Bruijn:** fijar el número de nodos de antemano → el tamaño del grafo no depende de la cantidad de lecturas.

```
Definir k (longitud del k-mer)
Nodos del grafo = todos los (k-1)-mers posibles (exactamente 4^(k-1) nodos posibles)
Para cada lectura:
    Partir la lectura en k-mers consecutivos
    Cada k-mer = arista entre su prefijo (k-1) y su sufijo (k-1)
    Si la arista ya existe → aumentar su peso (más lecturas la soportan)

Ensamblado = camino Euleriano (recorre todas las aristas una vez)
```

**Ventaja:** agregar más lecturas solo agrega peso a aristas existentes, no nuevos nodos → tamaño del grafo controlado.

**Problemas del grafo (y cómo se originan):**
- **Bifurcaciones:** errores de secuenciación crean k-mers alternativos (ej: ...CG vs ...TG en el mismo lugar)
- **Burbujas:** errores o repeticiones generan caminos paralelos que convergen. Hay que discernir si son errores o variantes reales

| k chico | k grande |
|---------|---------|
| Más conexiones, más sensible a errores | Menos conexiones, menos errores |
| Más rutas posibles (ambigüedades) | Más gaps (requiere mayor cobertura para conectar todo) |

**N50:** longitud L tal que el 50% del ensamble está en contigs de longitud ≥ L. A mayor N50, mejor ensamblado.

**SPAdes:** algoritmo multi-k (usa varios valores de k) — el más usado actualmente.

---

### Proyecto Genoma Humano — Historia y dos estrategias

#### Contexto histórico

| Año | Evento |
|-----|--------|
| 1988-1990 | Congreso de EEUU aprueba el HGP; Watson es director inicial; objetivo: 95% del genoma en 15 años |
| 1990-1994 | Primera etapa: generar **mapa físico** del genoma (marcadores moleculares en posiciones conocidas de cromosomas) |
| 1993 | Se suma Sanger Center (UK); también Japón, Francia, Italia → consorcio internacional |
| ~1993 | Venter (NIH) inventa secuenciación de cDNAs (ESTs): obtiene secuencias de 2000 genes de cerebro en días |
| ~1994 | Escándalo de patentes: NIH intenta patentar los 2000 genes sin avisar a Venter → Watson renuncia; lo reemplaza Francis Collins |
| 1995 | Venter (ya fuera del NIH) publica **primer genoma bacteriano** usando Whole Genome Shotgun → prueba de concepto |
| 1998 | Venter funda **Celera Genomics**; anuncia que secuenciará el genoma humano más rápido que el consorcio |
| 2000 | Venter publica **genoma de Drosophila melanogaster** (primer pluricelular por WGS); el consorcio se asusta |
| 2001 | **Anuncio conjunto** (mediado por Clinton): consorcio público publica en **Nature**; Celera publica en **Science** |

> Watson era editor de Nature y según la leyenda había prohibido publicar trabajos de Venter → Venter publicó en Science.

#### Dos estrategias de ensamblado

**Consorcio público — BAC-by-BAC (jerárquico):**

```
Mapa físico
    ↓
Fragmentar genoma en BACs (cromosomas artificiales bacterianos, ~500kb – 1Mb)
    ↓
Secuenciar cada BAC por separado (en labs distribuidos, ~1000 labs)
    ↓
Ensamblar cada BAC individualmente (grafo pequeño → más fácil de resolver)
    ↓
Ordenar BACs usando el mapa físico → genoma completo
```

Cobertura: **1.5×** → 92% del genoma cubierto. Publicado en Nature.

**Celera — Whole Genome Shotgun (WGS):**

```
Fragmentar TODO el genoma en pedazos de ~500-1000 bases
    ↓
Secuenciar 27 millones de lecturas de ~500 pb
    (+ 3 tamaños de insertos: 2kb, 10kb, 50kb → para armar scaffolds)
    ↓
Filtrar lecturas basura (repeticiones, contaminación)
    ↓
Overlapper: alinear todas vs todas, umbral ≥40 pb con ≤6% diferencias
    ↓
Unitigger: grafo de solapamiento → unitigs (73.6% del genoma)
    ↓
Scaffolder: usar pares de lecturas (2kb, 10kb, 50kb) para vincular contigs
    ↓
Gap-filling: meter lecturas sueltas en los agujeros entre contigs
```

Además Celera tomó las secuencias públicas del consorcio, las fragmentó simulando lecturas de 500 pb → obtuvo 16M de lecturas adicionales → cobertura total: **8×** (vs 1.5× del consorcio).

**¿Por qué Celera obtuvo mejor resultado?** Por mayor **cobertura**, no por mejor algoritmo. A mayor cobertura, más lecturas por posición → más fácil resolver ambigüedades → menos gaps.

**Ventaja del enfoque BAC:** cada BAC genera un grafo pequeño (~10,000 lecturas vs 27 millones). Grafo más pequeño = más fácil de resolver = mejor calidad por BAC. Pero la menor cobertura global de 1.5× limitó el resultado final.

#### ¿Qué es el genoma humano de referencia?

El genoma de referencia **no es el genoma de ninguna persona particular** — es un **consenso** construido a partir de múltiples donantes (consorcio) o de un pequeño grupo (Celera: los 5 fundadores).

El genoma de referencia tiene dos componentes:
1. **La secuencia**: las ~3.000 millones de letras de A, T, C, G ordenadas por cromosoma
2. **El mapa de anotación**: para cada posición, qué función biológica tiene (exón, intrón, promotor, intergénica, ARN no codificante, etc.)

El genoma de referencia sirve como punto de comparación para determinar las **variantes** de cada individuo.

#### Estadísticas del genoma humano

| Dato | Valor |
|------|-------|
| Tamaño haploide | ~3 Gb (3.000 millones de bases) |
| Tamaño diploide | ~6 Gb |
| Genes codificantes | ~20,000–25,000 |
| Exones (% del genoma) | ~1.5% |
| Repeticiones | ~45–50% |
| Genes sin función conocida (al publicarse en 2001) | ~41–45% |
| Diferencia entre dos individuos cualesquiera | ~1 por mil → **1–5 millones de variantes** |

> La complejidad humana no viene del número de genes (el gusano C. elegans tiene ~20,000 genes también), sino del splicing alternativo, modificaciones post-traduccionales y redes de regulación.

---

### BWT (Burrows-Wheeler Transform) — Algoritmo de mapeo

**Problema:** tengo el genoma de referencia (3.000 millones de bases) y una lectura de ~100 bases. ¿Cómo encuentro en qué posición del genoma cae esa lectura de forma eficiente?

**Solución naive:** alinear la lectura contra cada posición del genoma → O(N×m) = inviable para millones de lecturas.

**Solución BWT:** **indexar el genoma una vez** → buscar cualquier sub-secuencia en O(m). Base de BWA y Bowtie2 (10-20× más rápidos que otros alineadores).

#### Construcción del índice BWT (ejemplo con "GOGOL")

```
Paso 1: Agregar símbolo $ al final y circularizar
Secuencia: G O G O L $  (se lee en círculo: después de $ vuelve G)

Paso 2: Generar todas las rotaciones (una por posición)
  0: G O G O L $
  1: O G O L $ G
  2: G O L $ G O
  3: O L $ G O G
  4: L $ G O G O
  5: $ G O G O L

Paso 3: Ordenar las rotaciones alfabéticamente
($ < G < L < O en este ejemplo)
  $ G O G O L     ← fila 5
  G O G O L $     ← fila 0
  G O L $ G O     ← fila 2
  L $ G O G O     ← fila 4
  O G O L $ G     ← fila 1
  O L $ G O G     ← fila 3

Paso 4: El arreglo de sufijos SA = posición original de cada fila = [5, 0, 2, 4, 1, 3]
Paso 5: La BWT = última columna = L $ O O G G
```

El **arreglo de sufijos (SA)** es el índice: mapea posición en la tabla ordenada → posición en la secuencia original.

#### Búsqueda de una sub-secuencia con BWT

Para buscar una sub-secuencia (lectura) en el índice:

```
Ejemplo: buscar "GO" en "GOGOL"

1. Recorrer el índice buscando la primera fila que tiene "GO" como prefijo → fila 1 (índice 1)
2. Buscar la última fila que tiene "GO" como prefijo → fila 2 (índice 2)
3. El intervalo [1, 2] en el SA contiene las posiciones: SA[1]=0, SA[2]=2
   → "GO" aparece en posiciones 0 y 2 de "GOGOL" ✓
```

Para una lectura de largo m: se agrega una letra a la vez, el intervalo [Rmin, Rmax] se va achicando. Cuando Rmin = Rmax → la lectura aparece en **un solo lugar** del genoma → se mapea ahí.

**Tolerancia a mismatches:** si al agregar una letra el intervalo se vuelve vacío (la letra no calza), se ignora esa letra (como si no estuviera) y se continúa. Al final, la lectura se mapea con ese mismatch tolerado. Típicamente se toleran 2-5 mismatches.

**¿Qué pasa si Rmin ≠ Rmax al final?** La lectura está en múltiples posiciones (región repetitiva). Los programas descartan esas lecturas o las ponen tentativamente con un flag de calidad de mapeo bajo.

#### Resumen BWT

| Aspecto | Detalle |
|---------|---------|
| Complejidad de indexar | O(N) — se hace una sola vez |
| Complejidad de buscar una lectura | **O(m)** — proporcional al largo del read |
| Programa que lo usa | BWA-MEM, Bowtie2 |
| Qué indexa | El genoma de referencia |

---

### Pipeline Variant Calling (Genómica Humana)

```
FASTQ (reads NGS)
    ↓  BWA-MEM (mapeo contra genoma referencia con BWT)
BAM (reads mapeados — cada lectura tiene posición en el genoma)
    ↓  Picard (marcar duplicados de PCR, métricas de cobertura)
BAM deduplicado
    ↓  GATK HaplotypeCaller (variant calling)
VCF (variantes crudas)
    ↓  Filtrado por calidad (VQSR o hard filters)
    ↓  Anotación estructural + funcional (VEP / SnpEff / ANNOVAR)
        usando: ClinVar, gnomAD, OMIM, PolyPhen-2, SIFT
Variantes anotadas → filtrado progresivo
```

#### Llamado de variantes — la lógica

Una vez mapeadas las lecturas, en cada posición del genoma se cuenta:
- Cuántas lecturas tienen la misma base que la referencia
- Cuántas lecturas tienen una base distinta (variante)

Para que una variante sea "llamada":
- Debe aparecer en al menos **X lecturas** (umbral de profundidad de variante)
- Para un heterocigota en organismo diploide: esperar ~50% de lecturas con la variante; en la práctica se pide ≥25% para tolerar errores de secuenciación
- El resultado es un archivo **VCF** (Variant Call Format): cada línea = una variante con su posición, bases REF/ALT, calidad, y profundidad

#### Tipos de variantes

| Tipo | Descripción | Impacto |
|------|-------------|---------|
| **SNP / SNV** | Cambio de 1 base | Depende del cambio de aa |
| **Missense** | SNP que cambia el aminoácido | Puede alterar función |
| **Nonsense** | SNP que genera stop prematuro | Trunca la proteína |
| **Sinónima (silent)** | SNP que no cambia aa | Generalmente neutro |
| **Indel** | Inserción o deleción 1+ bases | Frameshift si no es múltiplo de 3 |
| **Frameshift** | Indel no múltiplo de 3 → cambia todo el marco de lectura | Severo |
| **SV** | Segmentos >50 pb (duplicaciones, inversiones, translocaciones) | Variable |

> Un genoma humano difiere del de referencia en **1–5 millones de variantes** (1 en cada 1000 bases). De esas, solo unas pocas son relevantes para una condición particular.

#### Filtrado progresivo de variantes

```
Todas las variantes del genoma         ~1.000.000 – 5.000.000
    ↓ con efecto en proteína (missense, nonsense, frameshift...)
                                        ~10.000 – 50.000
    ↓ raras en población (MAF < 1% en gnomAD)
                                        ~1.000 – 5.000
    ↓ segregan con la enfermedad (análisis familiar)
                                        ~10 – 100
    ↓ patogénicas según ClinVar / literatura
                                        1 – 5  ← candidatas causales
```

#### Anotación de variantes

**Anotación estructural:** dónde cae la variante en la estructura del gen
- Exón codificante → ¿cambia el aminoácido? ¿genera stop? ¿frameshift?
- Intrón → ¿afecta sitio de splicing?
- Promotor → ¿afecta expresión?
- Intergénica → efecto regulatorio posible

**Anotación funcional:** bases de datos relevantes
| BD | Qué aporta |
|----|-----------|
| **ClinVar** | ¿Ya se reportó esta variante? ¿Es patogénica/VUS/benigna? |
| **gnomAD** | ¿Qué tan frecuente es en la población sana? (variante patogénica para enfermedad rara → debe ser muy rara) |
| **OMIM** | ¿El gen está asociado a enfermedades mendelianas? |
| **PolyPhen-2 / SIFT** | Predicción computacional del efecto de la sustitución de aminoácido |
| **COSMIC** | ¿Es una mutación somática recurrente en cáncer? |

---

### Evolución del costo de secuenciación

| Año | Costo por genoma | Tecnología |
|-----|-----------------|-----------|
| 2001 (HGP) | ~2.700 millones USD | Sanger, consorcio de 1000 labs |
| 2001 (Celera) | ~100 millones USD | WGS + Sanger |
| 2007 (Watson) | ~1 millón USD | 454 pyrosequencing |
| 2008 | ~60.000 USD | Primeras plataformas NGS (ABI) |
| 2010 | ~5.000 USD | Illumina HiSeq |
| 2014 | **~1.000 USD** | Illumina HiSeq X |
| ~2020 | ~500 USD | Illumina NovaSeq |

> La caída supera ampliamente la Ley de Moore. Consecuencia directa: más genomas secuenciados → más datos → mejores bases de datos de variantes → diagnóstico más preciso.

---

### Bioinformática Traslacional y Medicina Personalizada

**Definición:** uso de herramientas bioinformáticas (análisis de genomas, variantes, expresión génica) directamente en el contexto clínico para diagnóstico, prevención o tratamiento.

**Concepto de medicina personalizada:**
- La medicina siempre fue "personalizada" (el médico considera el contexto del paciente)
- Lo nuevo es que ahora se puede agregar el **genoma individual** como fuente de información adicional, a costo razonable
- El diagnóstico de una enfermedad = genotipo + fenotipo + ambiente

**Aplicaciones:**
| Área | Ejemplo |
|------|---------|
| **Diagnóstico** | WES/WGS para enfermedades raras; diagnóstico prenatal no invasivo |
| **Prevención** | Riesgo genético de cáncer (BRCA1/2); estado portador de enfermedades recesivas |
| **Tratamiento** | Oncología de precisión (mutaciones del tumor → elección de quimio); farmacogenómica (cómo el paciente metaboliza una droga) |

#### Caso real: Bainbridge et al. 2011

**Pacientes:** dos hermanos con enfermedad neuro-muscular severa progresiva. No respondían a L-Dopa (tratamiento estándar para deficiencia dopaminérgica).

**Estrategia bioinformática:**
1. Modelo genético → herencia **autosómica recesiva** (padres sanos, hijos enfermos)
2. WGS de los hermanos y los padres
3. Buscar variantes raras que segreguen: los padres heterocigotos, los hijos homocigotos o heterocigotos compuestos
4. Filtrar por efecto en proteína + rareza en gnomAD/dbSNP

**Diagnóstico:** **heterocigotas compuestos para mutaciones en el gen SPR** (Sepiapterin Reductase):
- alelo 1: p.Arg150Gly (missense)
- alelo 2: p.Lys251X (nonsense — stop prematuro)

Ninguno de los dos alelos produce SPR funcional → **sin BH4** (tetrahidrobiopterina) → sin Dopamina **ni Serotonina** (BH4 es cofactor de la Tirosina Hidroxilasa y la Triptófano Hidroxilasa).

**Por qué no respondían a L-Dopa:** L-Dopa restaura la dopamina pero no la serotonina. Sin BH4, ambas vías están bloqueadas.

**Tratamiento correcto:** L-Dopa + 5-OH-Triptófano (precursor de serotonina) + BH4.

> **Heterocigota compuesto:** heterocigota a nivel de secuencia (dos mutaciones distintas, una en cada alelo) pero **homocigota funcional** (ningún alelo produce proteína funcional). Es la forma más común de herencia recesiva en poblaciones no consanguíneas.

---

## 10. Parcial 2018 — Resolución Completa

---

### Problema 1 — Programación Dinámica NW

**Secuencias:** CTGGCT (horizontal) y ATGCTG (vertical)
**Scoring:** Match = +2, Mismatch = 0, Gap = −1
**Algoritmo:** Needleman-Wunsch (alineamiento **global**)

#### a) Matriz de score completa

|   | − | C | T | G | G | C | T |
|---|---|---|---|---|---|---|---|
| **−** | 0 | −1 | −2 | −3 | −4 | −5 | −6 |
| **A** | −1 | 0 | −1 | −2 | −3 | −4 | −5 |
| **T** | −2 | −1 | 2 | 1 | 0 | −1 | −2 |
| **G** | −3 | −2 | 1 | 4 | 3 | 2 | 1 |
| **C** | −4 | −1 | 0 | 3 | 4 | 5 | 4 |
| **T** | −5 | −2 | 1 | 2 | 3 | **4** | **7** |
| **G** | −6 | −3 | 0 | 3 | 4 | **3** | **6** |

**Score final = 6** (celda esquina inferior derecha)

**Justificación de celdas clave:**

F(T,T) = F(5,6): T vs T → **match**
- Diagonal: F(4,5) + 2 = 5 + 2 = **7** ← elegido
- Arriba: F(4,6) − 1 = 4 − 1 = 3
- Izquierda: F(5,5) − 1 = 4 − 1 = 3

F(G,T) = F(6,6): G vs T → **mismatch**
- Diagonal: F(5,5) + 0 = 4
- Arriba: F(5,6) − 1 = 7 − 1 = **6** ← elegido
- Izquierda: F(6,5) − 1 = 3 − 1 = 2

F(G,C) = F(6,5): empate triple (diagonal=3, arriba=3, izquierda=3)

#### b) ¿Existe más de un alineamiento óptimo?

**Sí, existen exactamente dos.** El empate está en F(3,4) (fila G, columna G, cuarto elemento de fila 3):
- Diagonal: F(2,3) + 2 = 1 + 2 = 3 (G=G match)
- Izquierda: F(3,3) − 1 = 4 − 1 = 3

Ambos caminos son válidos → dos alineamientos óptimos.

**Alineamiento A** (camino diagonal en F(3,4)):
```
seq1 (CTGGCT): C  T  G  G  C  T  −
seq2 (ATGCTG): A  T  −  G  C  T  G
```

**Alineamiento B** (camino izquierda en F(3,4)):
```
seq1 (CTGGCT): C  T  G  G  C  T  −
seq2 (ATGCTG): A  T  G  −  C  T  G
```

#### c) Score de los alineamientos

Ambos tienen score = **6**.

**Verificación Alineamiento A:**
| Col | seq1 | seq2 | Score |
|-----|------|------|-------|
| 1 | C | A | mismatch = 0 |
| 2 | T | T | match = +2 |
| 3 | G | − | gap = −1 |
| 4 | G | G | match = +2 |
| 5 | C | C | match = +2 |
| 6 | T | T | match = +2 |
| 7 | − | G | gap = −1 |
**Total = 0+2−1+2+2+2−1 = 6** ✓

**Verificación Alineamiento B:**
| Col | seq1 | seq2 | Score |
|-----|------|------|-------|
| 1 | C | A | mismatch = 0 |
| 2 | T | T | match = +2 |
| 3 | G | G | match = +2 |
| 4 | G | − | gap = −1 |
| 5 | C | C | match = +2 |
| 6 | T | T | match = +2 |
| 7 | − | G | gap = −1 |
**Total = 0+2+2−1+2+2−1 = 6** ✓

---

### Problema 2 — Bases de Datos Primarias

Para cada categoría: clasificar según curación, redundancia, tipo (primaria/secundaria).

#### a) BD de nucleótidos: **RefSeq** (NCBI)

| Criterio | Clasificación | Justificación |
|----------|--------------|---------------|
| Curación | **Curada** | NCBI revisa y anota cada registro manualmente |
| Redundancia | **No redundante** | Un registro por molécula por organismo |
| Tipo | **Primaria** | Almacena secuencias crudas experimentales |

> Contraste: GenBank sería no curada y redundante.

#### b) BD de proteínas: **Swiss-Prot / UniProtKB reviewed**

| Criterio | Clasificación | Justificación |
|----------|--------------|---------------|
| Curación | **Curada** | Revisión manual por expertos; cada entrada validada con literatura |
| Redundancia | **No redundante** | 1 registro = 1 gen = 1 proteína; todas las isoformas en el mismo registro |
| Tipo | **Primaria** | Almacena datos directamente obtenidos de experimentos bioquímicos |

> Contraste: TrEMBL sería no curada y más redundante.

#### c) BD de expresión génica: **GEO** (NCBI)

| Criterio | Clasificación | Justificación |
|----------|--------------|---------------|
| Curación | **No curada** | Datos depositados por autores sin revisión sistemática del contenido biológico |
| Redundancia | **Redundante** | El mismo gen puede aparecer en miles de experimentos distintos; no hay criterio de unicidad |
| Tipo | **Primaria** | Almacena datos crudos de experimentos de expresión génica |

#### d) BD de organismo: **SGD** (Saccharomyces cerevisiae)

| Criterio | Clasificación | Justificación |
|----------|--------------|---------------|
| Curación | **Curada** | Equipo SGD revisa y actualiza manualmente la anotación funcional de cada gen |
| Redundancia | **No redundante** | Un registro por gen; variantes alélicas integradas en el mismo registro |
| Tipo | **Secundaria** | Integra y organiza datos de múltiples BD primarias (GenBank, UniProt, PubMed) + anotación funcional derivada |

#### e) BD de compuestos químicos: **PubChem** (NCBI)

| Criterio | Clasificación | Justificación |
|----------|--------------|---------------|
| Curación | **No curada** | Estructuras depositadas automáticamente desde múltiples fuentes; sin revisión manual sistemática |
| Redundancia | **Redundante** | El mismo compuesto puede aparecer bajo distintos nombres; PubChem intenta deduplicar (SID→CID) pero no elimina toda redundancia |
| Tipo | **Primaria** | Almacena datos de estructura química directamente obtenidos o sintetizados |

> Alternativa curada: **DrugBank** — curada, no redundante, secundaria.

#### f) BD de literatura: **Europe PMC** (EBI)

| Criterio | Clasificación | Justificación |
|----------|--------------|---------------|
| Curación | **Curada** | Los artículos han pasado por revisión editorial; metadatos (PMID, DOI, MeSH) revisados |
| Redundancia | **No redundante** | Cada artículo tiene un identificador único |
| Tipo | **Secundaria** | Recopila publicaciones que son en sí mismas productos secundarios del proceso científico |

#### g) BD libre: **PDB** (Protein Data Bank)

| Criterio | Clasificación | Justificación |
|----------|--------------|---------------|
| Curación | **Curada** | Deposición incluye validación estructural obligatoria (R-factor, clashscore, Ramachandran) |
| Redundancia | **Redundante** | La misma proteína puede tener cientos de estructuras (distintos ligandos, mutantes, condiciones) |
| Tipo | **Primaria** | Almacena datos crudos de experimentos estructurales (mapas de densidad electrónica, coordenadas atómicas) |

---

### Problema 3 — MSA, SP score, E-values, Prolina

**Contexto:** Motivo de N-glicosilación **N−X(≠P)−S/T**. MSA dado:
```
Prot1: L  E  N  K  T  V  A
Prot2: V  −  N  E  S  Y  A
Prot3: I  D  N  Q  T  I  A
```
**Penalidades gap affine:** apertura = −8, extensión = −3 → gap L=1: −8+(−3×1) = **−11**

#### a) Scores por par y árbol de similitud

**Par Prot1–Prot2** (S12):

| Pos | Prot1 | Prot2 | Score BLOSUM62 |
|-----|-------|-------|----------------|
| 1 | L | V | +1 |
| 2 | E | − | −11 (gap L=1) |
| 3 | N | N | +6 |
| 4 | K | E | +1 |
| 5 | T | S | +1 |
| 6 | V | Y | −1 |
| 7 | A | A | +4 |

**S(1,2) = 1−11+6+1+1−1+4 = 1**

**Par Prot1–Prot3** (S13):

| Pos | Prot1 | Prot3 | Score BLOSUM62 |
|-----|-------|-------|----------------|
| 1 | L | I | +2 |
| 2 | E | D | +2 |
| 3 | N | N | +6 |
| 4 | K | Q | +1 |
| 5 | T | T | +5 |
| 6 | V | I | +3 |
| 7 | A | A | +4 |

**S(1,3) = 2+2+6+1+5+3+4 = 23**

**Par Prot2–Prot3** (S23):

| Pos | Prot2 | Prot3 | Score BLOSUM62 |
|-----|-------|-------|----------------|
| 1 | V | I | +3 |
| 2 | − | D | −11 (gap L=1) |
| 3 | N | N | +6 |
| 4 | E | Q | +2 |
| 5 | S | T | +1 |
| 6 | Y | I | −1 |
| 7 | A | A | +4 |

**S(2,3) = 3−11+6+2+1−1+4 = 4**

**Resumen:**

| Par | Score |
|-----|-------|
| Prot1–Prot3 | **23** (más similares → se unen primero) |
| Prot2–Prot3 | 4 |
| Prot1–Prot2 | 1 (menos similares) |

**Árbol guía:**
```
        ┌── Prot1
   ┌────┤
   │    └── Prot3
───┤
   └──────── Prot2
```

**Justificación:** Prot1 y Prot3 comparten sustituciones conservativas sin gaps: L↔I (aliphatic), E↔D (ácidos), V↔I (aliphatic), T↔T (idéntico). Prot2 se aleja por el gap en posición 2 (−11) y la sustitución Y↔V poco conservativa (−1).

#### b) SP score total

SP = S(1,2) + S(1,3) + S(2,3) = 1 + 23 + 4 = **28**

Verificación columna por columna:

| Pos | S(1,2) | S(1,3) | S(2,3) | Subtotal |
|-----|--------|--------|--------|----------|
| 1 | 1 | 2 | 3 | 6 |
| 2 | −11 | 2 | −11 | −20 |
| 3 | 6 | 6 | 6 | 18 |
| 4 | 1 | 1 | 2 | 4 |
| 5 | 1 | 5 | 1 | 7 |
| 6 | −1 | 3 | −1 | 1 |
| 7 | 4 | 4 | 4 | 12 |
| **Total** | | | | **28** ✓ |

#### c) Ordenamiento de E-values: E12, E13, E23

E-value es exponencialmente inverso al score: E = K·m·n·e^(−λS).  
A mayor score → **menor** E-value (más significativo).

Scores: S13 = 23 > S23 = 4 > S12 = 1

**E₁₃ < E₂₃ < E₁₂**

El par 1–3 tiene el score más alto → menos probable de verse por azar → E-value más bajo = más significativo. El par 1–2 tiene score más bajo → más probable por azar → E-value más alto.

#### d) ¿Por qué la Prolina no se acepta en posición X del motivo?

**Razón estructural:** La Prolina es el único aa cuya cadena lateral forma un anillo con el nitrógeno del backbone (anillo pirrolidínico). Esto:

1. Elimina el H del N amídico → **no puede formar puentes H** con el backbone.
2. Restringe severamente el ángulo φ → **fuerza un quiebre (kink) en la estructura**.
3. **Impide que el tripéptido N-P-S/T adopte la conformación de β-turn** que necesita la enzima oligosacariltransferasa (OST) para reconocer y glicosilar el residuo N.

**¿Qué dice BLOSUM62?** La fila P en BLOSUM62 muestra scores negativos con prácticamente todos los aminoácidos (P−L = −3, P−I = −3, P−F = −4, P−W = −4, solo P−P = +7). Esto indica que la **Prolina es evolutivamente muy poco intercambiable** — sus sustituciones son seleccionadas negativamente porque alteran drásticamente la estructura local. La matriz refleja la excepcionalidad estructural de la prolina, que es exactamente la razón por la que elimina el sitio de glicosilación.

---

### Problema 4 — Bases de Datos Secundarias, PFAM y HMMer

#### A) Construcción del HMM con HMMer

**Problema bioinformático:** Dado un conjunto de sustratos confirmados de KC, construir un **profile-HMM** que capture la variabilidad posicional del motivo de reconocimiento para detectar nuevos sustratos en el proteoma de Av.

**Pipeline:**

```
1. Set de entrenamiento
   Secuencias de sustratos confirmados de KC
         ↓
2. MSA (CLUSTAL / MAFFT)
   Alinear las regiones de reconocimiento
         ↓
3. hmmbuild → KC_motivo.hmm
   Estima: P(aa | posición), P(match→insert), P(match→delete)
   = PROBLEMA DE ENTRENAMIENTO
         ↓
4. Calibración (automática con hmmbuild)
   Genera parámetros para calcular E-values
         ↓
5. hmmsearch → busca en proteoma de Av
   Lista de candidatos con score y E-value
   = PROBLEMA DE PUNTAJE
         ↓
6. Validación experimental (100 candidatos)
   Ensayos de unión y fosforilación
         ↓
7. Evaluación: curva ROC, sensibilidad, especificidad
```

**Métricas:**
- **Sensibilidad** = VP / (VP + FN): fracción de sustratos reales detectados.
- **Especificidad** = VN / (VN + FP): fracción de no-sustratos correctamente rechazados.
- Curva ROC: barrer umbrales de E-value → tradeoff sensibilidad/especificidad.

#### B) Interpretación del Logo del HMM

Un **sequence logo** representa cada posición del HMM como una columna de letras donde:
- **Altura total** = bits de información = conservación (máx ≈ 4.32 bits = log₂20)
- **Altura de cada letra** = frecuencia relativa (probabilidad de emisión) de ese aa en esa posición

**5 propiedades que se pueden leer:**

**1. Conservación posicional:** columnas altas = posiciones muy conservadas = **determinantes de especificidad** del reconocimiento por KC (si se mutan, eliminan la unión o la fosforilación).

**2. Identidad fisicoquímica del residuo dominante:** la letra más alta revela qué tipo de residuo requiere KC:
- K o R → quinasa basofílica (como PKA/PKC)
- D o E → quinasa acidofílica (como CK2)
- F, W, Y → bolsillo hidrofóbico aromático

**3. Tolerancia a sustituciones conservativas:** cuando una columna muestra 2–4 letras de altura apreciable, el modelo acepta variación conservativa. Si aparecen S y T juntas → KC solo requiere un hidroxilo; si K y R → requiere carga positiva, acepta cualquier básico.

**4. Posiciones variables:** columnas bajas = KC no impone restricciones ahí = residuos equivalentes a la "X" del motivo N-glicosilación.

**5. Extensión del motivo:** número total de columnas = largo del motivo lineal reconocido (20–25 aa en el enunciado).

#### C) Búsqueda en el proteoma de Av y falsos positivos

```bash
hmmsearch -E 0.001 KC_motivo.hmm proteoma_Av.fasta > candidatos.tbl
```

**Resultado:** Al evaluar 100 candidatos experimentalmente, 5 no se unen ni son fosforilados por KC.

**Nombre del resultado:** **Falsos Positivos (FP)**

El HMM predijo 100 proteínas como sustratos (positivos), pero 5 resultaron no serlo.
**Precisión = 95/100 = 95%.**

**¿Por qué ocurren los FP?**
1. El motivo existe en la secuencia pero está **enterrado en la estructura 3D** → KC no puede acceder in vivo.
2. Set de entrenamiento pequeño → el modelo aprendió señales espurias (sobreajuste).
3. **Contexto celular:** KC requiere co-localización, andamiajes o fosforilaciones previas no capturables por secuencia lineal.
4. Umbral de E-value permisivo → hits marginalmente similares al modelo.
