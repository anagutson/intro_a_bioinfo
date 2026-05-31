# TP2 — Alineamiento de secuencias y BLAST
## Bioinformática 2026


# Objetivo 1: Uso de BLAST en NCBI-web


## Objetivo 1.1 — Comparación de 2 secuencias

**¿Qué opciones de programas hay y para qué tipo de comparaciones sirven?**

En el formulario de BLAST aparecen cinco opciones:

- **blastn**: busca una secuencia nucleotídica en una base de nucleótidos
- **blastp**: busca una proteína en una base de proteínas (la que usamos en este ejercicio)
- **blastx**: traduce el query nucleotídico y lo busca en una base de proteínas
- **tblastn**: busca una proteína en una base de nucleótidos que el programa traduce
- **tblastx**: traduce tanto el query como la base y compara proteínas

**¿Qué otras matrices hay disponibles? ¿Conoce sus diferencias?**

Además de BLOSUM62, NCBI BLAST ofrece PAM30, PAM70, PAM250, BLOSUM45, BLOSUM50, BLOSUM80 y BLOSUM90. Las matrices BLOSUM se construyen a partir de bloques de secuencias alineadas: el número indica el porcentaje de identidad promedio de las secuencias usadas para construirla. BLOSUM62 es la más general y la que mejor funciona para la mayoría de las búsquedas. BLOSUM80 sirve cuando las secuencias a comparar son muy similares entre sí (>80% de identidad), porque puntúa mejor las sustituciones conservativas y penaliza más los cambios drásticos. BLOSUM45 es la opción para secuencias más divergentes. Las matrices PAM son una familia distinta, basada en modelos de evolución de aminoácidos, y se usan menos en la práctica.

**¿Qué particularidad tienen los gap costs? ¿A qué se debe?**

Los gap costs tienen dos componentes separados: **existence** (abrir el gap) y **extension** (cada residuo adicional dentro del gap). Esto se conoce como modelo de gap afín y existe porque biológicamente no es lo mismo abrir una inserción que extenderla. Abrir un gap corresponde a un evento evolutivo único (un indel), mientras que agregar un residuo más dentro de ese mismo gap es menos costoso desde el punto de vista evolutivo. Por eso la penalidad de apertura es mucho mayor que la de extensión (11 vs 1 en el default de BLOSUM62). Además, los valores disponibles cambian según la matriz porque cada combinación tiene sus propios parámetros estadísticos.

### Ejercicio 1.1

Comparamos la mioglobina humana (P02144, 154 aa) con la de caballo (P68082, 154 aa) y cerdo (P02189, 154 aa). Para todas las corridas usamos blastp desde el formulario "Align two or more sequences".

**i) ¿El alineamiento abarca la secuencia completa de las 2 proteínas? ¿Le parece que es un alineamiento local o global? Justifique.**

Sí, en los dos casos el alineamiento abarcó la secuencia completa (cobertura 100%). BLAST es un algoritmo de alineamiento local, pero cuando dos proteínas son tan similares como estas mioglobinas de mamíferos (~88% de identidad), el único HSP significativo termina extendiéndose por toda la longitud de la secuencia porque no hay ninguna región de baja similitud que haga conveniente cortarlo antes. El resultado se ve prácticamente igual a un alineamiento global, aunque algorítmicamente sigue siendo local.

**ii) ¿Se observan gaps en el dot-plot? ¿Cuántos? ¿Cómo se identifican?**

No hay gaps en ninguna de las dos comparaciones. El dot-plot muestra una única línea diagonal recta y continua de la posición 1 a la 154, sin ningún quiebre ni desplazamiento lateral, lo que es consistente con el 0/154 (0%) de gaps del alineamiento. Si hubiera gaps, la diagonal mostraría saltos o segmentos desplazados en paralelo: un gap en la query desplazaría la línea hacia abajo, un gap en el hit hacia la derecha. En este caso la diagonal es perfecta.

**iii) Analice los valores de Score, E-value, Identities, Positives y Gaps. ¿Qué cambia al modificar la matriz y los gap costs?**

Con parámetros default (BLOSUM62, gap costs 11:1):

| Organismo | Score bits (raw) | E-value | Identities | Positives | Gaps |
|---|---|---|---|---|---|
| Caballo (P68082) | 278 (711) | 7e-103 | 136/154 (88%) | 142/154 (92%) | 0/154 (0%) |
| Cerdo (P02189) | 276 (705) | 6e-102 | 136/154 (88%) | 144/154 (93%) | 0/154 (0%) |

Los dos alineamientos muestran una homología muy alta: E-values cercanos a 0, identidad del 88% y ningún gap. Los Positives (92-93%) superan las Identities porque varias de las posiciones que difieren son sustituciones conservativas, es decir, cambios entre aminoácidos de propiedades similares. Algo que nos llamó la atención: caballo y cerdo tienen exactamente la misma identidad (88%) pero el cerdo tiene un Positive más (93% vs 92%), lo que indica que sus sustituciones son ligeramente más conservativas que las del caballo.

Para analizar el efecto de los parámetros usamos humana vs caballo como referencia:

| Configuración | Score bits (raw) | E-value | Identities | Positives | Gaps |
|---|---|---|---|---|---|
| BLOSUM62, 11:1 (default) | 278 (711) | 7e-103 | 136/154 (88%) | 142/154 (92%) | 0/154 (0%) |
| BLOSUM80, 10:1 | 291 (666) | 3e-102 | 136/154 (88%) | 144/154 (93%) | 0/154 (0%) |
| BLOSUM62, 11:2 | 308 (711) | 5e-111 | 136/154 (88%) | 142/154 (92%) | 0/154 (0%) |
| BLOSUM62, 9:1 | 217 (711) | 3e-84 | 136/154 (88%) | 142/154 (92%) | 0/154 (0%) |

**a) BLOSUM80:** el raw score bajó de 711 a 666 porque BLOSUM80 tiene valores numéricos distintos para cada par de aminoácidos. Pero el score en bits subió (291 vs 278) y los Positives aumentaron levemente (93% vs 92%). Esto tiene sentido porque BLOSUM80 fue construida con secuencias más similares entre sí, entonces es más sensible a las sustituciones conservativas. El E-value empeoró un poco (3e-102 vs 7e-103) por diferencias en los parámetros estadísticos lambda y K propios de esta matriz.

**b) Gap costs 11:2 (mayor penalidad de extensión):** el raw score se mantuvo en 711 porque el alineamiento ya tenía 0 gaps, entonces penalizar más la extensión no cambia nada en el alineamiento en sí. Lo llamativo es que el score en bits subió bastante (308) y el E-value mejoró mucho (5e-111 vs 7e-103). Esto pasa porque con gaps más costosos, la probabilidad de encontrar un alineamiento igual de bueno por azar es menor, haciendo el resultado estadísticamente más significativo.

**c) Gap costs 9:1 (menor penalidad de extensión):** mismo razonamiento — el raw score quedó en 711 porque el alineamiento no tiene gaps. Pero el score en bits cayó fuerte (217) y el E-value empeoró bastante (3e-84 vs 7e-103). Al bajar el costo de los gaps, las secuencias aleatorias pueden alcanzar puntajes similares más fácilmente, con lo cual el mismo resultado vale estadísticamente mucho menos.

## Objetivo 1.2 — Búsqueda de proteínas homólogas a la Mb en humanos

### Ejercicio 1.2

**i) ¿Cómo mejoraría la búsqueda para obtener solo las homólogas a la Mb una sola vez?**

Usando **RefSeq Select proteins** con el organismo restringido a *Homo sapiens*. RefSeq Select tiene una sola entrada representativa por gen, lo que elimina todas las entradas duplicadas que aparecen en nr: distintas estructuras cristalográficas de la misma proteína, isoformas con accession diferente, secuencias redundantes de distintos proyectos de secuenciación. Con nr sin filtro, la mioglobina sola puede aparecer docenas de veces.

**ii) Describa brevemente cada proteína homóloga encontrada. Compare E-values, scores y alineamientos.**

Con BLASTp contra RefSeq Select en *Homo sapiens* con parámetros default encontramos 3 proteínas:

| Proteína | Accesion | Score (bits) | E-value | Cobertura | Identidad |
|---|---|---|---|---|---|
| myoglobin isoform 1 | NP_005359.1 | 312 | 1e-111 | 100% | 100% |
| cytoglobin (CYGB) | NP_599030.1 | 73.6 | 5e-17 | 95% | 28.67% |
| hemoglobin subunit zeta (HBZ) | NP_005323.1 | 71.6 | 1e-16 | 96% | 28.19% |

La mioglobina misma aparece primero con 100% de identidad (hit trivial, es la misma proteína). La citoglobina y la hemoglobina zeta tienen identidades muy bajas (~28%), lo que refleja la divergencia evolutiva de las globinas a lo largo de cientos de millones de años. Aun así, los E-values son significativos (5e-17 y 1e-16), lo que confirma que hay homología real aunque sea lejana. La diferencia de score entre la Mb (312 bits) y las otras dos (73-71 bits) es enorme, y se ve también en los alineamientos: con la Mb el alineamiento es perfecto, mientras que con CYGB y HBZ hay muchas diferencias y algunos gaps.

---

**iii) ¿Cuál es la globina humana más parecida a la Mioglobina? ¿Le sorprende algún resultado?**

La más parecida es la **citoglobina (CYGB)**, con E-value 5e-17 y 28.67% de identidad. Tiene sentido porque CYGB es, como la Mb, una globina intracelular monómera — a diferencia de las hemoglobinas, que circulan en la sangre como tetrámeros. Lo que nos sorprendió un poco es que la hemoglobina zeta (HBZ) aparezca con una identidad tan similar a la de CYGB (28.19%). HBZ es una hemoglobina embrionaria muy antigua y divergente, con diferencias tan grandes respecto a las otras Hb que termina pareciéndose a la Mb casi tanto como la citoglobina.

---

**iv) ¿Creen que encontraron todas las proteínas de la familia de las globinas que existen en el genoma humano?**

No. Con BLASTp y RefSeq Select solo aparecieron 3 proteínas. Sabemos que faltan al menos la neuroglobina (NGB), las hemoglobinas alfa (HBA1), beta (HBB), delta (HBD), épsilon (HBE1), gamma-1 (HBG1), gamma-2 (HBG2), mu (HBM) y theta-1 (HBQ1). Algunas de estas probablemente quedan por debajo del umbral de E-value con los parámetros por defecto dado lo baja que es su identidad con la Mb.

---

**v) Intente encontrar la mayor cantidad de globinas modificando los parámetros. ¿Qué parámetros modificaron y qué proteínas nuevas encontraron?**

Con PSI-BLAST (ver punto vi) recuperamos prácticamente todas las globinas sin tocar los parámetros de scoring. Si quisiéramos hacerlo con BLASTp simple, los cambios que probaríamos son: aumentar el E-value threshold a 10 para no perder hits distantes, cambiar la matriz a BLOSUM45 (más adecuada para secuencias con ~45% de identidad como estas globinas) y reducir el word size a 2 para que BLAST encuentre más palabras semilla en secuencias divergentes.

---

**vi) PSI-BLAST: ¿encontraron alguna proteína nueva? ¿Cuántas iteraciones hicieron?**

Hicimos **2 iteraciones** de PSI-BLAST con la Mb humana como query contra RefSeq Select en *Homo sapiens*.

**Iteración 1** — 9 secuencias con E < 0.005:

| Proteína | E-value | Identidad | Nueva respecto a BLASTp |
|---|---|---|---|
| myoglobin isoform 1 | 1e-111 | 100% | no |
| cytoglobin (CYGB) | 5e-17 | 28.67% | no |
| hemoglobin subunit zeta (HBZ) | 1e-16 | 28.19% | no |
| hemoglobin subunit theta-1 (HBQ1) | 5e-10 | 25.00% | sí |
| hemoglobin subunit gamma-2 (HBG2) | 5e-10 | 25.55% | sí |
| hemoglobin subunit gamma-1 (HBG1) | 7e-10 | 25.55% | sí |
| hemoglobin subunit alpha (HBA1) | 6e-09 | 27.52% | sí |
| hemoglobin subunit epsilon (HBE1) | 4e-08 | 24.09% | sí |
| hemoglobin subunit mu (HBM) | 7e-08 | 27.61% | sí |

Ya en la primera iteración aparecieron 6 globinas nuevas que BLASTp no había encontrado.

**Iteración 2** — 12 secuencias en total, 3 nuevas:

| Proteína nueva | E-value | Identidad |
|---|---|---|
| hemoglobin subunit delta (HBD) | 6e-57 | 25.52% |
| hemoglobin subunit beta (HBB) | 9e-56 | 24.83% |
| neuroglobin (NGB) | 6e-14 | 17.76% |

El hallazgo más interesante es la **neuroglobina (NGB)** en la segunda iteración, con solo 17.76% de identidad respecto a la Mb. BLASTp nunca la hubiera encontrado con esos parámetros porque esa identidad es demasiado baja para que el hit sea estadísticamente significativo usando BLOSUM62. PSI-BLAST pudo recuperarla porque después de la primera iteración construye una PSSM (Position-Specific Scoring Matrix) con todas las globinas encontradas. La PSSM asigna un puntaje diferente a cada aminoácido en cada posición de la secuencia según cuánto está conservada esa posición en la familia — es mucho más informativa que una matriz genérica como BLOSUM para detectar homólogos lejanos.

---

# Objetivo 2: Uso de BLAST en línea de comando

---

### Ejercicio 2.1

> **Repita la búsqueda con -outfmt 6 y visualice las columnas. Investigue cómo incluir el alineamiento. Repita para las otras 2 proteínas. ¿Qué puede concluir del análisis?**

El comando básico con outfmt 6:

```bash
blastp -query drrA.faa -db VFDB_setB_pro.fas -evalue 1e-5 -outfmt 6
```

Las 12 columnas del formato tabular por defecto son:
`qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore`

Para incluir las secuencias alineadas hay que agregar `qseq` y `sseq` en el string de outfmt:

```bash
blastp -query drrA.faa -db VFDB_setB_pro.fas -evalue 1e-5 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"
```

Resultados para las tres proteínas:

| Query | Mejor hit en VFDB | Identidad (%) | E-value |
|---|---|---|---|
| drrA | VFG018236 — drrA/sidM, *Legionella pneumophila* Ph-1 | 100.0 | 0.0 |
| map2k7 | VFG019987 — ppkA, *Pseudomonas aeruginosa* LESB58 | 30.2 | 1.63e-18 |
| ndk | VFG031478 — ndk, *Mycobacterium marinum* M | 58.8 | 3.21e-52 |

**drrA** tiene un match exacto (100% de identidad) con un factor de virulencia de VFDB: es la misma proteína, anotada en VFDB como factor de virulencia de *Legionella pneumophila*. **ndk** tiene un hit fuerte (58.8%, E-value 3e-52) contra NDK de *Mycobacterium marinum*, que está anotada como factor de virulencia porque bloquea la maduración del fagosoma; el resultado es biológicamente coherente. **map2k7** es el caso más interesante: tiene hits significativos, pero con solo 30% de identidad y cubriendo únicamente una región parcial de la proteína (el dominio quinasa). map2k7 es una quinasa eucariota de la vía de estrés JNK; la similitud refleja conservación del dominio catalítico quinasa, que está presente en muchas proteínas de organismos muy distintos. No sería correcto anotarla como factor de virulencia solo por este resultado: el BLAST detecta similitud de dominio, no equivalencia funcional.

---

### Ejercicio 2.2

> **Buscar todos los factores de virulencia de VFDB sobre el proteoma de S. aureus N315. Generar una tabla con: query, hit, cobertura de la query, cobertura del hit, identidad y especie del hit.**

Corrimos BLAST del proteoma completo de N315 contra VFDB con el siguiente comando (incluimos `qlen` y `slen` para poder calcular las coberturas en el script):

```bash
blastp -query staph_aureus_N315_proteome.faa -db VFDB_setB_pro.fas \
  -evalue 1e-5 -num_threads 4 -max_target_seqs 10 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq sseq stitle" \
  -out N315_vs_VFDB_blast.tsv
```

Para analizar la salida usamos Bio.SearchIO de Biopython (`scripts/ejercicio_2_2_biopython.py`):

```python
import Bio.SearchIO as bpio

campos = ["qseqid","sseqid","pident","length","mismatch","gapopen",
          "qstart","qend","sstart","send","evalue","bitscore","qlen","slen"]

resultados = list(bpio.parse("N315_vs_VFDB_blast.tsv", "blast-tab", fields=campos))

for query_result in resultados:
    for hit in query_result:
        for hsp in hit:
            cob_query = (hsp.query_end - hsp.query_start) * 100 / query_result.seq_len
            cob_hit   = (hsp.hit_end   - hsp.hit_start)   * 100 / hit.seq_len
            # filtrar por cobertura del hit >= 70% y armar la tabla
```

Con el criterio de cobertura del hit >= 70% identificamos **503 proteínas** del proteoma de N315 con similitud significativa a factores de virulencia de VFDB. La tabla completa está en `resultados/N315_vs_VFDB_biopython.tsv` con columnas: query, hit, cobertura_query (%), cobertura_hit (%), identidad (%), evalue, bitscore, especie_hit.

---

### Ejercicio 2.3

> **a) BLAST del proteoma de N315 contra la base de enterotoxinas. ¿Cuáles están presentes?**
> **b) Programar un criterio de decisión por identidad/cobertura. Si más de un hit cumple las expectativas, indicar que se encontró más de una.**

Usamos el siguiente criterio de decisión: identidad >= 80%, cobertura de la query >= 80% y cobertura del hit >= 80%. El script (`scripts/ejercicio_2_3_biopython.py`) lee la salida tabular, aplica estos tres filtros y para cada proteína query lista las enterotoxinas que los pasan. Si más de una enterotoxina cumple para la misma proteína, la decisión queda como `multiple`.

Enterotoxinas detectadas en el proteoma de *S. aureus* N315:

| Proteína N315 | Enterotoxina(s) | Decisión | Identidad (%) | E-value |
|---|---|---|---|---|
| WP_000034846.1 | sep | única | 99.6 | ~0 |
| WP_000278088.1 | sec3; sec2; sec1 | multiple | 100.0 | ~0 |
| WP_000475325.1 | selx | única | 96.1 | 8.6e-146 |
| WP_000713847.1 | sei | única | 100.0 | 1.5e-180 |
| WP_000736712.1 | seg | única | 100.0 | ~0 |
| WP_000746599.1 | sel | única | 99.2 | 2.3e-176 |
| WP_000821658.1 | sem; selv | multiple | 92.5 | 6.0e-167 |
| WP_001035599.1 | tsst-1 | única | 99.1 | 3.6e-174 |
| WP_001236362.1 | sen | única | 100.0 | ~0 |
| WP_010922839.1 | seo | única | 99.6 | ~0 |

Encontramos 10 enterotoxinas en total. Los dos casos `multiple` tienen una explicación biológica clara: sec1, sec2 y sec3 son variantes muy parecidas de la enterotoxina estafilocócica tipo C, con secuencias tan similares entre sí que con nuestro criterio de 80% de identidad no es posible distinguirlas; lo mismo pasa con sem y selv. En esos casos preferimos reportar `multiple` en lugar de elegir una al azar, ya que forzar una clasificación única sería incorrecto.
