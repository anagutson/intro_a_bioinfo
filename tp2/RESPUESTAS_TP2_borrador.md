# TP2 - Respuestas borrador

---

## Ejercicio 1.1 — Comparacion de mioglobinas con BLAST web

**Proteinas usadas:** mioglobina humana (P02144) vs. mioglobina de caballo (P68082) y cerdo (P02189).
Ambas tienen 154 aminoacidos, igual que la humana.

**i) ¿El alineamiento abarca la secuencia completa? ¿Local o global?**

Si, el alineamiento abarca la secuencia completa de ambas proteinas (cobertura 100% en los dos casos). BLAST es un algoritmo de alineamiento **local**, pero cuando las proteinas son muy homologas (identidad ~88%), el unico HSP significativo se extiende por toda la longitud de la secuencia, comportandose en la practica como un alineamiento global. Esto ocurre porque no hay regiones de baja similitud que hagan conveniente terminar el HSP antes de llegar al extremo de la proteina.

**ii) ¿Se observan gaps en el dot-plot?**

No se observan gaps. El dot-plot de ambas comparaciones (humana vs. caballo y humana vs. cerdo) muestra una **unica linea diagonal perfectamente recta y continua** desde la posicion 1 hasta la 154, sin quiebres ni desplazamientos laterales. Esto es consistente con los 0/154 (0%) gaps reportados en el alineamiento. Si hubiera gaps, la diagonal mostraria saltos o segmentos paralelos desplazados; si hubiera regiones sin homologia, aparecerian interrupciones en la linea.

**iii) Valores del alineamiento (BLOSUM62, gap costs 11:1):**

| Organismo | Score bits (raw) | E-value | Identities | Positives | Gaps |
|---|---|---|---|---|---|
| Caballo (P68082) | 278 (711) | 7e-103 | 136/154 (88%) | 142/154 (92%) | 0/154 (0%) |
| Cerdo (P02189) | 276 (705) | 6e-102 | 136/154 (88%) | 144/154 (93%) | 0/154 (0%) |

Ambas mioglobinas muestran alta identidad con la humana (~88%) y E-values cercanos a 0, indicando homologia fuerte y estadisticamente muy significativa. Interesante: caballo y cerdo tienen la misma identidad (88%) pero el cerdo tiene mas Positives (93% vs 92%), lo que indica que las sustituciones del cerdo tienden a ser mas conservativas (aminoacidos de propiedades similares).

**Efecto de cambiar parametros (human vs caballo como referencia):**

| Parametro | Score bits (raw) | E-value | Identities | Positives | Gaps |
|---|---|---|---|---|---|
| BLOSUM62, 11:1 (default) | 278 (711) | 7e-103 | 136/154 (88%) | 142/154 (92%) | 0/154 (0%) |
| BLOSUM80, 10:1 | 291 (666) | 3e-102 | 136/154 (88%) | 144/154 (93%) | 0/154 (0%) |
| BLOSUM62, 11:2 | 308 (711) | 5e-111 | 136/154 (88%) | 142/154 (92%) | 0/154 (0%) |
| BLOSUM62, 9:1 | 217 (711) | 3e-84 | 136/154 (88%) | 142/154 (92%) | 0/154 (0%) |

**Analisis de los cambios:**

- **BLOSUM80 (10:1):** El raw score bajo (666 vs 711) porque BLOSUM80 asigna puntajes distintos a cada par de aminoacidos. Sin embargo, el score en bits subio (291 vs 278) porque los parametros estadisticos (lambda, K) de BLOSUM80 hacen que el mismo nivel de homologia sea mas informativo. Los Positives subieron a 93%: BLOSUM80 reconoce mejor las sustituciones conservativas entre secuencias cercanas.

- **Gap costs 11:2 (mayor penalidad de extension):** El raw score no cambio (711) porque el alineamiento ya tenia 0 gaps — penalizar mas la extension no modifica un alineamiento que no tiene gaps. El score en bits subio mucho (308) y el E-value mejoro notablemente (5e-111 vs 7e-103) porque los parametros estadisticos cambian con los gap costs: con gaps mas costosos, es menos probable encontrar un alineamiento tan bueno por azar, haciendo el resultado mas significativo.

- **Gap costs 9:1 (menor penalidad de extension):** El raw score tampoco cambio (711) por la misma razon. Pero el score en bits bajo drasticamente (217) y el E-value empeoro mucho (3e-84 vs 7e-103). Con gaps baratos, las secuencias aleatorias pueden lograr scores altos mas facilmente, reduciendo la significancia estadistica del resultado.

---

## Ejercicio 1.2 — Busqueda de proteinas homologas a la mioglobina en humanos

**i) ¿Como mejorar la busqueda para evitar duplicados?**

Usar la base **RefSeq Select proteins** y restringir a *Homo sapiens*. RefSeq Select incluye una sola entrada representativa por gen, eliminando las multiples entradas de estructuras cristalograficas, isoformas y secuencias redundantes que aparecen en nr. Con nr sin filtro, la mioglobina misma aparece decenas de veces por distintos cristales o isoformas.

**ii) Proteinas homologas encontradas (corrida con RefSeq Select, Homo sapiens):**

| Accesion | Nombre | Score (bits) | E-value | Cobertura | Identidad | Comentario |
|---|---|---|---|---|---|---|
| NP_005359.1 | myoglobin isoform 1 | 312 | 1e-111 | 100% | 100% | La misma proteina (hit trivial) |
| NP_599030.1 | cytoglobin (CYGB) | 73.6 | 5e-17 | 95% | 28.67% | Globina intracelular, expresada en fibroblastos |
| NP_005323.1 | hemoglobin subunit zeta (HBZ) | 71.6 | 1e-16 | 96% | 28.19% | Hemoglobina embrionaria, la mas divergente de las Hb |

**iii) ¿Cual es la globina humana mas parecida a la Mb?**

De los resultados obtenidos, la mas parecida (excluyendo la propia Mb) es la **citoglobina (CYGB)** con E-value 5e-17 y 28.67% de identidad, levemente mejor que la hemoglobina zeta. La citoglobina es biologicamente la mas cercana a la mioglobina: ambas son proteinas monomeras intracelulares que unen oxigeno, a diferencia de las hemoglobinas que forman tetrámeros circulantes.

**iv) ¿Se encontraron todas las globinas humanas?**

No. Con RefSeq Select solo aparecieron 3 secuencias. La base **RefSeq Select** es muy restrictiva (una secuencia representativa por gen), por eso faltan: neuroglobina (NGB), hemoglobinas alfa (HBA1, HBA2), beta (HBB), delta (HBD), epsilon (HBE1), gamma (HBG1, HBG2), mu (HBM) y theta (HBQ1). Ademas, con el E-value por defecto (0.05) y BLOSUM62, las globinas mas divergentes pueden quedar justo por debajo del umbral de significancia.

**v) Parametros modificados para recuperar mas globinas:**

Usando **RefSeq protein** (en lugar de RefSeq Select) con los mismos parametros ya se recuperan mas hemoglobinas. Adicionalmente:
- Aumentar el E-value threshold a **10** para capturar homologos distantes
- Cambiar la matriz a **BLOSUM45** (optimizada para secuencias con ~45% de identidad, mas adecuada para buscar globinas divergentes)
- Reducir el **Word size a 2** (permite encontrar mas palabras semilla en secuencias divergentes)

**vi) PSI-BLAST:**

Se corrio PSI-BLAST con la misma query (mioglobina humana) contra RefSeq Select, restringido a *Homo sapiens*.

**Iteracion 1** (9 secuencias encontradas):

| Proteina | E-value | Identidad | Nota |
|---|---|---|---|
| myoglobin isoform 1 | 1e-111 | 100% | query |
| cytoglobin (CYGB) | 5e-17 | 28.67% | |
| hemoglobin subunit zeta (HBZ) | 1e-16 | 28.19% | |
| hemoglobin subunit theta-1 (HBQ1) | 5e-10 | 25.00% | nueva vs BLASTp |
| hemoglobin subunit gamma-2 (HBG2) | 5e-10 | 25.55% | nueva vs BLASTp |
| hemoglobin subunit gamma-1 (HBG1) | 7e-10 | 25.55% | nueva vs BLASTp |
| hemoglobin subunit alpha (HBA1) | 6e-09 | 27.52% | nueva vs BLASTp |
| hemoglobin subunit epsilon (HBE1) | 4e-08 | 24.09% | nueva vs BLASTp |
| hemoglobin subunit mu (HBM) | 7e-08 | 27.61% | nueva vs BLASTp |

Con solo una iteracion ya se recuperaron 6 proteinas que BLASTp no habia encontrado.

**Iteracion 2** (12 secuencias encontradas — 3 nuevas):

| Proteina | E-value | Identidad | Nota |
|---|---|---|---|
| hemoglobin subunit delta (HBD) | 6e-57 | 25.52% | nueva en iteracion 2 |
| hemoglobin subunit beta (HBB) | 9e-56 | 24.83% | nueva en iteracion 2 |
| neuroglobin (NGB) | 6e-14 | **17.76%** | nueva en iteracion 2 |

Se realizaron **2 iteraciones**. En la segunda aparecio la **neuroglobina (NGB)**, que tiene solo 17.76% de identidad con la mioglobina — BLASTp simple con BLOSUM62 nunca la hubiera recuperado. PSI-BLAST pudo encontrarla porque despues de la primera iteracion construye una **PSSM** (Position-Specific Scoring Matrix) que captura el patron de conservacion especifico del dominio globina en cada posicion, siendo mucho mas sensible que los puntajes genericos de BLOSUM para detectar homologos distantes.

---

## Objetivo 2: BLAST en linea de comando

---

## Ejercicio 2.1 — drrA, map2k7 y ndk contra VFDB (outfmt 6)

**Comando usado:**
```bash
blastp -query drrA.faa -db VFDB_setB_pro.fas -evalue 1e-5 -outfmt 6
```

Las 12 columnas del outfmt 6 por defecto son:
`qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore`

Para incluir el alineamiento se agregan `qseq` y `sseq`:
```bash
blastp -query drrA.faa -db VFDB_setB_pro.fas -evalue 1e-5 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"
```

**Resultados principales:**

| Query | Mejor hit en VFDB | Identidad (%) | E-value | Interpretacion |
|---|---|---|---|---|
| drrA | VFG018236 (drrA/sidM, *Legionella pneumophila* Ph-1) | 100.0 | 0.0 | Coincidencia exacta con factor de virulencia de VFDB. La proteina ES un factor de virulencia conocido. |
| map2k7 | VFG019987 (ppkA, *Pseudomonas aeruginosa*) | 30.2 | 1.63e-18 | Hit significativo pero parcial (solo dominio quinasa). map2k7 es una quinasa eucariota; la similitud refleja la conservacion del dominio catalogo quinasa, NO homologia funcional con un factor de virulencia. No deberia anotarse como virulencia. |
| ndk | VFG031478 (ndk, *Mycobacterium marinum* M) | 58.8 | 3.21e-52 | Hit fuerte contra NDK bacteriana anotada como factor de virulencia. NDK bloquea la maduracion del fagosoma. |

**Conclusion:** drrA y ndk tienen hits claros y biologicamente esperables. map2k7, siendo una quinasa eucariota, solo muestra similitud de dominio con quinasas bacterianas; no deberia anotarse como factor de virulencia solo por este resultado.

---

## Ejercicio 2.2 — Factores de virulencia de N315 (Bio.SearchIO)

**Script:** `scripts/ejercicio_2_2_biopython.py`

**Comando BLAST:**
```bash
blastp -query staph_aureus_N315_proteome.faa -db VFDB_setB_pro.fas \
  -evalue 1e-5 -max_target_seqs 10 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
  -out N315_vs_VFDB_blast.tsv
```

**Parseo con Bio.SearchIO:**
```python
import Bio.SearchIO as bpio
campos = ["qseqid","sseqid","pident","length","mismatch","gapopen",
          "qstart","qend","sstart","send","evalue","bitscore","qlen","slen"]
resultados = list(bpio.parse("N315_vs_VFDB_blast.tsv", "blast-tab", fields=campos))
for query_result in resultados:
    for hit in query_result:
        for hsp in hit:
            # calcular coberturas y filtrar...
```

**Resultado:** 503 proteinas de N315 con al menos un hit en VFDB con cobertura del hit >= 70%.
La tabla completa esta en `resultados/N315_vs_VFDB_biopython.tsv`.
Columnas: query, hit, cobertura_query (%), cobertura_hit (%), identidad (%), evalue, bitscore, especie_hit.

---

## Ejercicio 2.3 — Enterotoxinas en N315

**Script:** `scripts/ejercicio_2_3_biopython.py`

**Criterio de decision:** identidad >= 80%, cobertura query >= 80%, cobertura hit >= 80%.

**Enterotoxinas detectadas en N315:**

| Proteina N315 | Enterotoxina(s) | Decision | Identidad (%) | E-value |
|---|---|---|---|---|
| WP_000034846.1 | sep | unica | 99.6 | ~0 |
| WP_000278088.1 | sec3;sec2;sec1 | multiple | 100.0 | ~0 |
| WP_000475325.1 | selx | unica | 96.1 | 8.6e-146 |
| WP_000713847.1 | sei | unica | 100.0 | 1.5e-180 |
| WP_000736712.1 | seg | unica | 100.0 | ~0 |
| WP_000746599.1 | sel | unica | 99.2 | 2.3e-176 |
| WP_000821658.1 | sem;selv | multiple | 92.5 | 6.0e-167 |
| WP_001035599.1 | tsst-1 | unica | 99.1 | 3.6e-174 |
| WP_001236362.1 | sen | unica | 100.0 | ~0 |
| WP_010922839.1 | seo | unica | 99.6 | ~0 |

**Interpretacion del criterio "multiple":**
- WP_000278088.1 matchea sec3, sec2 y sec1 con identidad ~100%: las sec enterotoxinas son muy similares entre si y dificil de distinguir sin un criterio mas estricto.
- WP_000821658.1 matchea sem y selv: misma situacion, las variantes de enterotoxinas SEM/SELV son muy parecidas.

La tabla completa esta en `resultados/enterotoxinas_biopython.tsv`.
