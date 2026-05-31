# TP2 - Alineamiento de secuencias y BLAST

Bioinformatica 2026.

Este directorio deja preparada la parte local del TP2 y una guia corta para completar la parte web con los valores que devuelve NCBI BLAST.

## Estructura

- `TP2 Alin-BLAST - 2026.pdf`: consigna original.
- `scripts/descargar_datos.sh`: descarga las proteinas `drrA`, `map2k7`, `ndk`, VFDB set B, proteoma de *Staphylococcus aureus* N315 y la base de enterotoxinas enlazada en la consigna.
- `scripts/run_blast.sh`: arma bases BLAST y ejecuta las busquedas.
- `scripts/analizar_blast.py`: resume archivos BLAST tabulares.
- `scripts/decidir_enterotoxinas.py`: aplica un criterio identidad/cobertura para enterotoxinas.
- `data/`: datos descargados.
- `resultados/`: salidas BLAST y tablas finales.

## Entorno

En macOS, BLAST+ se instala con:

```bash
brew install blast
```

Python:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

En esta maquina tambien funcionan los scripts con el entorno ya creado en TP1:

```bash
../tp1/.venv/bin/python scripts/analizar_blast.py resultados/N315_vs_VFDB.tsv --out resultados/N315_vs_VFDB_resumen.tsv --min-scov 70
```

## Ejecucion

Desde `tp2`:

```bash
bash scripts/descargar_datos.sh
.venv/bin/python scripts/preparar_enterotoxinas.py
bash scripts/run_blast.sh
.venv/bin/python scripts/analizar_blast.py resultados/N315_vs_VFDB.tsv --out resultados/N315_vs_VFDB_resumen.tsv --min-scov 70
.venv/bin/python scripts/decidir_enterotoxinas.py resultados/N315_vs_enterotoxinas.tsv
```

## Objetivo 1.1 - BLAST web entre mioglobinas

Secuencia humana usada en TP1: mioglobina humana `P02144 / MYG_HUMAN`.

Para comparar 2 o 3 mioglobinas de otros organismos, buscar en UniProt entradas revisadas de mioglobina, por ejemplo:

| Organismo | Entrada UniProt sugerida |
|---|---|
| Caballo | `P68082` |
| Cerdo | `P02189` |
| Foca/cetaceo u otro mamifero | elegir una entrada revisada de `myoglobin` |

Puntos para responder:

- Si el alineamiento cubre casi toda la longitud de ambas proteinas, BLAST lo informa como un alineamiento local pero, en la practica, se comporta como global para estas homologas cercanas porque el mejor HSP se extiende por casi toda la secuencia.
- Los gaps se ven como interrupciones/desplazamientos en el dot-plot y como guiones `-` en el alineamiento.
- `E value` cercano a 0, score alto, identidad/positivos altos y pocos gaps indican homologia fuerte.
- BLOSUM80 suele favorecer sustituciones conservadas entre secuencias mas parecidas; puede cambiar el score crudo y, si hay regiones divergentes, acortar o modificar el alineamiento.
- Aumentar el costo de apertura/extension de gaps penaliza inserciones/deleciones; bajar esos costos permite mas gaps o extensiones de gaps si mejoran el score.

## Objetivo 1.2 - Homologas humanas de mioglobina

Para evitar repeticiones en NCBI BLAST web conviene usar RefSeq, limitar a *Homo sapiens*, revisar solo proteinas no redundantes y agrupar isoformas/estructuras repetidas por gen.

Globinas humanas esperables al relajar parametros o usar PSI-BLAST: mioglobina `MB`, hemoglobinas alfa/beta y variantes (`HBA1`, `HBA2`, `HBB`, `HBD`, `HBE1`, `HBG1`, `HBG2`, `HBM`, `HBQ1`, `HBZ`), neuroglobina `NGB` y citoglobina `CYGB`. La globina mas parecida a mioglobina suele ser citoglobina o neuroglobina entre las globinas no hemoglobinicas, pero hay que reportar la que devuelva el BLAST con mayor score y menor E-value en la corrida realizada.

## Objetivo 2 - Scripts principales (Bio.SearchIO)

### Ejercicio 2.1 — drrA, map2k7, ndk vs VFDB

Resultados en `resultados/drrA_outfmt6.tsv`, `map2k7_outfmt6.tsv`, `ndk_outfmt6.tsv`.

Las 12 columnas de outfmt 6 por defecto: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore.
Para incluir el alineamiento agregar `qseq sseq` al final del string de outfmt.

### Ejercicio 2.2 — N315 vs VFDB con Bio.SearchIO

```bash
# 1. Correr BLAST con qlen y slen incluidos
blastp -query data/staph_aureus_N315_proteome.faa -db data/VFDB_db \
  -evalue 1e-5 -num_threads 4 -max_target_seqs 10 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
  -out resultados/N315_vs_VFDB_blast.tsv

# 2. Generar tabla final (usa Bio.SearchIO internamente)
python3 scripts/ejercicio_2_2_biopython.py \
  resultados/N315_vs_VFDB_blast.tsv \
  resultados/N315_vs_VFDB_biopython.tsv
```

Resultado: 503 proteinas de N315 con hits en VFDB (cobertura hit >= 70%).

### Ejercicio 2.3 — N315 vs enterotoxinas

```bash
# 1. Correr BLAST
blastp -query data/staph_aureus_N315_proteome.faa -db data/enterotoxinas_db \
  -evalue 1e-5 -num_threads 4 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
  -out resultados/N315_vs_enterotoxinas_blast.tsv

# 2. Clasificar con criterio identidad/cobertura >= 80%
python3 scripts/ejercicio_2_3_biopython.py \
  resultados/N315_vs_enterotoxinas_blast.tsv \
  resultados/enterotoxinas_biopython.tsv
```

Resultado: 10 enterotoxinas detectadas en N315 (sep, selx, seg, sel, sem/selv, sei, sen, seo, tsst-1, sec1/sec2/sec3).
