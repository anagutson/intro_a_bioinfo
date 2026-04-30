# TP1 - Bases de datos primarias e introducción a Biopython

Bioinformática 2026 - FCEyN UBA Exactas.

Este TP combina exploración manual de bases de datos primarias (NCBI, UniProt)
con análisis programático de secuencias usando Biopython. Las búsquedas web
quedan documentadas como tablas con fecha de consulta y los análisis
programáticos usan copias locales de las secuencias para que el notebook sea
reproducible sin internet.

## Estructura

- `notebooks/TP1_Bioinformatica_Entrega_ejecutado.ipynb`: notebook ejecutado, autocontenido.
- `scripts/`: programas reutilizables del TP.
- `data/`: secuencias locales de mioglobina y datos de E. coli.
- `resultados/`: gráficos, tablas y archivos generados al correr los scripts y el notebook.
- `requirements.txt`: dependencias del entorno.

## Orden del notebook

1. **Objetivo 1a** — Recorriendo NCBI: bases de datos, campos relevantes, búsquedas con tags y respuestas a los 4 ejercicios 1a.
2. **Objetivo 1b** — Registro de mioglobina humana NP_976311.1: longitud, mRNA NM_203377.1, variantes de splicing, ubicación en cromosoma 22q12.3, genes adyacentes.
3. **Objetivo 1c** — UniProt P02144 / MYG_HUMAN: PTMs, exones, PDB y fecha de registro.
4. **Objetivo 2.1** — Generación de secuencias al azar (ADN, proteína uniforme y proteína con frecuencias sesgadas usando `numpy.random.choice`).
5. **Objetivo 2.2** — Frecuencias en secuencias y gráficos de barras.
6. **Objetivo 3** — Biopython básico: `Seq`, `complement`, `reverse_complement`, `transcribe`, `translate`, `gc_fraction`.
7. **Ejercicio 3.1** — Mioglobina: traducción directa vs reversa complementaria, comparación con la proteína real.
8. **Ejercicio 3.2** — Lectura de FASTA y GenBank con `SeqIO.read`. Incluye una demo de descarga programática con `Entrez.efetch` (con fallback si no hay internet).
9. **Ejercicio 3.3** — Lectura del registro UniProt en formato SwissProt.
10. **Objetivo 4.1** — Distribución de aminoácidos en ADN al azar vs mioglobina. Convergencia de la distribución con el tamaño de la secuencia. Comparación con la frecuencia esperada por el código genético (n_codones / 61).
11. **Objetivo 4.2** — Histograma de tamaños de proteínas de E. coli K-12 (UniProt reviewed), media y rango central del 70 %.
12. **Objetivo 4.3** — Distribución de aminoácidos en proteína al azar vs mioglobina.

## Entorno

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

## Cómo ejecutar

Abrir el notebook:

```bash
.venv/bin/jupyter notebook notebooks/TP1_Bioinformatica_Entrega_ejecutado.ipynb
```

Correr los scripts desde la terminal:

```bash
.venv/bin/python scripts/secuencias_azar.py --largo 900
.venv/bin/python scripts/biopython_mioglobina.py
.venv/bin/python scripts/descargas_biopython.py --email tu_email@example.com
```

## Datos en `data/`

| Archivo | Origen |
|---|---|
| `mioglobina_proteina.fasta` | NCBI Protein, RefSeq **NP_976311.1** (mioglobina humana). |
| `mioglobina_gen.fasta` | Región codificante del gen **MB**, descargada de NCBI Nucleotide. |
| `mioglobina_gen.gb` | Mismo registro en formato GenBank (con anotaciones). |
| `mioglobina_mrna.fasta` | mRNA **NM_203377.1**, NCBI RefSeq. |
| `mioglobina_proteina_uniprot.txt` | Registro **P02144 / MYG_HUMAN** de UniProt en formato SwissProt (texto). |
| `ecoli_k12_uniprot_reviewed_lengths.tsv` | Tabla de UniProt para *Escherichia coli* K-12 substr. MG1655 (taxon 83333), proteínas reviewed, columnas `Entry` y `Length`. |

Las búsquedas web (NCBI y UniProt) están documentadas en el notebook con la
fecha de consulta. Los conteos pueden cambiar si se rehacen, porque las
bases se actualizan continuamente.

## Archivos generados en `resultados/`

Generados por los scripts:

- `adn_azar.fasta`, `proteina_azar.fasta` — secuencias al azar guardadas en FASTA.
- `frecuencias_adn_azar.{csv,png}`, `frecuencias_proteina_azar.{csv,png}`, `frecuencias_traduccion_adn_azar.{csv,png}` — frecuencias de ADN al azar, proteína al azar y aminoácidos al traducir ADN al azar.
- `mioglobina_biopython_resumen.txt` — resumen de los registros FASTA con Biopython.

Generados por el notebook (prefijo `nb_`):

- `nb_frecuencias_adn_azar.png`, `nb_frecuencias_proteina_azar.png`, `nb_frecuencias_traduccion_adn_azar.png`
- `nb_comparacion_adn_azar_mioglobina.png`, `nb_comparacion_proteina_azar_mioglobina.png`
- `nb_convergencia_aa_azar.png` — convergencia de la distribución al esperado.
- `nb_codones_vs_frecuencia.png` — frecuencia esperada (codones/61) vs observada.
- `nb_histograma_tamanios_proteicos.png` — tamaños de proteínas de E. coli.
