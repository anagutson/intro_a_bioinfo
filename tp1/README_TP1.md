# TP1 - Bases de datos primarias e introducción a Biopython

Entrega preparada para Bioinformática 2026.

## Estructura

- `notebooks/TP1_Bioinformatica_Entrega.ipynb`: notebook limpio.
- `notebooks/TP1_Bioinformatica_Entrega_ejecutado.ipynb`: notebook ejecutado, recomendado para entregar.
- `scripts/`: programas reutilizables del TP.
- `data/`: secuencias locales de mioglobina en FASTA, GenBank y SwissProt.
- `resultados/`: gráficos, tablas y archivos generados.
- `requirements.txt`: dependencias del entorno.

## Entorno

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

## Ejecución

```bash
.venv/bin/python scripts/secuencias_azar.py --largo 900
.venv/bin/python scripts/biopython_mioglobina.py
```

El script `scripts/descargas_biopython.py` usa internet y requiere indicar un email para NCBI:

```bash
.venv/bin/python scripts/descargas_biopython.py --email tu_email@example.com
```

El notebook está pensado para correr sin internet. Las consultas de NCBI/UniProt quedan documentadas con fecha de consulta y los análisis programáticos usan copias locales de las secuencias.
