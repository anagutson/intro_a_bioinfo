"""
Ejercicio 2.1: BLAST local con outfmt 6 para drrA, map2k7 y ndk contra VFDB.

Comandos de preparación (ejecutar en terminal desde tp2/data/):
    wget -O drrA.faa "https://rest.uniprot.org/uniprotkb/Q29ST3.fasta"
    # Bajar map2k7 y ndk de UniProt de la misma forma

    wget "http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz"
    gunzip VFDB_setB_pro.fas.gz
    makeblastdb -dbtype prot -in VFDB_setB_pro.fas

Búsqueda con outfmt 6 (tabular):
    blastp -query drrA.faa -db VFDB_setB_pro.fas -evalue 1e-5 -outfmt 6 -out ../resultados/drrA_vs_vfdb.tsv
    blastp -query map2k7.faa -db VFDB_setB_pro.fas -evalue 1e-5 -outfmt 6 -out ../resultados/map2k7_vs_vfdb.tsv
    blastp -query ndk.faa -db VFDB_setB_pro.fas -evalue 1e-5 -outfmt 6 -out ../resultados/ndk_vs_vfdb.tsv
"""

import pandas as pd
from pathlib import Path

RESULTADOS_DIR = Path("../resultados")

COLUMNAS_OUTFMT6 = [
    "qseqid", "sseqid", "pident", "length",
    "mismatch", "gapopen", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore",
]


def leer_blast_tsv(ruta):
    return pd.read_csv(ruta, sep="\t", names=COLUMNAS_OUTFMT6)


def analizar_proteina(nombre, ruta_tsv):
    df = leer_blast_tsv(ruta_tsv)
    print(f"\n=== {nombre} ===")
    print(df)
    # TODO: analizar e interpretar resultados


if __name__ == "__main__":
    analizar_proteina("drrA", RESULTADOS_DIR / "drrA_vs_vfdb.tsv")
    analizar_proteina("map2k7", RESULTADOS_DIR / "map2k7_vs_vfdb.tsv")
    analizar_proteina("ndk", RESULTADOS_DIR / "ndk_vs_vfdb.tsv")
