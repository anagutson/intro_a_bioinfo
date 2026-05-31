#!/usr/bin/env python3
"""
Ejercicio 2.2 - TP2 Bioinformatica 2026
Analiza BLAST del proteoma de S. aureus N315 contra VFDB usando Bio.SearchIO.

Uso:
    python ejercicio_2_2_biopython.py <blast_tab.tsv> <salida.tsv>

El archivo blast_tab.tsv debe generarse con:
    blastp -query N315.faa -db VFDB_db -evalue 1e-5 -max_target_seqs 10 \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen \
                    qstart qend sstart send evalue bitscore qlen slen \
                    qseq sseq stitle" \
           -out blast_tab.tsv
"""
import re
import sys
from pathlib import Path

import Bio.SearchIO as bpio
import pandas as pd

CAMPOS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen",
]

COLUMNAS_CON_TITULO = CAMPOS + ["qseq", "sseq", "stitle"]


def especie_del_titulo(titulo: str) -> str:
    """Extrae especie desde el titulo VFDB, que termina con [Especie]."""
    matches = re.findall(r"\[([^\[\]]+)\]", str(titulo))
    return matches[-1] if matches else ""


def main(blast_tab: str, out_tsv: str, min_cob_hit: float = 70.0) -> None:
    resultados = list(bpio.parse(blast_tab, "blast-tab", fields=COLUMNAS_CON_TITULO))
    tabla = pd.read_csv(blast_tab, sep="\t", names=COLUMNAS_CON_TITULO)
    titulos = {}
    if "stitle" in tabla.columns:
        for _, row in tabla.iterrows():
            titulos[(row["qseqid"], row["sseqid"])] = row.get("stitle", "")

    filas = []

    for query_result in resultados:
        query_id = query_result.id
        query_len = query_result.seq_len

        # Iterar hits
        for hit in query_result:
            hit_id = hit.id
            hit_len = hit.seq_len
            titulo = titulos.get((query_id, hit_id), "")
            especie = especie_del_titulo(titulo)

            # Iterar hsps
            for hsp in hit:
                cob_query = (hsp.query_end - hsp.query_start) * 100 / query_len if query_len > 0 else 0
                cob_hit = (hsp.hit_end - hsp.hit_start) * 100 / hit_len if hit_len > 0 else 0

                if cob_hit >= min_cob_hit:
                    filas.append({
                        "query": query_id,
                        "hit": hit_id,
                        "cobertura_query (%)": round(cob_query, 2),
                        "cobertura_hit (%)": round(cob_hit, 2),
                        "identidad (%)": round(hsp.ident_pct, 3),
                        "evalue": hsp.evalue,
                        "bitscore": hsp.bitscore,
                        "especie_hit": especie,
                    })
                    break  # mejor HSP por hit

    df = pd.DataFrame(filas)
    if df.empty:
        print("No se encontraron hits con los criterios especificados.")
        return

    # Conservar el mejor hit por query (menor e-value, mayor bitscore)
    df = df.sort_values(["query", "evalue", "bitscore"], ascending=[True, True, False])
    df = df.groupby("query", as_index=False).first()

    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"Proteinas de N315 con hits en VFDB (cobertura_hit >= {min_cob_hit}%): {len(df)}")
    print(f"Tabla guardada en: {out_tsv}")
    print()
    print(df[["query", "hit", "cobertura_query (%)", "cobertura_hit (%)", "identidad (%)", "evalue", "especie_hit"]].head(10).to_string(index=False))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Uso: python ejercicio_2_2_biopython.py <blast_tab.tsv> <salida.tsv> [min_cob_hit]")
        sys.exit(1)
    min_cob = float(sys.argv[3]) if len(sys.argv) > 3 else 70.0
    main(sys.argv[1], sys.argv[2], min_cob)
