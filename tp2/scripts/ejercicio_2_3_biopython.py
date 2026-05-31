#!/usr/bin/env python3
"""
Ejercicio 2.3 - TP2 Bioinformatica 2026
Detecta enterotoxinas en el proteoma de S. aureus N315.

Nota: la base de enterotoxinas contiene IDs duplicados (ej. Gr11/seh-2p),
por lo que se usa pandas para parsear la salida tabular de BLAST.

Criterio de decision: identidad >= 80%, cobertura query >= 80%, cobertura hit >= 80%.
Si mas de un hit cumple el criterio para la misma proteina query, se indica 'multiple'.

Uso:
    python ejercicio_2_3_biopython.py <blast_tab.tsv> <salida.tsv>

El archivo blast_tab.tsv debe generarse con:
    blastp -query N315.faa -db enterotoxinas_db -evalue 1e-5 \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen \
                    qstart qend sstart send evalue bitscore qlen slen" \
           -out blast_tab.tsv
"""
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

COLUMNAS = [
    "query", "hit", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen",
]

MIN_IDENTITY = 80.0
MIN_COB_QUERY = 80.0
MIN_COB_HIT = 80.0


def main(blast_tab: str, out_tsv: str) -> None:
    df = pd.read_csv(blast_tab, sep="\t", names=COLUMNAS)

    # Calcular coberturas
    df["cob_query (%)"] = ((df["qend"] - df["qstart"]).abs() + 1) * 100 / df["qlen"]
    df["cob_hit (%)"] = ((df["send"] - df["sstart"]).abs() + 1) * 100 / df["slen"]

    # Aplicar criterio de decision
    filtrado = df[
        (df["pident"] >= MIN_IDENTITY) &
        (df["cob_query (%)"] >= MIN_COB_QUERY) &
        (df["cob_hit (%)"] >= MIN_COB_HIT)
    ].copy()

    if filtrado.empty:
        print("No se encontraron enterotoxinas con los criterios especificados.")
        return

    # Ordenar por e-value para quedarse con el mejor primero
    filtrado = filtrado.sort_values(["query", "evalue", "bitscore"], ascending=[True, True, False])

    candidatos = defaultdict(list)
    for _, row in filtrado.iterrows():
        candidatos[row["query"]].append(row["hit"])

    filas = []
    for query_id, hits in candidatos.items():
        nombres = list(dict.fromkeys(hits))  # eliminar duplicados manteniendo orden
        decision = "unica" if len(nombres) == 1 else "multiple"
        mejor = filtrado[filtrado["query"] == query_id].iloc[0]
        filas.append({
            "proteina_N315": query_id,
            "enterotoxina(s)": ";".join(nombres),
            "decision": decision,
            "mejor_identidad (%)": round(mejor["pident"], 3),
            "mejor_cob_query (%)": round(mejor["cob_query (%)"], 2),
            "mejor_cob_hit (%)": round(mejor["cob_hit (%)"], 2),
            "mejor_evalue": mejor["evalue"],
        })

    resultado = pd.DataFrame(filas).sort_values("proteina_N315")
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    resultado.to_csv(out_tsv, sep="\t", index=False)
    print(f"Enterotoxinas detectadas en N315: {len(resultado)}")
    print()
    print(resultado.to_string(index=False))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Uso: python ejercicio_2_3_biopython.py <blast_tab.tsv> <salida.tsv>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
