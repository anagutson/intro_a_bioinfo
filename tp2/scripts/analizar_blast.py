#!/usr/bin/env python3
import argparse
import math
import re
from pathlib import Path

import pandas as pd


COLUMNS = [
    "query",
    "hit",
    "identidad",
    "largo_alineamiento",
    "mismatches",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
    "qseq",
    "sseq",
    "titulo_hit",
]


def species_from_title(title: str) -> str:
    if title is None or (isinstance(title, float) and math.isnan(title)):
        return ""
    matches = re.findall(r"\[([^\[\]]+)\]", str(title))
    return matches[-1] if matches else ""


def load_blast(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", names=COLUMNS)
    if df.empty:
        return df
    numeric_cols = [
        "identidad",
        "largo_alineamiento",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qlen",
        "slen",
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["cobertura_query"] = ((df["qend"] - df["qstart"]).abs() + 1) * 100 / df["qlen"]
    df["cobertura_hit"] = ((df["send"] - df["sstart"]).abs() + 1) * 100 / df["slen"]
    df["especie_hit"] = df["titulo_hit"].map(species_from_title)
    return df


def best_hits(df: pd.DataFrame, min_identity: float, min_qcov: float, min_scov: float) -> pd.DataFrame:
    if df.empty:
        return df
    filtered = df[
        (df["identidad"] >= min_identity)
        & (df["cobertura_query"] >= min_qcov)
        & (df["cobertura_hit"] >= min_scov)
    ].copy()
    filtered = filtered.sort_values(
        ["query", "evalue", "bitscore", "identidad", "cobertura_hit"],
        ascending=[True, True, False, False, False],
    )
    return filtered


def main() -> None:
    parser = argparse.ArgumentParser(description="Resume salidas BLAST outfmt 6 del TP2.")
    parser.add_argument("blast_tsv", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--min-identity", type=float, default=30.0)
    parser.add_argument("--min-qcov", type=float, default=0.0)
    parser.add_argument("--min-scov", type=float, default=70.0)
    parser.add_argument("--top-per-query", type=int, default=1)
    args = parser.parse_args()

    df = load_blast(args.blast_tsv)
    selected = best_hits(df, args.min_identity, args.min_qcov, args.min_scov)
    if args.top_per_query > 0 and not selected.empty:
        selected = selected.groupby("query", as_index=False, group_keys=False).head(args.top_per_query)

    cols = [
        "query",
        "hit",
        "cobertura_query",
        "cobertura_hit",
        "identidad",
        "evalue",
        "bitscore",
        "especie_hit",
        "titulo_hit",
    ]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    selected[cols].to_csv(args.out, sep="\t", index=False)
    print(f"Filas de entrada: {len(df)}")
    print(f"Filas seleccionadas: {len(selected)}")
    print(f"Tabla escrita: {args.out}")


if __name__ == "__main__":
    main()
