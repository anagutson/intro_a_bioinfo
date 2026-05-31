#!/usr/bin/env python3
import argparse
from pathlib import Path

from analizar_blast import load_blast, best_hits


def main() -> None:
    parser = argparse.ArgumentParser(description="Clasifica enterotoxinas presentes en N315.")
    parser.add_argument("blast_tsv", type=Path)
    parser.add_argument("--out", type=Path, default=Path("resultados/enterotoxinas_presentes.tsv"))
    parser.add_argument("--min-identity", type=float, default=80.0)
    parser.add_argument("--min-qcov", type=float, default=80.0)
    parser.add_argument("--min-scov", type=float, default=80.0)
    args = parser.parse_args()

    df = load_blast(args.blast_tsv)
    selected = best_hits(df, args.min_identity, args.min_qcov, args.min_scov)
    if selected.empty:
        args.out.write_text("query\tenterotoxina\tdecision\n", encoding="utf-8")
        print("No hubo enterotoxinas que cumplan el criterio.")
        return

    rows = []
    for query, group in selected.groupby("query", sort=True):
        best = group.sort_values(["evalue", "bitscore"], ascending=[True, False])
        hits = best["hit"].drop_duplicates().tolist()
        decision = "unica" if len(hits) == 1 else "multiple"
        rows.append(
            {
                "query": query,
                "enterotoxina": ";".join(hits),
                "decision": decision,
                "mejor_identidad": best.iloc[0]["identidad"],
                "mejor_cobertura_query": best.iloc[0]["cobertura_query"],
                "mejor_cobertura_hit": best.iloc[0]["cobertura_hit"],
                "mejor_evalue": best.iloc[0]["evalue"],
            }
        )

    import pandas as pd

    out = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"Enterotoxinas candidatas: {len(out)}")
    print(f"Tabla escrita: {args.out}")


if __name__ == "__main__":
    main()
