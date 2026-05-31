#!/usr/bin/env python3
from pathlib import Path


raw = Path("data/enterotoxinas_raw.txt")
out = Path("data/enterotoxinas.faa")

lines = raw.read_text(encoding="utf-8", errors="replace").splitlines()
fasta = []
keeping = False

for line in lines:
    stripped = line.strip()
    if not stripped:
        continue
    if stripped.startswith(">"):
        keeping = True
        fasta.append(stripped)
    elif keeping:
        seq = "".join(ch for ch in stripped if ch.isalpha() or ch == "*").upper()
        if seq:
            fasta.append(seq)

if not fasta:
    raise SystemExit(
        "No se detectaron secuencias FASTA en data/enterotoxinas_raw.txt. "
        "Revise el documento fuente y copie la base en data/enterotoxinas.faa."
    )

out.write_text("\n".join(fasta) + "\n", encoding="utf-8")
print(f"Escribi {out} con {sum(1 for line in fasta if line.startswith('>'))} secuencias.")
