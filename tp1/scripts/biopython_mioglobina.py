"""Ejemplos de Biopython aplicados a la mioglobina humana."""

import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


DATA_DIR = Path("data")
RESULTADOS_DIR = Path("resultados")


def leer_fasta(ruta):
    return SeqIO.read(ruta, "fasta")


def resumen_record(record):
    return {
        "id": record.id,
        "name": record.name,
        "description": record.description,
        "longitud": len(record.seq),
        "secuencia_inicio": str(record.seq[:40]),
    }


def analizar_mioglobina():
    gen = leer_fasta(DATA_DIR / "mioglobina_gen.fasta")
    mrna = leer_fasta(DATA_DIR / "mioglobina_mrna.fasta")
    proteina = leer_fasta(DATA_DIR / "mioglobina_proteina.fasta")

    directa = gen.seq.translate(to_stop=True)
    reversa = gen.seq.reverse_complement().translate(to_stop=True)

    RESULTADOS_DIR.mkdir(exist_ok=True)
    with open(RESULTADOS_DIR / "mioglobina_biopython_resumen.txt", "w", encoding="utf-8") as archivo:
        archivo.write("Resumen de registros FASTA\n")
        for nombre, record in [("gen", gen), ("mRNA", mrna), ("proteína", proteina)]:
            archivo.write(f"\n{nombre}\n")
            for clave, valor in resumen_record(record).items():
                archivo.write(f"{clave}: {valor}\n")
        archivo.write(f"\nGC del gen: {gc_fraction(gen.seq) * 100:.2f}%\n")
        archivo.write(f"Traducción directa inicia: {directa[:80]}\n")
        archivo.write(f"Traducción reversa inicia: {reversa[:80]}\n")
        archivo.write(f"Proteína real inicia: {proteina.seq[:80]}\n")
        archivo.write(f"Directa coincide con proteína real: {str(directa) == str(proteina.seq)}\n")

    return gen, mrna, proteina, directa, reversa


def main():
    parser = argparse.ArgumentParser(description="Analiza la mioglobina humana con Biopython.")
    parser.parse_args()

    gen, mrna, proteina, directa, reversa = analizar_mioglobina()
    print("Gen:", resumen_record(gen))
    print("mRNA:", resumen_record(mrna))
    print("Proteína:", resumen_record(proteina))
    print(f"GC del gen: {gc_fraction(gen.seq) * 100:.2f}%")
    print("Traducción directa coincide con proteína real:", str(directa) == str(proteina.seq))
    print("Resumen guardado en resultados/mioglobina_biopython_resumen.txt")


if __name__ == "__main__":
    main()
