"""Descargas opcionales con Biopython.

El notebook usa datos locales para ser reproducible. Este script queda como
demostración de cómo se descargarían registros reales usando Entrez y ExPASy.
"""

import argparse
from pathlib import Path

from Bio import Entrez, ExPASy, SeqIO, SwissProt


RESULTADOS_DIR = Path("resultados")


def descargar_refseq_fasta(accession, email):
    Entrez.email = email
    with Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text") as handle:
        return handle.read()


def descargar_uniprot_swissprot(uniprot_id):
    handle = ExPASy.get_sprot_raw(uniprot_id)
    try:
        return SwissProt.read(handle)
    finally:
        handle.close()


def main():
    parser = argparse.ArgumentParser(description="Descarga registros de mioglobina con Biopython.")
    parser.add_argument("--email", required=True, help="Email requerido por NCBI Entrez.")
    parser.add_argument("--refseq", default="NP_976311.1", help="Accession RefSeq.")
    parser.add_argument("--uniprot", default="P02144", help="Accession UniProt.")
    args = parser.parse_args()

    RESULTADOS_DIR.mkdir(exist_ok=True)
    fasta = descargar_refseq_fasta(args.refseq, args.email)
    fasta_path = RESULTADOS_DIR / f"{args.refseq}.fasta"
    fasta_path.write_text(fasta, encoding="utf-8")
    record = SeqIO.read(fasta_path, "fasta")

    swiss = descargar_uniprot_swissprot(args.uniprot)
    print("RefSeq:", record.id, len(record.seq), "aa")
    print("UniProt:", swiss.entry_name, len(swiss.sequence), "aa")
    print("FASTA guardado en", fasta_path)


if __name__ == "__main__":
    main()
