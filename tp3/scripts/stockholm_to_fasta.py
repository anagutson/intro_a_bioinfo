"""
Convierte el alineamiento seed de HNOB descargado de Pfam
(formato Stockholm) a formato FASTA para usar en PRATT.

Uso:
    python scripts/stockholm_to_fasta.py data/HNOB_seed.sto resultados/HNOB_seed.fasta

El archivo .sto se descarga desde:
    https://www.ebi.ac.uk/interpro/entry/pfam/PF07730/ -> Alignments -> Seed -> Download
"""

import sys
from Bio import AlignIO

def convertir(entrada, salida):
    aln = AlignIO.read(entrada, "stockholm")
    AlignIO.write(aln, salida, "fasta")
    print(f"Convertido: {len(aln)} secuencias → {salida}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python stockholm_to_fasta.py <entrada.sto> <salida.fasta>")
        sys.exit(1)
    convertir(sys.argv[1], sys.argv[2])
