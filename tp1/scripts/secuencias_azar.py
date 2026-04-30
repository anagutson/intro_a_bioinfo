"""Generación y análisis simple de secuencias al azar."""

import argparse
from collections import Counter
import os
from pathlib import Path
import random

os.environ.setdefault("MPLCONFIGDIR", ".matplotlib-cache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


RESULTADOS_DIR = Path("resultados")
AMINOACIDOS = list("ACDEFGHIKLMNPQRSTVWY")
NUCLEOTIDOS = list("ACGT")
CODIGO_GENETICO = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def generar_adn(largo, semilla=7):
    generador = random.Random(semilla)
    return "".join(generador.choice(NUCLEOTIDOS) for _ in range(largo))


def generar_proteina(largo, frecuencias=None, semilla=7):
    rng = np.random.default_rng(semilla)
    if frecuencias is None:
        probabilidades = np.repeat(1 / len(AMINOACIDOS), len(AMINOACIDOS))
    else:
        probabilidades = np.array([frecuencias[aa] for aa in AMINOACIDOS], dtype=float)
        probabilidades = probabilidades / probabilidades.sum()
    return "".join(rng.choice(AMINOACIDOS, size=largo, p=probabilidades))


def traducir_adn(adn):
    aminoacidos = []
    for posicion in range(0, len(adn) - 2, 3):
        codon = adn[posicion : posicion + 3].upper()
        aminoacidos.append(CODIGO_GENETICO.get(codon, "X"))
    return "".join(aminoacidos)


def frecuencias_simbolos(secuencia, alfabeto=None):
    conteos = Counter(secuencia)
    if alfabeto is None:
        alfabeto = sorted(conteos)
    total = sum(conteos.get(simbolo, 0) for simbolo in alfabeto)
    datos = []
    for simbolo in alfabeto:
        cuenta = conteos.get(simbolo, 0)
        frecuencia = cuenta / total if total else 0
        datos.append((simbolo, cuenta, frecuencia))
    return pd.DataFrame(datos, columns=["símbolo", "conteo", "frecuencia"])


def graficar_frecuencias(dataframe, ruta, titulo):
    plt.figure(figsize=(9, 4))
    plt.bar(dataframe["símbolo"], dataframe["frecuencia"], color="#3b73b9")
    plt.ylabel("Frecuencia")
    plt.xlabel("Símbolo")
    plt.title(titulo)
    plt.tight_layout()
    plt.savefig(ruta, dpi=180)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Genera secuencias al azar y grafica frecuencias.")
    parser.add_argument("--largo", type=int, default=900, help="Largo de la secuencia de ADN.")
    parser.add_argument("--semilla", type=int, default=7, help="Semilla para reproducibilidad.")
    args = parser.parse_args()

    RESULTADOS_DIR.mkdir(exist_ok=True)
    adn = generar_adn(args.largo, semilla=args.semilla)
    proteina = generar_proteina(args.largo // 3, semilla=args.semilla)
    traduccion = traducir_adn(adn).replace("*", "")

    (RESULTADOS_DIR / "adn_azar.fasta").write_text(f">adn_azar_{args.largo}\n{adn}\n", encoding="utf-8")
    (RESULTADOS_DIR / "proteina_azar.fasta").write_text(
        f">proteina_azar_{len(proteina)}\n{proteina}\n", encoding="utf-8"
    )

    freq_adn = frecuencias_simbolos(adn, NUCLEOTIDOS)
    freq_prot = frecuencias_simbolos(proteina, AMINOACIDOS)
    freq_trad = frecuencias_simbolos(traduccion, AMINOACIDOS)

    freq_adn.to_csv(RESULTADOS_DIR / "frecuencias_adn_azar.csv", index=False)
    freq_prot.to_csv(RESULTADOS_DIR / "frecuencias_proteina_azar.csv", index=False)
    freq_trad.to_csv(RESULTADOS_DIR / "frecuencias_traduccion_adn_azar.csv", index=False)

    graficar_frecuencias(freq_adn, RESULTADOS_DIR / "frecuencias_adn_azar.png", "Frecuencias en ADN al azar")
    graficar_frecuencias(
        freq_prot,
        RESULTADOS_DIR / "frecuencias_proteina_azar.png",
        "Frecuencias en proteína al azar",
    )
    graficar_frecuencias(
        freq_trad,
        RESULTADOS_DIR / "frecuencias_traduccion_adn_azar.png",
        "Aminoácidos obtenidos al traducir ADN al azar",
    )

    print("ADN al azar:", adn[:60] + "...")
    print("Proteína al azar:", proteina[:60] + "...")
    print("Resultados guardados en resultados/")


if __name__ == "__main__":
    main()
