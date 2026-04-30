"""Tarea 3: conteo de palabras, grafico de frecuencias y nube de palabras."""

import argparse
from collections import Counter
import math
import os
from pathlib import Path
import re

os.environ.setdefault("MPLCONFIGDIR", ".matplotlib-cache")

import pandas as pd


RESULTADOS_DIR = Path("resultados")
DATA_DIR = Path("data")


STOPWORDS = {
    "a",
    "al",
    "con",
    "de",
    "del",
    "el",
    "en",
    "es",
    "la",
    "las",
    "lo",
    "los",
    "para",
    "por",
    "que",
    "se",
    "un",
    "una",
    "y",
}


def resolver_ruta_archivo(ruta):
    ruta = Path(ruta)
    if ruta.exists():
        return ruta

    ruta_en_data = DATA_DIR / ruta
    if ruta_en_data.exists():
        return ruta_en_data

    raise FileNotFoundError(
        f"No encontre '{ruta}'. Si queres usar el texto de prueba, ingresa 'data/texto_prueba.txt' "
        "o simplemente 'texto_prueba.txt'."
    )


def leer_texto(ruta):
    ruta = resolver_ruta_archivo(ruta)
    with open(ruta, "r", encoding="utf-8") as archivo:
        return archivo.read()


def tokenizar(texto):
    return re.findall(r"[A-Za-zÁÉÍÓÚÜÑáéíóúüñ]+", texto.lower())


def contar_palabra(texto, palabra):
    palabra = palabra.lower()
    return tokenizar(texto).count(palabra)


def frecuencias(texto, excluir_stopwords=True):
    palabras = tokenizar(texto)
    if excluir_stopwords:
        palabras = [palabra for palabra in palabras if palabra not in STOPWORDS]
    return Counter(palabras)


def guardar_frecuencias_csv(conteos, ruta):
    dataframe = pd.DataFrame(conteos.most_common(), columns=["palabra", "frecuencia"])
    dataframe.to_csv(ruta, index=False)


def guardar_nube_svg(conteos, ruta, max_palabras=40):
    palabras = conteos.most_common(max_palabras)
    ancho = 900
    alto = 620
    centro_x = ancho / 2
    centro_y = alto / 2
    max_frecuencia = max((frecuencia for _, frecuencia in palabras), default=1)
    colores = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#386cb0"]

    elementos = []
    for indice, (palabra, frecuencia) in enumerate(palabras):
        angulo = indice * 2.4
        radio = 18 + indice * 9
        x = centro_x + math.cos(angulo) * radio
        y = centro_y + math.sin(angulo) * radio
        tamanio = 14 + int(42 * frecuencia / max_frecuencia)
        color = colores[indice % len(colores)]
        elementos.append(
            f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="middle" '
            f'font-family="Arial" font-size="{tamanio}" fill="{color}">{palabra}</text>'
        )

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{ancho}" height="{alto}" viewBox="0 0 {ancho} {alto}">
  <rect width="100%" height="100%" fill="white"/>
  <text x="{centro_x}" y="32" text-anchor="middle" font-family="Arial" font-size="20">Nube de palabras</text>
  {"".join(elementos)}
</svg>
"""
    with open(ruta, "w", encoding="utf-8") as archivo:
        archivo.write(svg)


def guardar_nube_wordcloud(conteos, ruta):
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from wordcloud import WordCloud

    nube = WordCloud(
        width=1200,
        height=800,
        background_color="white",
        colormap="viridis",
        max_words=80,
        random_state=7,
    ).generate_from_frequencies(conteos)

    plt.figure(figsize=(10, 6.5))
    plt.imshow(nube, interpolation="bilinear")
    plt.axis("off")
    plt.tight_layout(pad=0)
    plt.savefig(ruta, dpi=180)
    plt.close()


def guardar_barras(conteos, ruta, cantidad=15):
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    dataframe = pd.DataFrame(conteos.most_common(cantidad), columns=["palabra", "frecuencia"])
    plt.figure(figsize=(9, 5))
    plt.barh(dataframe["palabra"], dataframe["frecuencia"], color="#386cb0")
    plt.gca().invert_yaxis()
    plt.title("Palabras mas frecuentes")
    plt.xlabel("Frecuencia")
    plt.tight_layout()
    plt.savefig(ruta, dpi=180)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Cuenta palabras y genera una nube de palabras.")
    parser.add_argument("--archivo", help="Ruta del archivo de texto.")
    parser.add_argument("--palabra", help="Palabra a buscar.")
    args = parser.parse_args()

    ruta = args.archivo or input("Ingrese la ruta del archivo de texto: ")
    palabra = args.palabra or input("Ingrese la palabra a buscar: ")

    texto = leer_texto(ruta)
    apariciones = contar_palabra(texto, palabra)
    conteos = frecuencias(texto)

    print(f"La palabra '{palabra}' aparece {apariciones} veces.")
    print("Frecuencias mas comunes:")
    for palabra_frecuente, frecuencia in conteos.most_common(15):
        print(f"{palabra_frecuente}: {frecuencia}")

    RESULTADOS_DIR.mkdir(exist_ok=True)
    guardar_frecuencias_csv(conteos, RESULTADOS_DIR / "tarea3_frecuencias.csv")
    guardar_nube_svg(conteos, RESULTADOS_DIR / "tarea3_nube_palabras.svg")
    guardar_nube_wordcloud(conteos, RESULTADOS_DIR / "tarea3_nube_wordcloud.png")
    guardar_barras(conteos, RESULTADOS_DIR / "tarea3_frecuencias_top15.png")
    print("Se guardaron los resultados en la carpeta resultados/")
    print("Excluir palabras muy frecuentes como articulos/preposiciones ayuda a ver mejor los temas del texto.")


if __name__ == "__main__":
    main()
