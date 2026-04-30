"""Tarea 1: numeros primos y comparacion de performance."""

import argparse
import os
from pathlib import Path
from time import perf_counter

os.environ.setdefault("MPLCONFIGDIR", ".matplotlib-cache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


RESULTADOS_DIR = Path("resultados")


def es_primo_basico(numero):
    """Devuelve True si numero es primo, probando todos los divisores menores."""
    if numero < 2:
        return False

    for divisor in range(2, numero):
        if numero % divisor == 0:
            return False
    return True


def es_primo_mejorado(numero):
    """Devuelve True si numero es primo, probando divisores hasta la raiz."""
    if numero < 2:
        return False
    if numero == 2:
        return True
    if numero % 2 == 0:
        return False

    divisor = 3
    while divisor * divisor <= numero:
        if numero % divisor == 0:
            return False
        divisor += 2
    return True


def primos_hasta(limite, metodo):
    primos = []
    for numero in range(1, limite + 1):
        if metodo(numero):
            primos.append(numero)
    return primos


def medir_tiempos(limites):
    mediciones = []
    for limite in limites:
        inicio = perf_counter()
        primos_hasta(limite, es_primo_basico)
        tiempo_basico = perf_counter() - inicio

        inicio = perf_counter()
        primos_hasta(limite, es_primo_mejorado)
        tiempo_mejorado = perf_counter() - inicio

        mediciones.append((limite, tiempo_basico, tiempo_mejorado))
    return mediciones


def guardar_csv(mediciones, ruta):
    dataframe = pd.DataFrame(
        mediciones,
        columns=["limite", "tiempo_basico_segundos", "tiempo_mejorado_segundos"],
    )
    dataframe.to_csv(ruta, index=False)


def puntos_svg(mediciones, columna, ancho, alto, margen, max_limite, max_tiempo):
    puntos = []
    for limite, tiempo_basico, tiempo_mejorado in mediciones:
        tiempo = tiempo_basico if columna == "basico" else tiempo_mejorado
        x = margen + (limite / max_limite) * (ancho - 2 * margen)
        y = alto - margen - (tiempo / max_tiempo) * (alto - 2 * margen)
        puntos.append(f"{x:.2f},{y:.2f}")
    return " ".join(puntos)


def guardar_svg(mediciones, ruta):
    ancho = 760
    alto = 420
    margen = 55
    max_limite = max(limite for limite, _, _ in mediciones)
    max_tiempo = max(max(tiempo_basico, tiempo_mejorado) for _, tiempo_basico, tiempo_mejorado in mediciones)

    linea_basica = puntos_svg(mediciones, "basico", ancho, alto, margen, max_limite, max_tiempo)
    linea_mejorada = puntos_svg(mediciones, "mejorado", ancho, alto, margen, max_limite, max_tiempo)

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{ancho}" height="{alto}" viewBox="0 0 {ancho} {alto}">
  <rect width="100%" height="100%" fill="white"/>
  <text x="{ancho / 2}" y="26" text-anchor="middle" font-family="Arial" font-size="18">Tiempo para calcular primos</text>
  <line x1="{margen}" y1="{alto - margen}" x2="{ancho - margen}" y2="{alto - margen}" stroke="black"/>
  <line x1="{margen}" y1="{margen}" x2="{margen}" y2="{alto - margen}" stroke="black"/>
  <text x="{ancho / 2}" y="{alto - 12}" text-anchor="middle" font-family="Arial" font-size="13">Limite superior</text>
  <text x="18" y="{alto / 2}" transform="rotate(-90 18,{alto / 2})" text-anchor="middle" font-family="Arial" font-size="13">Tiempo (segundos)</text>
  <polyline points="{linea_basica}" fill="none" stroke="#d43f3a" stroke-width="3"/>
  <polyline points="{linea_mejorada}" fill="none" stroke="#2166ac" stroke-width="3"/>
  <text x="{ancho - 220}" y="70" font-family="Arial" font-size="14" fill="#d43f3a">Metodo basico</text>
  <text x="{ancho - 220}" y="94" font-family="Arial" font-size="14" fill="#2166ac">Metodo mejorado</text>
  <text x="{margen}" y="{alto - 36}" text-anchor="middle" font-family="Arial" font-size="11">0</text>
  <text x="{ancho - margen}" y="{alto - 36}" text-anchor="middle" font-family="Arial" font-size="11">{max_limite}</text>
  <text x="{margen - 8}" y="{margen + 4}" text-anchor="end" font-family="Arial" font-size="11">{max_tiempo:.4f}</text>
</svg>
"""
    with open(ruta, "w", encoding="utf-8") as archivo:
        archivo.write(svg)


def guardar_grafico_png(mediciones, ruta):
    dataframe = pd.DataFrame(
        mediciones,
        columns=["Limite", "Metodo basico", "Metodo mejorado"],
    )

    plt.figure(figsize=(9, 5))
    plt.plot(dataframe["Limite"], dataframe["Metodo basico"], marker="o", label="Metodo basico")
    plt.plot(dataframe["Limite"], dataframe["Metodo mejorado"], marker="o", label="Metodo mejorado")
    plt.title("Tiempo para calcular numeros primos")
    plt.xlabel("Limite superior")
    plt.ylabel("Tiempo (segundos)")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(ruta, dpi=180)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Calcula numeros primos y compara performance.")
    parser.add_argument("--limite", type=int, help="Limite superior definido por el usuario.")
    args = parser.parse_args()

    print("Primos del 1 al 20:")
    print(primos_hasta(20, es_primo_basico))

    if args.limite is None:
        entrada = input("Ingrese un limite para buscar primos: ")
        limite_usuario = int(entrada)
    else:
        limite_usuario = args.limite

    print(f"Primos del 1 al {limite_usuario}:")
    print(primos_hasta(limite_usuario, es_primo_mejorado))

    RESULTADOS_DIR.mkdir(exist_ok=True)
    limites = [100, 500, 1000, 2000, 5000, 10000]
    mediciones = medir_tiempos(limites)
    guardar_csv(mediciones, RESULTADOS_DIR / "tarea1_tiempos_primos.csv")
    guardar_svg(mediciones, RESULTADOS_DIR / "tarea1_tiempos_primos.svg")
    guardar_grafico_png(mediciones, RESULTADOS_DIR / "tarea1_tiempos_primos.png")

    print("Mediciones de tiempo:")
    for limite, tiempo_basico, tiempo_mejorado in mediciones:
        print(f"limite={limite:5d} basico={tiempo_basico:.6f}s mejorado={tiempo_mejorado:.6f}s")
    print("Se guardaron los resultados en la carpeta resultados/")


if __name__ == "__main__":
    main()
