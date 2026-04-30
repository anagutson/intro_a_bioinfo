"""Tarea 1: números primos y comparación de performance."""

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
    """Devuelve True si número es primo, probando todos los divisores menores."""
    if numero < 2:
        return False

    for divisor in range(2, numero):
        if numero % divisor == 0:
            return False
    return True


def es_primo_mejorado(numero):
    """Devuelve True si número es primo, probando divisores hasta la raíz."""
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
        columns=["límite", "tiempo_básico_segundos", "tiempo_mejorado_segundos"],
    )
    dataframe.to_csv(ruta, index=False)


def guardar_grafico_png(mediciones, ruta):
    dataframe = pd.DataFrame(
        mediciones,
        columns=["Límite", "Método básico", "Método mejorado"],
    )

    figura, (eje_lineal, eje_log) = plt.subplots(1, 2, figsize=(12, 4.5))
    for eje, escala in zip((eje_lineal, eje_log), ("lineal", "log")):
        eje.plot(dataframe["Límite"], dataframe["Método básico"], marker="o", label="Método básico")
        eje.plot(dataframe["Límite"], dataframe["Método mejorado"], marker="o", label="Método mejorado")
        eje.set_xlabel("Límite superior")
        eje.set_ylabel("Tiempo (segundos)")
        eje.set_title(f"Tiempo para calcular primos ({escala})")
        eje.grid(alpha=0.3)
        eje.legend()
    eje_log.set_yscale("log")

    figura.tight_layout()
    figura.savefig(ruta, dpi=180)
    plt.close(figura)


def main():
    parser = argparse.ArgumentParser(description="Calcula números primos y compara performance.")
    parser.add_argument("--limite", type=int, help="Límite superior definido por el usuario.")
    args = parser.parse_args()

    print("Primos del 1 al 20:")
    print(primos_hasta(20, es_primo_basico))

    if args.limite is None:
        entrada = input("Ingrese un límite para buscar primos: ")
        limite_usuario = int(entrada)
    else:
        limite_usuario = args.limite

    print(f"Primos del 1 al {limite_usuario}:")
    print(primos_hasta(limite_usuario, es_primo_mejorado))

    RESULTADOS_DIR.mkdir(exist_ok=True)
    limites = [100, 500, 1000, 2000, 5000, 10000]
    mediciones = medir_tiempos(limites)
    guardar_csv(mediciones, RESULTADOS_DIR / "tarea1_tiempos_primos.csv")
    guardar_grafico_png(mediciones, RESULTADOS_DIR / "tarea1_tiempos_primos.png")

    print("Mediciones de tiempo:")
    for limite, tiempo_basico, tiempo_mejorado in mediciones:
        print(f"límite={limite:5d} básico={tiempo_basico:.6f}s mejorado={tiempo_mejorado:.6f}s")
    print("Se guardaron los resultados en la carpeta resultados/")


if __name__ == "__main__":
    main()
