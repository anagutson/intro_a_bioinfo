"""Tarea 2: generador de passwords a partir de nombre, fecha y DNI.

El TP sugiere "investigar diccionarios, listas y permutaciones". Por eso
implementamos dos enfoques:

- `generar_password_aleatorio`: arma claves al azar a partir del pool de
  caracteres y se queda con la primera que cumple los requisitos. Es el camino
  practico cuando hay muchos caracteres disponibles.
- `generar_password_permutaciones`: recorre permutaciones de 8 caracteres con
  `itertools.permutations` y devuelve la primera que cumple los requisitos.
  Es más didáctico y directo respecto a la consigna.
"""

import argparse
import random
from itertools import permutations


REQUISITOS_TEXTO = [
    "8 caracteres",
    "al menos 2 letras",
    "al menos 2 números",
    "sin fragmentos consecutivos de 2 caracteres tomados de nombre, fecha o DNI",
]


def solo_alfanumericos(texto):
    return "".join(caracter for caracter in texto if caracter.isalnum())


def pares_prohibidos(nombre, fecha, dni):
    """Devuelve todos los fragmentos de 2 caracteres presentes en los datos."""
    prohibidos = set()
    for dato in [nombre, fecha, dni]:
        limpio = solo_alfanumericos(dato).lower()
        for indice in range(len(limpio) - 1):
            prohibidos.add(limpio[indice : indice + 2])
    return prohibidos


def cumple_requisitos(password, prohibidos):
    letras = sum(caracter.isalpha() for caracter in password)
    numeros = sum(caracter.isdigit() for caracter in password)
    if len(password) != 8 or letras < 2 or numeros < 2:
        return False

    password_min = password.lower()
    for indice in range(len(password_min) - 1):
        if password_min[indice : indice + 2] in prohibidos:
            return False
    return True


def _separar_caracteres(nombre, fecha, dni):
    letras = [caracter for caracter in solo_alfanumericos(nombre) if caracter.isalpha()]
    numeros = [caracter for caracter in solo_alfanumericos(fecha + dni) if caracter.isdigit()]
    if len(letras) < 2 or len(numeros) < 2:
        raise ValueError("Se necesitan al menos 2 letras y 2 números en los datos de entrada.")
    return letras, numeros


def generar_password_aleatorio(nombre, fecha, dni, intentos=10000, semilla=None):
    """Estrategia practica: armar passwords al azar y filtrar por requisitos."""
    letras, numeros = _separar_caracteres(nombre, fecha, dni)
    caracteres = letras + numeros
    prohibidos = pares_prohibidos(nombre, fecha, dni)
    generador = random.Random(semilla) if semilla is not None else random.SystemRandom()

    for _ in range(intentos):
        base = [
            generador.choice(letras),
            generador.choice(letras),
            generador.choice(numeros),
            generador.choice(numeros),
        ]
        while len(base) < 8:
            base.append(generador.choice(caracteres))
        generador.shuffle(base)
        password = "".join(base)
        if cumple_requisitos(password, prohibidos):
            return password

    raise RuntimeError("No se pudo generar una clave válida con esos datos. Pruebe datos más variados.")


def generar_password_permutaciones(nombre, fecha, dni, semilla=7, max_pool=12):
    """Estrategia didáctica: recorrer permutaciones (consigna del TP).

    Como las permutaciones de 8 elementos son `n!/(n-8)!`, tomamos un pool
    chico (max_pool caracteres) elegido al azar para que termine en tiempo
    razonable. Devolvemos la primera permutación que cumple los requisitos.
    """
    letras, numeros = _separar_caracteres(nombre, fecha, dni)
    prohibidos = pares_prohibidos(nombre, fecha, dni)

    generador = random.Random(semilla)
    pool = sorted(set(letras + numeros))
    if len(pool) > max_pool:
        pool = generador.sample(pool, max_pool)
    generador.shuffle(pool)

    for combinacion in permutations(pool, 8):
        password = "".join(combinacion)
        if cumple_requisitos(password, prohibidos):
            return password

    raise RuntimeError("No se encontró una permutación válida. Probar con otro pool o semilla.")


# alias historico para compatibilidad con el notebook anterior
def generar_password(nombre, fecha, dni, intentos=10000):
    return generar_password_aleatorio(nombre, fecha, dni, intentos=intentos)


def main():
    parser = argparse.ArgumentParser(description="Genera passwords a partir de nombre, fecha y DNI.")
    parser.add_argument("--nombre", help="Nombre de la persona.")
    parser.add_argument("--fecha", help="Fecha de nacimiento.")
    parser.add_argument("--dni", help="DNI.")
    parser.add_argument(
        "--metodo",
        choices=["aleatorio", "permutaciones"],
        default="aleatorio",
        help="Estrategia de generacion.",
    )
    args = parser.parse_args()

    nombre = args.nombre or input("Ingrese nombre: ")
    fecha = args.fecha or input("Ingrese fecha de nacimiento (por ejemplo 20/04/2005): ")
    dni = args.dni or input("Ingrese DNI: ")

    if args.metodo == "permutaciones":
        password = generar_password_permutaciones(nombre, fecha, dni)
    else:
        password = generar_password_aleatorio(nombre, fecha, dni)

    print("Password generada:", password)
    print("Requisitos:")
    for linea in REQUISITOS_TEXTO:
        print(f"- {linea}")


if __name__ == "__main__":
    main()
