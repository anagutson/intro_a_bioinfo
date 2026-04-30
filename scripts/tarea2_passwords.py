"""Tarea 2: generador de passwords a partir de nombre, fecha y DNI."""

import argparse
import random


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


def generar_password(nombre, fecha, dni, intentos=10000):
    letras = list(solo_alfanumericos(nombre))
    numeros = list(solo_alfanumericos(fecha + dni))
    numeros = [caracter for caracter in numeros if caracter.isdigit()]
    caracteres = letras + numeros
    prohibidos = pares_prohibidos(nombre, fecha, dni)
    generador = random.SystemRandom()

    if len(letras) < 2 or len(numeros) < 2:
        raise ValueError("Se necesitan al menos 2 letras y 2 numeros en los datos de entrada.")

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

    raise RuntimeError("No se pudo generar una clave valida con esos datos. Pruebe datos mas variados.")


def main():
    parser = argparse.ArgumentParser(description="Genera passwords a partir de nombre, fecha y DNI.")
    parser.add_argument("--nombre", help="Nombre de la persona.")
    parser.add_argument("--fecha", help="Fecha de nacimiento.")
    parser.add_argument("--dni", help="DNI.")
    args = parser.parse_args()

    nombre = args.nombre or input("Ingrese nombre: ")
    fecha = args.fecha or input("Ingrese fecha de nacimiento (por ejemplo 20/04/2005): ")
    dni = args.dni or input("Ingrese DNI: ")

    password = generar_password(nombre, fecha, dni)
    print("Password generada:", password)
    print("Requisitos:")
    print("- 8 caracteres")
    print("- al menos 2 letras")
    print("- al menos 2 numeros")
    print("- sin fragmentos consecutivos de 2 caracteres tomados de nombre, fecha o DNI")


if __name__ == "__main__":
    main()
