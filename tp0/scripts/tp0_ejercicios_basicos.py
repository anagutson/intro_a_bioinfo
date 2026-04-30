"""Ejercicios básicos del TP0 de Bioinformática.

Este archivo muestra variables, listas, condicionales, loops, archivos y
funciones con ejemplos cortos que se pueden ejecutar desde la terminal.
"""

from pathlib import Path


RESULTADOS_DIR = Path("resultados")


def separador(titulo):
    print("\n" + "=" * 60)
    print(titulo)
    print("=" * 60)


def subtitulo(texto):
    print(f"\n--- {texto} ---")


def variables_y_print():
    separador("Ejercicio 4 - Variables, tipos, print y comentarios")

    entero = 1
    decimal = 2.0
    mi_cadena = "Hola"

    subtitulo("Variables creadas")
    print("entero =", entero, "| tipo:", type(entero))
    print("decimal =", decimal, "| tipo:", type(decimal))
    print("mi_cadena =", mi_cadena, "| tipo:", type(mi_cadena))

    subtitulo("Operaciones sencillas")
    print("entero + decimal =", entero + decimal)
    print("mi_cadena + ' mundo' =", mi_cadena + " mundo")
    print("Intento de mezclar cadena y número con +:")
    try:
        print(mi_cadena + entero)
    except TypeError as error:
        print("Python no permite sumar directamente str + int:", error)

    subtitulo("Comentarios")
    print("Las líneas que empiezan con # son comentarios: Python no las ejecuta.")
    # Esta línea es un comentario y no modifica el programa.


def listas_y_for():
    separador("Ejercicio 5 - Listas y for")

    mi_lista = [4, 1, 9, 1, "ADN"]
    subtitulo("Lista inicial")
    print("Lista original:", mi_lista)

    subtitulo("Algunas operaciones con listas")
    mi_lista.append("ARN")
    print("Luego de append:", mi_lista)
    print("Cantidad de veces que aparece 1:", mi_lista.count(1))
    print("Largo de la lista:", len(mi_lista))
    print("Elemento en posición 0:", mi_lista[0])

    subtitulo("Recorrido con for")
    print("Recorrido con for:")
    for elemento in mi_lista:
        print("elemento toma el valor:", elemento)

    print("En el for, la variable elemento toma sucesivamente cada valor de la lista.")


def condicionales_y_booleanos():
    separador("Ejercicio 5 bis - Booleanos, if, elif y else")

    b = True
    x = 2
    subtitulo("Variable booleana")
    print("b =", b)

    if b:
        print("b es true")
    else:
        print("b es false")

    subtitulo("Comparaciones")
    print("x =", x)
    print("x == 2:", x == 2)
    print("x == 3:", x == 3)
    print("x < 3:", x < 3)

    subtitulo("Condicional if / elif / else")
    if x < 0:
        print("x es negativo")
    elif x == 0:
        print("x es cero")
    else:
        print("x es positivo")


def loops_break_continue():
    separador("Loops, break y continue")

    primos = [2, 3, 5, 7]
    subtitulo("For sobre una lista")
    print("For sobre una lista de primos:")
    for numero in primos:
        print(numero)

    subtitulo("For usando range")
    print("For con range(1, 5):")
    for numero in range(1, 5):
        print(numero)

    suma = 0
    subtitulo("While")
    print("While hasta que suma sea 5:")
    while suma < 5:
        print(suma)
        suma = suma + 1

    subtitulo("continue")
    print("Ejemplo con continue: saltea el final del bloque cuando x == 2")
    for x in range(1, 5):
        print("empezando bloque", x)
        if x == 2:
            continue
        print("termine bloque", x)

    subtitulo("break")
    print("Ejemplo con break: corta el loop cuando x == 2")
    for x in range(1, 5):
        print("empezando bloque", x)
        if x == 2:
            break
        print("termine bloque", x)


def mostrar_archivo(nombre_archivo):
    with open(nombre_archivo, "r", encoding="utf-8") as archivo:
        print(archivo.read())


def archivos():
    separador("Input-output - Lectura y escritura de archivos")

    RESULTADOS_DIR.mkdir(exist_ok=True)
    nombre_archivo = RESULTADOS_DIR / "demofile.txt"

    print("\n1) Creo/sobrescribo el archivo con modo 'w'.")
    with open(nombre_archivo, "w", encoding="utf-8") as archivo:
        archivo.write("Línea inicial del archivo.\n")

    print("Contenido actual:")
    mostrar_archivo(nombre_archivo)

    print("2) Abro el mismo archivo con modo 'a' y agrego texto al final.")
    with open(nombre_archivo, "a", encoding="utf-8") as archivo:
        archivo.write("Ahora el archivo tiene más contenido agregado al final.\n")

    print("Contenido actual:")
    mostrar_archivo(nombre_archivo)

    print("3) Uso modo 'w' otra vez. Esto borra lo anterior y escribe contenido nuevo.")
    with open(nombre_archivo, "w", encoding="utf-8") as archivo:
        archivo.write("Uhhh! Borré el contenido anterior!\n")
        archivo.writelines(["línea1\n", "línea2\n"])

    print("Contenido final:")
    mostrar_archivo(nombre_archivo)


def mi_funcion(var1, var2):
    output = var1 + var2
    return output


def restar(var1, var2):
    return var1 - var2


def multiplicar(var1, var2):
    return var1 * var2


def dividir(var1, var2):
    return var1 / var2


def funciones():
    separador("Ejercicio 6 - Funciones")

    a = 10
    b = 5
    subtitulo("Función original del TP")
    print("La función mi_funcion(var1, var2) devuelve var1 + var2.")
    print("a =", a)
    print("b =", b)
    print("mi_funcion(a, b) =", mi_funcion(a, b))

    subtitulo("Misma función con cadenas")
    print('mi_funcion("bio", "informatica") =', mi_funcion("bio", "informatica"))
    print("Con strings, el operador + concatena texto.")

    subtitulo("Funciones modificadas con otras operaciones")
    print("restar(a, b) =", restar(a, b))
    print("multiplicar(a, b) =", multiplicar(a, b))
    print("dividir(a, b) =", dividir(a, b))


def main():
    variables_y_print()
    listas_y_for()
    condicionales_y_booleanos()
    loops_break_continue()
    archivos()
    funciones()


if __name__ == "__main__":
    main()
