"""Crea un notebook prolijo de entrega a partir de los scripts del TP0."""

from pathlib import Path
from textwrap import dedent

import nbformat as nbf


def code(source):
    return nbf.v4.new_code_cell(dedent(source).strip())


def markdown(source):
    return nbf.v4.new_markdown_cell(dedent(source).strip())


nb = nbf.v4.new_notebook()
nb["metadata"] = {
    "kernelspec": {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    },
    "language_info": {
        "name": "python",
        "pygments_lexer": "ipython3",
    },
}

nb["cells"] = [
    markdown(
        """
        # TP0 - Entornos y primeros programas

        **Bioinformatica 2026**

        Este notebook integra ejercicios basicos de Python y las tres tareas de programacion del TP0.
        """
    ),
    markdown("## Ejercicio 1 - Primer programa"),
    code('print("Este es mi primer programa de Python")'),
    markdown(
        """
        ## Ejercicio 4 - Variables, tipos, print y comentarios

        Se definen tres objetos: un entero, un decimal/flotante y una cadena de texto.
        Luego se prueban `print`, operaciones sencillas y comentarios.
        """
    ),
    code(
        """
        entero = 1
        decimal = 2.0
        mi_cadena = "Hola"

        print("--- Variables creadas ---")
        print("entero =", entero, "| tipo:", type(entero))
        print("decimal =", decimal, "| tipo:", type(decimal))
        print("mi_cadena =", mi_cadena, "| tipo:", type(mi_cadena))

        print("\\n--- Operaciones sencillas ---")
        print("entero + decimal =", entero + decimal)
        print("mi_cadena + ' mundo' =", mi_cadena + " mundo")

        print("\\nIntento de mezclar cadena y numero con +:")
        try:
            print(mi_cadena + entero)
        except TypeError as error:
            print("Python no permite sumar directamente str + int:", error)

        print("\\n--- Comentarios ---")
        print("Las lineas que empiezan con # son comentarios: Python no las ejecuta.")
        # Esta linea es un comentario y no modifica el programa.
        """
    ),
    markdown(
        """
        ## Ejercicio 5 - Listas y for

        Se crea una lista con numeros y cadenas, se prueban operaciones simples y se recorre
        con un `for`. En cada iteracion, la variable del loop toma el valor de un elemento.
        """
    ),
    code(
        """
        mi_lista = [4, 1, 9, 1, "ADN"]

        print("--- Lista inicial ---")
        print("Lista original:", mi_lista)

        print("\\n--- Algunas operaciones con listas ---")
        mi_lista.append("ARN")
        print("Luego de append:", mi_lista)
        print("Cantidad de veces que aparece 1:", mi_lista.count(1))
        print("Largo de la lista:", len(mi_lista))
        print("Elemento en posicion 0:", mi_lista[0])

        print("\\n--- Recorrido con for ---")
        for elemento in mi_lista:
            print("elemento toma el valor:", elemento)

        print("\\nEn el for, la variable elemento toma sucesivamente cada valor de la lista.")
        """
    ),
    markdown(
        """
        ## Booleanos y condicionales

        Se prueban variables booleanas, comparaciones y la estructura `if / elif / else`.
        """
    ),
    code(
        """
        b = True
        x = 2

        print("--- Variable booleana ---")
        print("b =", b)

        if b:
            print("b es true")
        else:
            print("b es false")

        print("\\n--- Comparaciones ---")
        print("x =", x)
        print("x == 2:", x == 2)
        print("x == 3:", x == 3)
        print("x < 3:", x < 3)

        print("\\n--- Condicional if / elif / else ---")
        if x < 0:
            print("x es negativo")
        elif x == 0:
            print("x es cero")
        else:
            print("x es positivo")
        """
    ),
    markdown(
        """
        ## Loops: for, while, continue y break

        Se prueban bucles `for`, un bucle `while` y la diferencia entre `continue` y `break`.
        """
    ),
    code(
        """
        primos = [2, 3, 5, 7]

        print("--- For sobre una lista ---")
        for numero in primos:
            print(numero)

        print("\\n--- For usando range(1, 5) ---")
        for numero in range(1, 5):
            print(numero)

        print("\\n--- While ---")
        suma = 0
        while suma < 5:
            print(suma)
            suma = suma + 1

        print("\\n--- continue ---")
        print("continue saltea el final del bloque cuando x == 2")
        for x in range(1, 5):
            print("empezando bloque", x)
            if x == 2:
                continue
            print("termine bloque", x)

        print("\\n--- break ---")
        print("break corta el loop cuando x == 2")
        for x in range(1, 5):
            print("empezando bloque", x)
            if x == 2:
                break
            print("termine bloque", x)
        """
    ),
    markdown(
        """
        ## Input-output - Lectura y escritura de archivos

        Se muestra la diferencia entre abrir un archivo con modo `"w"` y modo `"a"`.
        El modo `"w"` crea o sobrescribe; el modo `"a"` agrega contenido al final.
        """
    ),
    code(
        """
        from pathlib import Path

        proyecto = Path.cwd()
        while not (proyecto / "requirements.txt").exists() and proyecto != proyecto.parent:
            proyecto = proyecto.parent

        resultados_dir = proyecto / "resultados"
        resultados_dir.mkdir(exist_ok=True)
        nombre_archivo = resultados_dir / "demofile.txt"


        def mostrar_archivo(ruta):
            with open(ruta, "r", encoding="utf-8") as archivo:
                print(archivo.read())


        print("1) Creo/sobrescribo el archivo con modo 'w'.")
        with open(nombre_archivo, "w", encoding="utf-8") as archivo:
            archivo.write("Linea inicial del archivo.\\n")

        print("Contenido actual:")
        mostrar_archivo(nombre_archivo)

        print("2) Abro el mismo archivo con modo 'a' y agrego texto al final.")
        with open(nombre_archivo, "a", encoding="utf-8") as archivo:
            archivo.write("Ahora el archivo tiene mas contenido agregado al final.\\n")

        print("Contenido actual:")
        mostrar_archivo(nombre_archivo)

        print("3) Uso modo 'w' otra vez. Esto borra lo anterior y escribe contenido nuevo.")
        with open(nombre_archivo, "w", encoding="utf-8") as archivo:
            archivo.write("Uhhh! Borre el contenido anterior!\\n")
            archivo.writelines(["linea1\\n", "linea2\\n"])

        print("Contenido final:")
        mostrar_archivo(nombre_archivo)
        """
    ),
    markdown(
        """
        ## Ejercicio 6 - Funciones

        Se define una funcion que recibe dos valores, los procesa y devuelve un resultado.
        Luego se prueba cambiando datos y operaciones.
        """
    ),
    code(
        """
        def mi_funcion(var1, var2):
            output = var1 + var2
            return output


        def restar(var1, var2):
            return var1 - var2


        def multiplicar(var1, var2):
            return var1 * var2


        def dividir(var1, var2):
            return var1 / var2


        a = 10
        b = 5

        print("--- Funcion original del TP ---")
        print("La funcion mi_funcion(var1, var2) devuelve var1 + var2.")
        print("a =", a)
        print("b =", b)
        print("mi_funcion(a, b) =", mi_funcion(a, b))

        print("\\n--- Misma funcion con cadenas ---")
        print('mi_funcion("bio", "informatica") =', mi_funcion("bio", "informatica"))
        print("Con strings, el operador + concatena texto.")

        print("\\n--- Funciones modificadas con otras operaciones ---")
        print("restar(a, b) =", restar(a, b))
        print("multiplicar(a, b) =", multiplicar(a, b))
        print("dividir(a, b) =", dividir(a, b))
        """
    ),
    markdown(
        """
        ## Tarea 1 - Numeros primos

        Un numero primo es un numero entero mayor que 1 que solo puede dividirse exactamente
        por 1 y por si mismo. Para detectar si un numero es primo, alcanza con buscar si existe
        algun divisor intermedio.

        Pseudocodigo de la version basica:

        ```text
        para cada numero entre 1 y limite:
            si numero < 2:
                no es primo
            si no:
                asumir que es primo
                para cada divisor entre 2 y numero - 1:
                    si numero % divisor == 0:
                        no es primo
                si ningun divisor lo dividio:
                    guardar/imprimir el numero como primo
        ```

        El operador `%` calcula el resto de una division. Si `numero % divisor == 0`,
        significa que el divisor entra justo y entonces el numero no es primo.

        Pseudocodigo de la version mejorada:

        ```text
        para cada numero entre 1 y limite:
            si numero < 2:
                no es primo
            si numero == 2:
                es primo
            si numero es par y mayor que 2:
                no es primo
            si no:
                probar divisores impares: 3, 5, 7, ...
                detenerse cuando divisor * divisor > numero
                si ningun divisor lo dividio:
                    es primo
        ```

        Por eso, para saber si `n` es primo, alcanza con probar divisores hasta `raiz(n)`.
        Si no aparece ningun divisor hasta ese punto, el numero es primo.
        """
    ),
    code(
        """
        from time import perf_counter

        import matplotlib.pyplot as plt
        import pandas as pd


        def es_primo_basico(numero):
            if numero < 2:
                return False
            for divisor in range(2, numero):
                if numero % divisor == 0:
                    return False
            return True


        def es_primo_mejorado(numero):
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
            return [numero for numero in range(1, limite + 1) if metodo(numero)]


        print(primos_hasta(20, es_primo_basico))
        """
    ),
    code(
        """
        limites = [100, 500, 1000, 2000, 5000, 10000]
        mediciones = []

        for limite in limites:
            inicio = perf_counter()
            primos_hasta(limite, es_primo_basico)
            tiempo_basico = perf_counter() - inicio

            inicio = perf_counter()
            primos_hasta(limite, es_primo_mejorado)
            tiempo_mejorado = perf_counter() - inicio

            mediciones.append((limite, tiempo_basico, tiempo_mejorado))

        tiempos = pd.DataFrame(mediciones, columns=["limite", "basico", "mejorado"])
        tiempos
        """
    ),
    code(
        """
        tiempos.plot(x="limite", y=["basico", "mejorado"], marker="o", figsize=(8, 4), grid=True)
        plt.title("Tiempo para calcular numeros primos")
        plt.xlabel("Limite superior")
        plt.ylabel("Tiempo (segundos)")
        plt.tight_layout()
        plt.show()
        """
    ),
    markdown("## Tarea 2 - Passwords"),
    code(
        """
        import sys
        from pathlib import Path

        proyecto = Path.cwd()
        while not (proyecto / "requirements.txt").exists() and proyecto != proyecto.parent:
            proyecto = proyecto.parent

        sys.path.append(str(proyecto / "scripts"))

        from tarea2_passwords import generar_password, cumple_requisitos, pares_prohibidos

        nombre = "Elena Gilbert"
        fecha = "20/04/2005"
        dni = "40123456"
        password = generar_password(nombre, fecha, dni)

        print(password)
        print(cumple_requisitos(password, pares_prohibidos(nombre, fecha, dni)))
        """
    ),
    markdown(
        """
        ## Tarea 3 - Conteo y nube de palabras

        Para que la nube sea informativa se excluyen palabras muy frecuentes como articulos,
        preposiciones y conectores.
        """
    ),
    code(
        """
        import sys
        from pathlib import Path

        proyecto = Path.cwd()
        while not (proyecto / "requirements.txt").exists() and proyecto != proyecto.parent:
            proyecto = proyecto.parent

        sys.path.append(str(proyecto / "scripts"))

        import pandas as pd
        from tarea3_palabras import frecuencias, leer_texto, contar_palabra

        texto = leer_texto(proyecto / "data" / "texto_prueba.txt")
        print("bioinformatica aparece:", contar_palabra(texto, "bioinformatica"))

        conteos = frecuencias(texto)
        pd.DataFrame(conteos.most_common(15), columns=["palabra", "frecuencia"])
        """
    ),
    code(
        """
        import matplotlib.pyplot as plt
        from wordcloud import WordCloud

        nube = WordCloud(
            width=900,
            height=550,
            background_color="white",
            random_state=7,
        ).generate_from_frequencies(conteos)

        plt.figure(figsize=(9, 5.5))
        plt.imshow(nube, interpolation="bilinear")
        plt.axis("off")
        plt.show()
        """
    ),
]

Path("notebooks").mkdir(exist_ok=True)

with open("notebooks/TP0_Bioinformatica_Entrega.ipynb", "w", encoding="utf-8") as archivo:
    nbf.write(nb, archivo)

print("Notebook creado: notebooks/TP0_Bioinformatica_Entrega.ipynb")
