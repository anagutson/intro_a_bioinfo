"""Genera el notebook de entrega del TP0.

Notas para entrega:
- Este script *no* forma parte de la entrega: solo arma el `.ipynb`.
- El notebook resultante es autocontenido: la docente puede correrlo sin la
  carpeta `scripts/`. La unica dependencia externa es `data/texto_prueba.txt`,
  que tambien existe inline como respaldo.
"""

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

        Este notebook es autocontenido: cada celda incluye el codigo que necesita
        para correr, sin depender de archivos externos del repositorio.
        """
    ),
    markdown("## Ejercicio 1 - Primer programa"),
    code('print("Este es mi primer programa de Python")'),
    markdown(
        """
        ## Ejercicio 4 - Variables, tipos, print y comentarios

        Se definen tres objetos pedidos por el TP (entero, decimal y cadena), se
        prueban con `print`, se mezclan tipos en operaciones y se introduce el
        uso de comentarios. Tambien se llama a `help` sobre los objetos para
        explorar la documentacion incorporada.
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
        print("mi_cadena * 3 =", mi_cadena * 3)

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
        ### Interrogar objetos con `help`

        El TP pide usar `help(objeto)`. La salida es larga, asi que mostramos
        solo las primeras lineas para que se vea de que se trata.
        """
    ),
    code(
        """
        import io
        import contextlib


        def help_resumen(objeto, lineas=12):
            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                help(objeto)
            salida = buffer.getvalue().splitlines()
            for linea in salida[:lineas]:
                print(linea)
            if len(salida) > lineas:
                print(f"... ({len(salida) - lineas} lineas mas)")


        print("=== help(entero) ===")
        help_resumen(entero)
        print("\\n=== help(mi_cadena.upper) ===")
        help_resumen(mi_cadena.upper)
        """
    ),
    markdown(
        """
        ## Ejercicio 5 - Listas y for

        Se crea una lista mixta (numeros y cadenas), se prueban operaciones
        comunes y se recorre con `for`. En cada iteracion la variable del loop
        toma el valor del elemento actual.
        """
    ),
    code(
        """
        mi_lista = [4, 1, 9, 1, "ADN"]

        print("--- Lista inicial ---")
        print("Lista original:", mi_lista)

        print("\\n--- Algunas operaciones con listas ---")
        mi_lista.append("ARN")
        print("Luego de append('ARN'):", mi_lista)
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

        Se prueban variables booleanas, comparaciones y la estructura
        `if / elif / else`.
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

        Se prueban bucles `for`, `while`, y la diferencia entre `continue`
        (saltea el resto de la iteracion) y `break` (corta el loop entero).
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

        Se muestra la diferencia entre abrir un archivo con modo `"w"` y
        modo `"a"`. El modo `"w"` crea o sobrescribe el archivo, mientras que
        `"a"` agrega contenido al final.

        Para que el notebook sea autocontenido, el archivo se crea en una
        carpeta temporal.
        """
    ),
    code(
        """
        import tempfile
        from pathlib import Path

        carpeta = Path(tempfile.mkdtemp(prefix="tp0_demo_"))
        nombre_archivo = carpeta / "demofile.txt"


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

        Se define una funcion que recibe dos valores, los procesa y devuelve un
        resultado. Luego se prueba cambiando datos (numeros y strings) y
        operaciones (suma, resta, multiplicacion, division).
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

        Un numero primo es un numero entero mayor que 1 que solo puede dividirse
        exactamente por 1 y por si mismo. Para detectar si un numero es primo,
        alcanza con buscar si existe algun divisor intermedio.

        **Pseudocodigo de la version basica:**

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

        El operador `%` calcula el resto de una division. Si
        `numero % divisor == 0`, el divisor entra justo y entonces el numero no
        es primo.

        **Pseudocodigo de la version mejorada:**

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

        Para saber si `n` es primo alcanza con probar divisores hasta `raiz(n)`.
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


        print("Primos del 1 al 20:")
        print(primos_hasta(20, es_primo_basico))
        """
    ),
    markdown(
        """
        ### 5a - Medicion de tiempos

        Comparamos el metodo basico contra el mejorado para distintos limites.
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
    markdown(
        """
        ### 5b - Comparacion grafica

        Mostramos los mismos datos en escala lineal y en escala logaritmica.
        En escala log se ve mejor la diferencia de complejidad: el metodo
        basico es aproximadamente $O(n^2)$ y el mejorado es aproximadamente
        $O(n \\sqrt{n})$.
        """
    ),
    code(
        """
        figura, (eje_lineal, eje_log) = plt.subplots(1, 2, figsize=(12, 4.5))
        for eje, escala in zip((eje_lineal, eje_log), ("lineal", "log")):
            eje.plot(tiempos["limite"], tiempos["basico"], marker="o", label="basico")
            eje.plot(tiempos["limite"], tiempos["mejorado"], marker="o", label="mejorado")
            eje.set_xlabel("Limite superior")
            eje.set_ylabel("Tiempo (segundos)")
            eje.set_title(f"Tiempo para calcular primos ({escala})")
            eje.grid(alpha=0.3)
            eje.legend()
        eje_log.set_yscale("log")
        figura.tight_layout()
        plt.show()
        """
    ),
    markdown(
        """
        ## Tarea 2 - Passwords

        El programa pide nombre, fecha de nacimiento y DNI, y genera una clave
        de 8 caracteres con al menos 2 letras y 2 numeros, evitando fragmentos
        consecutivos de 2 caracteres tomados de los datos del usuario.

        Implementamos dos enfoques (la consigna sugiere "investigar
        diccionarios, listas y permutaciones"):

        1. **Aleatorio**: armamos passwords al azar y nos quedamos con la
           primera que cumple los requisitos.
        2. **Permutaciones**: recorremos `itertools.permutations` sobre un pool
           de caracteres y devolvemos la primera permutacion valida.
        """
    ),
    code(
        """
        import random
        from itertools import permutations


        def solo_alfanumericos(texto):
            return "".join(caracter for caracter in texto if caracter.isalnum())


        def pares_prohibidos(nombre, fecha, dni):
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


        def _separar(nombre, fecha, dni):
            letras = [c for c in solo_alfanumericos(nombre) if c.isalpha()]
            numeros = [c for c in solo_alfanumericos(fecha + dni) if c.isdigit()]
            if len(letras) < 2 or len(numeros) < 2:
                raise ValueError("Necesito al menos 2 letras y 2 numeros en los datos.")
            return letras, numeros


        def generar_aleatorio(nombre, fecha, dni, intentos=10000, semilla=7):
            letras, numeros = _separar(nombre, fecha, dni)
            caracteres = letras + numeros
            prohibidos = pares_prohibidos(nombre, fecha, dni)
            generador = random.Random(semilla)

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
            raise RuntimeError("No se pudo generar una clave valida.")


        def generar_permutaciones(nombre, fecha, dni, semilla=7, max_pool=12):
            letras, numeros = _separar(nombre, fecha, dni)
            prohibidos = pares_prohibidos(nombre, fecha, dni)

            generador = random.Random(semilla)
            pool = list(set(letras + numeros))
            if len(pool) > max_pool:
                pool = generador.sample(pool, max_pool)
            generador.shuffle(pool)

            for combinacion in permutations(pool, 8):
                password = "".join(combinacion)
                if cumple_requisitos(password, prohibidos):
                    return password
            raise RuntimeError("No se encontro una permutacion valida.")
        """
    ),
    code(
        """
        nombre = "Elena Gilbert"
        fecha = "20/04/2005"
        dni = "40123456"

        prohibidos = pares_prohibidos(nombre, fecha, dni)
        print("Pares prohibidos detectados (muestra):", sorted(prohibidos)[:10], "...")

        password_aleatoria = generar_aleatorio(nombre, fecha, dni)
        password_perm = generar_permutaciones(nombre, fecha, dni)

        print("\\nPassword (metodo aleatorio):     ", password_aleatoria)
        print("Cumple requisitos:                ", cumple_requisitos(password_aleatoria, prohibidos))

        print("\\nPassword (metodo permutaciones): ", password_perm)
        print("Cumple requisitos:                ", cumple_requisitos(password_perm, prohibidos))
        """
    ),
    markdown(
        """
        ### Validacion: passwords que NO cumplen

        Para mostrar que el chequeo funciona, evaluamos algunos casos
        construidos a mano:
        """
    ),
    code(
        """
        ejemplos = [
            ("Elena123",      "Tiene fragmento 'el' y 'le' del nombre y 7 caracteres"),
            ("12345678",      "No tiene letras"),
            ("abcdefgh",      "No tiene numeros"),
            ("xy20q4r5",      "Contiene '20' (de la fecha) y '04' (de la fecha)"),
            ("xy93q7r5",      "Deberia cumplir todo"),
        ]

        for password, comentario in ejemplos:
            ok = cumple_requisitos(password, prohibidos)
            estado = "OK   " if ok else "FALLA"
            print(f"{estado}  '{password}'  -> {comentario}")
        """
    ),
    markdown(
        """
        ## Tarea 3 - Conteo y nube de palabras

        Leemos un texto, contamos cuantas veces aparece una palabra dada y
        generamos una nube con la frecuencia de las palabras del texto.

        Para que el notebook sea autocontenido, el texto esta inline en la
        proxima celda. En el repositorio el mismo texto vive tambien en
        `data/texto_prueba.txt`.
        """
    ),
    code(
        """
        from collections import Counter
        import re

        TEXTO = (
            "La bioinformatica combina biologia, informatica y estadistica. "
            "En bioinformatica se analizan secuencias, genes, proteinas y datos biologicos. "
            "La programacion permite automatizar tareas, contar palabras, buscar patrones "
            "y comparar secuencias. Este texto de prueba sirve para contar palabras y "
            "generar una nube de palabras. Las palabras mas frecuentes en el texto "
            "suelen ser articulos y preposiciones, por eso conviene filtrarlas para que "
            "la nube de palabras muestre los temas reales del texto. Bioinformatica, "
            "secuencias, genes, datos y programacion son palabras tematicas relevantes."
        )

        STOPWORDS = {
            "a", "al", "con", "de", "del", "el", "en", "es", "la", "las", "lo",
            "los", "para", "por", "que", "se", "son", "su", "sus", "un", "una",
            "y", "este", "esta", "estos", "estas", "ser", "mas", "muestre",
            "suelen", "eso",
        }


        def tokenizar(texto):
            return re.findall(r"[A-Za-zÁÉÍÓÚÜÑáéíóúüñ]+", texto.lower())


        def contar_palabra(texto, palabra):
            return tokenizar(texto).count(palabra.lower())


        def frecuencias(texto, excluir_stopwords=True):
            palabras = tokenizar(texto)
            if excluir_stopwords:
                palabras = [p for p in palabras if p not in STOPWORDS]
            return Counter(palabras)


        palabra_buscada = "bioinformatica"
        print(f"'{palabra_buscada}' aparece:", contar_palabra(TEXTO, palabra_buscada), "veces")

        conteos = frecuencias(TEXTO)
        pd.DataFrame(conteos.most_common(15), columns=["palabra", "frecuencia"])
        """
    ),
    code(
        """
        from wordcloud import WordCloud

        nube = WordCloud(
            width=900,
            height=550,
            background_color="white",
            colormap="viridis",
            random_state=7,
        ).generate_from_frequencies(conteos)

        plt.figure(figsize=(9, 5.5))
        plt.imshow(nube, interpolation="bilinear")
        plt.axis("off")
        plt.title("Nube de palabras")
        plt.show()
        """
    ),
    markdown(
        """
        ### Reflexion

        - **¿Para que podria usar esta representacion?** Para tener una idea
          rapida de los temas dominantes en un texto: resumenes de articulos,
          analisis de respuestas abiertas en una encuesta, exploracion inicial
          de un corpus antes de un analisis mas formal. En bioinformatica
          podria servir, por ejemplo, para resumir terminos en abstracts de
          papers.
        - **¿Es necesario excluir alguna palabra?** Si. Articulos,
          preposiciones, conectores y verbos auxiliares (`la`, `el`, `de`,
          `que`, `y`, `por`, `para`, `es`, ...) aparecen mucho en cualquier
          texto en castellano y opacarian a las palabras tematicas. Sin
          stopwords la nube parece la misma para textos muy distintos.
        - **¿Estoy sesgada en alguna muletilla?** Probablemente si: al revisar
          mis propios textos noto que repito conectores como "entonces",
          "basicamente", "en realidad" mas seguido de lo necesario. La nube
          ayuda a verlo, porque las muletillas saltan en tamano grande sin
          aportar contenido tematico.
        """
    ),
]

destino = Path("notebooks/TP0_Bioinformatica_Entrega.ipynb")
destino.parent.mkdir(exist_ok=True)
with open(destino, "w", encoding="utf-8") as archivo:
    nbf.write(nb, archivo)
print(f"Notebook creado: {destino}")
