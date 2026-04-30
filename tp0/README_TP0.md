# TP0 - Entornos y primeros programas

Entrega preparada para ejecutar con Python 3 desde la terminal o desde Jupyter Notebook.

## Entorno

Los paquetes están instalados en un entorno virtual local `.venv`.

Para activarlo:

```bash
source .venv/bin/activate
```

Para recrearlo desde cero en otra computadora:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

## Estructura

- `scripts/`: programas Python del TP.
- `notebooks/`: notebook de entrega.
- `data/`: archivos de entrada para probar los programas.
- `resultados/`: archivos generados por los programas, como gráficos y CSV.
- `requirements.txt`: paquetes usados para recrear el entorno.

## Archivos principales

- `scripts/inicio.py`: primer programa pedido en el ejercicio 1.
- `scripts/tp0_ejercicios_basicos.py`: ejercicios básicos agrupados: variables, listas, condicionales, loops, archivos y funciones.
- `scripts/tarea1_primos.py`: números primos, límite definido por usuario, medición de tiempos y gráfico comparativo.
- `scripts/tarea2_passwords.py`: generador de passwords con nombre, fecha de nacimiento y DNI. Incluye una versión aleatoria y otra con `itertools.permutations`.
- `scripts/tarea3_palabras.py`: conteo de una palabra, frecuencia de palabras, gráfico de barras y nube de palabras.
- `data/texto_prueba.txt`: archivo de texto para probar la tarea 3.
- `notebooks/TP0_Bioinformatica_Entrega.ipynb`: notebook limpio y autocontenido.
- `notebooks/TP0_Bioinformatica_Entrega_ejecutado.ipynb`: notebook ejecutado, recomendado para entregar.

## Cómo ejecutar

```bash
.venv/bin/python scripts/inicio.py
.venv/bin/python scripts/tp0_ejercicios_basicos.py
.venv/bin/python scripts/tarea1_primos.py --limite 50
.venv/bin/python scripts/tarea2_passwords.py --nombre "Elena Gilbert" --fecha "20/04/2005" --dni "40123456" --metodo permutaciones
.venv/bin/python scripts/tarea3_palabras.py --archivo data/texto_prueba.txt --palabra bioinformatica
```

Para abrir el notebook:

```bash
.venv/bin/python -m notebook notebooks/TP0_Bioinformatica_Entrega_ejecutado.ipynb
```

## Pseudocódigo de la tarea 1

Un número primo es un entero mayor que 1 que solo puede dividirse exactamente por 1 y por sí mismo.
Para detectarlo, el programa busca si existe algún divisor intermedio.

Versión básica:

```text
para cada número entre 1 y límite:
    si número < 2:
        no es primo
    si no:
        asumir que es primo
        para cada divisor entre 2 y número - 1:
            si número % divisor == 0:
                no es primo
        si ningún divisor lo dividió:
            guardar/imprimir el número como primo
```

El operador `%` calcula el resto de una división. Si `número % divisor == 0`,
significa que el divisor entra justo y entonces el número no es primo.

Versión mejorada:

```text
para cada número entre 1 y límite:
    si número < 2:
        no es primo
    si número == 2:
        es primo
    si número es par y mayor que 2:
        no es primo
    si no:
        probar divisores impares: 3, 5, 7, ...
        detenerse cuando divisor * divisor > número
        si ningún divisor lo dividió:
            es primo
```

Por eso, para saber si `n` es primo, alcanza con probar divisores hasta `raíz(n)`.
Si no aparece ningún divisor hasta ese punto, el número es primo.

## Respuestas breves

En un loop `for`, la variable elegida toma de a uno los valores de la secuencia. Por ejemplo, en
`for x in [2, 3, 5]`, primero `x` vale 2, luego 3 y luego 5.

En la tarea 2 se muestran dos estrategias: una aleatoria, que es práctica para generar claves,
y otra con permutaciones, que responde de forma más directa a la consigna del TP.

En el conteo de palabras conviene excluir artículos, preposiciones y conectores muy frecuentes
porque suelen ocupar mucho lugar en la nube sin aportar demasiado significado temático. La nube
puede servir para explorar rápidamente los temas dominantes de un texto, por ejemplo abstracts
de papers o respuestas abiertas de una encuesta. También ayuda a detectar muletillas: si una
palabra aparece muy grande y no aporta contenido, probablemente convenga revisar el texto.
