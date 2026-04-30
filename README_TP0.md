# TP0 - Entornos y primeros programas

Entrega preparada para ejecutar con Python 3 desde la terminal.

## Entorno

Los paquetes estan instalados en un entorno virtual local `.venv`.

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
- `resultados/`: archivos generados por los programas, como graficos y CSV.
- `requirements.txt`: paquetes usados para recrear el entorno.

## Archivos principales

- `scripts/inicio.py`: primer programa pedido en el ejercicio 1.
- `scripts/tp0_ejercicios_basicos.py`: ejercicios basicos agrupados: variables, listas, condicionales, loops, archivos y funciones.
- `scripts/tarea1_primos.py`: numeros primos, limite definido por usuario, medicion de tiempos y graficos.
- `scripts/tarea2_passwords.py`: generador de passwords con nombre, fecha de nacimiento y DNI.
- `scripts/tarea3_palabras.py`: conteo de una palabra, frecuencia de palabras, grafico de barras y nube de palabras.
- `data/texto_prueba.txt`: archivo de texto para probar la tarea 3.
- `notebooks/TP0_Bioinformatica_Entrega_ejecutado.ipynb`: notebook ejecutado, recomendado para entregar.

## Como ejecutar

```bash
.venv/bin/python scripts/inicio.py
.venv/bin/python scripts/tp0_ejercicios_basicos.py
.venv/bin/python scripts/tarea1_primos.py --limite 50
.venv/bin/python scripts/tarea2_passwords.py --nombre "Elena Gilbert" --fecha "20/04/2005" --dni "40123456"
.venv/bin/python scripts/tarea3_palabras.py --archivo data/texto_prueba.txt --palabra bioinformatica
```

Para abrir el notebook:

```bash
.venv/bin/jupyter notebook notebooks/TP0_Bioinformatica_Entrega_ejecutado.ipynb
```

## Pseudocodigo de la tarea 1

Un numero primo es un entero mayor que 1 que solo puede dividirse exactamente por 1 y por si mismo.
Para detectarlo, el programa busca si existe algun divisor intermedio.

Version basica:

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

Version mejorada:

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

## Respuestas breves

En un loop `for`, la variable elegida toma de a uno los valores de la secuencia. Por ejemplo, en
`for x in [2, 3, 5]`, primero `x` vale 2, luego 3 y luego 5.

En el conteo de palabras conviene excluir articulos, preposiciones y conectores muy frecuentes
porque suelen ocupar mucho lugar en la nube sin aportar demasiado significado tematico.
