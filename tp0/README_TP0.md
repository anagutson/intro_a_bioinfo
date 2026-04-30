# TP0 - Entornos y primeros programas

Bioinformática 2026 - FCEyN UBA Exactas.

Este TP integra los ejercicios básicos de Python (variables, listas,
condicionales, loops, archivos y funciones) y las tres tareas de programación:
números primos, generador de passwords, y conteo + nube de palabras.

## Estructura

- `notebooks/TP0_Bioinformatica_Entrega_ejecutado.ipynb`: notebook con todas
  las celdas ya ejecutadas y las salidas visibles. Es **autocontenido**: se
  puede correr sin el resto del repositorio.
- `scripts/`: versiones CLI de los mismos ejercicios, organizadas por archivo.
- `data/texto_prueba.txt`: archivo de texto que usa la Tarea 3.
- `resultados/`: gráficos y CSV generados por los scripts.
- `requirements.txt`: lista de paquetes y versiones utilizadas.

## Entorno

Los paquetes se instalan en un entorno virtual local `.venv`:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

Para activarlo:

```bash
source .venv/bin/activate
```

## Archivos principales

- `scripts/inicio.py`: primer programa pedido en el ejercicio 1.
- `scripts/tp0_ejercicios_basicos.py`: ejercicios básicos agrupados (variables,
  listas, condicionales, loops, archivos y funciones).
- `scripts/tarea1_primos.py`: números primos con dos versiones (básica y
  mejorada con `sqrt(n)`), medición de tiempos y gráficos en escala lineal y
  logarítmica.
- `scripts/tarea2_passwords.py`: generador de passwords con dos enfoques
  (aleatorio y por permutaciones, según la consigna del TP).
- `scripts/tarea3_palabras.py`: conteo de palabras, gráfico de barras y nube
  de palabras con `wordcloud`.

## Cómo ejecutar

Abrir el notebook:

```bash
.venv/bin/jupyter notebook notebooks/TP0_Bioinformatica_Entrega_ejecutado.ipynb
```

Correr los scripts desde la terminal:

```bash
.venv/bin/python scripts/inicio.py
.venv/bin/python scripts/tp0_ejercicios_basicos.py
.venv/bin/python scripts/tarea1_primos.py --limite 50
.venv/bin/python scripts/tarea2_passwords.py --nombre "Elena Gilbert" --fecha "20/04/2005" --dni "40123456"
.venv/bin/python scripts/tarea2_passwords.py --nombre "Elena Gilbert" --fecha "20/04/2005" --dni "40123456" --metodo permutaciones
.venv/bin/python scripts/tarea3_palabras.py --archivo data/texto_prueba.txt --palabra bioinformatica
```

## Pseudocódigo de la Tarea 1

Un número primo es un entero mayor que 1 que sólo puede dividirse exactamente
por 1 y por sí mismo. Para detectarlo, el programa busca si existe algún
divisor intermedio.

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

Para saber si `n` es primo alcanza con probar divisores hasta `raíz(n)`. Si no
aparece ningún divisor hasta ese punto, el número es primo.

## Respuestas breves

- **Rol de `x` en un `for`:** la variable elegida toma de a uno los valores de
  la secuencia. Por ejemplo, en `for x in [2, 3, 5]`, primero `x` vale 2,
  luego 3 y luego 5.
- **Stopwords en la nube:** conviene excluir artículos, preposiciones y
  conectores muy frecuentes (`la`, `el`, `de`, `que`, `y`, `por`, ...) porque
  ocupan mucho lugar en la nube sin aportar significado temático.
- **Para qué sirve la nube de palabras:** para tener una vista rápida de los
  temas dominantes en un texto (resúmenes, abstracts, respuestas abiertas).
