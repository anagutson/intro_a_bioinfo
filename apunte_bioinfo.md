# Apunte Bioinformática — Examen a Libro Abierto

---

# ÍNDICE POR CLASE

- [Clase 1 — Marco Conceptual y Bases de Datos Primarias](#clase-1)
- [Clase 2 — Alineamiento, Matrices y Aminoácidos](#clase-2)
- [Clase 3 — BLAST y Alineamiento Múltiple](#clase-3)
- [Clase 4 — Bases de Datos Secundarias](#clase-4)
- [Clase 5 — NGS y Ensamblado de Genomas](#clase-5)
- [Clase 6 — Genómica Humana y Medicina Personalizada](#clase-6)
- [Parcial 2018 — Resolución](#parcial-2018)

---

## Clase 1 — Marco Conceptual y Bases de Datos Primarias

### Qué es la bioinformática y por qué importa la pregunta biológica

En clase el profe arrancó con una advertencia: la bioinformática no es un fin en sí mismo. El error que ve todo el tiempo es el de alguien que dice "quiero hacer bioinformática" sin tener ni idea de cuál es su pregunta biológica. Eso conduce a la frustración pura: uno empieza a probar páginas web, programas, bases de datos, y no llega a ningún lado. La pregunta biológica es la zanahoria adelante del burro. Sin ella, todo lo demás no tiene sentido.

La bioinformática en sentido amplio es el uso de la tecnología de la información para gestionar y analizar datos biológicos. En el sentido más práctico de esta materia: el uso de programas y algoritmos aplicados a estudios biológicos. Tiene tres aspectos: (1) desarrollo de teoría hacia algoritmos (codificar la biología en pasos lógico-matemáticos), (2) aplicación de esos algoritmos a datos para generar predicciones, (3) organización de los datos para que los algoritmos puedan correr.

### Jerarquía: Teoría → Algoritmo → Programa → Dato

**La teoría** es una idea, un modelo de cómo funciona un fenómeno biológico. La teoría celular, las leyes de Mendel, el dogma central, son teorías. No son informática todavía: están en papel, en libros, en la cabeza de los investigadores.

**El algoritmo** es el paso bisagra: cuando uno toma esa teoría y la traduce en una serie de pasos lógico-matemáticos con reglas específicas. Dado un input, producen un output. El algoritmo es la informatización de la teoría. Es independiente del lenguaje de programación: puede estar escrito en papel.

El ejemplo del profe para que se entienda la diferencia entre teoría y algoritmo fue **la división**. Todos saben dividir conceptualmente (6 caramelos para 3 amigos = 2 cada uno). Pero si decís: dividí 12.432.324 por 752, nadie lo puede hacer a mano de manera rápida. Sin embargo, todos aprendieron en la primaria el **algoritmo** de dividir: tomá la primera cifra, buscala en la tabla del divisor, anotá el cociente, restá, bajá la siguiente cifra, repetí. El algoritmo transforma el problema imposible en pasos simples. El dato son las tablas de multiplicar. El programa sería implementar eso en Python.

Otro ejemplo: **Mendel**. Mendel observó plantas y anotó los resultados (datos). De las observaciones formuló una teoría (alelos, dominante/recesivo, segregación). De la teoría bajó a un algoritmo: si cruzás dos heterocigotas, las frecuencias son 1:2:1. La bioinformática aparece cuando la cantidad de datos se vuelve tan grande que no se puede hacer a mano.

| Nivel | Definición | Ejemplo |
|-------|-----------|---------|
| **Teoría** | Modelo que describe el fenómeno biológico | Evolución por sustitución |
| **Algoritmo** | Serie de pasos lógico-matemáticos, independiente del lenguaje | Needleman-Wunsch |
| **Programa** | Implementación en un lenguaje específico | BLAST, CLUSTALW, hmmbuild |
| **Dato** | Lo medido experimentalmente; input del algoritmo | Secuencia FASTA |

> En el parcial se puede responder con el nombre del **algoritmo** O con el nombre del **programa** — ambas respuestas son equivalentes. "Haría un alineamiento de pares" = "usaría NW/SW" = "usaría BLAST".

**Observación importante del profe sobre eficiencia:** la eficiencia no importa cuando se diseña el algoritmo, importa cuando se implementa el programa. Si dos alternativas algorítmicas son del mismo orden de magnitud, siempre conviene la más clara. La eficiencia importa cuando la diferencia es de órdenes de magnitud (O(n²) vs O(n·log n) vs O(n)).

### Por qué el ADN es el origen de toda la bioinformática

Watson y Crick terminaron su paper con la frase: "it has not escaped our notice that with specific pairing we have postulated a possible copying mechanism for the genetic material". Con eso describían el mecanismo de replicación. Pero lo más profundo es que la secuencia lineal del ADN es una forma de contener información. El ADN es un lenguaje de cuatro letras, igual que el lenguaje humano es un lenguaje de 26 letras. Y la bioinformática es, en el fondo, el manejo de esa información.

La analogía del profe: descifrar el ADN es como los arqueólogos intentando entender los jeroglíficos antes de la Piedra Roseta. La secuencia está ahí, pero no se sabe qué significa. Los biólogos fueron descubriendo el código genético usando secuencias que sabían lo que codificaban como tabla de equivalencias.

La relación genotipo-fenotipo no es como un plano y un avión (bijección: si cambio el plano, cambio el avión). Es como una receta y una torta. Si tenés dos tortas que saben distinto, no podés saber de inmediato qué cambió en la receta. La bioinformática codifica los algoritmos que relacionan genotipo y fenotipo.

### Ortólogos y parálogos

La homología tiene dos sabores:

- **Ortólogos**: genes homólogos en distintas especies que derivan del mismo gen ancestral por **especiación**. El gen de la hemoglobina alfa en humano y en ratón son ortólogos. Tienden a conservar la misma función. Son los más confiables para inferir función por homología.
- **Parálogos**: genes homólogos dentro de la misma especie que derivan de una **duplicación génica**. La hemoglobina alfa y la beta en humano son parálogas. Pueden tener funciones similares pero diferenciadas. Menos confiables para inferir función directamente.

### Principio de "culpa por asociación"

"Guilt by association" es el fundamento detrás de prácticamente toda búsqueda en bases de datos biológicas. Si una secuencia desconocida se parece mucho a un grupo de secuencias con función conocida, se infiere que tiene la misma función, **sin necesidad de experimentar directamente**.

¿Por qué funciona? Por la homología. Si dos genes tienen secuencias muy similares, la explicación más parsimoniosa es que descienden de un ancestro común. Y si descienden de un ancestro común, lo más probable es que hayan heredado también su función.

**Ejemplo:** si tengo una proteína nueva y búsqueda en BD muestra que su secuencia se parece mucho a kinasas conocidas, puedo predecir que también es una kinasa. Eso ahorra años de experimentos. No es certeza absoluta pero es una hipótesis sólida.

El profe lo ilustró a nivel de aminoácido individual: si querés hacer una mutante con presupuesto para solo dos mutaciones, ¿cuál elegís? Buscás los aminoácidos subrepresentados (cisteínas, triptofanos, histidinas) y los mutás primero. Están subrepresentados precisamente porque la evolución los conservó con mayor presión selectiva: cuando aparecen, suelen estar en el sitio activo.

### Estructura típica del parcial

| Problema | Tema |
|----------|------|
| **P1** | Alineamiento / Programación Dinámica (NW o SW) |
| **P2** | Bases de datos primarias (criterios: curación, redundancia, primaria/secundaria) |
| **P3** | MSA y árbol guía (CLUSTALW, SP score, E-value, PROSITE) |
| **P4** | Bases de datos secundarias (HMM, logos, sensibilidad/especificidad) |

### Tips del profe (de la clase de repaso)

- Leer el parcial **completo** antes de empezar.
- Responder de lo más **fácil a lo más difícil**.
- La respuesta **simple es la correcta** — el parcial no tiene trampa.
- Cuando veas un logo, fijarse si las letras son aminoácidos (A/R/N...) o nucleótidos (A/T/G/C). ¡Son distintos!
- Simple answer = correct answer. No complicar de más.

---


### Por qué existen: historia y necesidad real

La historia de las bases de datos biológicas no empieza con internet. El profe la contó así:

Todo empieza con las **proteínas**. Sanger (una de las tres personas en toda la historia en recibir dos premios Nobel por dos descubrimientos distintos, junto con Marie Curie y Linus Pauling) primero desarrolló el método para secuenciar proteínas y décadas más tarde para secuenciar ADN. En la década de 1960 ya se conocían varias proteínas y alguien juntó toda esa información en el "Atlas of Protein Sequence and Structure": un registro por proteína con todo lo que se sabía. Primero era un libro. En 1972 ese libro fue pasado a cinta magnética.

En la década de 1980 Amos Bairoch digitalizó ese Atlas y creó **Swiss-Prot**. Lo distribuyó a través de la red universitaria precursora de internet. Ese fue el puntapié inicial de la bioinformática moderna: una base de datos digital de proteínas, accesible remotamente.

Al mismo tiempo, con el método de Sanger para secuenciar ADN volviéndose masivo, había que guardar las secuencias. En 1988, el NCBI de EE.UU. creó GenBank, el EMBL europeo creó su BD de nucleótidos (hoy ENA), y Japón creó el DDBJ. Ese mismo año los tres formaron el **INSDC** (International Nucleotide Sequence Database Collaboration): un consorcio donde los tres se sincronizan diariamente y mantienen el mismo formato. Si depositás en cualquiera de los tres, en 24 horas aparece en los otros dos.

Al principio todo era simple. Pero con el tiempo la cantidad de datos creció explosivamente y vinieron los problemas que justifican el diseño actual: información duplicada, errores, miles de entradas para el mismo gen.

### Criterios de clasificación

| Dimensión | Opciones | Definición |
|-----------|----------|-----------|
| **Curación** | Curada vs No curada | Curada = revisión manual por expertos antes de aceptar el registro. No curada = el autor deposita sin revisión sistemática. |
| **Redundancia** | Redundante vs No redundante | Redundante = el mismo dato puede aparecer en múltiples registros. |
| **Origen** | Primaria vs Secundaria | Primaria = datos experimentales crudos. Secundaria = resultado del análisis bioinformático de datos primarios. |

> **Regla:** para clasificar una BD, responder las 3 preguntas. Se puede ser curada y primaria (RefSeq, Swiss-Prot), no curada y primaria (GenBank, TrEMBL), o curada y secundaria (SGD, PFAM, KEGG).

### El problema de la redundancia: por qué un mismo gen puede tener 50 entradas en GenBank y 1 en RefSeq

El profe lo explicó con un ejemplo concreto: el Dr. K secuenció el gen de la fibronectina en 1985 y lo subió a GenBank. Diez años después, el Dr. L secuencia el mismo gen de la misma especie y también lo sube. ¿Y si difieren en una base? ¿Es una mutación, un polimorfismo, un error de secuenciación?

GenBank toma la decisión de repositorio: cada autor sube lo que quiere y es responsable del contenido. Por eso, si buscás el gen de la insulina humana en GenBank, podés encontrar 50, 100, 200 entradas.

Esta situación se volvió inmanejable con la secuenciación masiva. Entonces se creó **RefSeq** (Reference Sequence Database). RefSeq no es un repositorio: es una BD curada, donde para cada organismo modelo, un equipo establece cuál es la secuencia de referencia de cada molécula. El criterio: **una molécula física = un registro**.

Si el gen de la insulina humana tiene tres variantes de splicing, RefSeq tiene tres registros (NM_xxxxxx.1, NM_xxxxxx.2, NM_xxxxxx.3), uno por cada ARNm. En cambio, GenBank puede tener miles de entradas del mismo gen.

La analogía del profe: **GenBank es como un paper** (el autor presenta sus datos, pasa una revisión mínima de formato, y se publica sin que nadie verifique la biología). **RefSeq es como un review** (un experto revisa toda la literatura, consensúa, y publica la secuencia de referencia).

### Curada vs no curada: qué pasa concretamente cuando hay un error

**Definición del profe:** BD curada = entre el autor que genera el dato y el sistema que lo guarda, hay un proceso de revisión que verifica si ese dato es correcto. No curada (repositorio) = el autor sube directamente y nadie revisa.

**En un repositorio** como GenBank: un lab sube una secuencia con un error de secuenciación (un artefacto, una base incorrecta). Esa secuencia queda indefinidamente. Si alguien la usa para un alineamiento, puede llegar a conclusiones incorrectas. Históricamente hubo casos de genes reportados como "nuevos" que eran artefactos de secuenciación o contaminaciones.

**En una BD curada** como Swiss-Prot: antes de que una entrada entre, un curador humano revisa la secuencia, la contrasta con la literatura, verifica la anotación funcional, y solo entonces la acepta. Eso hace que Swiss-Prot tenga muchas menos entradas que TrEMBL, pero cada entrada vale más en términos de confiabilidad.

La situación práctica: si estás trabajando con un organismo nuevo o poco estudiado, es probable que no esté en las BDs curadas. Tenés que ir al repositorio y asumís que los datos tienen ruido. Si estás trabajando con un organismo modelo (humano, ratón, Drosophila, levadura), hay BDs curadas específicas donde la calidad es mucho mayor.

### Primaria vs secundaria: el principio más importante para el parcial

**BD primaria**: almacena el dato crudo obtenido directamente de experimentos. La secuencia de un gen (medida con un secuenciador), la estructura de una proteína (resuelta por cristalografía o cryo-EM), el perfil de expresión génica (medido por RNA-seq). El dato viene de un experimento de laboratorio o de una observación directa. GenBank, RefSeq, PDB, Swiss-Prot, TrEMBL, GEO: todas son primarias.

**BD secundaria**: almacena el resultado del análisis bioinformático de datos primarios. No hay nuevo experimento de laboratorio: hay análisis computacional. Las BDs de familias de proteínas (PFAM, PROSITE) se construyen tomando miles de secuencias de BDs primarias, alineándolas, detectando patrones conservados, y agrupándolas en familias. Eso es bioinformática aplicada a datos primarios.

Otros ejemplos de secundarias: SGD (integra GenBank + UniProt + literatura para levadura), FlyBase (lo mismo para Drosophila), KEGG (integra genomas y los mapea a pathways metabólicos), Ensembl (anota genomas de vertebrados integrando datos de múltiples fuentes).

La diferencia práctica: si querés la secuencia de un gen, vas a una BD primaria. Si querés saber a qué familia pertenece tu proteína o a qué pathway está asociada, vas a una BD secundaria.

### Bases de nucleótidos

| BD | Curación | Redundancia | Tipo | Notas |
|----|----------|-------------|------|-------|
| **GenBank** | No curada | Redundante | Primaria | Repositorio NCBI; el autor deposita y es responsable; INSDC |
| **EMBL/ENA** | No curada | Redundante | Primaria | Equivalente europeo (EBI); mismo formato INSDC |
| **DDBJ** | No curada | Redundante | Primaria | Equivalente japonés; sincronizado con INSDC |
| **RefSeq** | Curada | No redundante | Primaria | 1 registro por molécula por organismo; NCBI lo mantiene |
| **Ensembl** | Curada | No redundante | Secundaria | Browser genómico; anotación de vertebrados |

> GenBank = **paper**; RefSeq = **review**. Ambas primarias, pero RefSeq curada.

**Prefijos RefSeq:**
| Prefijo | Contenido |
|---------|-----------|
| NM\_ | ARNm verificado experimentalmente |
| NP\_ | Proteína verificada experimentalmente |
| XM\_ | ARNm predicho in silico |
| XP\_ | Proteína predicha in silico |
| NT\_ | contig genómico |
| NC\_ | cromosoma completo |

Si en el parcial ves un accession que empieza con NM\_ o NP\_, es RefSeq verificado. XM\_ o XP\_ es RefSeq pero predicción computacional.

### Swiss-Prot vs TrEMBL: calidad vs cobertura

**Swiss-Prot** es la BD de proteínas curada manualmente. Cada entrada tiene nombre del gen, función anotada, dominios, sitios activos, modificaciones post-traduccionales, variantes conocidas, todo verificado por un curador humano que leyó la literatura. El criterio de redundancia: "una proteína funcional = un registro". Todas las isoformas de una misma proteína (por splicing alternativo) van en el mismo registro, con una secuencia canónica y las variantes anotadas.

**TrEMBL** (Translated EMBL) es la traducción automática de todos los registros codificantes de GenBank y EMBL. No hay curación manual. Es enormemente grande (cientos de millones de entradas) pero con muchos errores: anotaciones incorrectas, proteínas hipotéticas, artefactos de secuenciación traducidos.

**UniProtKB** es la unión de Swiss-Prot y TrEMBL bajo una misma interfaz. Cuando buscás en UniProt, podés filtrar por "Swiss-Prot" (curada) o "TrEMBL" (no curada).

**Distinción de criterio de redundancia**: **RefSeq organiza por molécula física** (un ARNm = un registro, aunque sea una isoforma del mismo gen). **UniProt/Swiss-Prot organiza por unidad funcional** (un gen/proteína = un registro, con todas las isoformas anotadas dentro). Para un gen con 20 isoformas, RefSeq tiene 20 registros de ARNm y 20 registros de proteína; Swiss-Prot tiene 1 registro de proteína con las 20 isoformas anotadas.

### Bases de proteínas

| BD | Curación | Redundancia | Tipo | Notas |
|----|----------|-------------|------|-------|
| **Swiss-Prot** | Curada | No redundante | Primaria | 1 registro = 1 gen; todas las isoformas juntas |
| **TrEMBL** | No curada | Redundante | Primaria | Traducción automática de GenBank/EMBL |
| **UniProtKB** | Mixta | Baja | Primaria | Swiss-Prot + TrEMBL integrados |
| **GenPept** | No curada | Redundante | Primaria | Traducción automática del INSDC (NCBI) |
| **PDB** | Curada | Redundante | Primaria | Estructuras 3D; misma proteína puede tener cientos de estructuras |

### PDB: por qué hay 300 estructuras de la misma proteína

La **Protein Data Bank (PDB)** almacena estructuras tridimensionales de macromoléculas biológicas. Es primaria (datos de experimentos), curada (hay revisión antes de depositar), y redundante (múltiples estructuras por proteína).

¿Por qué redundante? Para la misma proteína puede haber cientos de estructuras depositadas porque:
- La misma proteína resuelta en distintas condiciones (pH diferente, temperatura diferente).
- La proteína con distintos ligandos: apo (sin ligando), con sustrato, con inhibidor, con cofactor.
- La proteína con distintas mutaciones puntuales para estudiar el efecto de cada mutación.
- La proteína en distintos contextos de cristalización que dan resoluciones distintas.
- Distintos laboratorios que resuelven la misma proteína de manera independiente.

Cada estructura tiene un identificador de 4 caracteres alfanuméricos (ej: 1HHO es una estructura de hemoglobina). Para el parcial: si te preguntan sobre una BD de estructuras, la respuesta es PDB. Es primaria, curada, y redundante.

### INSDC: el acuerdo de intercambio entre GenBank, EMBL y DDBJ

El INSDC formado en 1988 explica por qué GenBank, EMBL y DDBJ son esencialmente lo mismo pero tres cosas distintas. Los tres son repositorios no curados y redundantes de secuencias de nucleótidos. Los tres se sincronizan diariamente. El formato de los registros es el mismo. Para el usuario práctico, es indiferente consultar cualquiera de los tres: el contenido es el mismo. Existen los tres por razones históricas y políticas: cada región geográfica quería tener su nodo.

### Otras bases importantes

| BD | Tipo | Curación | Notas |
|----|------|----------|-------|
| **GEO** | Primaria | No curada | Repositorio NCBI de microarrays y RNA-seq |
| **ArrayExpress** | Primaria | No curada | Equivalente europeo (EBI) |
| **GTEx** | Secundaria | Curada | Expresión por tejido en humanos sanos |
| **dbSNP** | Primaria | No curada | Repositorio de todos los SNPs; IDs `rs` |
| **ClinVar** | Primaria | Curada | Variantes con significado clínico (Patogénica/VUS/Benigna) |
| **OMIM** | Secundaria | Curada | Genes y enfermedades mendelianas |
| **gnomAD** | Secundaria | Curada | Frecuencias alélicas en ~125.000 genomas de referencia |
| **COSMIC** | Secundaria | Curada | Mutaciones somáticas en cáncer |
| **SGD** | Secundaria | Curada | S. cerevisiae; integra GenBank + UniProt + literatura |
| **FlyBase** | Secundaria | Curada | D. melanogaster |
| **KEGG** | Secundaria | Curada | Pathways metabólicos y de señalización |
| **IntAct** | Secundaria | Curada | Interacciones proteína-proteína |
| **DrugBank** | Secundaria | Curada | Fármacos con targets, mecanismo, farmacocinética |

### Lo que se pregunta en el parcial sobre bases de datos

**"Dar un ejemplo de BD curada":** RefSeq, Swiss-Prot, PDB, SGD, FlyBase, ClinVar, KEGG. Argumentar: hay revisión humana antes de que el dato entre.

**"Dar un ejemplo de BD no curada":** GenBank, TrEMBL, GEO, ArrayExpress, dbSNP. Argumentar: el autor deposita directamente, no hay verificación sistemática.

**"Dar un ejemplo de BD primaria":** GenBank, RefSeq, Swiss-Prot, TrEMBL, PDB, GEO. Argumentar: el dato viene de experimentos de laboratorio.

**"Dar un ejemplo de BD secundaria":** SGD, FlyBase, KEGG, Ensembl, PFAM, PROSITE, IntAct, DrugBank. Argumentar: el contenido se genera por análisis bioinformático de datos primarios, no por nuevos experimentos.

**"¿Por qué una búsqueda de fibronectina en GenBank devuelve 10.000 resultados?":** GenBank es no curada y redundante. Cada laboratorio que secuenció cualquier fragmento o variante subió su secuencia. Para la secuencia de referencia, hay que ir a RefSeq o Swiss-Prot.

**"¿Por qué hay 300 estructuras de la misma proteína en el PDB?":** el PDB es redundante porque la misma proteína puede cristalizarse en distintas condiciones, con distintos ligandos, con distintas mutaciones, en distintos estados conformacionales. Cada estructura aporta información biológica diferente.

**"¿Qué es el INSDC?":** el International Nucleotide Sequence Database Collaboration, formado en 1988 entre GenBank (NCBI, EE.UU.), EMBL/ENA (EBI, Europa) y DDBJ (Japón). Los tres se sincronizan diariamente y usan el mismo formato. Depositar en cualquiera equivale a depositar en los tres.

---

## Clase 2 — Alineamiento de Secuencias, Matrices y Aminoácidos

### Por qué importan las propiedades fisicoquímicas en bioinformática

En bioinformática, los aminoácidos no son solo letras de un alfabeto de 20 caracteres. Cada letra trae consigo una identidad química que determina tres cosas clave:

1. **Qué sustituciones son toleradas evolutivamente** → eso es exactamente lo que codifica BLOSUM62. Si dos aminoácidos comparten propiedades fisicoquímicas, es probable que uno pueda reemplazar al otro sin destruir la función → la sustitución se observa en la evolución → score positivo en la matriz.

2. **Qué residuos aparecen en qué contextos estructurales.** Los hidrofóbicos se entierran en el núcleo de la proteína, lejos del solvente. Los polares y cargados aparecen en la superficie. La Glicina aparece en giros porque es flexible. La Prolina aparece en kinks porque la fuerza.

3. **Qué motivos de secuencia tienen significado funcional.** El motivo de N-glicosilación es N-{P}-[S/T] (Asn, cualquier cosa excepto Pro, Ser o Thr). Si hubiera Pro en esa segunda posición, no funciona → porque la Pro distorsiona la estructura local que reconoce la enzima transferasa.

### Clasificación por grupos fisicoquímicos — con razonamiento

**Hidrofobicidad** importa porque el plegamiento de proteínas en solución acuosa está impulsado principalmente por el efecto hidrofóbico: los residuos apolares tienden a quedar en el interior, alejados del agua.

```
HIDROFÓBICOS: I  V  L  M  A  G  C  F  Y  W
```

**Carga eléctrica** importa para interacciones electrostáticas. Los residuos cargados aparecen en la superficie de la proteína. Forman puentes salinos, sitios de unión a DNA (Lys y Arg básicas interactúan con el fosfato), y sitios activos enzimáticos.

```
ÁCIDOS (carga −): D (Asp, pKa ≈3.9)   E (Glu, pKa ≈4.1)
BÁSICOS (carga +): K (Lys, pKa ≈10.5)  R (Arg, pKa ≈12.5)  H (His, pKa ≈6)
```

**Capacidad de formar puentes de hidrógeno** — los residuos polares (con grupos –OH, –NH₂, –CONH₂) forman puentes H con el solvente y entre sí. Esos puentes contribuyen a la estructura secundaria (hélices α y láminas β).

```
POLARES sin carga: S  T  C  N  Q
```

**Aromaticidad** — los aromáticos tienen anillos planos que apilan entre sí (interacciones π-π) y absorben UV a ~280 nm. Esa absorbancia se usa en laboratorio para cuantificar proteínas.

```
AROMÁTICOS: F  Y  W
```

### Los 20 aminoácidos — referencia completa

| 1L | 3L  | Nombre           | Grupo principal        | Notas clave para el parcial |
|:--:|:---:|------------------|------------------------|-----------------------------|
| A  | Ala | Alanina          | Pequeño, hidrofóbico   | Frecuente en α-hélices |
| R  | Arg | Arginina         | Básico (+)             | pKa ≈12.5 → carga + muy estable a pH fisiológico; interacciona con DNA |
| N  | Asn | Asparagina       | Polar sin carga        | Sitio de N-glicosilación: motivo **N-{P}-[S/T]** |
| D  | Asp | Ácido aspártico  | Ácido (−)              | pKa ≈3.9; cadena corta; score D↔E = +2 |
| C  | Cys | Cisteína         | Polar / hidrofóbico    | Puentes disulfuro; sitios activos nucleofílicos; 2do más raro |
| Q  | Gln | Glutamina        | Polar sin carga        | N↔Q y D↔E son sustituciones conservativas |
| E  | Glu | Ácido glutámico  | Ácido (−)              | pKa ≈4.1; E = D + un CH₂ |
| G  | Gly | Glicina          | Tiny / flexible        | **Sin cadena lateral → máxima flexibilidad → loops, giros, colágeno (Gly-X-Y)** |
| H  | His | Histidina        | Básico / aromático     | **Único pKa ≈6 → puede ceder o tomar H⁺ a pH 7 → catálisis ácido-base (tríada Ser-His-Asp)** |
| I  | Ile | Isoleucina       | Alifático, hidrofóbico | Isómero de Leu; L↔I es sustitución muy conservativa |
| L  | Leu | Leucina          | Alifático, hidrofóbico | El más frecuente en proteínas globulares |
| K  | Lys | Lisina           | Básico (+)             | pKa ≈10.5; sitio de ubiquitinación → degradación; interacciona con DNA |
| M  | Met | Metionina        | Hidrofóbico, azufre    | Codón de inicio AUG → toda proteína empieza con Met |
| F  | Phe | Fenilalanina     | Aromático              | Núcleo hidrofóbico; absorbe UV a 280nm levemente |
| P  | Pro | Prolina          | Especial / rígido      | **Anillo pirrolidina unido al N del backbone → no puede donar H en puentes H → rompe hélices α y láminas β → fuerza kinks → NO va en posición X del motivo N-glicosilación** |
| S  | Ser | Serina           | Polar pequeño (–OH)    | Fosforilación (Ser-kinasas); O-glicosilación; posición [S/T] en N-glicosilación |
| T  | Thr | Treonina         | Polar (–OH)            | Fosforilación; posición [S/T] en motivo N-glicosilación |
| W  | Trp | Triptófano       | Aromático, grande      | **El más raro (~1%) → casi siempre crítico funcionalmente → score W↔W = +11 → mayor absorbancia UV a 280nm** |
| Y  | Tyr | Tirosina         | Aromático + polar (–OH)| Fosforilación (Tyr-kinasas → señalización); absorbe fuerte a 280nm |
| V  | Val | Valina           | Alifático, hidrofóbico | Volumen: V < L = I; ramificado → desestabiliza hélices α |

### Diagrama de grupos fisicoquímicos

```
                      TINY
                    ┌──────────┐
                    │  G   A   │
                    └──────────┘
                    (subconjunto de SMALL)
┌─────────────────────────────────────────────────────┐
│                      SMALL                          │
│         G   A   C   S   N   D   T   P               │
└─────────────────────────────────────────────────────┘

ALIPHATIC:     I   V   L
AROMATIC:      F   Y   W

HYDROPHOBIC (se entierran en el núcleo):
   I   V   L   M   A   G   C   F   Y   W

POLAR (hacen puentes H, suelen estar en superficie):
   S   T   C   N   Q   D   E   K   R   H

CHARGED (carga neta a pH 7):
   D   E   → negativas / ácidas
   K   R   → positivas / básicas
   H       → básico solo a pH ≤6; neutro a pH 7
```

**Regla para BLOSUM:** dos aminoácidos que comparten grupo → sustitución conservativa → score positivo en BLOSUM62.

### Aminoácidos especiales con su rol mecanístico

**Glicina (G):** la glicina es el único aminoácido sin cadena lateral (solo un H). Eso le da libertad de rotación en el backbone que ningún otro aminoácido tiene. Aparece en loops y giros (turns) donde la cadena polipeptídica necesita doblar en un ángulo imposible para otro residuo. También en colágeno (repetición Gly-Pro-Hyp): el Gly va en el interior del triple hélice, donde no cabe ningún otro residuo.

**Prolina (P):** la prolina tiene su cadena lateral cíclica unida al nitrógeno del backbone, formando un anillo de pirrolidina de 5 miembros. Consecuencias: (1) no puede donar H al backbone → rompe la red de puentes H de la hélice α y la lámina β; (2) fuerza ángulos ψ/φ inusuales → introduce un kink (quiebre) en la cadena polipeptídica; (3) no va en la posición X del motivo N-glicosilación N-X-[S/T]: la oligosacaril-transferasa no reconoce el motivo si X es Pro. BLOSUM: scores negativos con prácticamente todos los demás aminoácidos.

**Cisteína (C):** el grupo tiol (–SH) permite: (1) puentes disulfuro — dos Cys pueden oxidarse para formar un puente covalente S–S, estabiliza proteínas extracelulares (anticuerpos, proteínas secretadas); (2) nucleofilia: el azufre es buen nucleófilo a pH fisiológico → sitios activos de proteasas de cisteína. Es el segundo aminoácido más raro (después de Trp).

**Histidina (H):** el único aminoácido cuyo pKa (≈6.0) cae en el rango fisiológico (pH 6–8). Puede actuar como donante o aceptor de protones (H⁺) a pH 7 → rol central en catálisis ácido-base. Aparece en la **tríada catalítica** de muchas hidrolasas: **Ser-His-Asp**. La His transfiere el H⁺ de la Ser al Asp, activando la Ser como nucleófilo para atacar el enlace peptídico. También coordina iones metálicos (Zn²⁺, Fe²⁺) en sitios activos metalo-enzimáticos.

**Triptófano (W):** el más raro de los 20 aminoácidos (~1% de todos los residuos). Casi siempre crítico funcionalmente: si un Trp está conservado en una familia proteica (aparece prominente en el logo de un HMM), es una señal fortísima de que es esencial. Score W↔W = +11 en BLOSUM62 (la señal evolutiva más fuerte de la matriz). Mayor absorbancia UV a 280 nm de todos los aminoácidos → domina la lectura de absorbancia para cuantificar proteína en laboratorio.

### Modificaciones post-traduccionales (PTMs)

Las PTMs son cambios químicos que ocurren luego de que la proteína fue sintetizada en el ribosoma. Regulan función, localización y vida media de la proteína.

| Modificación | Aminoácido(s) | Función biológica | Nota clave |
|---|---|---|---|
| **Fosforilación** | Ser (S), Thr (T), Tyr (Y) | Señalización celular: activa o inactiva enzimas, recluta proteínas. Reversible (fosfatasas lo quitan). | Las kinasas reconocen secuencias consenso. Ser-kinasas, Thr-kinasas, Tyr-kinasas son clases distintas. |
| **N-glicosilación** | Asn (N) en el motivo **N-{P}-[S/T]** | Plegamiento correcto en el RE, estabilidad térmica, reconocimiento celular. | La Pro en posición X impide la modificación → N-P-[S/T] NO se glicosilará. |
| **O-glicosilación** | Ser (S), Thr (T) | Regulación de transcripción; señalización alternativa a la fosforilación. | Compite con fosforilación en los mismos residuos Ser/Thr. |
| **Ubiquitinación** | Lys (K) | Marca proteínas para degradación por el proteasoma 26S. También regula transcripción (histonas). | La Lys tiene pKa ≈10.5 → cargada + a pH 7. |
| **Puentes disulfuro** | Cys–Cys | Estabilidad covalente en proteínas secretadas/extracelulares. Los anticuerpos (IgG) dependen de múltiples puentes S–S. | El par Cys–Cys en posiciones compatibles del alineamiento es un indicador estructural conservado. |
| **Acetilación** | Lys (histonas), Met N-terminal | En histonas: abre la cromatina, activa transcripción. En Met N-terminal: regula estabilidad de la proteína. | Las marcas epigenéticas (H3K27ac, etc.) son Lys acetiladas en colas de histonas. |

---


### Por qué el sistema match/mismatch no alcanza

Cuando uno arma un alineamiento con programación dinámica (NW o SW), necesita un número para cada par de aminoácidos alineados. La opción más cruda sería: match = +1, mismatch = −1. El problema es que eso trata todos los cambios como si fueran igualmente probables o igualmente dañinos, y la biología nos dice que eso es falso.

El razonamiento: pensá en dos proteínas que divergieron hace 50 millones de años. En ese tiempo hubo mutaciones, pero las que sobrevivieron son las que no destruyeron la función. Si una Leucina (hidrofóbica, apolar) muta a una Isoleucina (también hidrofóbica, también apolar), la proteína probablemente sigue funcionando igual. Esa sustitución va a aparecer muchas veces cuando uno mira alineamientos de proteínas reales. En cambio, si esa Leucina muta a un Ácido glutámico (cargado negativamente), el interior hidrofóbico se va a desestabilizar, la proteína probablemente se pierde, y esa mutación raramente se fija.

La idea entonces es: **observar qué sustituciones ocurren realmente en la evolución y codificarlas como puntaje**. La fórmula conceptual es:

```
score(a, b) = log [ P(a y b alineados en la naturaleza) / P(a y b por azar) ]
```

- Si el par aparece más de lo esperado al azar → logaritmo > 0 → score positivo.
- Si el par aparece menos de lo esperado al azar → logaritmo < 0 → score negativo.

Esto explica por qué Trp↔Trp tiene score 11: el Triptófano es rarísimo (~1% de todos los residuos), entonces encontrar dos Trp alineados es estadísticamente muy improbable al azar → señal evolutiva muy fuerte → score alto.

### Familia PAM — modelo evolutivo extrapolado

**PAM** significa *Point Accepted Mutation* (mutación puntual aceptada por la selección natural).

**Cómo se construye:** Margaret Dayhoff y colaboradores en los 70s agarraron alineamientos de proteínas que eran muy cercanamente relacionadas (más del 85% de identidad de secuencia). Con eso contaron, para cada par de aminoácidos (i, j), cuántas veces se observó que i fue reemplazado por j → matriz de probabilidades de sustitución para 1 PAM de distancia evolutiva.

**¿Qué es 1 PAM?** Es la cantidad de evolución necesaria para que, en promedio, 1 de cada 100 aminoácidos de una proteína cambie (1% de divergencia). Es una unidad de distancia evolutiva, no de tiempo.

**La extrapolación:** la matriz de 1 PAM se puede multiplicar por sí misma (multiplicación de matrices de Markov) para simular el efecto de múltiples pasos evolutivos. PAM250 = PAM1 multiplicada 250 veces. Eso modela proteínas que tienen 250 unidades de divergencia acumulada, equivalente a secuencias que comparten ~20% de identidad residual.

```
PAM1   → proteínas casi idénticas (>85% de identidad)
PAM40  → proteínas moderadamente relacionadas
PAM120 → proteínas con ~40% de identidad
PAM250 → proteínas muy divergentes (~20% de identidad)
```

### Familia BLOSUM — observación directa de bloques conservados

**BLOSUM** significa *BLOcks SUbstitution Matrix*. Henikoff & Henikoff (1992) tomaron la base de datos BLOCKS y contaron directamente las sustituciones observadas entre pares de secuencias. No hay extrapolación: se observan directamente las sustituciones que aparecen en bloques conservados reales.

El número en el nombre indica el **umbral de identidad** con el que se agruparon las secuencias:

```
BLOSUM80 → bloques con secuencias ≥80% idénticas → para proteínas muy similares
BLOSUM62 → bloques con ≥62% de identidad → el estándar de facto (default de BLAST)
BLOSUM45 → bloques con ≥45% de identidad → para proteínas más divergentes
```

**La confusión de los números (¡importante para el parcial!):** En PAM, número grande = más divergencia (PAM250 = muy lejano). En BLOSUM, número grande = más similitud de las secuencias usadas para construirla (BLOSUM80 = construida con secuencias muy parecidas). Son escalas opuestas.

```
Para proteínas CERCANAS:   BLOSUM80  o  PAM40
Para proteínas ALEJADAS:   BLOSUM45  o  PAM250
Default (la mayoría):      BLOSUM62
```

### La diferencia conceptual central entre PAM y BLOSUM

| | PAM | BLOSUM |
|---|---|---|
| Origen de los datos | Proteínas muy cercanas (>85% id) | Bloques conservados en familias diversas |
| Mecanismo | Modelo de Markov: multiplicación iterativa de PAM1 | Conteo estadístico directo |
| Escala del número | Número alto = más divergencia | Número alto = más similitud en los datos |
| Limitación | La extrapolación asume proceso evolutivo estacionario | Requiere muchas secuencias con bloques bien definidos |

En clase se enfatizó que BLOSUM es generalmente preferida porque no extrapola, sino que observa directamente qué sustituciones aparecen en familias proteicas reales.

### Cómo leer la matriz — ejemplos concretos con BLOSUM62

```
Par         Score   Razón biológica
───────────────────────────────────────────────────────────────────
W ↔ W       +11    Trp rarísimo, encontrar dos = señal evolutiva fuerte
P ↔ P       +7     Pro casi nunca cambia → su auto-match vale mucho
F ↔ Y       +3     Ambos aromáticos, Y tiene solo un OH extra
V ↔ I       +3     Ambos alifáticos, tamaños similares
L ↔ I       +2     Ambos alifáticos + hidrofóbicos
D ↔ E       +2     Ambos ácidos (carga −), E es D + un CH₂
K ↔ R       +2     Ambos básicos (carga +)
S ↔ T       +1     Ambos polares pequeños con grupo –OH
K ↔ E        0     Cargas opuestas, sin conservación de función
L ↔ D       −4     Hidrofóbico vs. ácido; grupos completamente distintos
P ↔ L       −3     Pro rompe estructura; sin similitud funcional
```

**El caso especial de Prolina:** En el parcial preguntaron "¿por qué Prolina no puede estar en el motivo de N-glicosilación?" y "¿qué dice la matriz BLOSUM de los cambios a Prolina?" La respuesta es que la Prolina tiene su cadena lateral unida al nitrógeno del backbone, formando un anillo rígido. Eso le impide donar hidrógeno en la red de puentes H de las estructuras secundarias, y fuerza un ángulo de torsión distinto al de todos los demás aminoácidos. Insertar una Pro donde había cualquier otro aminoácido es casi siempre estructuralmente disruptivo → sustitución rarísima en la evolución → **todos los scores de Pro contra otros aminoácidos son negativos**.

### Regla de oro para elegir matriz

```
Situación                                    Matriz recomendada
──────────────────────────────────────────────────────────────────
Comparación estándar / búsqueda en BLAST     BLOSUM62  (default)
Proteínas muy similares (>70% id)            BLOSUM80 o PAM40
Detección de homólogos lejanos (<30% id)     BLOSUM45 o PAM250
PSI-BLAST (iterativo)                        BLOSUM62 para la primera ronda;
                                             después PSSM sitio-específica
```

---


### ¿Qué es un alineamiento y para qué sirve?

Antes de meterse con los algoritmos, hay que entender *por qué* alinear secuencias en primer lugar.

El profe arrancó la clase con un experimento mental: imaginá que sos un biólogo marciano que llega a la Tierra y quiere estudiar la fauna. Tenés una muestra: un perro, un gato, un cocodrilo, una mosca, un mosquito. Un amigo te dice "tengo un tigre que no conozco, ¿qué me podés decir de él?". Sin saber nada del tigre, mirás su tamaño, contás sus patas, ves si tiene pelo, y automáticamente lo agrupás con el perro y el gato. Y entonces podés predecir: tiene cuatro patas, es mamífero, da de mamar a sus crías, tiene sangre caliente.

Eso es **inferencia por asociación** (el profe dice "dime con quién andas y te diré quién eres"). Lo mismo aplica a las secuencias biológicas. Si una secuencia nueva es muy parecida a secuencias cuya función ya conocés, podés inferir que tiene la misma función. El alineamiento es la herramienta que te permite medir ese parecido de manera rigurosa.

**El alineamiento es una comparación de secuencias.** Consiste en poner las secuencias en filas y tratar de que, en cada columna, queden los caracteres que se "corresponden" o son equivalentes entre sí. Así como un lingüista que alinea la palabra "sal" en castellano, francés, italiano y portugués nota que la primera letra y la última son iguales y que difieren solo en el medio (y puede separar claramente ese grupo de la palabra "azúcar" en los mismos idiomas), un biólogo que alinea secuencias de proteínas puede notar qué posiciones están conservadas y cuáles cambiaron.

---

### Dot plot: el alineamiento visual

El dot plot es el método más simple para visualizar la similitud entre dos secuencias. Es un método **gráfico**: no te da directamente el alineamiento escrito, sino que te muestra las regiones de similitud como una imagen.

#### Cómo se construye

1. Ponés una secuencia en el eje horizontal (columnas).
2. Ponés la otra secuencia en el eje vertical (filas).
3. En cada celda (i, j), ponés un punto **si el carácter i de la secuencia vertical coincide con el carácter j de la secuencia horizontal**.
4. Dejás vacío donde no hay coincidencia.

```
Ejemplo: alinear ATCG (horizontal) vs ACGT (vertical)

     A  T  C  G
A  [ .           ]
C  [    .        ]
G  [       .     ]
T  [  .          ]
```

#### Cómo se lee el dot plot

Lo importante es lo que pasa con la diagonal. Si las dos secuencias son idénticas, todos los puntos van a estar exactamente en la diagonal principal (de arriba-izquierda a abajo-derecha), formando una línea perfecta.

Cuando las secuencias son similares pero no idénticas:
- **Seguir por la diagonal** = hay un match.
- **Saltar la diagonal un paso hacia la derecha** = estás metiendo un **gap en la secuencia vertical**. Tuviste que saltear un carácter de la horizontal para poder seguir encontrando matches.
- **Saltar la diagonal un paso hacia abajo** = gap en la secuencia horizontal.

```
Ejemplo:
Sec. 1 (horizontal): A T C G A
Sec. 2 (vertical):   A C G A

     A  T  C  G  A
A  [ .           . ]
C  [    .          ]
G  [       .       ]
A  [           .   ]

Alineamiento que surge: seguir la diagonal con un gap en la T.

Sec1:  A T C G A
Sec2:  A - C G A   (gap donde estaba la T)
```

#### Cómo interpretar distintos patrones en el dot plot

| Patrón en el dot plot | Qué significa biológicamente |
|---|---|
| Diagonal perfecta de extremo a extremo | Las dos secuencias son idénticas |
| Diagonal continua con algunos huecos | Secuencias muy similares, pequeñas variaciones |
| Diagonal con saltos (quiebres) | Hay gaps / inserciones / deleciones |
| Diagonal solo en una región | Solo parte de la secuencia alinea → candidato a alineamiento local |
| Diagonal en la anti-diagonal | Una secuencia es el reverso de la otra (palíndromo) |
| El mismo segmento aparece dos veces en la diagonal | Hubo una duplicación interna |
| Sin diagonal reconocible | Las secuencias no se parecen globalmente |

#### El filtro del dot plot (window size)

Si las secuencias son largas y tienen muchos matches aislados por azar, el dot plot queda lleno de puntos y no se puede leer nada. Por eso se aplica un **filtro de ventana (window size)**:
- Solo se dibuja un punto si hay, dentro de una ventana de `w` posiciones, al menos `k` matches consecutivos.
- Ejemplo: ventana = 5, mínimo 4 matches → solo se conservan segmentos donde hay 4 o más matches en 5 posiciones consecutivas.

Un filtro muy estricto puede ocultar similitud real; uno muy laxo deja pasar ruido.

---

### Sistema de puntuación: ¿cómo sé cuál alineamiento es mejor?

Cuando el dot plot muestra varias diagonales posibles, necesitás una manera de decidir cuál es mejor. Eso se hace con una **función de puntuación (scoring function)**.

#### El problema con "match = 1, mismatch = 0"

Mirá este ejemplo: tenés dos alineamientos alternativos de las mismas secuencias proteicas. Uno alinea la tirosina con el triptófano (dos aminoácidos raros y funcionalmente similares). El otro alinea serina con alanina (aminoácidos comunes). Con una función simple (match = +1, mismatch = 0), los dos dan el **mismo puntaje**. Pero biológicamente, alinear triptófano con tirosina es mucho más significativo: son aminoácidos raros, y si aparecen en la misma posición en dos proteínas distintas, probablemente cumplan un rol funcional conservado.

Conclusión: para proteínas, necesitás una función que refleje el **grado de parecido entre aminoácidos**, no solo si son iguales o no. Ver Clase 2 para las matrices PAM y BLOSUM.

**Para DNA:** en general basta con match = +1, mismatch = -1. Las cuatro bases no tienen tanto "grado de parecido" entre sí como los 20 aminoácidos.

---

### Penalización por gaps: lineal vs affine

#### Por qué los gaps son especiales

Un gap en el alineamiento representa una inserción o deleción ocurrida a lo largo de la evolución. Lo importante biológicamente es que **una sola deleción puede eliminar varios nucleótidos o aminoácidos contiguos en un solo evento mutacional**.

> Una deleción de 10 bases es **un evento** evolutivo, no diez eventos independientes.

Entonces, tiene más sentido biológico permitir **un gap largo** que forzar **muchos gaps chicos**. Un gap largo (de 10 posiciones) probablemente refleja una sola deleción real.

#### Gap lineal (el más simple)

```
W(k) = d × k       (d = costo por posición, k = largo del gap)

Ejemplo: d = -2, gap de longitud 3  → costo = -6
         d = -2, tres gaps de largo 1 → costo total = -6

Con gap lineal, da lo mismo 1 gap de 3 que 3 gaps de 1.
No refleja la biología.
```

#### Gap affine (el realista)

```
W(k) = a + b × k

donde:
  a = penalización por ABRIR el gap    (costo fijo, se paga solo la primera vez)
  b = penalización por EXTENDER el gap (costo por cada posición adicional)

Biológicamente: a >> b
  → abrir un gap es caro (implica un evento de inserción/deleción)
  → extender el gap es barato (agregar una posición más al mismo evento)
```

**Ejemplo típico para DNA:** a = -5 (abrir), b = -1 (extender)

```
Un gap de 1 posición:     -5 + (-1×1) = -6
Un gap de 5 posiciones:   -5 + (-1×5) = -10
Cinco gaps de 1 posición: 5 × (-6)    = -30

→ Un gap de 5 es mucho más barato que cinco gaps de 1.
→ El algoritmo va a preferir los gaps largos concentrados.
```

Con gap affine: venir desde la izquierda (misma fila) = extender gap (costo b). Venir desde la diagonal = abrir gap nuevo (costo a).

---

### Programación dinámica: NW y SW

#### El problema computacional

Para encontrar el mejor alineamiento entre dos secuencias, en principio habría que probar todos los alineamientos posibles y quedarse con el de mayor puntaje. La cantidad de alineamientos posibles crece de manera explosiva → inviable para secuencias reales.

La solución es **programación dinámica**: en vez de resolver el problema grande directamente, lo dividís en problemas más chicos y vas construyendo la solución de a poco.

La clave es el **principio de optimalidad de Bellman**: si el alineamiento óptimo de dos secuencias completas pasa por alinear los primeros `i` caracteres de la primera con los primeros `j` de la segunda de cierta manera, entonces esa sub-alineación también tiene que ser óptima para esos prefijos. Eso garantiza que se puede construir la solución óptima global a partir de sub-soluciones óptimas locales.

---

### Needleman-Wunsch (NW): alineamiento global

#### Cuándo usarlo

Cuando las dos secuencias son de longitud similar y querés alinearlas **de punta a punta**. Por ejemplo: comparar el gen de una proteína en dos especies donde ambas secuencias corresponden al gen completo.

#### Lógica de la recurrencia

Construís una matriz `F` de tamaño `(n+1) × (m+1)` donde `n` = largo de la secuencia X (vertical) y `m` = largo de la secuencia Y (horizontal).

Para llegar a la celda `F(i,j)`, solo tenés tres opciones de dónde venir:

```
1. Diagonal ↖: venías alineando caracteres. Alineaste xᵢ con yⱼ.
   Puntaje tentativo: F(i-1, j-1) + s(xᵢ, yⱼ)   [s = score de la sustitución]

2. Desde arriba ↑: pusiste un gap en la secuencia Y (horizontal).
   Puntaje tentativo: F(i-1, j) + gap_penalty

3. Desde la izquierda ←: pusiste un gap en la secuencia X (vertical).
   Puntaje tentativo: F(i, j-1) + gap_penalty
```

Te quedás con el máximo de las tres:

```
F(i,j) = max {
    F(i-1, j-1) + s(xᵢ, yⱼ)    ← diagonal (match/mismatch)
    F(i-1, j)   + gap_penalty    ← arriba    (gap en Y)
    F(i, j-1)   + gap_penalty    ← izquierda (gap en X)
}
```

#### Inicialización

La primera fila y la primera columna representan alinear contra nada (puro gap):

```
F(0, 0) = 0
F(i, 0) = i × gap_penalty    para todo i > 0
F(0, j) = j × gap_penalty    para todo j > 0
```

Con gap lineal = -1: la primera fila es 0, -1, -2, -3, ... y la primera columna igual.

#### El traceback

Al llenar la matriz, en cada celda anotás también **de dónde viniste** (diagonal, arriba o izquierda). Esto forma la **matriz de traceback**.

Una vez llenada la matriz:
1. Empezás en la celda `F(n, m)` (esquina inferior derecha).
2. Seguís las flechas hacia atrás hasta llegar a `F(0, 0)`.
3. Cada paso del traceback construye el alineamiento al revés:
   - Flecha diagonal ↖: alineaste xᵢ con yⱼ (match o mismatch).
   - Flecha arriba ↑: gap en Y (ponés `-` en la secuencia horizontal, avanzás solo en la vertical).
   - Flecha izquierda ←: gap en X (ponés `-` en la secuencia vertical, avanzás solo en la horizontal).

---

### Ejemplo resuelto paso a paso: NW

**Secuencias:**
```
X (vertical):   A C G T
Y (horizontal): A G T
```

**Parámetros:** match = +1, mismatch = -1, gap lineal = -1

**Paso 1: Inicialización**

```
        ""   A    G    T
""    [  0   -1   -2   -3 ]
A     [ -1    .    .    . ]
C     [ -2    .    .    . ]
G     [ -3    .    .    . ]
T     [ -4    .    .    . ]
```

**Paso 2: Llenado celda a celda**

`F(1,1)`: alinear A (de X) con A (de Y).
- Diagonal: F(0,0) + s(A,A) = 0 + 1 = **+1**
- Arriba: F(0,1) + gap = -1 + (-1) = -2
- Izquierda: F(1,0) + gap = -1 + (-1) = -2
→ `F(1,1) = +1` (viene de diagonal ↖)

`F(1,2)`: alinear A (de X) con G (de Y).
- Diagonal: F(0,1) + s(A,G) = -1 + (-1) = -2
- Arriba: F(0,2) + gap = -2 + (-1) = -3
- Izquierda: F(1,1) + gap = +1 + (-1) = **0**
→ `F(1,2) = 0` (viene de izquierda ←)

`F(2,1)`: alinear C (de X) con A (de Y).
- Diagonal: F(1,0) + s(C,A) = -1 + (-1) = -2
- Arriba: F(1,1) + gap = +1 + (-1) = **0**
- Izquierda: F(2,0) + gap = -2 + (-1) = -3
→ `F(2,1) = 0` (viene de arriba ↑)

`F(2,2)`: alinear C (de X) con G (de Y).
- Diagonal: F(1,1) + s(C,G) = +1 + (-1) = **0**
- Arriba: F(1,2) + gap = 0 + (-1) = -1
- Izquierda: F(2,1) + gap = 0 + (-1) = -1
→ `F(2,2) = 0` (viene de diagonal ↖)

`F(3,2)`: alinear G (de X) con G (de Y).
- Diagonal: F(2,1) + s(G,G) = 0 + 1 = **+1**
- Arriba: F(2,2) + gap = 0 + (-1) = -1
- Izquierda: F(3,1) + gap = -1 + (-1) = -2
→ `F(3,2) = +1` (viene de diagonal ↖)

`F(4,3)`: alinear T (de X) con T (de Y).
- Diagonal: F(3,2) + s(T,T) = +1 + 1 = **+2**
- Arriba: F(3,3) + gap = 0 + (-1) = -1
- Izquierda: F(4,2) + gap = 0 + (-1) = -1
→ `F(4,3) = +2` (viene de diagonal ↖)

**Matriz completa:**

```
        ""   A    G    T
""    [  0   -1   -2   -3 ]
A     [ -1   +1    0   -1 ]
C     [ -2    0    0   -1 ]
G     [ -3   -1   +1    0 ]
T     [ -4   -2    0   +2 ]
```

**Paso 3: Traceback** (empezando en F(4,3) = +2)

1. `F(4,3)` → diagonal ↖ → alineamos **T** (X) con **T** (Y) → match
2. `F(3,2)` → diagonal ↖ → alineamos **G** (X) con **G** (Y) → match
3. `F(2,1)` → arriba ↑ → gap en Y → alineamos **C** (X) contra **-**
4. `F(1,1)` → diagonal ↖ → alineamos **A** (X) con **A** (Y) → match

**El alineamiento resultante:**

```
X:   A  C  G  T
Y:   A  -  G  T
     |     |  |
  match  gap  match  match

Score = +1 (A-A) + (-1) (gap en C) + +1 (G-G) + +1 (T-T) = +2  ✓
```

---

### Smith-Waterman (SW): alineamiento local

#### Cuándo usarlo

Cuando las secuencias son de longitud diferente, o cuando solo comparten una región (un dominio). Por ejemplo: una proteína que tiene un dominio conocido incrustado en una cadena mucho más larga.

Hoy en día, el alineamiento local es el más usado en la práctica. Si hay regiones que no se parecen, el algoritmo global las fuerza a alinearse de todos modos (con puntaje muy negativo), ensuciando el resultado. Con SW, esas regiones directamente no aparecen.

#### La única diferencia con NW: el piso de cero

```
F(i,j) = max {
    0                             ← NUNCA baja de cero
    F(i-1, j-1) + s(xᵢ, yⱼ)
    F(i-1, j)   + gap_penalty
    F(i, j-1)   + gap_penalty
}
```

El cero actúa como un "reset": cuando el alineamiento acumulado se vuelve tan malo que baja de cero, es mejor olvidar todo y empezar un alineamiento nuevo desde cero que seguir arrastrando un puntaje negativo.

#### Comparación directa NW vs SW

```
                  NW (global)               SW (local)
────────────────────────────────────────────────────────────
Inicialización    F(i,0) = i × gap           F(i,0) = 0
primera fila/col  F(0,j) = j × gap           F(0,j) = 0

Recurrencia       max(↖, ↑, ←)              max(0, ↖, ↑, ←)

Inicio traceback  Esquina (n,m)              Celda con valor máximo

Fin traceback     Celda (0,0)                Celda con valor 0

Resultado         Toda X alineada            Mejor región en común
                  contra toda Y              entre cualquier parte de X e Y
```

---

### Múltiples alineamientos óptimos

A veces al llenar la matriz, dos o más opciones de la recurrencia dan exactamente el mismo valor en la misma celda. Ejemplo:

```
F(i-1, j-1) + s(xᵢ, yⱼ) = F(i-1, j) + gap_penalty
```

En ese caso, hay un **empate**: el algoritmo anota **ambas flechas** en la matriz de traceback.

Cuando hacés el traceback y llegás a una celda con dos flechas, el camino se **bifurca**: podés seguir por cualquiera de los dos caminos, y cada uno va a darte un alineamiento diferente, pero **con el mismo score total**.

**Ejemplo conceptual:**
```
X:  A A T
Y:  A T

Con gap = -1, match = +1, mismatch = -1:

Alineamiento 1:   A A T      Score: +1 (A-A) + (-1)(gap) + +1(T-T) = +1
                  A - T

Alineamiento 2:   A A T      Score: (-1)(gap) + +1(A-A) + +1(T-T) = +1
                  - A T
```

Para encontrarlos todos:
1. Llenás la matriz de scores normalmente.
2. En el traceback, cada vez que hay empate → bifurcás el camino.
3. Cada camino completo es un alineamiento óptimo distinto.

### Cómo leer el alineamiento desde el traceback

| Movimiento en el traceback | Columna del alineamiento |
|---|---|
| Diagonal ↖ | xᵢ alineado con yⱼ (match si son iguales, mismatch si difieren) |
| Arriba ↑ | xᵢ alineado con `-` (gap en la secuencia horizontal Y) |
| Izquierda ← | `-` alineado con yⱼ (gap en la secuencia vertical X) |

### Complejidad computacional

Tanto NW como SW tienen:
- **Tiempo:** O(n × m) — se calcula cada celda de la matriz una sola vez.
- **Espacio:** O(n × m) — hay que guardar toda la matriz para el traceback.

Para secuencias de n = m = 1000, eso es ~1.000.000 operaciones. Para secuencias genómicas de millones de bases, se usan variantes heurísticas como BLAST.

### Resumen visual del flujo completo

```
Tengo dos secuencias
        ↓
¿Son similares en largo y quiero alinear todo?
    ├─ SÍ → Needleman-Wunsch (global)
    └─ NO → Smith-Waterman (local)
        ↓
Elijo función de scoring:
  · DNA: match/mismatch simples + penalización por gap
  · Proteínas: matriz de sustitución (BLOSUM 62 por default)
        ↓
Elijo tipo de gap:
  · Gap lineal: W(k) = d × k  (simple, no diferencia 1 gap largo de muchos cortos)
  · Gap affine: W(k) = a + b×k  (realista: abrir es caro, extender es barato)
        ↓
Lleno la matriz F de programación dinámica
  → guardo también la matriz de traceback (de dónde viene cada celda)
        ↓
Traceback:
  NW: empiezo en (n,m), termino en (0,0)
  SW: empiezo en el máximo global, termino en una celda con valor 0
        ↓
Leo el alineamiento:
  diagonal ↖ → match / mismatch
  arriba   ↑ → gap en secuencia horizontal
  izquierda← → gap en secuencia vertical
        ↓
¿Hay empates en el traceback?
  SÍ → múltiples alineamientos óptimos (mismo score, distinta forma)
  NO → un único alineamiento óptimo
```

---

## Clase 3 — BLAST y Alineamiento Múltiple

### La intuición fundamental: buscar una secuencia es como buscar en un libro gigante

Antes de ver el algoritmo, hay que entender *por qué* existe BLAST y qué problema conceptual resuelve.

Supongamos que tenemos un fragmento de 50–100 aminoácidos de una proteína que acabamos de obtener experimentalmente y queremos saber qué proteína es. La base de datos tiene 10 millones de proteínas. No sabemos el nombre, solo tenemos la secuencia. Entonces necesitamos buscar en el **campo secuencia**.

Cuando buscamos "hemoglobin" en el campo nombre de una BD, la base de datos consulta un **índice** (una lista ordenada, como el índice de un libro) y encuentra las entradas que coinciden. Eso es rápido. El problema con las secuencias es que no puedo indexarlas como texto, porque lo que me interesa no es la igualdad exacta sino la **similitud**. Necesito comparar mi secuencia contra cada secuencia del campo secuencia usando **alineamiento**, y eso es O(n·m) por par → inviable para 10M de registros.

> **Analogía:** es como si en Google quisieras buscar una *frase* en todos los sitios web, pero la frase puede tener pequeñas variaciones. Google no puede revisar cada página; tiene que tener pre-indexado algo más fino.

La solución de BLAST: indexar las secuencias de la BD **en trozos pequeños (words)**, buscar esos trozos exactos en el índice, y solo desde esos hits exactos extender el alineamiento. Transforma O(n·m) por registro a casi O(m) en total.

---

### Los 4 pasos de BLAST con visualización detallada

#### Paso 1 — Construcción de la word list (fragmentar la query)

La query se fragmenta en **words** de tamaño fijo (word size w=3 para proteínas, w=11 para nucleótidos).

```
Query de 13 aminoácidos:
  A  L  V  G  T  T  Y  H  H  V  D  R  R

Words generadas (w=3):
  pos 1: ALV
  pos 2: LVG
  pos 3: VGT
  pos 4: GTT
  pos 5: TTY
  pos 6: TYH
  pos 7: YHH
  pos 8: HHV
  pos 9: HVD
  pos 10: VDR
  pos 11: DRR
  → 11 words para una query de largo 13 (m−w+1)
```

Para proteínas, BLAST además genera **vecinas** de cada word: todas las secuencias de largo w cuyo score BLOSUM62 contra la word original sea mayor que un umbral T (típicamente T=11).

```
Vecinas de "ALV" (score ≥ 11 con BLOSUM62):
  ALV, ALI, ALA, ALM, AIV, AVI, ...
  → unas ~50 vecinas para w=3, T=11
```

**Por qué tiene sentido:** si la proteína buscada tiene LEU donde la query tiene VAL (cambio conservativo), la word "ALI" pega en el índice donde "ALV" no hubiera pegado. El umbral T controla el balance: T alto → lista corta → más rápido pero menos sensible; T bajo → más lento pero detecta más homólogos distantes.

#### Paso 2 — Búsqueda de hits en el índice pre-construido

La BD tiene pre-indexadas todas sus secuencias en words. El índice es una tabla que, para cada posible word de largo w (20³ = 8.000 entradas para w=3), lista todos los registros que la contienen y en qué posiciones.

```
ÍNDICE PRE-CONSTRUIDO (fragmento):

Word  | Registro | Posición en registro
------|----------|---------------------
ALI   | reg_0472 | pos 34, 118
ALI   | reg_1203 | pos 7
ALV   | reg_0472 | pos 35
ALV   | reg_5503 | pos 12
```

BLAST busca cada word de la word list en ese índice → obtiene hits **(registro, posición_en_BD, posición_en_query)**.

> **Clave:** este índice se calcula **una sola vez** (cuando se actualiza la BD) y se reutiliza en cada búsqueda. Por eso una búsqueda de 300 aa contra 1M de registros tarda solo 10–30 segundos.

Los registros que aparecen repetidamente para muchas words son candidatos a contener una región similar a la query.

#### Paso 3 — Extensión sin gaps (Ungapped Extension)

Los hits se posicionan usando sus coordenadas para armar un dot-plot interno. Los hits en la misma diagonal (posición_BD − posición_query = constante) pertenecen al mismo alineamiento potencial. BLAST extiende a izquierda y derecha desde cada hit **sin introducir gaps**, sumando scores BLOSUM62:

```
← extiende       HIT      extiende →
query: ...  L  V [A  L  V] G  T  T  ...
BD:    ...  L  V [A  L  L] G  T  T  ...
           acumula score mientras los residuos sean parecidos
```

La extensión se **detiene** cuando el score acumulado cae por debajo de un umbral X₀. La región resultante con score ≥ S es un **HSP (High-Scoring Segment Pair)**.

```
Score acumulado durante la extensión:
  4  8  12  16  14  12  10  8  3  [cae bajo umbral: CORTE]
              ↑ punto máximo del HSP
```

**Por qué sin gaps primero:** hacer Smith-Waterman completo para cada hit sería demasiado costoso. La extensión sin gaps es O(largo del HSP) y filtra la enorme mayoría del ruido.

#### Paso 4 — Extensión con gaps y cálculo de E-value

Solo para los HSPs que superaron el umbral S, BLAST hace un **alineamiento local tipo Smith-Waterman** (ahora con gaps) centrado en el HSP.

```
ANTES (extensión sin gaps):
  query: A L V G T T Y H
  BD:    A L I G T T Y H

DESPUÉS (con gaps):
  query: A L V G T - T Y H
         | | . | |   | | |
  BD:    A L I G T K T Y H
         ↑ gap insertado en la query
```

**Nota sobre múltiples alineamientos por registro:** BLAST es *local*, puede reportar **más de un alineamiento por el mismo registro** si hay dos regiones de similitud separadas por una región mala. BLAST = Basic **Local** Alignment Search Tool.

---

### El E-value: qué significa realmente

**El recorrido conceptual:**

1. **Score crudo (raw score):** suma de valores BLOSUM62 + penalidades de gap. No está normalizado.

2. **Bit-score:** normalización del score crudo con los parámetros λ y K de la distribución estadística. Permite comparar alineamientos de distinto largo.

3. **P-value:** probabilidad de obtener ese score o mayor por azar, dado el modelo estadístico. Problema: cuando busco en una BD con millones de registros, hago millones de comparaciones. Si el p-value es 0.01 y hay 1M de registros, obtengo 10.000 hits por azar.

4. **E-value:** corrección del p-value por el tamaño de la BD. Es el **número esperado de hits con ese score o mayor que obtendrías por azar** en una BD de ese tamaño.

```
Fórmula: E = K · m · N · e^(−λS*)

  m  = largo de la query (residuos)
  N  = tamaño total de la BD (número total de residuos)
  S* = bit-score
  K, λ = constantes dependientes de la matriz y parámetros de gap
```

```
Interpretación práctica:

  E = 1e-50  →  imposible por azar → hit muy significativo
  E = 1e-5   →  1 hit por azar cada 100.000 búsquedas → bueno
  E = 0.01   →  1 hit por azar cada 100 búsquedas → aceptable para homólogos distantes
  E = 1      →  1 hit por azar en ESTA búsqueda → dudoso
  E = 10     →  10 hits por azar → ruido, descartar
```

**Por qué el E-value depende del tamaño de la BD:**

El mismo alineamiento con el mismo score tendrá distinto E-value en BDs distintas:

```
  BD pequeña (N = 10^6 residuos):  E = K · m · 10^6 · e^(-λS)
  BD grande  (N = 10^9 residuos):  E = K · m · 10^9 · e^(-λS)  ← 1000 veces mayor (peor)
```

**Consecuencia práctica:** si restringís la búsqueda de BLAST a un organismo o taxón, la BD efectiva es más chica → el mismo hit tendrá E-value más favorable (menos ruido). Restricción de BD ahorra tiempo (menos registros) y mejora estadística (menos ruido).

**Cuándo usar cada métrica:**
- Score crudo: solo para analizar la matriz de sustitución.
- Bitscore: para comparar alineamientos entre sí (sin búsqueda en BD).
- E-value: siempre cuando se hace una búsqueda en BD y querés saber si el hit es significativo.

---

### Blast2Sequences vs búsqueda en base de datos

| | Búsqueda en BD | Blast2Sequences |
|---|---|---|
| **Qué hace** | Compara 1 query vs millones de secuencias | Compara 2 secuencias específicas entre sí |
| **Para qué** | Descubrir homólogos desconocidos | Comparar proteínas ya elegidas |
| **E-value** | Tiene sentido estadístico pleno | No tiene el mismo sentido (BD de 1 secuencia) |
| **Análogo** | Búsqueda en Google | Comparar dos documentos específicos |

Blast2Sequences es equivalente a hacer un alineamiento tipo Smith-Waterman entre dos secuencias que ya elegiste.

---

### PSI-BLAST: BLAST iterativo con matriz posición-específica

#### El problema que resuelve

BLAST estándar usa BLOSUM62 igual en todas las posiciones. Pero en una familia de proteínas hay posiciones muy conservadas (sitio activo, contacto con ligando) y posiciones variables (superficies expuestas). PSI-BLAST aprende esa heterogeneidad y es más sensible.

```
Posición 44 del MSA (catálisis):
  Seq1: D     Seq2: D     Seq3: D     Seq4: D
  → 100% conservada. BLOSUM62 da el mismo score a D→N aquí que en otra posición.
  → PSI-BLAST hace ese cambio MUY costoso en esta posición.

Posición 97 del MSA (loop expuesto):
  Seq1: A     Seq2: S     Seq3: T     Seq4: G
  → Muy variable. PSI-BLAST hace cualquier cambio BARATO aquí.
```

#### El ciclo iterativo ronda 1 → PSSM → ronda 2

```
RONDA 1:
  query → BLAST con BLOSUM62 genérica
        → hits con E < umbral (default 0.001)
        → conjunto inicial de homólogos (alta confianza)

CONSTRUCCIÓN DEL MSA:
  Alinear todas las secuencias homólogas entre sí

       pos1  pos2  pos3  pos4  ...
  seq1:  A     C     D     K
  seq2:  A     -     D     K
  seq3:  S     C     D     R
  seq4:  A     C     E     K

CONSTRUCCIÓN DE LA PSSM:
  Para cada posición i y cada aminoácido a:
    PSSM[i][a] = log(frec_observada(a en columna i) / frec_esperada_azar(a))

  → Posición conservada: score alto para el aa correcto, muy negativo para otros
  → Posición variable:   scores similares para muchos aminoácidos

  Resultado: UNA MATRIZ DIFERENTE por posición (no una sola BLOSUM62 global)

            A    C    D    E    F    G    ...
  pos 1  [ +2   -3   -2   -2   -4   +1   ...]  ← A y S comunes
  pos 2  [ -2  +4   -3   -3   -5   -3   ...]  ← C muy conservada
  pos 3  [ -3   -4   +5   -1   -5   -4   ...]  ← D muy conservada

RONDA 2:
  query → buscar en BD con PSSM[i][a] en lugar de BLOSUM62
        → detecta homólogos remotos que ronda 1 perdió
        → agregar nuevos hits al MSA (si E < umbral)
        → recalcular PSSM → RONDA 3, RONDA 4, ... hasta convergencia
```

#### El peligro del PSI-BLAST drift

Si un **falso positivo** entra al MSA en alguna ronda, contamina la PSSM → en la siguiente ronda la PSSM distorsionada encuentra más proteínas parecidas al FP → entran más FP → la PSSM se desvía progresivamente de la familia original.

```
Ronda 1: query ─PSSM₁─→ {seqs familia + 1 FP}
                                              ↓
Ronda 2:       ─PSSM₂─→ {seqs familia + 5 FP}   ← PSSM contaminada
                                              ↓
Ronda 3:       ─PSSM₃─→ {seqs alejadas + muchos FP} ← "drift"
```

**Cómo evitarlo:** umbral más estricto (E < 1e-5), revisión manual de cada ronda, no hacer más rondas de las necesarias (2–3 suele ser suficiente).

> **Lo que importa retener:** la PSSM es una **matriz diferente por posición** (no una sola BLOSUM62 global). Eso es exactamente lo que se llama **perfil** y es el punto de partida para los HMMs de la clase siguiente.

---

### Tipos de BLAST: razonamiento biológico de cuándo usar cada uno

```
              QUERY                    BD
blastp:       proteína         →       proteínas
blastn:       nucleótidos      →       nucleótidos
blastx:       nucl (6 marcos)  →       proteínas
tblastn:      proteína         →       nucl (6 marcos)
tblastx:      nucl (6 marcos)  →       nucl (6 marcos)
```

**blastp:** el uso más común. Buscás homólogos de una proteína conocida. El más sensible para detectar homología proteica porque compara directamente residuos.

**blastn:** útil para buscar secuencias genómicas muy parecidas (mismo gen en especie muy cercana). Menos sensible para homólogos distantes porque el alfabeto de solo 4 letras produce muchas coincidencias por azar, y los nucleótidos divergen más rápido que las proteínas.

**blastx:** tenés una secuencia de ADN nueva (fragmento de PCR, contig de novo) y querés saber si codifica alguna proteína conocida. BLAST traduce la query en **6 marcos de lectura** y busca cada traducción en la BD de proteínas.

```
Marco +1:  ATG GCT TTA → M  A  L
Marco +2:   TGG CTT TAA → W  L  *
Marco +3:    GGC TTT AA → G  F  ...
Marco −1, −2, −3: hebra complementaria inversa en 3 marcos
```

**tblastn:** tenés una proteína conocida y querés buscarla en un genoma sin anotar. La BD de nucleótidos se traduce on-the-fly en 6 marcos y se compara con la query proteica. Imprescindible cuando el genoma target no tiene proteínas anotadas todavía.

**tblastx:** compara dos genomas no anotados entre sí. Ambos, query y BD, se traducen en 6 marcos. Muy costoso computacionalmente, se usa raramente.

**Regla de sensibilidad (detectar homólogos remotos):**
```
blastp > blastx > tblastn
```
Por eso cuando tenés secuencia de ADN pero querés encontrar proteínas homólogas, usás blastx: la comparación a nivel proteico (20 letras) es mucho más informativa que la de nucleótidos (4 letras), y los cambios sinónimos del ADN no oscurecen la señal.

---


### Por qué el MSA es fundamental

Un alineamiento de pares nos dice si dos secuencias son similares. Para entender una **familia de proteínas** necesitamos ver muchas secuencias juntas. El MSA permite:

1. **Identificar posiciones conservadas:** si en una columna del MSA todas las secuencias tienen el mismo aminoácido, esa posición está bajo presión selectiva → probablemente funcional (sitio activo, contacto con ligando).

2. **Construir perfiles (PSSM):** el MSA es la base para construir matrices sitio-específicas usadas en PSI-BLAST y HMMs.

3. **Construir árboles filogenéticos:** las distancias evolutivas entre pares se calculan a partir del MSA y permiten reconstruir relaciones entre organismos.

4. **Transferir anotación funcional:** si una posición está conservada en la familia y en una proteína conocida es el sitio catalítico, podés inferir que esa posición es funcional en todas las demás.

5. **Entrada para HMMs:** los Hidden Markov Models de secuencia (Clase 4) se construyen a partir de un MSA y permiten búsquedas mucho más sensibles que BLAST.

**El problema de escala:** el MSA óptimo por programación dinámica tiene complejidad O(L^N) donde L es el largo y N el número de secuencias. Para N=10, L=300: ~10^24 operaciones → inviable. Se necesitan heurísticas.

---

### CLUSTALW: el algoritmo progresivo paso a paso con razonamiento

CLUSTALW es el algoritmo de MSA más clásico. Es **progresivo**: en lugar de alinear todo junto, agrega secuencias de a una (o grupos), siguiendo un árbol de similitud.

**Razonamiento detrás:** si dos secuencias son muy similares, su alineamiento va a ser casi perfecto desde el inicio. Tiene sentido empezar por lo que sabemos bien (pares más similares) y construir hacia lo más incierto (pares más distantes).

```
Analogía: es como ordenar cronológicamente fotos → empezás agrupando
las del mismo año y luego decidís cómo encajan los grupos entre sí.
```

#### Paso 1 — Alineamiento todos-contra-todos y matriz de similitud

Se alinean todos los pares de secuencias (Smith-Waterman o similar) y se calcula el bitscore de cada par. Con N secuencias hay N(N−1)/2 pares.

```
Para 4 secuencias: 4×3/2 = 6 comparaciones
Para 8 secuencias: 8×7/2 = 28 comparaciones
```

Resultado: matriz de scores (ejemplo):

```
MATRIZ DE SCORES (bitscore):
        Prot1   Prot2   Prot3   Prot4
Prot1     —      18      72      15
Prot2    18       —      14      81
Prot3    72      14       —      11
Prot4    15      81      11       —

Pares más similares: (P2, P4) = 81 y (P1, P3) = 72
Pares más distantes: (P3, P4) = 11
```

#### Paso 2 — Construcción del árbol guía (Neighbor-Joining)

A partir de la matriz de distancias se construye un árbol que une los pares más cercanos primero:

```
ÁRBOL GUÍA:

        ┌─── Prot2
   ┌────┤                  ← nodo A: los más similares (P2+P4, score 81)
   │    └─── Prot4
───┤
   │    ┌─── Prot1
   └────┤                  ← nodo B: segundo par más similar (P1+P3, score 72)
        └─── Prot3

Los dos nodos se unen por el mayor score inter-bloque = score(P2,P1) = 18.
```

#### Paso 3 — Alineamiento progresivo siguiendo el árbol

**Regla fundamental: una vez hecho un alineamiento parcial, NO SE TOCA.**

```
Iteración 1: Alinear Prot2 con Prot4 (par más similar, score=81):
  Prot2: M A I K - T L
  Prot4: M A I K G T L   ← gap insertado en Prot2 pos 5. FIJO.

Iteración 2: Alinear Prot1 con Prot3 (siguiente par, score=72):
  Prot1: M - V K T L
  Prot3: M C V K T L     ← gap en Prot1 pos 2. FIJO.

Iteración 3: Alinear bloque {Prot2, Prot4} contra bloque {Prot1, Prot3}
  → Se alinea el representante de mayor score inter-bloque: P2 vs P1 (score=18)
  → Los gaps que recibe P2 se propagan a P4 (misma fila del bloque A)
  → Los gaps que recibe P1 se propagan a P3 (misma fila del bloque B)

MSA FINAL:
  Prot2: M  A  I  K  -  T  L
  Prot4: M  A  I  K  G  T  L
  Prot1: M  -  -  K  -  T  L
  Prot3: M  C  -  K  -  T  L
         c1 c2 c3 c4 c5 c6 c7
```

**¿Por qué el alineamiento no se toca?** El algoritmo es greedy (avaro): toma la mejor decisión local en cada paso pero no tiene visión global. Una vez tomada una decisión, no la revisa. Esto hace posible terminar en tiempo razonable, pero los errores tempranos se propagan.

---

### Score de un bloque vs otro: promedio de pares

Para decidir si una secuencia nueva se pega a un bloque existente o forma una nueva rama, se compara el score de la secuencia contra cada miembro del bloque y se usa el **promedio**:

```
Ejemplo: tengo bloque {1, 2, 3, 4} y secuencia 7. ¿La 7 se une al bloque o a la 9?

  Score(7 vs 9) = 14

  Score(7 vs bloque {1,2,3,4}):
    Score(7,1) = 3
    Score(7,2) = 5
    Score(7,3) = 2
    Score(7,4) = 14
    Promedio = (3+5+2+14)/4 = 6

  6 < 14  →  la 7 se une a la 9 (mayor score = más similar = une primero)
```

---

### SP score (Sum of Pairs): cálculo detallado con ejemplo

El SP score es la métrica estándar para evaluar la calidad de un MSA. Para cada columna del MSA, suma los scores entre todos los pares de caracteres en esa columna. Luego suma todas las columnas.

```
SP = Σ_columnas Σ_(pares i<j) score(seq_i[col], seq_j[col])
```

**Reglas:**
- gap con residuo → penalidad normal (gap penalty)
- **gap con gap → contribución = 0** (no penaliza, por convención)
- dos residuos iguales → match score
- dos residuos distintos → mismatch score

**Ejemplo con scoring simple (match=+1, mismatch=0, gap−residuo=−1, gap−gap=0), 3 secuencias, 4 columnas:**

```
       col1  col2  col3  col4
Seq1:   A     C     -     T
Seq2:   A     C     G     T
Seq3:   A     T     G     T

Pares: (1,2), (1,3), (2,3) → 3 pares por columna

COLUMNA 1: A, A, A
  (1,2): score(A,A) = +1
  (1,3): score(A,A) = +1
  (2,3): score(A,A) = +1
  Subtotal = 3

COLUMNA 2: C, C, T
  (1,2): score(C,C) = +1
  (1,3): score(C,T) = mismatch = 0
  (2,3): score(C,T) = 0
  Subtotal = 1

COLUMNA 3: -, G, G
  (1,2): score(-,G) = gap−residuo = −1
  (1,3): score(-,G) = −1
  (2,3): score(G,G) = match = +1
  Subtotal = −1

COLUMNA 4: T, T, T
  (1,2): score(T,T) = +1
  (1,3): score(T,T) = +1
  (2,3): score(T,T) = +1
  Subtotal = 3

SP TOTAL = 3 + 1 + (−1) + 3 = 6
```

**Para 4 secuencias**: SP_total = suma de los 6 pares posibles = S(1,2) + S(1,3) + S(1,4) + S(2,3) + S(2,4) + S(3,4). Se calcula columna a columna sumando todos los pares.

---

### La limitación crítica del enfoque progresivo

**El problema:**

```
SITUACIÓN:
  Tenés 5 secuencias: A, B, C, D, E
  A y B son muy similares → se alinean primero → FIJO

  Supongamos que cuando se considera C, D, E, resulta que había que
  poner un gap en posición 7 del alineamiento (A,B).
  → Pero ese alineamiento ya no se puede modificar.
  → El error se propaga a todo el MSA final.
```

**Consecuencias:**
- Si las secuencias más similares tienen alguna peculiaridad no representativa de la familia, ese ruido se propaga.
- Para secuencias muy divergentes (identidad < 20%) CLUSTALW da resultados pobres.

**Soluciones parciales:**
- Refinamiento por SP score.
- Métodos más avanzados: MUSCLE, MAFFT, T-Coffee.
- Métodos probabilísticos con HMM: consideran toda la información global (Clase 4).

---

### El hilo conductor completo (como lo presentó el profe)

El profe enmarcó esto como una clase de **transición e integración**. El flujo conceptual es:

```
Bases de datos (registros, campos, índices)
    ↓ ¿cómo busco por secuencia en vez de por nombre?
BLAST (indexar el campo secuencia con words → búsqueda rápida)
    ↓ los hits de BLAST son secuencias homólogas
MSA (alinear múltiples secuencias con CLUSTALW)
    ↓ el MSA revela conservación posición por posición
PSSM (matriz sitio-específica = una matriz por posición del MSA)
    ↓ usar la PSSM para buscar en vez de BLOSUM62 genérica
PSI-BLAST (BLAST iterativo con PSSM → detecta homólogos remotos)
    ↓
Perfiles / HMM (Clase 4: modelo probabilístico del MSA)
```

Cada paso es la *entrada* del siguiente. Entender este flujo es esencial para no ver estos temas como capítulos aislados.

---

## Clase 4 — Bases de Datos Secundarias — HMM, PROSITE, InterPro

### El problema: una BD primaria no es un clasificador

Las BDs primarias (Clase 1) guardan registros individuales: una secuencia, un experimento, una estructura. Cuando yo tengo una proteína nueva y quiero saber a qué familia pertenece, puedo hacer BLAST contra Swiss-Prot y ver qué proteínas conocidas se parecen. Pero hay un límite: dos proteínas de la misma familia pueden tener menos del 30% de identidad de secuencia (zona twilight) y BLAST las pierde.

La diferencia conceptual fundamental es esta:

- **BD primaria**: el registro = una secuencia (resultado de un experimento). La pregunta que responde: "¿hay una secuencia parecida a la mía?"
- **BD secundaria**: el registro = un **modelo matemático** que describe una familia entera. La pregunta que responde: "¿pertenece mi secuencia a esta familia?"

En una BD secundaria el dato no viene de un experimento, viene del **análisis de datos de otras BDs**. El registro es una abstracción: captura la esencia matemática de qué tienen en común todos los miembros de la familia, incluyendo los que tienen solo 20% de identidad de secuencia entre sí.

### ¿Qué es un dominio? ¿Qué es una familia?

**Familia**: grupo de proteínas evolutivamente relacionadas (mismo ancestro común) que comparten función. El término viene del árbol filogenético.

**Dominio**: segmento de secuencia proteica con tres propiedades simultáneas:
1. **Función propia** — el dominio tiene actividad biológica característica por sí solo. Ej: dominio kinasa, dominio de unión ATP, dominio de unión al ADN.
2. **Estructura independiente** — el dominio se pliega correctamente incluso si se lo separa del resto de la proteína. Esto permite estudiarlo en forma aislada.
3. **Evolución independiente** — los dominios se "mezclan y combinan" evolutivamente (domain shuffling). Las proteínas de organismos superiores son más complejas precisamente porque tienen más dominios.

Ejemplo: un anticuerpo tiene dominios variables (reconocen el antígeno) y dominio Fc (interacciona con células inmunes). Cambiando el dominio Fc se obtienen las diferentes clases de inmunoglobulinas (IgG, IgM, etc.) que tienen la misma especificidad pero diferente función efectora. Los exones suelen corresponderse con dominios — el splicing alternativo es en parte un mecanismo de ensamblar combinaciones de dominios.

### El desafío: construir el modelo matemático

Para que una BD secundaria sea útil, el modelo matemático de cada registro tiene que ser capaz de **clasificar**: dada una secuencia nueva, decir si pertenece o no. Para hacer eso hay que resolver tres problemas:

| Problema | Pregunta | Descripción |
|----------|----------|-------------|
| **Entrenamiento** | ¿Cómo estimo los parámetros del modelo? | A partir de las secuencias conocidas de la familia, determinar los parámetros que definen el modelo de ese registro |
| **Puntaje (scoring)** | ¿Qué tan bien encaja una secuencia nueva en este modelo? | Dado el modelo ya entrenado, asignar un número a la secuencia nueva |
| **Alineamiento** | ¿Cómo se posiciona la secuencia nueva dentro del modelo? | Encontrar cuál es la correspondencia óptima entre posiciones de la secuencia y posiciones del modelo |

Estos tres problemas están **íntimamente ligados**: para calcular el puntaje hay que hacer el alineamiento; para hacer el alineamiento hay que tener el modelo entrenado.

**La lógica del umbral:**
Una vez construido el modelo, se evalúan dos grupos de secuencias:
- Secuencias que **sí** deben pertenecer a la familia → distribución de scores "positivos"
- Secuencias que **no** pertenecen → distribución de scores "negativos"

Lo ideal: dos histogramas bien separados, sin solapamiento. En la práctica hay solapamiento → se elige un umbral de corte. Las secuencias que caen del lado equivocado son los falsos positivos (FP) y falsos negativos (FN).

---

### Solución 1: PROSITE — expresiones regulares (alta especificidad, baja sensibilidad)

La forma más simple de modelar un motivo es con un lenguaje formal que describa exactamente qué aminoácidos se toleran en cada posición. PROSITE usa una sintaxis inspirada en expresiones regulares:

| Constructor | Significado | Ejemplo |
|-------------|-------------|---------|
| Letra | Aminoácido fijo obligatorio | `N` = solo Asn |
| `x` | Cualquier aminoácido | `N-x-T` |
| `[AB]` | A **o** B en esa posición | `[ST]` = Ser o Thr |
| `{AB}` | Cualquiera **excepto** A o B | `{P}` = cualquiera menos Pro |
| `A(n)` | Repetición exacta n veces | `x(3)` = tres cualquieras |
| `x(n,m)` | Entre n y m repeticiones | `x(2,4)` |
| `x(0,n)` | Elemento opcional (0 a n veces) | `x(0,1)` = puede o no estar |

> Los guiones `-` son separadores opcionales para legibilidad.

**Ejemplos:**

| Patrón | Cómo leerlo |
|--------|-------------|
| `[AC]-x-V-x(4)-{ED}` | (A o C) — cualquiera — V — 4×cualquiera — (no E ni D) |
| `A-x(0,1)-{V}` | A — (un aa opcional) — (no V) |
| `[ST]-C-x(1,2)-M` | (S o T) — C — (1 o 2 aa) — M |

**Motivo de N-glicosilación:** `N-{P}-[ST]`
- Asn − (cualquier aa **excepto** Pro) − (Ser **o** Thr)
- **¿Por qué no Prolina?** La Pro tiene su nitrógeno backbone atrapado en un anillo pirrolidínico → no puede rotar → fuerza un quiebre (kink) rígido en la cadena peptídica → el tripéptido N-P-S/T no puede adoptar la conformación de β-turn que necesita la enzima OST para realizar la glicosilación. Por eso Pro está **excluida** de esa posición.

**Limitación de PROSITE:** solo captura variabilidad posición a posición mediante un lenguaje discreto. No aprende qué tan frecuente es cada aminoácido en cada posición, ni modela inserciones/deleciones con probabilidades. Resultado: alta especificidad (pocas falsas alarmas) pero baja sensibilidad (se pierde muchos miembros divergentes de la familia).

---

### Solución 2: PSSM (Position-Specific Scoring Matrix)

Una PSSM es una tabla de 20 columnas (una por aminoácido) × L filas (una por posición del motivo). Cada celda contiene `score(posición i, aa a) = log[ frecuencia observada / frecuencia de fondo ]`.

Se construye a partir del MSA de las secuencias conocidas de la familia:
- Posición muy conservada (ej: siempre Asp) → score muy alto para Asp, muy negativo para el resto
- Posición variable → scores similares para muchos aminoácidos

**Ventaja sobre PROSITE:** cuantifica la variabilidad en lugar de usar "sí/no". Es lo que usa PSI-BLAST para refinar su búsqueda iterativamente.

**Limitación de la PSSM:** trata todas las posiciones independientemente. No modela explícitamente inserciones ni deleciones — solo funciona bien si las secuencias que comparo tienen exactamente el mismo largo que el motivo.

---

### Solución 3: HMM de perfil — el modelo más potente

En 1996 se demostró que un perfil de secuencias (como la PSSM) puede describirse como un **Modelo Oculto de Markov (Hidden Markov Model)**. Esto importa porque los HMMs tienen toda una teoría matemática detrás que permite resolver los 3 problemas de forma rigurosa y eficiente.

#### ¿Qué es un HMM?

Un HMM es una caja negra que **genera secuencias**. Por dentro tiene:
- **Estados** (ocultos — no los observamos directamente)
- **Probabilidades de transición** entre estados: P(ir del estado A al estado B)
- **Probabilidades de emisión**: para cada estado, la probabilidad de que "escupa" cada símbolo (aminoácido o nucleótido)

La clave: los estados están **ocultos**. Solo vemos la secuencia de símbolos emitidos, nunca la secuencia de estados que la generó. Eso explica el nombre "oculto".

**Ejemplo simplificado** (para intuir el concepto):
Supongamos un HMM con 2 estados: Estado 1 emite A y T con frecuencia (zona AT-rica), Estado 2 emite G y C con frecuencia (zona GC-rica). Si el modelo está el 98% del tiempo en Estado 1 y 2% en Estado 2, las secuencias generadas tendrán mayoría de A y T. Si veo una región de la secuencia con muchas G y C, puedo inferir que el modelo estaba probablemente en Estado 2 ahí — aunque no puedo saberlo con certeza. Los estados permanecen ocultos.

#### Los 3 tipos de estados en un HMM de perfil de proteínas

Para modelar una familia de proteínas, el HMM usa tres tipos de estados:

```
          ┌─────────────────── Match ──────────────────┐
          │  Estado M(i): corresponde a la columna i   │
          │  del MSA. Emite los 20 aa con probabilidad │
          │  = frecuencia observada en esa columna.    │
          └────────────────────────────────────────────┘

          ┌─────────────────── Insert ─────────────────┐
          │  Estado I(i): puede aparecer entre dos      │
          │  estados Match. Modela una inserción        │
          │  (aminoácido extra que no corresponde a     │
          │  ninguna columna del MSA). Emite con        │
          │  probabilidad de fondo (uniforme).          │
          └────────────────────────────────────────────┘

          ┌─────────────────── Delete ─────────────────┐
          │  Estado D(i): estado silente — no emite     │
          │  ningún símbolo. Modela una deleción        │
          │  (un gap en esa posición). Se "salta"       │
          │  el estado Match correspondiente.           │
          └────────────────────────────────────────────┘
```

La arquitectura del modelo para una familia con L columnas en el MSA es:
```
Inicio → M(1) → M(2) → M(3) → … → M(L) → Fin
           ↕      ↕      ↕
          I(1)  I(2)  I(3)      ← estados de inserción
           ↕      ↕      ↕
          D(1)  D(2)  D(3)      ← estados de deleción (silentes)
```

Las probabilidades de transición (P(M→M), P(M→I), P(M→D)) se estiman a partir de la frecuencia de gaps observados en el MSA de entrenamiento.

#### Problema 1: Entrenamiento — construir el modelo desde el MSA

```
Paso 1: Tomar el MSA de las secuencias conocidas de la familia
Paso 2: Para cada columna i del MSA:
        → Contar la frecuencia de cada aminoácido → probabilidades de emisión de M(i)
        → Contar la frecuencia de gaps → probabilidades de ir a D(i) o I(i)
Paso 3: Suavizar las probabilidades con frecuencias de fondo
        (porque nunca tengo infinitas secuencias para cada columna)
Resultado: el perfil HMM = archivo .hmm con todos los parámetros
Herramienta: hmmbuild
```

**¿Por qué suavizar con frecuencias de fondo?** Si en 10 secuencias de entrenamiento la posición 5 siempre tiene Phe, la probabilidad estimada de Ala en esa posición sería 0. Pero probablemente existan miembros de la familia donde esa posición tiene Tyr (similar a Phe). El suavizado mezcla la observación con la frecuencia general de cada aminoácido, evitando asignar probabilidad 0 a algo solo porque no apareció en el set de entrenamiento.

#### Problema 2: Scoring — ¿qué puntaje le doy a una secuencia nueva?

El puntaje de una secuencia S contra el HMM es la **probabilidad de que el HMM haya generado esa secuencia S**. Matemáticamente: P(S | modelo).

Intuitivamente: si corriera el HMM millones de veces y contara cuántas veces genera exactamente la secuencia S, esa fracción es el puntaje. En la práctica se calcula analíticamente usando programación dinámica.

Después del scoring se compara contra el umbral:
- Score > umbral → la secuencia probablemente pertenece a la familia
- Score < umbral → probablemente no pertenece

#### Problema 3: Alineamiento — ¿cómo se "encaja" la secuencia en el modelo?

Alinear la secuencia contra el HMM significa determinar **qué estado corresponde a cada posición de la secuencia**: ¿este aminoácido corresponde al estado M(3), o al estado I(2), o al D(4)? Esa asignación es el alineamiento.

La solución óptima es la asignación de estados que maximiza la probabilidad P(S, secuencia_de_estados | modelo). Se calcula con programación dinámica igual que NW/SW, pero en lugar de llenar una matriz de scores se llena una matriz de probabilidades.

**Conexión con inserciones y deleciones:**
- Si la secuencia tiene un aminoácido "de más" en la posición 5 → ese aminoácido se alinea contra el estado I(4) o I(5) → es una inserción
- Si la secuencia tiene un gap en la posición 3 → el estado M(3) se mapea a un D(3) → es una deleción

Esta es la gran ventaja del HMM sobre la PSSM: **modela inserciones y deleciones de forma explícita y probabilística**, por eso es mucho más potente para detectar homólogos divergentes.

---

### Pipeline HMMer completo

```
1. Secuencias conocidas de la familia (experimentalmente validadas)
        ↓
2. MSA de esas secuencias (con CLUSTALW, MAFFT, etc.)
        ↓
3. hmmbuild   →   construye el perfil HMM (.hmm)
   = PROBLEMA DE ENTRENAMIENTO
        ↓
4. Calibración del umbral:
   - Evaluar las secuencias de entrenamiento con el modelo → distribución de scores "positivos"
   - Evaluar secuencias al azar → distribución de scores "negativos"
   - Elegir el umbral que separa las dos distribuciones (minimiza FP + FN)
        ↓
5. hmmsearch  →   buscar en el proteoma o en una BD con ese umbral
   = PROBLEMA DE PUNTAJE (+ ALINEAMIENTO para cada hit)
        ↓
6. Candidatos con score > umbral → validación experimental
        ↓
7. Evaluar calidad: sensibilidad / especificidad / curva ROC
```

---

### Sequence Logo — leer el modelo visualmente

El logo es la visualización del HMM. Cada columna del logo corresponde a un estado Match del HMM.

- **Altura total de la columna** = información (bits) = qué tan conservada está esa posición
  - Máximo teórico ≈ 4.32 bits (= log₂20 — posición completamente conservada en 1 aminoácido)
  - Mínimo = 0 bits — posición completamente variable (todos los aa igualmente probables)
- **Tamaño relativo de cada letra** = probabilidad de ese aminoácido en esa posición
- **Letra más grande** = aminoácido más frecuente (más probable que el HMM emita en ese estado Match)

**Qué se puede leer de un logo:**
1. **Posiciones críticas** — columnas altas: ese residuo es casi imposible de cambiar sin perder función
2. **Identidad fisicoquímica del residuo dominante** — ¿es cargado, hidrófobo, polar? Dice qué tipo de interacción necesita la proteína ahí
3. **Tolerancia a sustituciones conservativas** — si la columna tiene altura media y 2-3 letras similares (ej: L, I, V), acepta variación dentro del grupo hidrófobo
4. **Posiciones variables** — columnas bajas: no involucradas en reconocimiento específico
5. **Largo del motivo** — número total de columnas

> **¡OJO en el parcial!** Verificar si el logo es de proteínas (letras de aminoácidos: A, R, N, D...) o de ADN (A, T, G, C). Si tiene T y U con frecuencias similares es ARN. Si tiene principalmente H, I, L, V, F... es proteína.

---

### Métricas de evaluación del modelo

Una vez que se busca en el proteoma, los resultados se clasifican en:

```
                      Predicción del modelo
                    ┌──────────┬──────────┐
                    │ POSITIVO │ NEGATIVO │
       ┌────────────┼──────────┼──────────┤
Verdad │ POSITIVO   │    VP    │    FN    │
real   │ NEGATIVO   │    FP    │    VN    │
       └────────────┴──────────┴──────────┘

Sensibilidad  = VP / (VP + FN)   ← fracción de miembros reales detectados
Especificidad = VN / (VN + FP)   ← fracción de no-miembros correctamente rechazados
```

**Falso positivo:** secuencia que supera el umbral pero no es miembro real de la familia. Causas posibles:
- El motivo existe en la secuencia pero está enterrado/inaccesible estructuralmente
- El umbral está demasiado bajo (muy permisivo)
- Set de entrenamiento pequeño → el modelo aprendió ruido estadístico

**Falso negativo:** miembro real de la familia que no supera el umbral. Causas posibles:
- Secuencia muy divergente (zona twilight)
- El umbral está demasiado alto (muy estricto)
- La región de la familia está poco representada en el set de entrenamiento

**Curva ROC:** se barre el umbral de score y en cada valor se calcula sensibilidad y especificidad → eje Y = sensibilidad, eje X = 1 − especificidad. Un modelo perfecto pasa por la esquina superior izquierda (sensibilidad=1, especificidad=1). Diagonal = clasificación aleatoria.

**Tensión entre sensibilidad y especificidad:** bajar el umbral aumenta sensibilidad (se detectan más miembros) pero baja especificidad (entran más FP). Subir el umbral hace lo contrario. El umbral óptimo depende del objetivo biológico.

---

### ¿Por qué HMM > BLAST para familias?

| Criterio | BLAST | HMM de perfil |
|----------|-------|---------------|
| Modelo | Secuencia individual | Toda la familia a la vez |
| Inserciones/deleciones | Solo penaliza con gap penalty uniforme | Modela con probabilidades posición-específicas |
| Posiciones conservadas | No las distingue de las variables | Les asigna mayor peso automáticamente |
| Detección de homólogos remotos | Pierde en zona twilight (<30% id) | Detecta hasta 10-20% de identidad |
| Construcción del modelo | No necesita entrenamiento | Requiere MSA de entrenamiento curado |

---

### Bases de datos secundarias principales

| BD | Método | Contenido | Notas |
|----|--------|-----------|-------|
| **Pfam** | HMM (hmmbuild/hmmsearch) | Familias y dominios proteicos | Seeds curados manualmente; la más usada |
| **PROSITE** | Regex + PSSM | Patrones cortos y perfiles | Alta especificidad; motivos funcionales |
| **SMART** | HMM | Dominios de señalización | Incluye contexto arquitectónico (qué dominios aparecen juntos) |
| **TIGRFAMs** | HMM | Familias de función específica | Útil en genómica microbiana |
| **CDD** (NCBI) | RPS-BLAST | Dominios conservados | Integrado en BLAST de NCBI |

### InterPro — la BD de BDs secundarias

InterPro (EBI) consolida todas las BDs de dominios/familias en una sola entrada unificada. Un registro InterPro agrupa entradas equivalentes de Pfam, PROSITE, SMART, etc. que describen el mismo dominio o familia.

- **Herramienta:** **InterProScan** — busca una secuencia contra todas las bases miembro en una sola corrida
- **Output:** lista de dominios encontrados con posición, E-value, BD de origen y términos GO (Gene Ontology)

**Tipos de entradas InterPro:**

| Tipo | Descripción |
|------|-------------|
| **Family** | Proteínas evolutivamente relacionadas que comparten función |
| **Domain** | Unidad estructural/funcional que aparece en distintos contextos proteicos |
| **Homologous superfamily** | Comparten estructura 3D sin similitud de secuencia detectable |
| **Repeat** | Secuencia corta repetida en tandem dentro de la proteína |
| **Site** | Sitio activo, de unión a ligando, o de modificación post-traduccional |

---

## Clase 5 — NGS y Ensamblado de Genomas

### Introducción: el problema de leer el ADN

Desde que Watson y Crick describieron la doble hélice, los biólogos entendieron que la información biológica está en la secuencia de bases del ADN. El desafío fue desarrollar métodos para leer esa secuencia. La historia de la secuenciación es la historia de dos problemas que se fueron resolviendo progresivamente:

1. **¿Cómo identificar cada base?** → necesito un mecanismo que me diga si es A, T, C o G
2. **¿Cómo escalar?** → un genoma humano tiene 3.000 millones de bases; no puedo leerlas de a una

Las soluciones a estos dos problemas evolucionaron en generaciones tecnológicas. Cada generación resolvió las limitaciones de la anterior y abrió nuevas preguntas.

---

### Solución 1: Sanger (1ª generación) — leer un fragmento a la vez

Sanger (Premio Nobel ×2) inventó el primer método práctico basado en **didesoxinucleótidos terminadores (ddNTPs)**. La idea: si la polimerasa copia una hebra molde pero de vez en cuando incorpora un nucleótido "roto" que le impide seguir, va a generar fragmentos de todos los largos posibles (1, 2, 3… N bases). Si además marco **en qué base termina cada fragmento**, puedo ordenarlos por tamaño y leer la secuencia.

**¿Qué es un ddNTP?** Un nucleótido normal tiene un -OH en el carbono 3' que permite que la polimerasa enganche el siguiente nucleótido. Un ddNTP (didesoxinucleótido) no tiene ese -OH → la polimerasa lo incorpora, pero ahí se detiene: no puede agregar nada más. Es un terminador.

**Paso a paso:**

```
Hebra molde (lo que quiero secuenciar): 5'- G A T C A G -3'

Hago 4 tubos en paralelo. En cada uno pongo:
  • polimerasa + todos los dNTPs normales (para copiar)
  • 10% de un ddNTP distinto (para terminar al azar en esa base)

  Tubo G: 10% ddGTP  → la polimerasa copia hasta que incorpora una ddG → para
  Tubo A: 10% ddATP  → para cuando incorpora una ddA
  Tubo T: 10% ddTTP  → para cuando incorpora una ddT
  Tubo C: 10% ddCTP  → para cuando incorpora una ddC

Al final de cada tubo tengo una mezcla de fragmentos de todos los largos
posibles que terminan en esa base. Por ejemplo, el Tubo G contiene:
  • el fragmento de largo 1 (termina en la primera G de la copia)
  • el fragmento de largo 4 (termina en la segunda G)
  • etc.
```

**El gel — leer la secuencia:**

Un gel de electroforesis es como una carrera de fragmentos: los más cortos migran más rápido y llegan más abajo. El gel tiene **4 calles** (carriles verticales), una por tubo:

```
               Calle G   Calle A   Calle T   Calle C
  ARRIBA        |         |         |         |
  (fragmentos   |         |         |    ●    |   ← largo 5 → termina en T → 5ª letra = T
  largos)       |         |    ●    |         |   ← largo 4 → termina en A → 4ª letra = A
                |    ●    |         |         |   ← largo 3 → termina en G → 3ª letra = G
                |         |    ●    |         |   ← largo 2 → termina en A → 2ª letra = A
  ABAJO         |    ●    |         |         |   ← largo 1 → termina en G → 1ª letra = G
  (fragmentos
  cortos)

  Leyendo de abajo hacia arriba la calle donde aparece cada banda:
  G → A → G → A → T → ...   = la secuencia de la copia (complemento de la molde)
```

**Versión moderna — Sanger fluorescente + capilar:**

El gel con 4 calles es incómodo. La solución fue darle a cada ddNTP un **color fluorescente distinto** → ahora todo cabe en 1 sola reacción (no 4 tubos). Los fragmentos se separan por tamaño en un capilar (un tubo finísimo) y un detector lee el color de cada fragmento al pasar:

```
Fragmentos pasando por el capilar (del más corto al más largo):

  tiempo →  [●rojo] [●verde] [●rojo] [●azul] [●amarillo] ...
               G        A       G       A        T        ...
              
  El detector registra: rojo=G, verde=A, rojo=G, azul=A, amarillo=T → GAGAT...
  Eso es el electroferograma: un gráfico de picos de color = la secuencia
```

**Limitación que motiva la siguiente generación:** todo esto secuencia **un solo fragmento** por reacción. Para un genoma humano (3.000 millones de bases), necesitaría millones de reacciones → décadas de trabajo. Solución: paralelizar masivamente.

---

### Herramienta clave: PCR (Mullis, Premio Nobel)

La PCR permite:
1. **Amplificación exponencial** de un fragmento específico a partir de una sola molécula
2. **Selectividad**: en una mezcla compleja (ej: genoma entero), amplifica solo el fragmento flanqueado por los primers
3. Amplificación **clonal**: todas las copias resultantes son idénticas

La PCR es la base del paso de amplificación clonal en las tecnologías de 2ª generación. Sin PCR (o algún equivalente), no hay señal detectable.

---

### Preparación de librería — paso previo a todos los NGS

El ADN genómico crudo no puede entrar directamente al secuenciador: es demasiado largo, no tiene puntos de anclaje y no se puede distinguir un fragmento de otro. La preparación de librería lo transforma en una mezcla de fragmentos cortos, etiquetados y listos para ser leídos.

---

#### Paso 1 — Fragmentación

El genoma humano tiene ~3 000 millones de pares de bases. Las tecnologías de 2ª generación leen ~100–600 pb por lectura, así que primero hay que romper el ADN en pedacitos de ese tamaño.

**Fragmentación mecánica (sonicación):**
```
Genoma completo:
═══════════════════════════════════════════════════════
           ↓ ultrasonido (ondas de presión)
    ════   ════════   ═══   ════════   ═══   ════
    ~200pb   ~350pb   ~150pb  ~400pb   ~180pb  ~300pb
```
Las ondas de sonido rompen el ADN al azar → fragmentos de tamaño variable. Luego se selecciona el rango de tamaño deseado con geles o beads (se descartan los demasiado cortos o largos).

**Fragmentación enzimática:** enzimas de restricción o transposasas cortan en sitios específicos o semi-aleatorios. Más reproducible que la sonicación pero introduce un sesgo hacia los sitios de corte de la enzima.

---

#### Paso 2 — Agregar adaptadores

Los fragmentos recién rotos tienen extremos irregulares y no tienen nada que los identifique. Se les ligan adaptadores: **secuencias cortas de ADN conocido** que se pegan a ambos extremos de cada fragmento.

```
Fragmento crudo:
        ATCGGCTATGCATTAGC...GCTATGCA
                  ↓ ligación de adaptadores
   [ADAPT]ATCGGCTATGCATTAGC...GCTATGCA[ADAPT]
    ~20pb                               ~20pb
```

¿Para qué sirven los adaptadores?
- **Anclaje:** en Illumina, los adaptadores se hibridan a los primers del flow cell para que el fragmento quede pegado a la superficie.
- **Primer de inicio:** la polimerasa arranca desde el adaptador para leer el fragmento.
- **Identificador (índice/barcode):** si se mezclan muestras de varios pacientes en un mismo flow cell, se usa un adaptador con un código único por muestra → al final se puede separar qué lectura viene de quién.

```
Librería completa (todos los fragmentos tienen los mismos adaptadores en los extremos):

   [ADAPT]─── fragmento 1 ───[ADAPT]
   [ADAPT]─── fragmento 2 ───[ADAPT]
   [ADAPT]─── fragmento 3 ───[ADAPT]
   [ADAPT]─── fragmento 4 ───[ADAPT]
      ↑
   (mismo adapt. en todos → el secuenciador puede anclarlos a todos)
```

---

#### Paso 3 — [Opcional] Enriquecimiento en regiones de interés

Secuenciar todo el genoma humano (~3Gb) es caro. Si solo me interesa el **exoma** (las regiones codificantes, ~1% del genoma) o un **panel de genes de enfermedad**, puedo enriquecer esas regiones antes de secuenciar. Hay dos estrategias:

---

**Estrategia A — Captura por hibridación (exoma o panel)**

La idea: preparar sondas complementarias a las regiones que me interesan y usarlas como imanes moleculares.

```
Paso A1 — librería completa en solución:
   [AD]─exón 1─[AD]      [AD]─intrón─[AD]     [AD]─exón 2─[AD]
   [AD]─exón 3─[AD]      [AD]─basura─[AD]      [AD]─exón 4─[AD]
   (todos mezclados, la mayoría son intrones o regiones no deseadas)

Paso A2 — agregar sondas biotiniladas complementarias a los exones:
   Sonda: ████████████ (complementaria al exón, con biotina en el extremo)
              ||||||||||||
   [AD]─exón 1─[AD]   ← se hibrida con la sonda
   
   Los intrones y regiones no deseadas no hibridan → quedan sueltos.

Paso A3 — agregar bolitas con estreptavidina y un imán:
   Estreptavidina ─── bolita magnética
       │
   biotina ─── sonda ─── [AD]─exón─[AD]   ← pegado

   Imán: ▓▓▓▓▓▓▓▓▓▓▓▓
                         ← bolitas (con exones pegados) van al imán
                         ← intrones y basura quedan en solución → se descartan

Paso A4 — lavar y eluir:
   Resultado: solo los fragmentos de exones, listos para secuenciar
```

Ejemplo real: secuenciación de exoma completo (WES). Hay ~20 000 genes codificantes → el panel de sondas cubre todas las regiones exónicas del genoma. El resultado es ~50× menos datos que el genoma completo, a la misma profundidad de lectura.

---

**Estrategia B — Amplificación selectiva por PCR**

En lugar de capturar con sondas, simplemente se hacen primers que flanquean solo la región de interés y se amplifica solo eso.

```
Genoma/mezcla de fragmentos:
   ───────────────[REGIÓN DE INTERÉS]───────────────
                  ↑                 ↑
               Primer F          Primer R
               (forward)         (reverse)

PCR: solo amplifica lo que está entre los dos primers.
Todo el resto del genoma → queda a concentración bajísima.

Resultado:
   [REGIÓN DE INTERÉS]  × 10⁶ copias  ← lista para secuenciar
   resto del genoma     × 1 copia     ← ignorado
```

Ejemplo real del profesor: diagnóstico de COVID-19. La muestra nasal tiene ADN humano + posiblemente ARN viral. Con RT-PCR se hace:
1. Retrotranscripción: ARN viral → ADN copia
2. PCR con primers específicos para genes del virus (ej. gen N, gen E)
3. Solo el ADN viral se amplifica → detectable; el genoma humano no interfiere

---

**Resumen del pipeline completo:**

```
ADN genómico
    ↓ sonicación → fragmentos ~150–1000 pb con extremos irregulares
    ↓ ligation de adaptadores → [ADAPT]─fragmento─[ADAPT]
    ↓ [opcional] enriquecimiento:
        ┌─ Captura: sondas + imán → quedan solo los fragmentos de interés
        └─ PCR: primers flanqueantes → amplifican solo la región deseada
    ↓
  LIBRERÍA LISTA → entra al secuenciador (454, Illumina, PacBio, etc.)
```

Los adaptadores son esenciales: el secuenciador los usa como sitio de anclaje y como primer para empezar a leer.

---

### Pair-End Reads — concepto fundamental

Un fragmento de ADN doble cadena tiene **dos extremos** (ends). Los secuenciadores de 2ª generación leen ambos → generan **dos lecturas apareadas**:

```
5'───[Lectura 1]──────────────────[Lectura 2]───3'
      ←~100 pb→   (interior no leído)  ←~100 pb→
      │←──────────── inserto ~1–5 kb ───────────→│
```

- Las dos lecturas vienen del **mismo fragmento** → tamaño del inserto conocido de antemano
- Si sé dónde mapea una lectura, sé a qué distancia debe estar la otra
- **Fundamental para el scaffolding** (Paso 2 del ensamblado): permiten ordenar contigs y estimar el tamaño del gap entre ellos

---

### Profundidad (Depth) vs Cobertura (Coverage)

Son dos métricas distintas que se confunden fácilmente:

**Profundidad (depth):** número de veces que una base concreta aparece en las lecturas.
- Es un valor por posición, varía a lo largo del genoma
- Se resume como **profundidad media** = promedio sobre todas las posiciones

**Cobertura (coverage):** porcentaje del genoma cubierto por al menos X lecturas (X = umbral elegido antes de calcularlo).
- Ejemplo: "cobertura del 90% a 10X" → el 90% de las bases se leyó al menos 10 veces

```
Profundidad target (teórica) = (N_reads × L_read) / L_genoma  ← se fija ANTES del experimento

Profundidad media real  ← se calcula DESPUÉS del ensamblado
Cobertura real          ← se calcula DESPUÉS (necesita definir el umbral primero)
```

Calidad del experimento: si profundidad media real ≈ profundidad target → pocas lecturas se perdieron en ruido.

**Estándar clínico:** 30X para variant calling confiable. "Billetera mata galán": más cobertura → mejor resultado.

### Phred Quality Score

```
Q = −10 · log₁₀(P_error)

Q20 → P = 0.01  (1 error cada 100 bases)   → 99%
Q30 → P = 0.001 (1 error cada 1000 bases)  → 99.9%  ← estándar mínimo aceptable
```

### Formato FASTQ

4 líneas por read:
```
@ID            ← identificador
ACGTACGT...    ← secuencia
+              ← separador
IIHHGF...      ← calidades Phred en ASCII (Phred+33)
```

---

### Solución 2: 2ª Generación — paralelizar millones de reacciones

**Idea central:** en lugar de hacer una reacción de Sanger por fragmento, hacer millones de reacciones en paralelo en el mismo chip. El costo por letra baja porque los reactivos se reparten entre millones de lecturas simultáneas.

**Restricción compartida de todas las plataformas 2ª gen:** la señal de una sola molécula de ADN es indetectable. Para amplificar la señal hay que amplificar la cantidad de moléculas → **amplificación clonal** previa a la secuenciación. Cada fragmento se transforma en un cluster/colonia de miles de copias idénticas. Esto impone lecturas cortas (~100–600 pb) porque la PCR de amplificación clonal no funciona bien con fragmentos muy largos.

#### 2ª Gen — 454 Pyrosequencing (Roche, históricamente primera)

**Analogía general:** imaginate miles de tubitos de ensayo en miniatura (pocillos), cada uno haciendo su propia reacción de Sanger en paralelo. En lugar de leer el color de los fragmentos en un gel, detectás luz cuando se incorpora un nucleótido.

**Paso 1 — Amplificación clonal en emulsión (PCR en gotitas):**

```
Tengo mis fragmentos de ADN con adaptadores en los extremos.
Los mezclo con bolitas (beads) cubiertas de primers que se pegan a esos adaptadores.
Proporciones calibradas → en promedio 1 fragmento se pega a 1 bolita.

Ahora agito esa mezcla con aceite → se forma una EMULSIÓN:
el agua queda en gotitas microscópicas suspendidas en aceite.

Aceite  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        (  gotita 1  )   (  gotita 2  )   (  gotita 3  )
        [ bolita+frag1]  [ bolita+frag2]  [ bolita+frag3]
        (  PCR aquí  )   (  PCR aquí  )   (  PCR aquí  )
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Como cada gotita es un reactor aislado (el aceite la separa de las demás),
la PCR dentro de cada gotita amplifica SOLO el fragmento que tiene esa bolita.
Al terminar: cada bolita está rodeada de miles de copias del MISMO fragmento.
```

**Paso 2 — Cargar en la PicoTiterPlate:**

```
La PicoTiterPlate es una placa con miles de pocillos, cada uno del tamaño
exacto de 1 bolita. Como 1 bolita = 1 secuencia, cada pocillo tiene 1 secuencia distinta.

 ┌──┬──┬──┬──┬──┐
 │● │● │● │● │● │  ← bolitas con fragmentos amplificados
 ├──┼──┼──┼──┼──┤
 │● │● │● │● │● │  (en realidad hay cientos de miles de pocillos)
 └──┴──┴──┴──┴──┘
```

**Paso 3 — Secuenciación por pirofosfato (detección de luz):**

```
Agrego a todos los pocillos: polimerasa + sulfurilasa + luciferasa + luciferina.

Luego ciclo nucleótidos de a UNO:

Ciclo 1: agrego T a toda la placa
  → en los pocillos cuya siguiente base a copiar es A (complemento de T),
    la polimerasa incorpora la T y libera PIROFOSFATO (PPi)
  PPi → sulfurilasa → ATP → luciferasa → LUZ  💡
  
  Foto de la placa:  pocillo 1=luz💡  pocillo 2=oscuro  pocillo 3=luz💡 ...
  → registro: pocillo 1 incorporó T, pocillo 3 incorporó T

Lavo todo el T no incorporado.

Ciclo 2: agrego G
  → nueva foto → sé cuáles incorporaron G

Ciclo 3: agrego A ... y así sucesivamente.

Leyendo la secuencia de ciclos en los que cada pocillo emitió luz → secuencia de ese fragmento.
```

**⚠ Problema de los homopolímeros:**

```
¿Qué pasa si la secuencia tiene TTTT (4 T's seguidas)?

Cuando agrego T en el ciclo correspondiente, la polimerasa pone las 4T de GOLPE
(no para entre medio) → se liberan 4× más pirofosfato → 4× más luz.

El equipo calibra: 1 luz = 1T, doble luz = 2T, triple = 3T...
Pero con 5 o más T iguales seguidas, la escala se vuelve imprecisa → errores de lectura.

AAAAA → ¿5 A's? ¿6 A's? El equipo no lo sabe bien. Eso es el problema del homopolímero.
```

#### 2ª Gen — Ion Torrent (misma arquitectura, distinta detección)

La reacción de la luciferasa (454) requiere agregar enzimas caras que se degradan. Ion Torrent dice: "¿para qué luz si puedo medir directamente la química?" Cuando la polimerasa incorpora un nucleótido, libera un H⁺. Ion Torrent pone un **sensor de pH debajo de cada pocillo** — básicamente un pHímetro microscópico.

```
Arquitectura: IDÉNTICA a 454 hasta el paso 2
  (emulsión PCR → bolitas → pocillos de la PicoTiterPlate)

Paso 3: DETECCIÓN por cambio de pH

 ┌─────────────────┐
 │  pocillo con    │
 │  bolita+poli-   │  ← agrego nucleótido T
 │  merasa         │
 ├─────────────────┤
 │  sensor de pH   │  ← semiconductor de silicio
 └────────┬────────┘
          │
     señal eléctrica

Si la polimerasa incorpora la T → libera H⁺ → el pH cae → el sensor detecta la caída
Si no había que incorporar T (no es complementario) → no hay H⁺ → pH no cambia → silencio

Ciclo de nucleótidos: igual que 454 (T → lavar → G → lavar → A → lavar → C → lavar → ...)
Mido el cambio de pH en cada ciclo → sé qué base se incorporó en cada pocillo.
```

Ventaja: más barato y más rápido (no necesita luciferasa ni luciferina).
⚠ **Mismo problema de homopolímeros**: TTTT libera 4 H⁺ → mayor caída de pH, pero difícil distinguir 4T de 5T.

#### 2ª Gen — Illumina (dominante, ~70% del mercado)

Illumina resolvió la amplificación clonal de forma más elegante: en lugar de emulsionar en gotitas de aceite, **ancla los fragmentos directamente a una superficie sólida** (la flow cell) y los amplifica ahí mismo formando colonias locales llamadas clusters.

**Paso 1 — Bridge Amplification (amplificación en puente):**

```
La flow cell es una placa de vidrio con millones de primers cortos pegados a su superficie.

1. El fragmento (con adaptadores) se pega a uno de esos primers y se copia:

   superficie  ──────────────────────────────────────
                ↑primer  ↑primer  ↑primer  ↑primer
               [frag1]

2. El fragmento copiado se dobla y pesca el primer de al lado → forma un "puente":

   superficie  ──────────────────────────────────────
               (──────)   ← puente entre los dos primers
              original    copia

3. La polimerasa copia el puente → 2 copias. Se repite → 4 → 8 → 16 → ...
   Todas las copias quedan en el MISMO punto de la superficie = CLUSTER.

   superficie  ──────────────────────────────────────
               [cluster 1]   [cluster 2]   [cluster 3]
               (frag A ×1000) (frag B ×1000) (frag C ×1000)

1 flow cell tiene ~1.600 millones de clusters = 1.600 millones de reacciones en paralelo.
```

**Paso 2 — Secuenciación por síntesis con terminadores reversibles:**

Los 4 nucleótidos tienen DOS modificaciones: un fluoroforo de color distinto por base, y un bloqueo en el -OH 3' que impide incorporar más de uno por ciclo.

```
Ciclo 1: agrego los 4 dNTPs modificados a la flow cell.
  Cada cluster incorpora exactamente 1 nucleótido (el bloqueo frena ahí).
  Saco UNA FOTO de toda la flow cell:

  [cluster1: 🔵]  [cluster2: 🔴]  [cluster3: 🟢]  [cluster4: 🟡] ...
      A                T                C                G

  Paso de química: corto el bloqueo y el fluoroforo → el nucleótido queda libre.

Ciclo 2: agrego de nuevo los 4 dNTPs → cada cluster incorpora su siguiente base → nueva foto.
  [cluster1: 🔴]  [cluster2: 🔵]  [cluster3: 🔴]  [cluster4: 🟢] ...
      T                A                T                C

Ciclo 3: foto → ciclo 4: foto → ... → N ciclos = N bases leídas por cluster.

Leyendo el color del cluster 1 en cada foto: 🔵🔴🟢... → A, T, C... = su secuencia.
```

Ventaja: altísima precisión (~0.1% de error), sin problema de homopolímeros (se incorpora de a 1 aunque haya TTTT).
Desventaja: lecturas cortas (75–300 pb) — la señal se degrada después de muchos ciclos.

---

### Solución 3: 3ª Generación — eliminar la amplificación clonal

Las plataformas de 2ª generación tienen dos costos inherentes: (1) reactivos para la amplificación clonal y (2) lecturas cortas forzadas por ella. La 3ª generación elimina ambos detectando **molécula única** — para lo cual se necesita amplificar la señal por medios físicos en lugar de multiplicar las moléculas.

#### 3ª Gen — PacBio SMRT (secuenciación por síntesis, molécula única)

**El problema a resolver:** la señal de fluorescencia de UN SOLO fluoroforo es demasiado débil para detectarla — por eso la 2ª gen amplificaba la cantidad de moléculas. PacBio resuelve esto de otra manera: en lugar de multiplicar las moléculas, **concentra la luz** en un volumen nanoscópico usando los Zero-Mode Waveguides (ZMW).

```
¿Qué es un ZMW? Un pozo extremadamente pequeño (diámetro de ~70 nm,
mucho menor que la longitud de onda de la luz visible ~500 nm).
Por las leyes de la óptica, la luz no puede propagarse dentro de algo más
pequeño que su longitud de onda → queda atrapada y concentrada en el fondo.

Dentro de cada ZMW:
  • 1 polimerasa fija en el fondo (no hay amplificación clonal)
  • 1 molécula de ADN que se introduce → la polimerasa empieza a copiarla
  • Los 4 dNTPs con fluoroforos distintos circulan en la solución encima

Cuando la polimerasa incorpora, por ejemplo, una T:
  • La T con su fluoroforo queda retenida en el fondo del ZMW durante ~ms
  • La luz concentrada en ese volumen hace que el fluoroforo emita con señal
    ~1000× más intensa que en solución → DETECTABLE a molécula única
  • Después del ms, el fluoroforo se libera y la polimerasa sigue con la siguiente base

  tiempo →  T💡  A💡  G💡  T💡  C💡 ...
              ↑    ↑    ↑    ↑    ↑
            color distinto por base → leés la secuencia en tiempo real
```

Sin necesidad de amplificar → el fragmento puede ser TAN LARGO COMO QUIERA la polimerasa copiarlo → lecturas de 1–50 kb o más.

Desventaja: la polimerasa comete más errores (~1% vs 0.1% de Illumina). Para compensar se secuencia con mayor profundidad o se hace circular el mismo fragmento varias veces.

---

#### 3ª Gen — Oxford Nanopore (sin síntesis, lectura directa)

PacBio todavía necesita polimerasa y nucleótidos. Nanopore hace algo conceptualmente más radical: **leer el ADN directamente**, sin sintetizar nada, igual que leer un libro sin tener que copiarlo.

```
Analogía: imaginá una manguera de agua. Si la tapás con el dedo, el flujo cae.
Si la tapás con algo más grande, cae más. Nanopore hace lo mismo pero con iones
y moléculas de ADN.

Configuración:
  • Una membrana con un PORO PROTEICO (derivado de una proteína de bacteria,
    modificado por ingeniería genética para que le encajen bien las 4 bases)
  • A un lado: voltaje positivo → los iones fluyen a través del poro → corriente medible
  • El ADN simple cadena, cargado negativamente, es "jalado" hacia el voltaje positivo
    → pasa a través del poro de a una base por vez

Lo que pasa mientras el ADN pasa:

  corriente
  sin ADN:  ──────────────────────────────────  (corriente constante = poro abierto)

  con ADN:  ───────┐    ┌──┐  ┌──────┐  ┌───   (cada base tapa el poro diferente)
            base A ↓T   ↓G  ↓G       ↓C  ↓A
                   (cada caída tiene amplitud característica de esa base)

  El cuello del poro tiene aminoácidos que forman puentes de hidrógeno con cada base
  de forma diferencial → A tapa diferente que T, que G, que C → caída distinta.
  El detector registra la corriente ms a ms → identifica cada base directamente.
```

Sin síntesis, sin amplificación → sin límite de largo de lectura. El mismo poro puede leer fragmentos de megabases.

Ventajas: lecturas larguísimas, portátil (tamaño USB), no gasta reactivos de síntesis.
Desventaja: los poros proteicos son descartables (el cartucho con ~5000 poros se usa una vez). Tasa de error más alta (~5–10%). Futuro: reemplazar los poros proteicos por poros sintéticos de grafeno o semiconductores → más estables, más baratos, reutilizables.

---

### Resumen comparativo de plataformas

| Plataforma | Gen | Ampl. clonal | Long. lectura | Error | Detección |
|---|---|---|---|---|---|
| 454 (Roche) | 2ª (discontinuado) | Emulsión PCR | ~400 pb | ~1% | Bioluminiscencia (PPi → luz) |
| Ion Torrent | 2ª | Emulsión PCR | 200–600 pb | ~1% | pH (H⁺) |
| Illumina | 2ª | Bridge PCR | 75–300 pb | ~0.1% | Fluorescencia + terminador reversible |
| PacBio SMRT | 3ª | No (mol. única) | 1–50 kb | ~1% | Fluorescencia en ZMW |
| Oxford Nanopore | 3ª | No (mol. única) | >100 kb | ~5–10% | Corriente iónica |

---

### Ensamblado de genomas

```
LECTURAS (millones de fragmentos cortos de 100–300 pb)
    ↓ Paso 1: solapamiento → Contigs
CONTIGS (secuencias continuas, secuencia base a base conocida)
    ↓ Paso 2: pair-end reads → Scaffolds
SCAFFOLDS (contigs ordenados, huecos de tamaño conocido)
    ↓ Paso 3: mapa físico (STS) → Cromosomas
GENOMA COMPLETO
```

Las bases de datos almacenan los contigs. En RefSeq tienen prefijo `NT_`. Al buscar un gen, siempre se obtiene la secuencia del contig correspondiente.

#### Paso 1 — Generar contigs (OLC vs De Bruijn)

**Razonamiento biológico:** el genoma que quiero ensamblar es desconocido. Todo lo que tengo son millones de lecturas. Si varias lecturas se solapan con alta identidad en sus extremos, deben venir de la misma región del genoma y puedo unirlas.

**Método OLC (Overlap-Layout-Consensus) — el clásico:**
1. **Overlap:** comparar todas las lecturas contra todas → para N lecturas, N² comparaciones → con millones de reads, extremadamente costoso.
2. **Layout:** construir grafo donde nodos = lecturas y aristas = solapamientos. Buscar camino **hamiltoniano** (pasa por cada nodo una vez). El hamiltoniano es NP-completo → intratable con millones de reads.
3. **Consensus:** con el orden definido, calcular la secuencia consenso.

**Método De Bruijn — el estándar actual para NGS:**

El grafo de De Bruijn resuelve el problema de OLC: la complejidad del grafo **no aumenta con el número de lecturas**.

- Los **nodos** son los (k-1)-mers. El número máximo es 4^(k-1) — fijo para un k dado, sin importar cuántas lecturas haya.
- Cada k-mer de una lectura genera **una arista** entre su prefijo (primeras k-1 bases) y su sufijo (últimas k-1 bases).
- Agregar más lecturas solo añade peso a aristas existentes, no nuevos nodos.

```
Lectura: AACGG (k=4)
k-mers: AACG, ACGG
Aristas: (AAC → ACG) y (ACG → CGG)

Otra lectura también tiene AACG:
    → la arista AAC→ACG aumenta su peso
    → NO se agrega ningún nodo nuevo
```

**Ensamblado = camino Euleriano** (recorre cada arista una vez). El euleriano tiene algoritmo eficiente O(E), mientras que el hamiltoniano es NP-completo → ventaja enorme.

**Filtrado de errores:** un k-mer real aparece en ~100 aristas (cobertura). Un k-mer generado por error de secuenciación tiene peso 1. Los k-mers con peso anormalmente bajo se descartan → simplifica el grafo.

**Problemas del grafo:**

| Estructura | Causa | Solución |
|---|---|---|
| **Bifurcación** | Error en extremo de lectura | El camino de menor peso es el error |
| **Burbuja** | Error en el medio o polimorfismo real | Si ≤2 diferencias: colapsar (error). Si más: variante real |
| **Ciclo/bucle** | Región repetitiva | Duplicar el nodo; cada copia queda en un contig separado |

**¿Cómo elegir k?**
- k pequeño: más conexiones, más sensible a errores, más ambigüedades.
- k grande: menos ambigüedades, pero necesita mayor cobertura.
- **SPAdes** (programa más usado): usa múltiples k iterativamente.

**N50:** longitud L tal que el 50% del genoma ensamblado está en contigs de largo ≥ L. A mayor N50, mejor.

**Regla fundamental del profe:** "villeteras mata galán" — un algoritmo brillante con poca cobertura siempre pierde frente a un algoritmo mediocre con mucha cobertura.

#### Paso 2 — Generar scaffolds (pair-end reads)

**El problema después del paso 1:** los contigs son islas. Sé la secuencia interna de cada isla, pero no sé cuál va antes ni a qué distancia.

**¿Qué es un pair-end read?** Al preparar la biblioteca, el ADN se fragmenta en trozos de tamaño controlado (2 kb, 10 kb, 50 kb). De cada fragmento se secuencian solo los dos extremos. El resultado es un par de lecturas que vienen del mismo fragmento original.

**Diagrama de scaffolding:**
```
Contig A                          Contig B
|============================|...hueco...|============================|
        ---->                                    <----
      lectura L                               lectura R
        └──────────── inserto ~2 kb ──────────────────┘

Si lectura L cae al final del Contig A
y su par (lectura R) cae al principio del Contig B
→ A y B están separados por ≈ (2000 - longL - longR) bases
→ se puede estimar el tamaño del hueco
```

**¿Por qué 3 tamaños de inserto?** (ej.: 2 kb, 10 kb, 50 kb)
- Insertos cortos (2 kb): unen contigs cercanos, scaffolding fino.
- Insertos medios (10 kb): saltan sobre repeticiones cortas que romperían los contigs.
- Insertos largos (50 kb): unen scaffolds grandes, ordenan regiones muy separadas.

**Gap-filling:** los huecos se intentan rellenar con lecturas sueltas pair-ends de lecturas ya presentes en los contigs adyacentes al hueco.

#### Contig vs scaffold vs cromosoma

| | Contig | Scaffold | Cromosoma |
|---|---|---|---|
| ¿Secuencia completa conocida? | Sí, base a base | Solo en partes de contigs | Sí (genoma finalizado) |
| ¿Hay huecos? | No | Sí, de **tamaño conocido** | No |
| ¿Se sabe posición en genoma? | No | No | Sí |
| Información que lo genera | Solapamiento de lecturas | Pair-end reads | Mapa físico (STS) |

---

## Clase 6 — Genómica Humana y Medicina Personalizada

### Proyecto Genoma Humano — Historia y dos estrategias

#### Contexto histórico

| Año | Evento |
|-----|--------|
| 1988-1990 | Congreso de EEUU aprueba el HGP; Watson es director inicial; objetivo: 95% del genoma en 15 años |
| 1990-1994 | Primera etapa: generar **mapa físico** del genoma (marcadores moleculares en posiciones conocidas de cromosomas) |
| 1993 | Se suma Sanger Center (UK); también Japón, Francia, Italia → consorcio internacional |
| ~1993 | Venter (NIH) inventa secuenciación de cDNAs (ESTs): obtiene secuencias de 2000 genes de cerebro en días |
| ~1994 | Escándalo de patentes: NIH intenta patentar los 2000 genes sin avisar a Venter → Watson renuncia; lo reemplaza Francis Collins |
| 1995 | Venter (ya fuera del NIH) publica **primer genoma bacteriano** usando Whole Genome Shotgun → prueba de concepto |
| 1998 | Venter funda **Celera Genomics**; anuncia que secuenciará el genoma humano más rápido que el consorcio |
| 2000 | Venter publica **genoma de Drosophila melanogaster** (primer pluricelular por WGS); el consorcio se asusta |
| 2001 | **Anuncio conjunto** (mediado por Clinton): consorcio público publica en **Nature**; Celera publica en **Science** |

> Watson era editor de Nature y según la leyenda había prohibido publicar trabajos de Venter → Venter publicó en Science.

#### Dos estrategias de ensamblado

**Consorcio público — BAC-by-BAC (jerárquico):**

```
Mapa físico
    ↓
Fragmentar genoma en BACs (cromosomas artificiales bacterianos, ~500kb – 1Mb)
    ↓
Secuenciar cada BAC por separado (en labs distribuidos, ~1000 labs)
    ↓
Ensamblar cada BAC individualmente (grafo pequeño → más fácil de resolver)
    ↓
Ordenar BACs usando el mapa físico → genoma completo
```

Cobertura: **1.5×** → 92% del genoma cubierto. Publicado en Nature.

**Celera — Whole Genome Shotgun (WGS):**

```
Fragmentar TODO el genoma en pedazos de ~500-1000 bases
    ↓
Secuenciar 27 millones de lecturas de ~500 pb
    (+ 3 tamaños de insertos: 2kb, 10kb, 50kb → para armar scaffolds)
    ↓
Filtrar lecturas basura (repeticiones, contaminación)
    ↓
Overlapper: alinear todas vs todas, umbral ≥40 pb con ≤6% diferencias
    ↓
Unitigger: grafo de solapamiento → unitigs (73.6% del genoma)
    ↓
Scaffolder: usar pares de lecturas (2kb, 10kb, 50kb) para vincular contigs
    ↓
Gap-filling: meter lecturas sueltas en los agujeros entre contigs
```

Además Celera tomó las secuencias públicas del consorcio, las fragmentó simulando lecturas de 500 pb → obtuvo 16M de lecturas adicionales → cobertura total: **8×** (vs 1.5× del consorcio).

**¿Por qué Celera obtuvo mejor resultado?** Por mayor **cobertura**, no por mejor algoritmo. A mayor cobertura, más lecturas por posición → más fácil resolver ambigüedades → menos gaps.

**Ventaja del enfoque BAC:** cada BAC genera un grafo pequeño (~10,000 lecturas vs 27 millones). Grafo más pequeño = más fácil de resolver = mejor calidad por BAC. Pero la menor cobertura global de 1.5× limitó el resultado final.

#### ¿Qué es el genoma humano de referencia?

El genoma de referencia **no es el genoma de ninguna persona particular** — es un **consenso** construido a partir de múltiples donantes (consorcio) o de un pequeño grupo (Celera: los 5 fundadores).

El genoma de referencia tiene dos componentes:
1. **La secuencia**: las ~3.000 millones de letras de A, T, C, G ordenadas por cromosoma
2. **El mapa de anotación**: para cada posición, qué función biológica tiene (exón, intrón, promotor, intergénica, ARN no codificante, etc.)

El genoma de referencia sirve como punto de comparación para determinar las **variantes** de cada individuo.

#### Estadísticas del genoma humano

| Dato | Valor |
|------|-------|
| Tamaño haploide | ~3 Gb (3.000 millones de bases) |
| Tamaño diploide | ~6 Gb |
| Genes codificantes | ~20,000–25,000 |
| Exones (% del genoma) | ~1.5% |
| Repeticiones | ~45–50% |
| Genes sin función conocida (al publicarse en 2001) | ~41–45% |
| Diferencia entre dos individuos cualesquiera | ~1 por mil → **1–5 millones de variantes** |

> La complejidad humana no viene del número de genes (el gusano C. elegans tiene ~20,000 genes también), sino del splicing alternativo, modificaciones post-traduccionales y redes de regulación.

---

### GWAS — Genome Wide Association Studies

**¿Qué es?** Un GWAS es un estudio epidemiológico a escala genómica. En lugar de secuenciar el genoma completo de cada individuo (costoso), se genotipifican ~500,000 a varios millones de SNPs en un chip (microarray de genotipado). Se compara la frecuencia de cada SNP entre casos (personas con la enfermedad) y controles (personas sin la enfermedad).

Para cada SNP, se calcula si su frecuencia es significativamente diferente entre casos y controles. Un SNP con p < 5×10⁻⁸ (umbral corregido por ~1 millón de comparaciones múltiples) se considera **asociado** con la enfermedad.

**Intuición estadística:** si en 10,000 personas con diabetes tipo 2 el alelo A en una posición aparece en el 40% de los cromosomas, pero en 10,000 controles sanos aparece solo en el 30%, esa diferencia es estadísticamente significativa → indica que en esa región del genoma hay algo relevante.

**¿Qué permite descubrir?**
- Regiones genómicas (loci) asociadas con enfermedades complejas (diabetes, esquizofrenia, cáncer, hipertensión, autismo).
- Cuantificar la heredabilidad de un rasgo: ¿cuánta variación poblacional se explica por variantes genéticas comunes?
- Construir **Polygenic Risk Scores (PRS)**: un puntaje de riesgo individual que combina la información de miles de SNPs asociados.

**Limitaciones:**
- El GWAS identifica regiones, no variantes causales directas.
- Detecta variantes comunes (MAF >1–5%) con efectos generalmente pequeños.
- Los primeros GWAS se hicieron casi exclusivamente en poblaciones europeas; las asociaciones pueden no transferirse a otras poblaciones.

---

### BWT (Burrows-Wheeler Transform) — Algoritmo de mapeo

**Problema:** tengo el genoma de referencia (3.000 millones de bases) y una lectura de ~100 bases. ¿Cómo encuentro en qué posición del genoma cae esa lectura de forma eficiente?

**Solución naive:** alinear la lectura contra cada posición del genoma → O(N×m) = inviable para millones de lecturas.

**Solución BWT:** **indexar el genoma una vez** → buscar cualquier sub-secuencia en O(m). Base de BWA y Bowtie2 (10-20× más rápidos que otros alineadores).

#### Construcción del índice BWT (ejemplo con "GOGOL")

```
Paso 1: Agregar símbolo $ al final y circularizar
Secuencia: G O G O L $  (se lee en círculo: después de $ vuelve G)

Paso 2: Generar todas las rotaciones (una por posición)
  0: G O G O L $
  1: O G O L $ G
  2: G O L $ G O
  3: O L $ G O G
  4: L $ G O G O
  5: $ G O G O L

Paso 3: Ordenar las rotaciones alfabéticamente
($ < G < L < O en este ejemplo)
  $ G O G O L     ← fila 5
  G O G O L $     ← fila 0
  G O L $ G O     ← fila 2
  L $ G O G O     ← fila 4
  O G O L $ G     ← fila 1
  O L $ G O G     ← fila 3

Paso 4: El arreglo de sufijos SA = posición original de cada fila = [5, 0, 2, 4, 1, 3]
Paso 5: La BWT = última columna = L $ O O G G
```

El **arreglo de sufijos (SA)** es el índice: mapea posición en la tabla ordenada → posición en la secuencia original.

#### Búsqueda de una sub-secuencia con BWT

Para buscar una sub-secuencia (lectura) en el índice:

```
Ejemplo: buscar "GO" en "GOGOL"

1. Recorrer el índice buscando la primera fila que tiene "GO" como prefijo → fila 1 (índice 1)
2. Buscar la última fila que tiene "GO" como prefijo → fila 2 (índice 2)
3. El intervalo [1, 2] en el SA contiene las posiciones: SA[1]=0, SA[2]=2
   → "GO" aparece en posiciones 0 y 2 de "GOGOL" ✓
```

Para una lectura de largo m: se agrega una letra a la vez, el intervalo [Rmin, Rmax] se va achicando. Cuando Rmin = Rmax → la lectura aparece en **un solo lugar** del genoma → se mapea ahí.

**Tolerancia a mismatches:** si al agregar una letra el intervalo se vuelve vacío (la letra no calza), se ignora esa letra (como si no estuviera) y se continúa. Al final, la lectura se mapea con ese mismatch tolerado. Típicamente se toleran 2-5 mismatches.

**¿Qué pasa si Rmin ≠ Rmax al final?** La lectura está en múltiples posiciones (región repetitiva). Los programas descartan esas lecturas o las ponen tentativamente con un flag de calidad de mapeo bajo.

#### Resumen BWT

| Aspecto | Detalle |
|---------|---------|
| Complejidad de indexar | O(N) — se hace una sola vez |
| Complejidad de buscar una lectura | **O(m)** — proporcional al largo del read |
| Programa que lo usa | BWA-MEM, Bowtie2 |
| Qué indexa | El genoma de referencia |

---

### Pipeline Variant Calling (Genómica Humana)

```
FASTQ (reads NGS)
    ↓  BWA-MEM (mapeo contra genoma referencia con BWT)
BAM (reads mapeados — cada lectura tiene posición en el genoma)
    ↓  Picard (marcar duplicados de PCR, métricas de cobertura)
BAM deduplicado
    ↓  GATK HaplotypeCaller (variant calling)
VCF (variantes crudas)
    ↓  Filtrado por calidad (VQSR o hard filters)
    ↓  Anotación estructural + funcional (VEP / SnpEff / ANNOVAR)
        usando: ClinVar, gnomAD, OMIM, PolyPhen-2, SIFT
Variantes anotadas → filtrado progresivo
```

#### Llamado de variantes — la lógica

Una vez mapeadas las lecturas, en cada posición del genoma se cuenta:
- Cuántas lecturas tienen la misma base que la referencia
- Cuántas lecturas tienen una base distinta (variante)

Para que una variante sea "llamada":
- Debe aparecer en al menos **X lecturas** (umbral de profundidad de variante)
- Para un heterocigota en organismo diploide: esperar ~50% de lecturas con la variante; en la práctica se pide ≥25% para tolerar errores de secuenciación
- El resultado es un archivo **VCF** (Variant Call Format): cada línea = una variante con su posición, bases REF/ALT, calidad, y profundidad

#### Tipos de variantes

| Tipo | Descripción | Impacto |
|------|-------------|---------|
| **SNP / SNV** | Cambio de 1 base | Depende del cambio de aa |
| **Missense** | SNP que cambia el aminoácido | Puede alterar función |
| **Nonsense** | SNP que genera stop prematuro | Trunca la proteína |
| **Sinónima (silent)** | SNP que no cambia aa | Generalmente neutro |
| **Indel** | Inserción o deleción 1+ bases | Frameshift si no es múltiplo de 3 |
| **Frameshift** | Indel no múltiplo de 3 → cambia todo el marco de lectura | Severo |
| **SV** | Segmentos >50 pb (duplicaciones, inversiones, translocaciones) | Variable |

> Un genoma humano difiere del de referencia en **1–5 millones de variantes** (1 en cada 1000 bases). De esas, solo unas pocas son relevantes para una condición particular.

#### Filtrado progresivo de variantes

```
Todas las variantes del genoma         ~1.000.000 – 5.000.000
    ↓ con efecto en proteína (missense, nonsense, frameshift...)
                                        ~10.000 – 50.000
    ↓ raras en población (MAF < 1% en gnomAD)
                                        ~1.000 – 5.000
    ↓ segregan con la enfermedad (análisis familiar)
                                        ~10 – 100
    ↓ patogénicas según ClinVar / literatura
                                        1 – 5  ← candidatas causales
```

#### Anotación de variantes

**Anotación estructural:** dónde cae la variante en la estructura del gen
- Exón codificante → ¿cambia el aminoácido? ¿genera stop? ¿frameshift?
- Intrón → ¿afecta sitio de splicing?
- Promotor → ¿afecta expresión?
- Intergénica → efecto regulatorio posible

**Anotación funcional:** bases de datos relevantes
| BD | Qué aporta |
|----|-----------|
| **ClinVar** | ¿Ya se reportó esta variante? ¿Es patogénica/VUS/benigna? |
| **gnomAD** | ¿Qué tan frecuente es en la población sana? (variante patogénica para enfermedad rara → debe ser muy rara) |
| **OMIM** | ¿El gen está asociado a enfermedades mendelianas? |
| **PolyPhen-2 / SIFT** | Predicción computacional del efecto de la sustitución de aminoácido |
| **COSMIC** | ¿Es una mutación somática recurrente en cáncer? |

---

### Evolución del costo de secuenciación

| Año | Costo por genoma | Tecnología |
|-----|-----------------|-----------|
| 2001 (HGP) | ~2.700 millones USD | Sanger, consorcio de 1000 labs |
| 2001 (Celera) | ~100 millones USD | WGS + Sanger |
| 2007 (Watson) | ~1 millón USD | 454 pyrosequencing |
| 2008 | ~60.000 USD | Primeras plataformas NGS (ABI) |
| 2010 | ~5.000 USD | Illumina HiSeq |
| 2014 | **~1.000 USD** | Illumina HiSeq X |
| ~2020 | ~500 USD | Illumina NovaSeq |

> La caída supera ampliamente la Ley de Moore. Consecuencia directa: más genomas secuenciados → más datos → mejores bases de datos de variantes → diagnóstico más preciso.

---

### Bioinformática Traslacional y Medicina Personalizada

**Concepto del profe:** la medicina personalizada no significa "solo mirar el genoma". Significa tener más información específica de ese individuo para personalizar diagnóstico, prevención y tratamiento: genotipo individual + fenotipo + ambiente + historia clínica. El genoma individual se puede agregar ahora como fuente de información adicional porque el costo bajó a niveles razonables.

| Área | Ejemplo |
|------|---------|
| **Diagnóstico** | WES/WGS para enfermedades raras; diagnóstico prenatal no invasivo |
| **Prevención** | Riesgo genético de cáncer (BRCA1/2); estado portador de enfermedades recesivas |
| **Tratamiento** | Oncología de precisión (mutaciones del tumor → elección de quimio); farmacogenómica |

#### Ejemplo 1 — HER2 y Herceptin (oncología de precisión)

El gen **HER2** codifica un receptor de membrana que promueve la proliferación celular. En ~15–20% de los cánceres de mama, el gen HER2 está amplificado en las células tumorales (muchas copias → sobreexpresión masiva → proliferación descontrolada).

**Herceptin (trastuzumab)** es un anticuerpo monoclonal que bloquea específicamente el receptor HER2. Es altamente efectivo en pacientes **HER2+**, pero completamente inútil en pacientes **HER2-**. Por eso, antes de prescribir Herceptin se realiza obligatoriamente genotipado del tumor (IHC + FISH) para determinar el estado HER2.

```
Biopsia tumoral → Determinación HER2 (IHC / FISH)
    ↓
HER2+  → Herceptin + quimioterapia → respuesta mucho mejor
HER2-  → Herceptin no indicado → protocolo diferente
```

#### Ejemplo 2 — BRCA1/2 y cáncer de mama/ovario

Las mutaciones de pérdida de función en **BRCA1** y **BRCA2** (genes de reparación del ADN) aumentan drásticamente el riesgo:

| | Población general | Portadora BRCA1 | Portadora BRCA2 |
|---|---|---|---|
| Cáncer de mama (lifetime) | ~12% | **55–72%** | **45–69%** |
| Cáncer de ovario | ~1.3% | **44%** | **17%** |

Las portadoras pueden optar por: vigilancia intensiva (mamografías + RMN anuales desde los 25 años), quimioprevención, o cirugía preventiva (mastectomía bilateral profiláctica, que reduce el riesgo de cáncer de mama en ~90%).

> **Caso público:** Angelina Jolie (2013) anunció su mastectomía preventiva al descubrir que era portadora de una mutación BRCA1.

#### Ejemplo 3 — Farmacogenómica (CYP2C19 y clopidogrel)

Clopidogrel es un anticoagulante para prevenir trombosis en stents coronarios. Es un pro-fármaco: necesita ser metabolizado por **CYP2C19** para volverse activo. El gen *CYP2C19* tiene variantes de pérdida de función relativamente comunes que reducen o eliminan la actividad de la enzima.

Un paciente con dos alelos de pérdida de función es un **metabolizador pobre**: no convierte el clopidogrel a su forma activa → el fármaco no funciona → mayor riesgo de trombosis. La FDA recomienda genotipado de CYP2C19 antes de prescribir clopidogrel.

#### Ejemplo 4 — Diagnóstico de enfermedades raras: caso Bainbridge et al. 2011

**Pacientes:** dos hermanos con enfermedad neuro-muscular severa progresiva. No respondían a L-Dopa (tratamiento estándar para deficiencia dopaminérgica).

**Estrategia bioinformática (secuenciación de trio: paciente + padres):**
1. Modelo genético → herencia **autosómica recesiva** (padres sanos, hijos enfermos)
2. WGS de los hermanos y los padres
3. Buscar variantes raras que segreguen: los padres heterocigotos, los hijos **heterocigotos compuestos**
4. Filtrar por efecto en proteína + rareza en gnomAD/dbSNP

**Diagnóstico:** heterocigotos compuestos para mutaciones en el gen **SPR** (Sepiapterin Reductase):
- alelo 1: p.Arg150Gly (missense)
- alelo 2: p.Lys251X (nonsense — stop prematuro)

Ninguno de los dos alelos produce SPR funcional → **sin BH4** (tetrahidrobiopterina) → sin Dopamina **ni Serotonina** (BH4 es cofactor de la Tirosina Hidroxilasa y la Triptófano Hidroxilasa).

**Por qué no respondían a L-Dopa:** L-Dopa restaura la dopamina pero no la serotonina. Sin BH4, ambas vías están bloqueadas.

**Tratamiento correcto:** L-Dopa + 5-OH-Triptófano (precursor de serotonina) + BH4.

> **Heterocigota compuesto:** el paciente tiene dos mutaciones distintas en el mismo gen, una en cada alelo (una heredada del padre, otra de la madre). A nivel de secuencia es heterocigoto en cada posición mutada; a nivel funcional es equivalente a homocigoto recesivo: ninguno de los dos alelos produce proteína funcional. Es la forma más común de herencia recesiva en poblaciones no consanguíneas.

**Filtrado progresivo para enfermedades raras:**
```
Todas las variantes del genoma:          ~4–7 millones
    ↓ Con efecto en proteína (missense, nonsense, frameshift, splicing)
                                         ~10,000–50,000
    ↓ Raras en población (MAF <0.1% en gnomAD para enfermedad rara)
                                         ~1,000–5,000
    ↓ En genes biológicamente plausibles para la enfermedad
                                         ~100–500
    ↓ Que segregan con la enfermedad (análisis familiar o de novo)
                                         ~5–50
    ↓ Clasificadas patogénicas en ClinVar o predicción deletérea fuerte
                                         1–5 candidatas causales
```

---

## Parcial 2018 — Resolución Completa

---

### Problema 1 — Programación Dinámica NW

**Secuencias:** CTGGCT (horizontal) y ATGCTG (vertical)
**Scoring:** Match = +2, Mismatch = 0, Gap = −1
**Algoritmo:** Needleman-Wunsch (alineamiento **global**)

#### a) Matriz de score completa

|   | − | C | T | G | G | C | T |
|---|---|---|---|---|---|---|---|
| **−** | 0 | −1 | −2 | −3 | −4 | −5 | −6 |
| **A** | −1 | 0 | −1 | −2 | −3 | −4 | −5 |
| **T** | −2 | −1 | 2 | 1 | 0 | −1 | −2 |
| **G** | −3 | −2 | 1 | 4 | 3 | 2 | 1 |
| **C** | −4 | −1 | 0 | 3 | 4 | 5 | 4 |
| **T** | −5 | −2 | 1 | 2 | 3 | **4** | **7** |
| **G** | −6 | −3 | 0 | 3 | 4 | **3** | **6** |

**Score final = 6** (celda esquina inferior derecha)

**Justificación de celdas clave:**

F(T,T) = F(5,6): T vs T → **match**
- Diagonal: F(4,5) + 2 = 5 + 2 = **7** ← elegido
- Arriba: F(4,6) − 1 = 4 − 1 = 3
- Izquierda: F(5,5) − 1 = 4 − 1 = 3

F(G,T) = F(6,6): G vs T → **mismatch**
- Diagonal: F(5,5) + 0 = 4
- Arriba: F(5,6) − 1 = 7 − 1 = **6** ← elegido
- Izquierda: F(6,5) − 1 = 3 − 1 = 2

F(G,C) = F(6,5): empate triple (diagonal=3, arriba=3, izquierda=3)

#### b) ¿Existe más de un alineamiento óptimo?

**Sí, existen exactamente dos.** El empate está en F(3,4) (fila G, columna G, cuarto elemento de fila 3):
- Diagonal: F(2,3) + 2 = 1 + 2 = 3 (G=G match)
- Izquierda: F(3,3) − 1 = 4 − 1 = 3

Ambos caminos son válidos → dos alineamientos óptimos.

**Alineamiento A** (camino diagonal en F(3,4)):
```
seq1 (CTGGCT): C  T  G  G  C  T  −
seq2 (ATGCTG): A  T  −  G  C  T  G
```

**Alineamiento B** (camino izquierda en F(3,4)):
```
seq1 (CTGGCT): C  T  G  G  C  T  −
seq2 (ATGCTG): A  T  G  −  C  T  G
```

#### c) Score de los alineamientos

Ambos tienen score = **6**.

**Verificación Alineamiento A:**
| Col | seq1 | seq2 | Score |
|-----|------|------|-------|
| 1 | C | A | mismatch = 0 |
| 2 | T | T | match = +2 |
| 3 | G | − | gap = −1 |
| 4 | G | G | match = +2 |
| 5 | C | C | match = +2 |
| 6 | T | T | match = +2 |
| 7 | − | G | gap = −1 |
**Total = 0+2−1+2+2+2−1 = 6** ✓

**Verificación Alineamiento B:**
| Col | seq1 | seq2 | Score |
|-----|------|------|-------|
| 1 | C | A | mismatch = 0 |
| 2 | T | T | match = +2 |
| 3 | G | G | match = +2 |
| 4 | G | − | gap = −1 |
| 5 | C | C | match = +2 |
| 6 | T | T | match = +2 |
| 7 | − | G | gap = −1 |
**Total = 0+2+2−1+2+2−1 = 6** ✓

---

### Problema 2 — Bases de Datos Primarias

Para cada categoría: clasificar según curación, redundancia, tipo (primaria/secundaria).

#### a) BD de nucleótidos: **RefSeq** (NCBI)

| Criterio | Clasificación | Justificación |
|----------|--------------|---------------|
| Curación | **Curada** | NCBI revisa y anota cada registro manualmente |
| Redundancia | **No redundante** | Un registro por molécula por organismo |
| Tipo | **Primaria** | Almacena secuencias crudas experimentales |

> Contraste: GenBank sería no curada y redundante.

#### b) BD de proteínas: **Swiss-Prot / UniProtKB reviewed**

| Criterio | Clasificación | Justificación |
|----------|--------------|---------------|
| Curación | **Curada** | Revisión manual por expertos; cada entrada validada con literatura |
| Redundancia | **No redundante** | 1 registro = 1 gen = 1 proteína; todas las isoformas en el mismo registro |
| Tipo | **Primaria** | Almacena datos directamente obtenidos de experimentos bioquímicos |

> Contraste: TrEMBL sería no curada y más redundante.

#### c) BD de expresión génica: **GEO** (NCBI)

| Criterio | Clasificación | Justificación |
|----------|--------------|---------------|
| Curación | **No curada** | Datos depositados por autores sin revisión sistemática del contenido biológico |
| Redundancia | **Redundante** | El mismo gen puede aparecer en miles de experimentos distintos; no hay criterio de unicidad |
| Tipo | **Primaria** | Almacena datos crudos de experimentos de expresión génica |

#### d) BD de organismo: **SGD** (Saccharomyces cerevisiae)

| Criterio | Clasificación | Justificación |
|----------|--------------|---------------|
| Curación | **Curada** | Equipo SGD revisa y actualiza manualmente la anotación funcional de cada gen |
| Redundancia | **No redundante** | Un registro por gen; variantes alélicas integradas en el mismo registro |
| Tipo | **Secundaria** | Integra y organiza datos de múltiples BD primarias (GenBank, UniProt, PubMed) + anotación funcional derivada |

#### e) BD de compuestos químicos: **PubChem** (NCBI)

| Criterio | Clasificación | Justificación |
|----------|--------------|---------------|
| Curación | **No curada** | Estructuras depositadas automáticamente desde múltiples fuentes; sin revisión manual sistemática |
| Redundancia | **Redundante** | El mismo compuesto puede aparecer bajo distintos nombres; PubChem intenta deduplicar (SID→CID) pero no elimina toda redundancia |
| Tipo | **Primaria** | Almacena datos de estructura química directamente obtenidos o sintetizados |

> Alternativa curada: **DrugBank** — curada, no redundante, secundaria.

#### f) BD de literatura: **Europe PMC** (EBI)

| Criterio | Clasificación | Justificación |
|----------|--------------|---------------|
| Curación | **Curada** | Los artículos han pasado por revisión editorial; metadatos (PMID, DOI, MeSH) revisados |
| Redundancia | **No redundante** | Cada artículo tiene un identificador único |
| Tipo | **Secundaria** | Recopila publicaciones que son en sí mismas productos secundarios del proceso científico |

#### g) BD libre: **PDB** (Protein Data Bank)

| Criterio | Clasificación | Justificación |
|----------|--------------|---------------|
| Curación | **Curada** | Deposición incluye validación estructural obligatoria (R-factor, clashscore, Ramachandran) |
| Redundancia | **Redundante** | La misma proteína puede tener cientos de estructuras (distintos ligandos, mutantes, condiciones) |
| Tipo | **Primaria** | Almacena datos crudos de experimentos estructurales (mapas de densidad electrónica, coordenadas atómicas) |

---

### Problema 3 — MSA, SP score, E-values, Prolina

**Contexto:** Motivo de N-glicosilación **N−X(≠P)−S/T**. MSA dado:
```
Prot1: L  E  N  K  T  V  A
Prot2: V  −  N  E  S  Y  A
Prot3: I  D  N  Q  T  I  A
```
**Penalidades gap affine:** apertura = −8, extensión = −3 → gap L=1: −8+(−3×1) = **−11**

#### a) Scores por par y árbol de similitud

**Par Prot1–Prot2** (S12):

| Pos | Prot1 | Prot2 | Score BLOSUM62 |
|-----|-------|-------|----------------|
| 1 | L | V | +1 |
| 2 | E | − | −11 (gap L=1) |
| 3 | N | N | +6 |
| 4 | K | E | +1 |
| 5 | T | S | +1 |
| 6 | V | Y | −1 |
| 7 | A | A | +4 |

**S(1,2) = 1−11+6+1+1−1+4 = 1**

**Par Prot1–Prot3** (S13):

| Pos | Prot1 | Prot3 | Score BLOSUM62 |
|-----|-------|-------|----------------|
| 1 | L | I | +2 |
| 2 | E | D | +2 |
| 3 | N | N | +6 |
| 4 | K | Q | +1 |
| 5 | T | T | +5 |
| 6 | V | I | +3 |
| 7 | A | A | +4 |

**S(1,3) = 2+2+6+1+5+3+4 = 23**

**Par Prot2–Prot3** (S23):

| Pos | Prot2 | Prot3 | Score BLOSUM62 |
|-----|-------|-------|----------------|
| 1 | V | I | +3 |
| 2 | − | D | −11 (gap L=1) |
| 3 | N | N | +6 |
| 4 | E | Q | +2 |
| 5 | S | T | +1 |
| 6 | Y | I | −1 |
| 7 | A | A | +4 |

**S(2,3) = 3−11+6+2+1−1+4 = 4**

**Resumen:**

| Par | Score |
|-----|-------|
| Prot1–Prot3 | **23** (más similares → se unen primero) |
| Prot2–Prot3 | 4 |
| Prot1–Prot2 | 1 (menos similares) |

**Árbol guía:**
```
        ┌── Prot1
   ┌────┤
   │    └── Prot3
───┤
   └──────── Prot2
```

**Justificación:** Prot1 y Prot3 comparten sustituciones conservativas sin gaps: L↔I (aliphatic), E↔D (ácidos), V↔I (aliphatic), T↔T (idéntico). Prot2 se aleja por el gap en posición 2 (−11) y la sustitución Y↔V poco conservativa (−1).

#### b) SP score total

SP = S(1,2) + S(1,3) + S(2,3) = 1 + 23 + 4 = **28**

Verificación columna por columna:

| Pos | S(1,2) | S(1,3) | S(2,3) | Subtotal |
|-----|--------|--------|--------|----------|
| 1 | 1 | 2 | 3 | 6 |
| 2 | −11 | 2 | −11 | −20 |
| 3 | 6 | 6 | 6 | 18 |
| 4 | 1 | 1 | 2 | 4 |
| 5 | 1 | 5 | 1 | 7 |
| 6 | −1 | 3 | −1 | 1 |
| 7 | 4 | 4 | 4 | 12 |
| **Total** | | | | **28** ✓ |

#### c) Ordenamiento de E-values: E12, E13, E23

E-value es exponencialmente inverso al score: E = K·m·n·e^(−λS).  
A mayor score → **menor** E-value (más significativo).

Scores: S13 = 23 > S23 = 4 > S12 = 1

**E₁₃ < E₂₃ < E₁₂**

El par 1–3 tiene el score más alto → menos probable de verse por azar → E-value más bajo = más significativo. El par 1–2 tiene score más bajo → más probable por azar → E-value más alto.

#### d) ¿Por qué la Prolina no se acepta en posición X del motivo?

**Razón estructural:** La Prolina es el único aa cuya cadena lateral forma un anillo con el nitrógeno del backbone (anillo pirrolidínico). Esto:

1. Elimina el H del N amídico → **no puede formar puentes H** con el backbone.
2. Restringe severamente el ángulo φ → **fuerza un quiebre (kink) en la estructura**.
3. **Impide que el tripéptido N-P-S/T adopte la conformación de β-turn** que necesita la enzima oligosacariltransferasa (OST) para reconocer y glicosilar el residuo N.

**¿Qué dice BLOSUM62?** La fila P en BLOSUM62 muestra scores negativos con prácticamente todos los aminoácidos (P−L = −3, P−I = −3, P−F = −4, P−W = −4, solo P−P = +7). Esto indica que la **Prolina es evolutivamente muy poco intercambiable** — sus sustituciones son seleccionadas negativamente porque alteran drásticamente la estructura local. La matriz refleja la excepcionalidad estructural de la prolina, que es exactamente la razón por la que elimina el sitio de glicosilación.

---

### Problema 4 — Bases de Datos Secundarias, PFAM y HMMer

#### A) Construcción del HMM con HMMer

**Problema bioinformático:** Dado un conjunto de sustratos confirmados de KC, construir un **profile-HMM** que capture la variabilidad posicional del motivo de reconocimiento para detectar nuevos sustratos en el proteoma de Av.

**Pipeline:**

```
1. Set de entrenamiento
   Secuencias de sustratos confirmados de KC
         ↓
2. MSA (CLUSTAL / MAFFT)
   Alinear las regiones de reconocimiento
         ↓
3. hmmbuild → KC_motivo.hmm
   Estima: P(aa | posición), P(match→insert), P(match→delete)
   = PROBLEMA DE ENTRENAMIENTO
         ↓
4. Calibración (automática con hmmbuild)
   Genera parámetros para calcular E-values
         ↓
5. hmmsearch → busca en proteoma de Av
   Lista de candidatos con score y E-value
   = PROBLEMA DE PUNTAJE
         ↓
6. Validación experimental (100 candidatos)
   Ensayos de unión y fosforilación
         ↓
7. Evaluación: curva ROC, sensibilidad, especificidad
```

**Métricas:**
- **Sensibilidad** = VP / (VP + FN): fracción de sustratos reales detectados.
- **Especificidad** = VN / (VN + FP): fracción de no-sustratos correctamente rechazados.
- Curva ROC: barrer umbrales de E-value → tradeoff sensibilidad/especificidad.

#### B) Interpretación del Logo del HMM

Un **sequence logo** representa cada posición del HMM como una columna de letras donde:
- **Altura total** = bits de información = conservación (máx ≈ 4.32 bits = log₂20)
- **Altura de cada letra** = frecuencia relativa (probabilidad de emisión) de ese aa en esa posición

**5 propiedades que se pueden leer:**

**1. Conservación posicional:** columnas altas = posiciones muy conservadas = **determinantes de especificidad** del reconocimiento por KC (si se mutan, eliminan la unión o la fosforilación).

**2. Identidad fisicoquímica del residuo dominante:** la letra más alta revela qué tipo de residuo requiere KC:
- K o R → quinasa basofílica (como PKA/PKC)
- D o E → quinasa acidofílica (como CK2)
- F, W, Y → bolsillo hidrofóbico aromático

**3. Tolerancia a sustituciones conservativas:** cuando una columna muestra 2–4 letras de altura apreciable, el modelo acepta variación conservativa. Si aparecen S y T juntas → KC solo requiere un hidroxilo; si K y R → requiere carga positiva, acepta cualquier básico.

**4. Posiciones variables:** columnas bajas = KC no impone restricciones ahí = residuos equivalentes a la "X" del motivo N-glicosilación.

**5. Extensión del motivo:** número total de columnas = largo del motivo lineal reconocido (20–25 aa en el enunciado).

#### C) Búsqueda en el proteoma de Av y falsos positivos

```bash
hmmsearch -E 0.001 KC_motivo.hmm proteoma_Av.fasta > candidatos.tbl
```

**Resultado:** Al evaluar 100 candidatos experimentalmente, 5 no se unen ni son fosforilados por KC.

**Nombre del resultado:** **Falsos Positivos (FP)**

El HMM predijo 100 proteínas como sustratos (positivos), pero 5 resultaron no serlo.
**Precisión = 95/100 = 95%.**

**¿Por qué ocurren los FP?**
1. El motivo existe en la secuencia pero está **enterrado en la estructura 3D** → KC no puede acceder in vivo.
2. Set de entrenamiento pequeño → el modelo aprendió señales espurias (sobreajuste).
3. **Contexto celular:** KC requiere co-localización, andamiajes o fosforilaciones previas no capturables por secuencia lineal.
4. Umbral de E-value permisivo → hits marginalmente similares al modelo.
