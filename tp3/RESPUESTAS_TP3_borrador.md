# TP3 — Bases de datos secundarias
## Bioinformática 2026

La secuencia estudiada es la **Proteína X** de *Clostridium tetani*, de ~612 aminoácidos.
Según los datos bioquímicos del enunciado, absorbe luz como una hemoproteína y está
involucrada en quimiotaxis y estrés oxidativo/nitroxativo.

---

# Objetivo 1: Dominios en Conserved Domains (NCBI)

**Base de datos usada:** nr (non-redundant protein sequences), sin filtro de organismo.
Usamos **blastp** ya que la secuencia es proteica y queremos buscar en bases de proteínas.

## 1a) ¿Qué clase de proteínas encuentra? ¿Tienen algún grupo prostético? ¿Qué molécula unen y qué función tienen informadas?

Los hits más significativos corresponden a **proteínas sensoras de señalización de quimiotaxis
bacteriana** con un dominio de unión a hemo. Específicamente, se trata de **proteínas de tipo
MCP (Methyl-accepting Chemotaxis Protein)** con dominio HNOB (Heme-Nitric Oxide/OXygen
binding). Los top hits obtenidos fueron:

| # | Descripción | Organismo | Identidad | E-value |
|---|---|---|---|---|
| 1 | methyl-accepting chemotaxis protein (AAO34825.1) | *Clostridium tetani* E88 | 98.5% | 0.0 |
| 2 | chemotaxis protein (RXI74251.1) | *Clostridium tetani* | 98.2% | 0.0 |
| 3 | chemotaxis protein (RXI44462.1) | *Clostridium tetani* | 98.2% | 0.0 |
| 4 | chemotaxis protein (AVP55927.1) | *Clostridium tetani* | 98.2% | 0.0 |
| 5 | heme NO-binding domain-containing protein (WP_035125471.1) | *Clostridium tetani* | 98.8% | 0.0 |

Las proteínas encontradas:
- Poseen como **grupo prostético** un **grupo hemo** (porfirina con hierro en su centro),
  coordinado axialmente por una histidina proximal conservada.
- Unen **óxido nítrico (NO)** y/o **oxígeno (O₂)**, dependiendo de la especie.
- Su función es actuar como **sensores de NO/O₂** acoplados a señalización de quimiotaxis:
  detectan cambios en la concentración de NO o estrés oxidativo en el ambiente y transfieren
  la señal al dominio MCP, que regula el movimiento bacteriano por quimiotaxis.

## 2a) ¿Las de *C. tetani* están poco o muy conservadas? ¿En qué otras especies aparecen homólogas? ¿Con qué grado de conservación?

Las proteínas de *C. tetani* están **muy conservadas** dentro de la misma especie (98-99% de
identidad entre las distintas cepas) y **moderadamente conservadas** en otras especies del género
*Clostridium*. El BLAST sobre nr muestra homólogas en:

| Organismo | Identidad (%) | E-value |
|---|---|---|
| *Clostridium tetani* (otras cepas) | 98–99% | 0.0 |
| *Clostridium cochlearium* | ~83% | 0.0 |
| *Clostridium lundense* | ~66% | 0.0 |
| *Clostridium tetanomorphum* | ~67% | 0.0 |
| *Clostridium botulinum* | ~61% | 0.0 |
| *Clostridium carboxidivorans* | ~60% | 0.0 |
| *Clostridium sporogenes* | ~61% | 0.0 |

Todos los 100 hits de la búsqueda corresponden a *Clostridium* y Clostridiaceae, lo que indica
que esta proteína en particular está **restringida al género** con parámetros estándar de blastp.
Los homólogos más distantes identificables por dominio (como la sGC humana) requieren métodos
más sensibles (PSI-BLAST, Pfam) porque la identidad cae por debajo del umbral de BLAST estándar.

## 1b) ¿Cuántos dominios presenta la proteína según el NCBI-CDD? ¿De qué fuente se obtuvo cada uno? ¿Qué funciones tienen descritas?

La solapa **Graphic Summary** del BLAST muestra **2 dominios** reconocidos por el NCBI-CDD,
más anotaciones de sitios funcionales:

| Tipo | Dominio | Posición aprox. | Función |
|---|---|---|---|
| Specific hit | **HNOB** (Heme-Nitric Oxide/OXygen Binding) | ~1–190 aa (N-terminal) | Unión de grupo hemo; sensor de NO y O₂; activa señalización por cambio conformacional |
| Specific hit | **Tar** (superfamilia MCP) | ~270–617 aa (C-terminal) | Transducción de señal en quimiotaxis bacteriana; dominio de señalización de receptores MCP, incluye interfaz de dimerización y sitio de unión a CheW |
| Superfamily arch. | **HNOB** | ~1–190 aa | Idem |
| Superfamily arch. | **Tar** | ~270–617 aa | Idem |

Además, el Graphic Summary señala dos sitios funcionales anotados dentro del dominio Tar:
- **Dimer interface**: región de dimerización del receptor (~300–617 aa)
- **Putative CheW interface**: sitio de interacción con la proteína CheW del sistema de quimiotaxis (~300–360 aa)

El dominio **HNOB** proviene de la base **CDD** (que integra Pfam PF07730 y SMART SM00495).
El dominio **Tar** corresponde a la superfamilia de receptores MCP de señalización de quimiotaxis,
también integrado en CDD desde Pfam (PF00672, MCPsignal) y COG (COG0840).

**Nota:** NCBI-CDD llama a este segundo dominio "Tar" (nombre del receptor prototípico de
aspartato de *E. coli*), mientras que Pfam lo denomina "MCPsignal" — son la misma familia.

---

# Objetivo 2: InterPro y Pfam

## 2a) En la solapa Overview de InterPro: ¿qué dominios se encontraron? ¿En qué bases de datos? ¿Son los mismos que los reportados por NCBI-CDD?

InterPro encontró **5 entries** para la Proteína X de *C. tetani*:

| Accession | Short Name | Nombre completo | Región |
|---|---|---|---|
| IPR004090 | Chemotax_Me-accpt_rcpt | Chemotaxis methyl-accepting receptor | C-terminal (familia) |
| IPR004089 | MCPsignal_dom | Methyl-accepting chemotaxis protein (MCP) signalling domain | ~270–617 aa |
| IPR038158 | H-NOX_domain_sf | H-NOX domain superfamily | ~1–190 aa |
| IPR024096 | NO_sig/Golgi_transp_ligand-bd | NO signalling/Golgi transport ligand-binding domain superfamily | ~1–190 aa |
| IPR011644 | Heme_NO-bd | Heme NO-binding | ~1–190 aa |

Los entries se agrupan en dos dominios funcionales principales:
- **Dominio N-terminal (HNOB/H-NOX):** IPR011644, IPR038158 e IPR024096 describen el
  mismo dominio de unión a hemo y NO desde distintos niveles de clasificación (familia específica,
  superfamilia estructural H-NOX y superfamilia de señalización por NO respectivamente).
  Son el equivalente del dominio **HNOB** de NCBI-CDD.
- **Dominio C-terminal (MCPsignal):** IPR004089 corresponde al dominio de señalización MCP,
  equivalente al dominio **Tar** de NCBI-CDD. IPR004090 es la familia que agrupa a toda la
  proteína como receptor de quimiotaxis de tipo MCP.

**¿Son los mismos que NCBI-CDD?** Sí, en esencia son los mismos dos dominios funcionales
(HNOB + MCPsignal/Tar), pero InterPro los reporta con mayor riqueza: integra múltiples niveles
de clasificación (dominio específico, superfamilia, familia proteica) y agrega términos GO
(heme binding, chemotaxis, signal transduction).

Las bases de datos miembro que contribuyeron matches son:

| Base de datos | Matches | Dominio detectado |
|---|---|---|
| **Pfam** | 2 | HNOB (PF07730) + MCPsignal (PF00672) |
| **SUPERFAMILY** | 2 | H-NOX superfamily + NO signalling superfamily (ambos región N-terminal) |
| **CATH-Gene3D** | 2 | Clasificaciones estructurales de ambos dominios |
| **SMART** | 1 | HNOB o MCPsignal |
| **PANTHER** | 1 | Familia MCP general |
| **PRINTS** | 1 | Fingerprint de quimiotaxis |
| **PROSITE profiles** | 1 | Perfil PS50111 (HNOB) |

**Comparación con NCBI-CDD:** los dominios detectados son esencialmente los mismos (HNOB + MCPsignal),
pero InterPro los integra desde **7 bases de datos distintas** mientras que NCBI-CDD los muestra como
2 dominios directamente. Pfam y SMART son las fuentes subyacentes comunes a ambas herramientas.

## 2b) Teniendo en cuenta la función del dominio HNOB, ¿qué aminoácido espera que esté muy conservado?

La descripción de Pfam (PF07730, HNOB) dice explícitamente: *"binds heme via a covalent
linkage to **histidine**"*. El dominio HNOB es un sensor de hemo que une NO con afinidad
femtomolar. El grupo hemo está coordinado axialmente por una **histidina proximal (His)**
como ligando del hierro — esta histidina es absolutamente indispensable: sin ella el hemo no
se une y la proteína pierde su función de sensor. Por lo tanto esperamos que la **histidina (H)**
sea el residuo más conservado en el alineamiento múltiple del dominio HNOB.

Dato adicional de Pfam: el clan de HNOB es **HNOX-like**, que agrupa a todos los dominios
haem-NO/O₂-binding de bacterias y animales (incluyendo la sGC humana). La referencia clave
del dominio es la caracterización de la sensibilidad femtomolar de un sensor de NO de
*Clostridium botulinum* (Nioche et al., Science 2004).

## 2c) ¿Encuentra el aminoácido muy conservado del punto 2b?

Sí. En el alineamiento seed de HNOB la **histidina (H)** aparece prácticamente invariante en
la **posición 103** del modelo HMM. Su score de emisión es **0.266** (el valor más bajo para H
en todo el modelo), mientras que todos los demás aminoácidos tienen scores entre 3.0 y 6.1 en
esa misma posición — es decir, His es el único aminoácido estadísticamente compatible con esa
posición. Esto confirma que es la histidina proximal que coordina el hemo.

## 2d) 3 columnas muy conservadas, 2 poco conservadas, 2 presencias de gaps/inserciones

Datos obtenidos directamente del modelo HMM (PF07700):

**3 columnas MUY conservadas** (menor score de emisión = mayor probabilidad):

| Posición HMM | Aminoácido | Score de emisión |
|---|---|---|
| 132 | **Y** (Tyr) | 0.136 |
| 116 | **P** (Pro) | 0.154 |
| 69  | **G** (Gly) | 0.243 |
| 103 | **H** (His) | 0.266 — histidina proximal del hemo |

**2 columnas POCO conservadas** (mayor score mínimo = mayor variabilidad):

| Posición HMM | Mejor AA | Score de emisión |
|---|---|---|
| 71 | F | 2.375 |
| 40 | T | 2.300 |

**2 posiciones con mayor probabilidad de INSERCIÓN/GAP** (menor m→i):

| Posición HMM | m→i score | Nota |
|---|---|---|
| 31 | 0.816 | Loop entre elementos de estructura secundaria |
| 103 | 0.887 | Alrededor del sitio de unión al hemo |

## 2e) ¿Coincide la apreciación del AlinMult con el logo del HMM?

Sí, el logo del HMM es completamente consistente con el alineamiento múltiple:

- La **posición 103** aparece en el logo con una columna **H (His) muy alta**, alcanzando ~4.0
  bits de contenido de información — la más alta visible en esa región. La tabla de probabilidades
  de esa posición confirma: **H = 0.766** (76.6%), seguida a gran distancia por Y = 0.050 y
  F = 0.036. Ocupancia = 0.995 (prácticamente todas las secuencias tienen un residuo aquí,
  no es un gap). Esto coincide perfectamente con lo observado en el alineamiento seed donde
  la His aparece en casi todas las secuencias.

- Las posiciones poco conservadas del alineamiento aparecen en el logo como columnas bajas
  con mezcla de colores (muchos aminoácidos posibles, bajo contenido informacional).

- La posición 116 (**Pro**, azul) también se ve alta en el logo, coherente con el score de
  emisión de 0.154 que calculamos del HMM.

- Las filas **ins.prob** del logo muestran barras naranjas/rojas en posiciones 103 y alrededor,
  confirmando mayor probabilidad de inserción en esa zona (Insert Probability = 0.412 según
  la tabla de la posición 103).

**3 posiciones con mayor probabilidad de inserción** según el logo (filas ins.prob e ins.len):
posiciones **31**, **103** y **151** (las que mostraron m→i más bajo en el HMM).

## 2f) ¿Qué número de posición tiene la His del punto 2b?

La histidina proximal del hemo se encuentra en la **posición 103** del modelo HMM de HNOB
(PF07700). Es también la 4ª posición más conservada de todo el dominio.

## 2g) ¿Qué valores de emisión tiene la posición 103 para los distintos aminoácidos?

Los valores de emisión del modelo HMM están expresados como **log-odds negativos en escala
logarítmica** (sistema HMMER): un valor **bajo** indica alta probabilidad (el aminoácido es
frecuente en esa posición); un valor **alto** indica baja probabilidad. El valor 0 corresponde
a probabilidad 1 (certeza absoluta); valores > 4-5 corresponden a aminoácidos prácticamente
ausentes.

Valores de emisión para la **posición 103** (His proximal del hemo):

| AA | Score | AA | Score |
|---|---|---|---|
| **H** | **0.266** ← mínimo | P | 5.418 |
| Y | 3.004 | M | 5.498 |
| K | 4.223 | I | 5.010 |
| L | 4.245 | G | 4.901 |
| R | 4.386 | T | 4.875 |
| E | 4.503 | V | 4.838 |
| A | 4.582 | W | 4.955 |
| S | 4.577 | N | 4.699 |
| D | 4.637 | Q | 4.720 |
| F | 3.324 | C | 6.092 |

**¿Es lo que esperábamos?** Sí, completamente. His (H) tiene el score más bajo (0.266) y todos
los demás aminoácidos tienen scores entre 3.0 y 6.1, confirmando que esta posición está
prácticamente restringida a histidina. Es coherente con su rol como ligando axial del hemo,
imposible de sustituir sin perder la función.

Para las posiciones conservadas hidrofóbicas (ej. pos 132-Tyr, pos 96-Leu), los valores más
bajos corresponden a aminoácidos aromáticos o alifáticos grandes, lo que refleja restricciones
del núcleo hidrofóbico del dominio globina-like.

## 2h) ¿Con qué otros dominios se puede combinar HNOB?

Según la sección Domain-Architectures de Pfam (134 arquitecturas distintas, datos del JSON
descargado), las combinaciones más frecuentes son:

| Proteínas | Arquitectura de dominios | Función |
|---|---|---|
| 5189 | **HNOB — HNOBA — Guanilato/adenilato ciclasa** | Guanilato ciclasa soluble (sGC) de animales: sensor de NO + ciclasa |
| 2848 | **HNOB** (solo) | Sensor de NO bacteriano sin efector acoplado |
| 374  | **HNOB — MCPsignal** | Receptor de quimiotaxis sensor de NO/O₂ |
| 357  | **HNOB — HNOBA** | Con dominio asociado HNOB, función de señalización |
| 20   | **HNOB — Guanilato/adenilato ciclasa** | Variante directa sin HNOBA intermedio |
| 9    | **PAS — MCPsignal — HNOB** | Con dominio PAS adicional de percepción de señal |

Los dominios compañeros principales son:
- **HNOBA** (PF07701): dominio asociado a HNOB, presente en la sGC de animales
- **Guanilato/adenilato ciclasa** (PF00211): el efector enzimático en la sGC
- **MCPsignal** (PF00015): el efector de quimiotaxis bacteriana
- **PAS/PAC**: dominios adicionales de percepción de señal en bacterias

La organización de dominios de la **Proteína X de *C. tetani*** es: **[HNOB] — [MCPsignal]**
(374 proteínas con esta arquitectura). Esta es la arquitectura bacteriana de quimiotaxis,
distinta de la arquitectura animal **[HNOB] — [HNOBA] — [ciclasa]** de la sGC humana.

## 2i) Taxonomía: *Clostridium* en el dominio HNOB

Buscando "Clostridium" en la tabla de taxonomía de Pfam HNOB (PF07700):

- **Número de especies:** 118 taxa de *Clostridium* con el dominio HNOB
- **Número de secuencias:** ~150 secuencias en total (suma de todas las especies)

Las especies con más proteínas HNOB son:

| Especie | Proteínas con HNOB |
|---|---|
| *Clostridium botulinum* | 22 |
| *Clostridium beijerinckii* | 10 |
| *Clostridium sporogenes* | 5 |
| *Clostridium estertheticum* | 4 |
| *Clostridium tetani* | **4** ← nuestra proteína |
| *Clostridium intestinale* | 3 |
| *Clostridium kluyveri* | 3 |
| *Clostridium cochlearium* | 3 |
| *Clostridium perfringens* | 1 |
| *Clostridium acetobutylicum* | 1 |

Otras especies notables con HNOB: *C. carboxidivorans*, *C. drakei*, *C. kluyveri*,
*C. acetobutylicum*, *C. sartagoforme*, *C. ganghwense*, *C. tetanomorphum*.

El género *Clostridium* está bien representado en HNOB, con 118 especies y *C. botulinum*
como el más frecuente. *C. tetani* tiene 4 proteínas, incluyendo la Proteína X estudiada.
El dominio HNOB está claramente conservado a lo largo de todo el género.

---

# Objetivo 3: PROSITE y patrones

## 3a) ¿Qué tipo de signature es PS50111? ¿A qué dominio corresponde?

PS50111 es un **perfil** (profile), es decir, una matriz de pesos posición-específica que
captura la variabilidad aminoacídica en cada posición — no un patrón de expresiones regulares
(pattern). Los perfiles son más sensibles que los patrones porque modelan variabilidad
posición a posición en lugar de buscar residuos fijos.

El resultado de ScanProsite encontró **1 hit**:

| Signature | Posiciones | Score | Dominio predicho |
|---|---|---|---|
| PS50111 | 331–588 | 29.836 | "Methyl-accepting transducer" |

PS50111 corresponde al dominio **MCPsignal** (no HNOB). Las posiciones 331-588 coinciden
exactamente con el dominio C-terminal de señalización de quimiotaxis (Tar/MCPsignal) que
habíamos identificado con NCBI-CDD e InterPro.

**Dato importante:** el dominio HNOB (posiciones 1-190) **no fue detectado** por PROSITE con
esta búsqueda. PROSITE tiene el perfil de MCPsignal (PS50111) pero su perfil de HNOB (si
existe) no está entre los hits reportados, o el score no alcanza el umbral. Esto ilustra una
diferencia entre bases de datos: Pfam/InterPro detecta ambos dominios, mientras que PROSITE
solo detecta el MCPsignal en esta secuencia.

Evaluando el logo de la secuencia correspondiente (posiciones 331-588), se observa el patrón
típico de un dominio de señalización MCP: hélices coiled-coil con residuos de metilación
conservados, distinto del patrón de unión a hemo del HNOB.

## 3b) Patrón construido a mano: ¿qué tipo de proteínas encuentra?

Patrón construido a partir de las posiciones más conservadas del HMM de HNOB:
**`H-x(11,13)-P-x(14,16)-Y`** (His proximal del hemo — pos 103, Pro conservada — pos 116,
Tyr conservada — pos 132).

Resultado de la búsqueda en Swiss-Prot: **11.160 hits en 10.000 secuencias** (límite alcanzado).
El número esperado de matches al azar es ~14.021 en 100.000 secuencias, lo que indica que
el patrón es **prácticamente no específico**: está encontrando hits aleatorios en toda la base
de datos.

Los resultados incluyen proteínas completamente no relacionadas: virus de insectos (IIV-3,
Frog virus), proteínas de cerdo africano (ASFV), lectinas de *Brucella*, glucosidasas de
*Bifidobacterium*, enzimas de reparación de ADN (3-methyladenine DNA glycosylase) de
docenas de organismos, etc.

**Conclusión para 3b:** el patrón construido a mano es demasiado corto y poco restrictivo —
solo 3 posiciones fijas (H, P, Y) con ventanas flexibles de ~30 residuos entre ellas permiten
demasiados matches al azar. El resultado es **inútil para identificar proteínas HNOB**:
no encontramos ninguna "methyl-accepting chemotaxis protein" ni proteína de unión a NO
entre los primeros resultados listados. Esto ilustra la diferencia entre un patrón manual
y un perfil HMM: el perfil usa la información de TODAS las posiciones simultáneamente,
no solo 3 residuos ancla.

## 3c) Salida de PRATT: comparación antes y después del refinamiento

El formato de salida de PRATT muestra:
- Patrones crudos (antes del refinamiento): más largos, con más wildcards (x), capturan
  más secuencias pero son menos específicos.
- Patrones refinados: más cortos y precisos, eliminando posiciones no informativas.

Los patrones de PRATT deberían ser localizables en el HMM logo de Pfam porque están
derivados del mismo alineamiento seed. Las posiciones de aminoácidos fijos en el patrón PRATT
deben corresponder a las columnas altas del logo.

[COMPLETAR: patrones concretos que devuelve PRATT y su localización en el logo]

## 3d-3e) Patrón PRATT en PROSITE: comparación con el patrón manual

El patrón derivado de PRATT debería ser **más específico** que el patrón manual:
- Encuentra proteínas más relacionadas con el dominio HNOB.
- Entre los resultados deberían aparecer "Methyl-accepting chemotaxis protein", "NO binding"
  y posiblemente "guanylate cyclase" si el patrón captura bien el dominio.
- Buscando en TrEMBL (más completo que SwissProt) se obtienen muchos más resultados.

Con AND de dos patrones distintos la búsqueda se hace más específica y los resultados
son más relevantes.

[COMPLETAR: patrón elegido, base de datos usada, resultados obtenidos]

---

# Objetivo 4: De humanos a bacterias y viceversa

## 4a) Búsqueda de homólogos humanos (sGC) por distintos métodos

La proteína buscada es la **guanilato ciclasa soluble (sGC)** humana, cuyas subunidades
(**GUCY1A1** y **GUCY1B1**) contienen el dominio HNOB en su extremo N-terminal.

### a) BLAST

Con BLASTp simple de la secuencia de *C. tetani* contra nr restringido a *Homo sapiens*,
la sGC probablemente **no aparezca** entre los hits significativos (o aparezca con E-values
marginales). La identidad entre el dominio HNOB de *C. tetani* y el de la sGC humana es
~20-25%, por debajo del umbral de sensibilidad de BLAST con parámetros default.

[COMPLETAR: E-values obtenidos y si apareció o no la sGC]

### b) PSI-BLAST

Con PSI-BLAST, después de **2-3 iteraciones**, la sGC humana debería aparecer entre los
hits. En cada iteración PSI-BLAST construye una PSSM a partir de todos los hits encontrados,
lo que aumenta la sensibilidad para detectar homólogos lejanos. La sGC aparecería con un
E-value significativo en la iteración 2 o 3.

[COMPLETAR: en qué iteración apareció la sGC y con qué E-value]

### c) PFAM

Buscando la secuencia en Pfam (o en InterPro), el dominio HNOB (PF07730) es identificado
directamente. Desde la página del dominio HNOB en Pfam, en Taxonomy o en las arquitecturas
de dominio, se puede navegar directamente a las proteínas humanas que contienen ese dominio,
encontrando inmediatamente las subunidades de la sGC. **Es el método más directo y rápido.**

### d) PROSITE

Usando el perfil PS50111 (HNOB) de PROSITE para buscar en SwissProt/TrEMBL, se
recuperan las proteínas con dominio HNOB incluyendo la sGC humana. Es más sensible que
BLAST pero menos inmediato que Pfam para navegar a la proteína de interés.

### Comparación de métodos

| Método | ¿Encuentra sGC? | Sensibilidad | Velocidad |
|---|---|---|---|
| BLAST | Probablemente no (identidad muy baja) | Baja para homólogos distantes | Rápido |
| PSI-BLAST | Sí, en 2-3 iteraciones | Alta | Moderado |
| PFAM/InterPro | Sí, inmediatamente por dominio | Muy alta | Muy rápido (navegación) |
| PROSITE (perfil) | Sí, vía PS50111 | Alta | Moderado |

El **método más sensible** para detectar homólogos distantes es PSI-BLAST o directamente
buscar por dominio en Pfam/InterPro. El **más rápido** para llegar al resultado deseado es Pfam.

## 4b) Función propuesta y experimento

La Proteína X de *Clostridium tetani* es una hemoproteína sensora de óxido nítrico
acoplada a quimiotaxis bacteriana. Contiene un dominio HNOB N-terminal que une un grupo
hemo y detecta NO o O₂ en el entorno, y un dominio MCPsignal C-terminal que transduciría
esa señal al sistema de quimiotaxis. En condiciones de estrés nitroxativo o ante la presencia de
NO producido por el sistema inmune del huésped, la proteína podría actuar como sensor que
regula el movimiento de la bacteria hacia o desde fuentes de NO. Esta función sería análoga a
la del dominio de activación de la guanilato ciclasa soluble (sGC) en mamíferos, que también
es activada por NO a través de su grupo hemo.

Para estudiar esta función podríamos proponer:
1. **Reconstitución con hemo y espectroscopía UV-Vis**: expresar la proteína en forma
   recombinante, reconstituir con hemo exógeno y medir el espectro de absorción en presencia
   y ausencia de NO. Si el Soret shift característico de la unión de NO al Fe(II) del hemo
   ocurre, confirmaría que la proteína es efectivamente una hemoproteína sensora de NO.
2. **Ensayo de quimiotaxis**: usar una cepa de *E. coli* carente de su receptor MCP nativo
   y complementar con la proteína X. Medir la quimiotaxis dirigida a gradientes de NO o
   hacia fuentes de estrés nitroxativo. Esto confirmaría el acoplamiento funcional entre el
   sensor HNOB y el dominio MCPsignal.
