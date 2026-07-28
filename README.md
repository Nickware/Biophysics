# Biofísica

La **biofísica** es una ciencia interdisciplinaria que aplica los principios y métodos de la física para estudiar los sistemas y procesos biológicos en todos los niveles, desde lo molecular hasta lo ecológico. Busca explicar cómo funcionan los organismos vivos utilizando leyes físicas, modelos matemáticos y técnicas experimentales precisas. Entre sus áreas de estudio principales se encuentran la estructura y función de las proteínas, membranas biológicas, dinámica molecular, bioenergética, biomecánica y bioelectromagnetismo, entre otras.[^1][^2]

Este repositorio reúne herramientas computacionales usadas en biofísica molecular. Este documento es el índice: describe brevemente cada herramienta y su rol en el flujo de trabajo, y enlaza al README dedicado de cada una, donde vive el detalle técnico (instalación, API verificada, ejemplos).

## Biofísica y la simulación computacional

La biofísica moderna saca enorme provecho de la simulación y el análisis computacional para comprender mecanismos complejos, como el plegamiento de proteínas, interacciones entre biomoléculas, dinámica de membranas y mecanismos enzimáticos. Las herramientas de este repositorio complementan ese objetivo desde tres roles distintos: **simulación**, **análisis** y **visualización**.

---

## Herramientas del repositorio

| Herramienta | Descripción breve | Relación con la biofísica | Documentación |
|---|---|---|---|
| **ESPResSo** (`espressomd`) | Motor de dinámica molecular (MD) para sistemas de materia blanda: polímeros, coloides y biomoléculas coarse-grained. | Permite modelar sistemas biológicos "de grano grueso" (membranas, ADN, proteínas), dando información sobre procesos físicos y propiedades emergentes en biología molecular.[^3] | [`/espressomd/README.md`](./espressomd/README.md) |
| **MDAnalysis** | Librería en Python para analizar trayectorias de dinámica molecular de múltiples motores de simulación (GROMACS, NAMD, AMBER, entre otros). | Extrae, cuantifica y visualiza propiedades físicas de biomoléculas simuladas: conformación, difusión, interacciones intermoleculares.[^4] | [`/mdanalysis/README.md`](./mdanalysis/README.md) |
| **PyMOL** | Software de visualización y análisis 3D de estructuras moleculares. | Indispensable en biofísica estructural para visualizar proteínas, ácidos nucleicos y complejos moleculares, y para producir figuras de resultados.[^5] | [`/pymol/README.md`](./pymol/README.md) |
| **MMTK** (Molecular Modelling Toolkit) | Biblioteca en Python para simulación molecular: dinámica molecular, minimización de energía, análisis de modos normales, con varios campos de fuerza para biomoléculas (Amber 94/99, modelos de red elástica). | Simula el comportamiento dinámico de proteínas, ácidos nucleicos y membranas para responder preguntas de dinámica, energía y estructura.[^6] | *(sin README dedicado en este repositorio todavía)* |
| **pDynamo** | Biblioteca open-source para simulaciones con funciones de energía de química cuántica, mecánica molecular y combinadas (QC/MM). | Modela reacciones químicas enzimáticas a nivel atómico, combinando técnicas cuánticas y clásicas: catálisis, mutaciones, rutas de reacción.[^7] | *(sin README dedicado en este repositorio todavía)* |

---

## Cómo se combinan en una investigación

- **Simulación:** ESPResSo, MMTK y pDynamo generan trayectorias y predicen movimientos y cambios de estado de las moléculas.
- **Análisis:** MDAnalysis procesa esas trayectorias y extrae información cuantitativa (variaciones estructurales, contactos atómicos, fluctuaciones).
- **Visualización:** PyMOL traduce los resultados en representaciones 3D para el análisis visual y la comunicación científica.

Para ejemplos de código verificados sobre cómo conectar estas herramientas entre sí (por ejemplo, el puente correcto `espressomd.MDA_ESP` entre ESPResSo y MDAnalysis, o el flujo MDAnalysis → PyMOL para figuras), ver la sección de integración dentro de cada README dedicado — no se duplica aquí para evitar que el código quede desactualizado en dos lugares a la vez.

---

## Nota sobre esta revisión

Una versión anterior de este README citaba, bajo la fila de MMTK, dos referencias (`mmtk.io` y `github.com/mmtk/mmtk-core`) que en realidad corresponden a un proyecto homónimo no relacionado: **"MMTk" (Memory Management Toolkit)**, un framework de investigación en recolección de basura para runtimes de lenguajes de programación (usado, por ejemplo, en el compilador de Rust), sin relación con simulación molecular. Se corrigió para apuntar únicamente a la fuente correcta del MMTK de biofísica. También se retiraron dos referencias sin relación real con el software (dominios que solo coincidían en la palabra "espresso"), y se corrigió la ortografía de `espressomd` (no "expressomd") en toda la tabla.

## Referencias

[^1]: Biofísica - Wikipedia, la enciclopedia libre. https://es.wikipedia.org/wiki/Biofísica
[^2]: Biofísica: qué es, áreas de estudio y aplicaciones - Ferrovial. https://www.ferrovial.com/es/stem/biofisica/
[^3]: ESPResSo » Extensible Simulation Package for the Research on Soft Matter. https://espressomd.org/wordpress/
[^4]: MDAnalysis. https://www.mdanalysis.org
[^5]: PyMOL - Wikipedia. https://en.wikipedia.org/wiki/PyMOL
[^6]: The Molecular Modelling Toolkit. http://dirac.cnrs-orleans.fr/MMTK.html
[^7]: pDynamo. https://www.pdynamo.org
