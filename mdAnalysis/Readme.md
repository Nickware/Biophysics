# MDAnalysis

MDAnalysis es una librería de Python de código abierto diseñada específicamente para analizar trayectorias de dinámica molecular y estructuras biomoleculares, como proteínas, ácidos nucleicos o membranas.[1][3]

## Qué es MDAnalysis

MDAnalysis proporciona una capa de alto nivel para cargar, manipular y analizar datos estructurales y de trayectorias generados por programas de simulación como GROMACS, NAMD, AMBER, CHARMM, LAMMPS o DL_POLY, además de archivos PDB y otros formatos estándar.  Está escrita en Python con partes críticas en C y se integra estrechamente con NumPy para tratar coordenadas atómicas como arreglos eficientes en memoria.[3][4][6]

## Capacidades principales

Permite hacer tareas típicas de análisis de dinámica molecular como cálculo de RMSD, RMSF, radios de giro, funciones de distribución radial, distancias, ángulos, dihedros, densidades y contactos intermoleculares, entre muchas otras.  Incluye un potente lenguaje de selección de átomos similar al de paquetes clásicos (por ejemplo, estilo CHARMM), lo que facilita trabajar con subconjuntos específicos de la molécula, como cadenas, residuos activos o tipos de átomos.[6][3]

## Uso en biofísica y simulación

En biofísica computacional se utiliza para estudiar la flexibilidad de proteínas, cambios conformacionales, estabilidad de complejos, interacción proteína-ligando y organización de membranas, todo a partir de trayectorias de dinámica molecular.  Su diseño orientado a objetos hace que sea sencillo escribir scripts personalizados de análisis, lo que permite a grupos de investigación implementar métricas y protocolos específicos sin reescribir manejadores de formatos o utilidades básicas.[5][3][6]

## Ecosistema y extensibilidad

La librería sirve como base para herramientas adicionales (“kits”) de análisis más especializados, que aprovechan su infraestructura para implementar métodos avanzados sin duplicar el manejo de datos.  Gracias a su naturaleza en Python y su integración con el ecosistema científico (NumPy, SciPy, matplotlib), se adapta bien a flujos de trabajo interactivos en cuadernos Jupyter y al procesamiento de sistemas con millones de partículas.[2][1][3]

[1](https://www.mdanalysis.org)
[2](https://www.mdanalysis.org/pages/mdakits/)
[3](https://pmc.ncbi.nlm.nih.gov/articles/PMC3144279/)
[4](https://bio.tools/mdanalysis)
[5](https://academic.oup.com/bioinformatics/article/33/17/2768/3859177)
[6](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.21787)
[7](https://www.youtube.com/watch?v=p3OUUnHXQjU)
