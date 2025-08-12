## pDynamo

**pDynamo** es una biblioteca open-source diseñada para la simulación de sistemas moleculares, integrando métodos de química cuántica (QC), mecánica molecular (MM) y enfoques híbridos QC/MM que permiten estudiar procesos a nivel atomístico en bioquímica y química computacional.[^1][^2][^3]

### Características principales

- **Métodos QC:** Permite simulaciones con teoría del funcional de la densidad (*DFT*), *Hartree-Fock* y métodos semiempíricos como AM1, MNDO, PM3 y RM1, usando bases gaussianas.
- **Mecánica Molecular:** Soporta varios campos de fuerza estándar como AMBER, CHARMM y OPLS-AA.[^2]
- **Simulación híbrida QC/MM:** Puede combinar (e incluso acoplar a software externo) métodos cuánticos y clásicos, ideal para estudiar reacciones químicas enzimáticas, catálisis o mutaciones.
- **Simulaciones dinámicas:** Incorpora algoritmos avanzados para dinámica molecular y simulaciones Monte Carlo, permitiendo explorar el comportamiento temporal de biomoléculas, sus interacciones y rutas de reacción.[^4][^1]
- **Análisis estructurales:** Realiza optimizaciones geométricas, búsqueda de estados de transición, cálculos de propiedades (cargas, dipolos), análisis de modos normales, entre otros.
- **Manipulación y conversión de archivos:** Lee y escribe en múltiples formatos (XYZ, MOL, PDB, SMILES), facilitando su integración con otras herramientas y el manejo de datos estructurales.[^5]
- **Extensible y modular:** Está escrito principalmente en Python y Cython, lo que facilita la personalización mediante scripts Python; además, ofrece add-ons para cálculos de protonación y búsqueda de estados de transición.[^2]
- **Uso y aprendizaje:** Cuenta con tutoriales que abarcan desde simulación de dinámica molecular (por ejemplo, para el bALA) hasta la comparación de diferentes algoritmos de simulación y análisis de datos.[^6][^4]


### Contexto en biofísica y química computacional

pDynamo se utiliza en investigaciones de biofísica para modelar y simular mecanismos moleculares complejos, como la dinámica de proteínas, reacciones enzimáticas y reconocimiento molecular. Su flexibilidad para integrar métodos cuánticos y clásicos lo hace ideal para estudios donde es necesario un equilibrio entre precisión (zona reactiva tratada cuánticamente) y eficiencia computacional (resto del sistema tratado clásicamente).[^7][^1]

### Ejemplo de flujo de trabajo con pDynamo

- **Preparación:** Importación de la estructura molecular desde archivos XYZ, PDB, MOL o cadenas SMILES.
- **Definición de modelos:** Asignación de modelos de energía QC, MM o híbridos.
- **Simulación:** Ejecución de dinámica molecular, optimización geométrica o cálculos de trayectorias de reacción.
- **Análisis:** Extracción de propiedades estructurales y energéticas, visualización y post-procesamiento con software complementario.


### Integración con otras herramientas

pDynamo puede ser integrado en flujos de trabajo junto a herramientas como PyMOL (visualización), MDAnalysis (análisis de trayectorias), y otras plataformas computacionales, aportando capacidades avanzadas para simulaciones y análisis en biofísica computacional.[^1][^5]

***

En resumen, **pDynamo** es una potente y flexible plataforma para la simulación, análisis y modelado en química teórica y biofísica computacional, apta para investigaciones académicas y científicas avanzadas.

<div style="text-align: center">⁂</div>

[^1]: https://www.pdynamo.org

[^2]: https://github.com/pdynamo/pDynamo3

[^3]: https://pubmed.ncbi.nlm.nih.gov/26636368/

[^4]: https://www.pdynamo.org/tutorials/molecular-dynamics

[^5]: https://www.pdynamo.org/tutorials/molecular-systems

[^6]: https://sites.google.com/site/pdynamowiki/how-tos

[^7]: https://pubmed.ncbi.nlm.nih.gov/36449463/

[^8]: https://pubs.acs.org/doi/abs/10.1021/acs.jcim.2c01239

[^9]: https://pubs.acs.org/doi/10.1021/ct800092p

[^10]: https://mybiosoftware.com/tag/pdynamo

[^11]: https://dynamomd.com/index.php/features.html

[^12]: https://mybiosoftware.com/pdynamo-simulation-molecular-systems.html

[^13]: https://pubs.acs.org/doi/abs/10.1021/ct800092p

[^14]: https://rm1.sparkle.pro.br/rm1-software/pdynamo

[^15]: https://www.pdynamo.org/documentation/publications

