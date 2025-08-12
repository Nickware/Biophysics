# Análisis de Conformación de Glicina

## Descripción

Este script está diseñado para manipular y analizar las coordenadas atómicas de diferentes confórmeros de la glicina. Su objetivo principal es calcular la **desviación cuadrática media (RMSD)** entre dos estructuras moleculares. El script toma como referencia una estructura principal, la traslada a su centro de masa y la alinea con sus ejes principales de inercia. Luego, compara otras estructuras con esta referencia para determinar qué tan similares son antes y después de la superposición.

## Funcionalidad

El script realiza los siguientes pasos:

1.  **Carga de estructuras**: Lee los archivos de coordenadas `.xyz` de los confórmeros de la glicina (`gly-iin.xyz` y `gly-ip.xyz`).
2.  **Referencia de estructura**: Selecciona una de las estructuras como referencia, la traslada a su centro de masa y la orienta según sus **ejes principales de inercia**.
3.  **Cálculo de la matriz de inercia**: Imprime la matriz de inercia de la estructura de referencia antes y después de la reorientación para demostrar el efecto de la alineación.
4.  **Comparación de estructuras**: Itera sobre las estructuras restantes para calcular el RMSD con respecto a la estructura de referencia.
      * Calcula el **RMSD inicial** (antes de la superposición).
      * Superpone la estructura actual a la de referencia para encontrar la mejor alineación posible.
      * Calcula el **RMSD final** (después de la superposición), lo que indica la similitud estructural intrínseca.
5.  **Generación de tabla**: Muestra los resultados del RMSD en una tabla clara para comparar las estructuras.

## Dependencias

El script utiliza un módulo llamado `Definitions` y funciones como `XYZFile_ToSystem` y `XYZFile_ToCoordinates3`, que no están incluidas en el código proporcionado. Se debe tener instaladas estas funciones, que forman parte de una biblioteca local que maneja las operaciones de E/S de archivos y las manipulaciones de coordenadas.

  * `os`: Módulo estándar de Python para interactuar con el sistema operativo.
  * `Definitions`: Módulo local que contiene funciones para la definición de sistemas moleculares.
  * `XYZFile_ToSystem`: Función para leer un archivo `.xyz` y construir un objeto de sistema molecular.
  * `XYZFile_ToCoordinates3`: Función para leer un archivo `.xyz` y obtener un objeto de coordenadas 3D.

## Uso

Para ejecutar este script, se requiere tener los archivos `.xyz` de las estructuras de la glicina en la ruta especificada por `xyzPath` dentro del script.

1.  Asegúrarse de que las dependencias estén correctamente instaladas y configuradas.
2.  Ejecuta el script desde la línea de comandos:

<!-- end list -->

```bash
python nombre_del_script.py
```

El resultado será impreso en el `logFile`, mostrando la matriz de inercia y la tabla de desviaciones RMS entre las estructuras.

## Archivos de Ejemplo

El script espera encontrar los siguientes archivos en la ruta de las coordenadas:

  * `gly-iin.xyz`: Coordenadas del confórmero Gly-IIn de la glicina.
  * `gly-ip.xyz`: Coordenadas del confórmero Gly-Ip de la glicina.

## Requisitos

- Python 3.x
- pDynamo instalado y configurado correctamente.
- Archivos `.xyz` de entrada.
