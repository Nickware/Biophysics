## Análisis Visual de Fluctuaciones Térmicas en Proteínas con PyMOL

Este archivo `README.md` explica cómo utilizar un script en PyMOL para analizar visualmente las fluctuaciones térmicas o el desorden estructural en proteínas, empleando los valores de los factores B (B-factors).

### Requisitos

- **PyMOL** instalado (versión recomendada: 2.x o superior).
- Archivo de estructura de proteína en formato **PDB**.
- Script `etiquetar_b_factors.py` (ver fuente).

### Ejemplos de Archivos PDB de Proteínas Conocidas

A continuación, se presentan ejemplos de proteínas ampliamente estudiadas y sus identificadores PDB (PDB ID), ideales para analizar fluctuaciones térmicas o desorden estructural con PyMOL:


| Proteína | Organismo | PDB ID | Descripción breve |
| :-- | :-- | :-- | :-- |
| Hemoglobina | Humano | 1A3N | Forma desoxigenada de la hemoglobina humana. |
| Mioglobina | Ballena/Vertebrados | 1MBN | Primera proteína cuya estructura se resolvió en 3D. |
| Lisozima | Humano | 1LZ1 | Enzima antimicrobiana presente en secreciones humanas. |
| Lisozima | Pollo | 2CDS | Variante clásica de lisozima de Gallus gallus. |

### Detalles adicionales

- **Hemoglobina (1A3N):** Es la proteína responsable del transporte de oxígeno en la sangre. Su estructura ha sido clave para entender la cooperatividad y los cambios conformacionales asociados a la unión de oxígeno.
- **Mioglobina (1MBN):** Almacena oxígeno en el músculo y fue la primera proteína cuya estructura tridimensional se resolvió, marcando un hito en la biología estructural.
- **Lisozima (1LZ1, 2CDS):** Enzima que degrada la pared bacteriana, ampliamente utilizada como modelo en estudios de estabilidad y dinámica proteica.

Se pueden descargar estos archivos directamente desde el sitio oficial del Protein Data Bank (PDB) usando sus identificadores.

## Instrucciones de Uso

1. **Descarga un archivo PDB**
Ejemplos de proteínas conocidas:
    - Hemoglobina humana: `1A3N`
    - Mioglobina: `1MBN`
    - Lisozima humana: `1LZ1`
    - Lisozima de pollo: `2CDS`
2. **Carga el archivo PDB en PyMOL**
En la consola de PyMOL:

```
fetch 1A3N
```

3. **Carga el script en PyMOL**
En la consola de PyMOL:

```
run etiquetar_b_factors.py
```

4. **Ejecuta la función sobre el objeto cargado**
Sustituye `nombre_del_objeto` por el nombre del objeto (por ejemplo, `1A3N`):

```
etiquetar_b_factors nombre_del_objeto
```


## ¿Qué hace el script?

- **Elimina moléculas de agua** para limpiar la visualización.
- **Colorea la proteína** según los valores de B-factor, usando una escala de colores tipo arcoíris.
- **Etiqueta cada carbono alfa (CA)** con su valor de B-factor.
- **Ajusta el fondo** a blanco para mejorar la presentación.
- **Imprime en consola** los valores mínimo y máximo de B-factor encontrados.


## Ejemplo de Visualización

Al ejecutar el script, obtendrá una representación donde las regiones más flexibles o desordenadas (mayor B-factor) se distinguen fácilmente por el color y las etiquetas numéricas sobre cada residuo.

## Notas adicionales

- Es posible adaptar el script para otras selecciones o para exportar imágenes.
- El análisis visual de B-factors es útil para identificar regiones flexibles, desordenadas o con alta movilidad en proteínas.


## Referencias

- Ejemplo y adaptación del script: <a href="https://shaker.umh.es/computing/oe-11.html">shaker.umh.es/computing/oe-11.html</a>
- Más información sobre scripting en PyMOL: <a href="https://wiki.pymol.org">wiki.pymol.org</a>
