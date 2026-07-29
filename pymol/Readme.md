# PyMOL

PyMOL es un sistema de visualización molecular ampliamente utilizado en biología estructural y química computacional, capaz de generar imágenes 3D de alta calidad de moléculas pequeñas y macromoléculas biológicas (proteínas, ácidos nucleicos, complejos proteína-ligando).

Fue creado en el año 2000 por **Warren Lyford DeLano**, quien fundó DeLano Scientific LLC para su distribución comercial. Tras su fallecimiento en 2009, el desarrollo y mantenimiento del proyecto pasó a **Schrödinger, Inc.**, que lo sigue manteniendo hoy.

> **Nota sobre licenciamiento:** PyMOL no es simplemente "código abierto" sin matices. Existe una versión *open-source* (repositorio `schrodinger/pymol-open-source`, licencia estilo BSD) mantenida por la comunidad y por Schrödinger, pero la versión completa ("incentive build"), con todas las funciones —incluyendo `morph` para animaciones de interpolación entre conformaciones— requiere una licencia paga tras un periodo de prueba de 30 días. Vale la pena tenerlo claro antes de diseñar un flujo de trabajo que dependa de funciones exclusivas de la versión con licencia.

## Tabla de contenidos

- [Instalación](#instalación)
- [API Python nativa](#api-python-nativa)
- [Lenguaje de selección](#lenguaje-de-selección)
- [Integración con software de simulación](#integración-con-software-de-simulación)
- [Pipeline: simulación → análisis → visualización](#pipeline-simulación--análisis--visualización)
- [Plugins oficiales](#plugins-oficiales)
- [Automatización de películas](#automatización-de-películas)
- [Integración con Jupyter](#integración-con-jupyter)
- [Ventajas y limitaciones](#ventajas-y-limitaciones)
- [Referencias](#referencias)

## Instalación

**Versión open-source, vía conda (recomendada para uso académico sin licencia):**

`pymol-bundle` solo publica builds para Python 3.9–3.12. Si tu entorno `base` de conda está fijado a otra versión (por ejemplo Python 3.13, algo común en instalaciones recientes de Anaconda/Miniforge), el solver falla con un `LibMambaUnsatisfiableError` porque no puede satisfacer esa versión de Python junto con los requisitos de `pymol-bundle` — no es un problema del paquete, es un choque de versión de Python en el entorno donde intentas instalarlo.

La solución es **no instalarlo en `base`**, sino crear un entorno dedicado con una versión de Python compatible:

```bash
conda create -n pymol-env -c conda-forge -c schrodinger python=3.12 pymol-bundle
conda activate pymol-env
```

Si el solver por defecto (`classic`) es muy lento, usa `mamba` (normalmente ya viene con `libmamba` como solver por defecto en conda ≥23.10):

```bash
conda install -n base conda-libmamba-solver
conda create -n pymol-env -c conda-forge -c schrodinger python=3.12 pymol-bundle --solver=libmamba
```

**Versión con licencia (incentive build):** se descarga desde `pymol.org`, requiere una clave de licencia de Schrödinger tras el periodo de prueba.

**Verificar la instalación (con el entorno `pymol-env` activo):**

```bash
pymol -cq -d "print(cmd.get_version())"
```

### Alternativa: Distrobox con Debian (vía `apt`)

Si el canal de `schrodinger` en conda te da paquetes con dependencias de librerías compartidas mal declaradas (por ejemplo, un `ImportError` de una `.so` que el propio paquete debería haber traído), la ruta más simple es instalar PyMOL desde los repositorios de Debian dentro de un contenedor Distrobox — el mismo patrón que usamos para ESPResSo:

```bash
distrobox create --name pymol-dev --image debian:trixie
distrobox enter pymol-dev

sudo apt update
sudo apt install -y pymol
```

**Verificar:**

```bash
pymol -cq -d "print(cmd.get_version())"
```

Esta build viene de los repositorios de Debian, no del canal de Schrödinger, así que corresponde a la versión **open-source** de PyMOL (sin `morph` ni otras funciones de la versión con licencia), pero evita por completo los problemas de resolución de dependencias del canal de conda.

## API Python nativa

```python
from pymol import cmd

# Descargar un complejo real desde la RCSB PDB: 3PTB = tripsina bovina
# con el inhibidor benzamidina (código de residuo BEN) ya presente en el
# mismo archivo — no hay que inventar un .mol2 aparte para el ligando.
cmd.fetch("3ptb", "complejo", type="pdb")

# Separar receptor y ligando en objetos distintos a partir del complejo descargado
cmd.create("prot", "complejo and polymer.protein")
cmd.create("lig", "complejo and resn BEN")

# Selección y representación
cmd.select("active_site", "prot within 5 of lig")
cmd.show("cartoon", "prot")
cmd.show("sticks", "lig")
cmd.color("cyan", "prot")

# Renderizado y exportación
cmd.bg_color("white")
cmd.ray(1200, 900)
cmd.png("figura.png", dpi=300)
```

`polymer.protein` es un selector de PyMOL que toma solo la cadena polipeptídica (excluye agua, iones y el propio ligando); `resn BEN` aísla la benzamidina por su código de residuo de tres letras tal como aparece en el archivo de la PDB.

## Lenguaje de selección

PyMOL tiene su propio lenguaje de selección (distinto del de MDAnalysis, aunque conceptualmente similar):

```python
# Por cadena y rango de residuos
cmd.select("dominio", "chain A and resi 100-150")

# Por proximidad a un ligando (equivalente al "around" de MDAnalysis)
cmd.select("sitio_activo", "prot within 5 of resn LIG")

# Combinación booleana
cmd.select("polares_superficie", "resn SER+THR+ASN+GLN and not (name C+N+O+CA)")

# Selección por tipo de átomo/elemento
cmd.select("metales", "elem Zn+Fe+Mg")
```

Palabras clave frecuentes: `chain`, `resi`, `resn`, `name`, `elem`, `within X of`, `byres`, y operadores booleanos (`and`, `or`, `not`). A diferencia de las selecciones dinámicas de MDAnalysis (`updating=True`), en PyMOL una selección creada con `cmd.select` es estática: si el objeto cambia (por ejemplo, al cargar otro estado o trayectoria), hay que volver a ejecutar el `select` para actualizarla.

## Integración con software de simulación

| Software | Formato típico de entrada a PyMOL | Uso habitual |
|---|---|---|
| GROMACS | PDB/GRO tras `gmx trjconv` | Estructuras y frames de trayectorias MD |
| NAMD | DCD (vía VMD o carga directa de coordenadas) | Biomoléculas grandes |
| LAMMPS | Dump convertido a PDB/XYZ | Sistemas coloidales o de materia blanda |
| AMBER/cpptraj | PDB/NetCDF | Estructuras y trayectorias |
| pDynamo | PDB/XYZ | Estructuras de cálculos QM/MM |

PyMOL no lee directamente todos los formatos de trayectoria de cada motor de simulación; en la práctica, el paso intermedio más común es exportar a PDB (para una estructura o un frame) o usar una librería como MDAnalysis para leer la trayectoria nativa y escribir los frames de interés en un formato que PyMOL sí soporta bien.

## Pipeline: simulación → análisis → visualización

Un patrón real y verificable para combinar MDAnalysis (análisis numérico) con PyMOL (visualización), sin inventar módulos que no existen en la API pública de MDAnalysis:

```python
from pymol import cmd

# 0. Obtener una topología real como punto de partida (en vez de asumir
#    que ya existe un prot.pdb local): descargar 3PTB y quedarnos solo
#    con la proteína, sin agua ni el ligando de cristalización.
cmd.fetch("3ptb", "complejo", type="pdb")
cmd.save("prot.pdb", "complejo and polymer.protein")
```

```python
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

# 1. Cargar trayectoria y calcular RMSD respecto al primer frame
#    (prot.pdb = topología real obtenida en el paso 0; traj.dcd = tu propia
#    trayectoria de dinámica molecular, por ejemplo generada con GROMACS/NAMD)
u = mda.Universe("prot.pdb", "traj.dcd")
backbone = u.select_atoms("backbone")

R = rms.RMSD(u, select="backbone").run()
rmsd_series = R.results.rmsd[:, 2]  # columna con el RMSD por frame

# 2. Clustering simple sobre una matriz de distancias RMSD par-a-par
#    (aquí se asume que 'rmsd_matrix' ya fue calculada por separado,
#    p. ej. con MDAnalysis.analysis.encore o comparando frames dos a dos)
# distancias_condensadas = squareform(rmsd_matrix, checks=False)
# Z = linkage(distancias_condensadas, method="average")
# clusters = fcluster(Z, t=3, criterion="maxclust")

# 3. Guardar el frame representativo de un cluster de interés como PDB
frame_representativo = 0  # índice de ejemplo
u.trajectory[frame_representativo]
backbone.write("cluster_rep.pdb")
```

```python
# 4. Visualizar en PyMOL el frame representativo
from pymol import cmd

cmd.load("cluster_rep.pdb", "cluster_rep")
cmd.show("cartoon", "cluster_rep")
cmd.color("marine", "cluster_rep")
cmd.zoom("cluster_rep")
cmd.png("cluster_representativo.png", dpi=300)
```

> El clustering de trayectorias en MDAnalysis no tiene una clase pública `analysis.cluster.ClusterAnalysis`; para clustering conformacional real conviene usar `MDAnalysis.analysis.encore` (que sí es parte de la librería) o construir la matriz de distancias a mano y pasarla a `scipy.cluster.hierarchy`, como se esboza arriba.

## Plugins oficiales

| Plugin | Función |
|---|---|
| APBS Electrostatics | Cálculo y visualización de potencial electrostático |
| CEalign | Superposición estructural basada en alineamiento de contactos |
| Caver | Identificación de túneles y cavidades en proteínas |

Estos plugins se instalan y gestionan desde el menú `Plugin` de la interfaz de PyMOL o vía el gestor de plugins (`Plugin → Plugin Manager`).

## Automatización de películas

La interpolación entre conformaciones (`morph`) es una función **exclusiva de la versión con licencia** (incentive build), no está disponible en la versión open-source:

```python
from pymol import cmd

cmd.fetch("1akeA", "conf1", async_=0)
cmd.fetch("4akeA", "conf2", async_=0)
cmd.align("conf1", "conf2")

# morph(nombre_objeto, selección1, selección2)
cmd.morph("transicion", "conf1", "conf2")

cmd.mplay()          # reproducir en el visor
cmd.mpng("frame_")   # exportar frames como PNG numerados
```

Si se trabaja con la versión open-source, la alternativa es una interpolación lineal simple entre coordenadas (mucho más tosca), o herramientas de terceros como el módulo `psico` (`morpheasy_linear`).

## Integración con Jupyter

```python
from pymol import cmd
import ipywidgets as widgets

def mostrar_estructura(pdb_id):
    cmd.reinitialize()
    cmd.fetch(pdb_id, async_=0)
    cmd.show("cartoon")
    cmd.zoom()

widgets.interact(mostrar_estructura, pdb_id=["1abc", "2xyz"])
```

Para usar PyMOL de forma interactiva en Jupyter se requiere una build de PyMOL con soporte para el modo `-cq`/API embebida (por ejemplo `pymol2` o el paquete `pymol-open-source` instalado en el mismo entorno que el kernel).

## Ventajas y limitaciones

| Aspecto | Detalle |
|---|---|
| API Python completa | Scripts reproducibles para figuras y análisis |
| Renderizado con ray-tracing | Figuras de calidad de publicación |
| Selección expresiva | Sintaxis compacta por cadena, residuo, elemento, proximidad |
| Limitación de licencia | Funciones como `morph` requieren la versión con licencia |
| Formatos de trayectoria | Lectura nativa limitada; conviene apoyarse en MDAnalysis/VMD para formatos de motores de simulación |

## Referencias

- Página oficial: https://pymol.org
- Wiki de comandos (PyMOL Wiki): https://pymolwiki.org
- Repositorio open-source: https://github.com/schrodinger/pymol-open-source
- Documentación de MDAnalysis (integración de análisis): https://docs.mdanalysis.org
- Wikipedia, "PyMOL": https://en.wikipedia.org/wiki/PyMOL
- Wikipedia, "Warren Lyford DeLano": https://en.wikipedia.org/wiki/Warren_Lyford_DeLano
