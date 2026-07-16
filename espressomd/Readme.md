# ESPResSo (espressomd)

ESPResSo (*Extensible Simulation Package for Research on Soft matter*) es un motor de dinámica molecular orientado a **materia blanda** y coloides, usado en física, química y biofísica computacional para estudiar polímeros, suspensiones coloidales, membranas y biopolímeros mediante representaciones coarse-grained.

`espressomd` es el módulo de **Python** (nota: es "espressomd", una sola palabra — no "expressomd") que envuelve el núcleo en C++ del motor, permitiendo definir sistemas, parámetros, fuerzas y protocolos de simulación mediante scripts de alto nivel, integrados con NumPy y SciPy.

## Tabla de contenidos

- [Instalación](#instalación)
- [Capacidades principales](#capacidades-principales)
- [Construyendo un sistema: partículas y polímeros](#construyendo-un-sistema-partículas-y-polímeros)
- [Integración con MDAnalysis (flujo correcto)](#integración-con-mdanalysis-flujo-correcto)
- [Exportar a VMD/PyMOL (formato VTF)](#exportar-a-vmdpymol-formato-vtf)
- [Integración con otro software de simulación](#integración-con-otro-software-de-simulación)
- [Pipeline: simulación → análisis → visualización](#pipeline-simulación--análisis--visualización)
- [Ventajas y limitaciones](#ventajas-y-limitaciones)
- [Referencias](#referencias)

## Instalación

A diferencia de MDAnalysis o PyMOL, **ESPResSo no se instala con un simple `pip install`**: es un motor en C++ que se compila con CMake, porque muchas características (electrostática, GPU, métodos hidrodinámicos) se activan o desactivan en tiempo de compilación según el archivo `myconfig.hpp`.

### Paso a paso en un contenedor Distrobox con Debian

Esta ruta aísla toda la toolchain de compilación (Boost, FFTW, HDF5, etc.) dentro de un contenedor Debian, sin tocar el sistema anfitrión. Requiere tener **Distrobox** y un motor de contenedores (`podman` o `docker`) ya instalados en el host — eso sí corre en el sistema anfitrión, no dentro del contenedor:

```bash
# En el HOST: instalar distrobox si no lo tienes (ejemplo con el script oficial)
curl -s https://raw.githubusercontent.com/89luca89/distrobox/main/install | sh -s -- --prefix ~/.local
```

**1. Crear y entrar al contenedor Debian (trixie, la estable actual):**

```bash
distrobox create --name espresso-dev --image debian:trixie
distrobox enter espresso-dev
```

A partir de aquí, todos los comandos se ejecutan **dentro** del contenedor `espresso-dev`.

**2. Actualizar el índice de paquetes e instalar las dependencias de compilación:**

```bash
sudo apt update
sudo apt install -y build-essential cmake cmake-curses-gui git \
    python3-dev python3-venv python3-pip openmpi-bin \
    libboost-all-dev libfftw3-dev libfftw3-mpi-dev \
    libhdf5-dev libhdf5-openmpi-dev libgsl-dev freeglut3-dev
```

**3. Crear el entorno virtual de Python y las dependencias del lado Python:**

```bash
python3 -m venv ~/espresso_env
source ~/espresso_env/bin/activate
cd ~/espresso_env/
. espresso_env/bin/activate
python3 -m pip install cython numpy scipy packaging setuptools h5py
```

**4. Clonar y compilar ESPResSo:**

```bash
git clone --depth=1 https://github.com/espressomd/espresso.git
cd espresso
mkdir build && cd build
cmake ..
make -j"$(nproc)"
```

**5. Verificar la instalación (todavía dentro del contenedor, dentro de `build/`):**

```bash
./pypresso -c "import espressomd; print(espressomd.features())"
```

Si ves una lista de features en lugar de un `ImportError`, la compilación fue exitosa. El script `pypresso` es el intérprete de Python que ESPResSo genera durante la compilación con las rutas del núcleo C++ ya configuradas; úsalo en vez de un `python3` genérico para correr tus scripts (`./pypresso mi_script.py`).

**6. Para volver a entrar al contenedor en sesiones futuras (desde el host, sin recompilar):**

```bash
distrobox enter espresso-dev
source ~/espresso_env/bin/activate
cd espresso/build
```

**Dependencias de compilación (ejemplo Ubuntu 24.04, fuera de un contenedor):**

```bash
sudo apt install build-essential cmake cmake-curses-gui python3-dev openmpi-bin \
    libboost-all-dev libfftw3-dev libfftw3-mpi-dev libhdf5-dev libhdf5-openmpi-dev \
    python3-pip libgsl-dev freeglut3-dev
```

**Dependencias de Python (idealmente en un entorno virtual):**

```bash
python3 -m venv espresso_env
source espresso_env/bin/activate
python3 -m pip install cython numpy scipy packaging setuptools h5py
```

**Compilación:**

```bash
git clone --depth=1 https://github.com/espressomd/espresso.git
cd espresso
mkdir build && cd build
cmake ..
make -j$(nproc)
```

**Alternativas más rápidas para no compilar localmente:**
- **Docker**: `ghcr.io/espressomd/docker/ubuntu:<tag>` trae una build ya compilada.
- **GitHub Codespaces**: el repositorio oficial soporta compilación automática en la nube para quien tenga cuenta de GitHub.

**Verificar la instalación** (desde el directorio `build`):

```bash
./pypresso -c "import espressomd; print(espressomd.features())"
```

## Capacidades principales

- Modelos coarse-grained para polímeros, membranas, surfactantes, coloides cargados y nanopartículas.
- Interacciones electrostáticas (Ewald, P3M), potenciales de Lennard-Jones, enlaces, ángulos, restricciones y campos externos (cizalla, flujos, campos eléctricos).
- Métodos hidrodinámicos como Lattice-Boltzmann o dinámica browniana para incluir efectos de solvente efectivo.

## Construyendo un sistema: partículas y polímeros

En ESPResSo **no existen clases de polímero que se añadan como "actor"** al sistema. Las partículas se agregan una a una (o mediante bucles/NumPy) con `system.part.add()`, y los enlaces se definen explícitamente:

```python
import espressomd
import numpy as np

system = espressomd.System(box_l=[32.0, 32.0, 100.0])
system.time_step = 0.01
system.cell_system.skin = 0.4

# Definir un tipo de enlace armónico para la cadena
from espressomd.interactions import HarmonicBond
bond = HarmonicBond(k=10.0, r_0=1.0)
system.bonded_inter.add(bond)

# Construir una cadena lineal de 20 beads
# Nota: desde ESPResSo 4.2.0 se eliminó el operador [] sobre system.part
# (para evitar ambigüedad con listas de partículas no contiguas).
# La forma soportada es guardar el handle que devuelve .add(),
# o usar system.part.by_id(i) si se necesita recuperarlo por id.
n_monomeros = 20
particulas = []
for i in range(n_monomeros):
    p = system.part.add(pos=[16.0, 16.0, i * 1.0], type=0)
    particulas.append(p)
    if i > 0:
        p.add_bond((bond, particulas[i - 1]))

# Interacción no enlazada tipo Lennard-Jones
system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=1.0, sigma=1.0, cutoff=2.5, shift="auto"
)

# Termostato de Langevin — 'seed' es obligatorio desde ESPResSo 4.x
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)

system.integrator.run(10000)
```

`system.actors` queda reservado para solvers físicos (P3M para electrostática, Lattice-Boltzmann, magnetostática), no para topología molecular.

## Integración con MDAnalysis (flujo correcto)

El README original planteaba escribir un archivo y luego leerlo con `mda.Universe(...)`. La documentación oficial de ESPResSo describe algo más directo: la clase **`espressomd.MDA_ESP`** expone los datos de partículas de ESPResSo como un *stream* que MDAnalysis puede leer **en memoria**, sin pasar por un archivo intermedio:

```python
import espressomd
from espressomd import MDA_ESP
import MDAnalysis as mda

system = espressomd.System(box_l=[32.0, 32.0, 100.0])
# ... construir el sistema y correr la simulación ...

eos = MDA_ESP.Stream(system)
u = mda.Universe(eos.topology, eos.trajectory)

# Análisis normal de MDAnalysis sobre el sistema de ESPResSo
backbone = u.select_atoms("type 0")
rg = backbone.radius_of_gyration()
print(f"Radio de giro: {rg:.3f}")
```

Este camino evita el problema de fondo del README original: no hay ningún "formato VMD" que MDAnalysis pueda leer directamente (VMD es el nombre del *programa* de visualización, no un formato de archivo).

## Exportar a VMD/PyMOL (formato VTF)

Si el objetivo es visualizar en VMD (o convertir para PyMOL), el formato correcto es **VTF** (VTF/VSF/VCF — *VTF Structure/Coordinate Format*), escrito a un descriptor de archivo abierto, no con un método `writevmd()` sobre `system.part` (que no existe):

```python
from espressomd.io.writer import vtf

fp = open("trayectoria.vtf", mode="w+t")
vtf.writevsf(system, fp)   # bloque de topología, una sola vez
for i in range(100):
    system.integrator.run(100)
    vtf.writevcf(system, fp)  # bloque de coordenadas, cada frame
fp.close()
```

ESPResSo tampoco tiene un escritor nativo de PDB colgado de `system.part`; para obtener PDB, el camino soportado es pasar por MDAnalysis (vía `MDA_ESP`) y usar sus métodos de escritura (`AtomGroup.write("archivo.pdb")`).

## Integración con otro software de simulación

| Software | Vía de integración | Uso típico |
|---|---|---|
| GROMACS | Exportar coordenadas/topología vía MDAnalysis | Validación cruzada, mapeo CG↔atomístico |
| LAMMPS | Formatos de coordenadas compartidos (XYZ) | Benchmarking, comparación de coloides |
| OpenMM | Sin puente directo; se comparan resultados a nivel de trayectoria/observables | Validación cross-engine |
| HOOMD-blue | Sin puente directo; ambos son motores GPU para materia blanda | Benchmarking de rendimiento |

A diferencia de lo que sugiere un README que promete "puentes nativos" hacia todo, en la práctica la mayoría de estas integraciones pasan por **MDAnalysis como capa común**, no por exportadores directos entre motores.

## Pipeline: simulación → análisis → visualización

```python
import espressomd
from espressomd import MDA_ESP
import MDAnalysis as mda
from MDAnalysis.analysis import rms

# 1. Simulación en ESPResSo (ver sección anterior para construir el sistema)
system = espressomd.System(box_l=[32.0, 32.0, 100.0])
# ... construir sistema, correr integrator ...

# 2. Puente a MDAnalysis en memoria
eos = MDA_ESP.Stream(system)
u = mda.Universe(eos.topology, eos.trajectory)

# 3. Análisis
R = rms.RMSD(u, select="type 0").run()

# 4. Exportar el frame final como PDB para visualizar en PyMOL/VMD
u.atoms.write("frame_final.pdb")
```

```python
# 5. Visualización en PyMOL
from pymol import cmd
cmd.load("frame_final.pdb", "sistema")
cmd.show("spheres", "sistema")
cmd.zoom()
cmd.png("figura.png", dpi=300)
```

## Ventajas y limitaciones

| Aspecto | Detalle |
|---|---|
| Interfaz Python | Scripts reproducibles, sin archivos de entrada rígidos |
| Rendimiento | Paralelización MPI y soporte GPU (CUDA) para ciertos métodos |
| Instalación | Requiere compilación con CMake; no hay `pip install espressomd` simple |
| Integración con análisis | El puente documentado y soportado es `MDA_ESP` hacia MDAnalysis, no exportadores directos a otros motores |
| Exportación para visualización | Formato nativo es VTF (VMD); PDB se obtiene indirectamente vía MDAnalysis |

## Referencias

- Documentación oficial: https://espressomd.github.io/doc/
- Instalación: https://espressomd.github.io/doc/installation.html
- Entrada/Salida (formatos VTF, puente `MDA_ESP`): https://espressomd.github.io/doc/io.html
- Repositorio: https://github.com/espressomd/espresso
- Documentación de MDAnalysis: https://docs.mdanalysis.org
