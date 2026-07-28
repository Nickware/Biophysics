# MDAnalysis

MDAnalysis es una librería de Python de código abierto diseñada para analizar trayectorias de dinámica molecular (DM) y estructuras biomoleculares, como proteínas, ácidos nucleicos o membranas.

Proporciona una capa de alto nivel para cargar, manipular y analizar datos estructurales y de trayectorias generados por motores de simulación como GROMACS, NAMD, AMBER, CHARMM, LAMMPS o DL_POLY, además de archivos PDB y otros formatos estándar. Está escrita en Python con partes críticas en C, y se integra estrechamente con NumPy para tratar las coordenadas atómicas como arreglos eficientes en memoria.

## Tabla de contenidos

- [Dependencias](#dependencias)
- [Instalación](#instalación)
- [Verificación con pruebas (tests)](#verificación-con-pruebas-tests)
- [Lenguaje de selección de átomos](#lenguaje-de-selección-de-átomos)
- [Comparación con MDTraj y pytraj](#comparación-con-mdtraj-y-pytraj)
- [Capacidades principales](#capacidades-principales)
- [Ecosistema y extensibilidad](#ecosistema-y-extensibilidad)

## Dependencias

MDAnalysis tiene un núcleo mínimo de dependencias y un conjunto extendido, opcional, requerido por módulos específicos de análisis:

**Núcleo (instalación mínima):**
- Python ≥ 3.9 (revisar la matriz de versiones soportadas en cada release)
- NumPy
- SciPy
- GSD, MMTF-python y otros parsers de formato (instalados automáticamente según el paquete)

**Análisis completo (tag `[analysis]`):**
- matplotlib
- scikit-learn
- seaborn
- networkx
- tqdm
- joblib

**Para pruebas y ejemplos:**
- `MDAnalysisTests` (paquete separado, ~90 MB, con archivos de ejemplo)
- `pytest`

**Opcional (para OpenMP / paralelismo):**
- Un compilador con soporte OpenMP instalado en el sistema. Los canales de conda solo distribuyen builds seriales; para paralelismo hace falta instalar vía `pip` con un toolchain OpenMP funcional.

## Instalación

### Con pip (mínima)

```bash
pip install --upgrade MDAnalysis
```

### Con pip (con todas las dependencias de análisis)

```bash
pip install --upgrade "MDAnalysis[analysis]"
```

### Con conda / mamba

```bash
mamba install -c conda-forge mdanalysis
```

> Nota: las builds de conda/mamba solo soportan cálculo serial. Si se necesita paralelismo vía OpenMP, instalar con `pip` y contar con un compilador OpenMP disponible en el sistema.

### Paquete de pruebas y archivos de ejemplo

```bash
pip install --upgrade MDAnalysisTests
# o, con conda:
conda install -c conda-forge mdanalysistests
```

### Verificar la instalación

```bash
python -c "import MDAnalysis as mda; print(mda.__version__)"
```

## Verificación con pruebas (tests)

Un ejemplo mínimo de test con `pytest`, usando los archivos de ejemplo que trae `MDAnalysisTests`, para validar que la instalación carga correctamente una topología y una trayectoria:

```python
# test_mdanalysis_basico.py
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD

def test_carga_universo():
    """Verifica que un Universe se construye con el número esperado de átomos."""
    u = mda.Universe(PSF, DCD)
    assert len(u.atoms) > 0
    assert u.trajectory.n_frames > 0

def test_seleccion_basica():
    """Verifica que una selección de átomos de backbone no esté vacía."""
    u = mda.Universe(PSF, DCD)
    backbone = u.select_atoms("backbone")
    assert len(backbone) > 0
```

Ejecución:

```bash
pytest test_mdanalysis_basico.py -v
```

Esto sirve como humo (smoke test) tras cualquier instalación o actualización del entorno, antes de correr análisis más complejos.

## Lenguaje de selección de átomos

Uno de los puntos fuertes de MDAnalysis es su lenguaje de selección, inspirado en el estilo CHARMM/VMD, que permite aislar subconjuntos moleculares con sintaxis compacta sobre un objeto `Universe`:

```python
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD

u = mda.Universe(PSF, DCD)

# Selección por tipo de residuo
proteina = u.select_atoms("protein")

# Selección por nombre de átomo y rango de residuos
ca_activos = u.select_atoms("name CA and resid 10-50")

# Selección geométrica: átomos dentro de un radio de un ligando
sitio_activo = u.select_atoms("protein and around 5.0 resname LIG")

# Combinación booleana de criterios
hidrofobicos_superficie = u.select_atoms(
    "resname ALA VAL LEU ILE PHE and not backbone"
)

# Selecciones dinámicas (se recalculan en cada frame de la trayectoria)
contactos_dinamicos = u.select_atoms("around 4.0 resname LIG", updating=True)
for ts in u.trajectory:
    print(ts.frame, len(contactos_dinamicos))
```

Las palabras clave más comunes incluyen `protein`, `backbone`, `name`, `resname`, `resid`, `segid`, `around`, `same residue as`, y operadores booleanos (`and`, `or`, `not`). El parámetro `updating=True` es particularmente útil para análisis de contactos o sitios de unión que cambian a lo largo de la trayectoria.

## Comparación con MDTraj y pytraj

Las tres librerías cubren necesidades similares en el análisis de DM, pero con énfasis distintos:

| Aspecto | MDAnalysis | MDTraj | pytraj |
|---|---|---|---|
| Enfoque principal | Análisis general, orientado a objetos (`Universe`, `AtomGroup`) | Manejo eficiente de trayectorias como arreglos NumPy | Interfaz Python sobre el motor de análisis de AmberTools (`cpptraj`) |
| Lenguaje de selección | Selección estilo CHARMM/VMD, muy expresiva y con selecciones dinámicas | Selección más limitada, basada en expresiones tipo SQL/DSL propio | Máscaras estilo Amber (similares a sintaxis de `cpptraj`) |
| Formatos soportados | Amplia variedad (GROMACS, NAMD, AMBER, CHARMM, LAMMPS, DL_POLY, PDB, etc.) | Amplia variedad, fuerte en formatos comprimidos propios (HDF5, XTC) | Fuerte integración nativa con formatos de Amber; también soporta otros vía cpptraj |
| Rendimiento | Bueno, con partes críticas en C; paralelismo OpenMP opcional vía pip | Muy eficiente en I/O y en operaciones vectorizadas sobre NumPy | Muy eficiente para flujos de trabajo puramente Amber, aprovecha cpptraj en C++ |
| Ecosistema | MDAKits (extensiones especializadas), fuerte comunidad multi-motor | Buena integración con OpenMM y el stack científico de Python | Pensado para usuarios ya inmersos en el ecosistema AmberTools |
| Curva de aprendizaje | Moderada; API orientada a objetos intuitiva para quien viene de scripting | Moderada; requiere pensar en términos de arreglos NumPy | Más específica; requiere familiaridad con convenciones de Amber/cpptraj |

En términos de resultados numéricos (RMSD, RMSF, radios de giro, distancias), las tres librerías deberían converger a valores equivalentes para el mismo sistema y la misma definición de selección, ya que implementan los mismos algoritmos estándar de la literatura de DM. Las diferencias prácticas suelen aparecer en:

- **Precisión numérica y manejo de PBC** (condiciones periódicas de contorno): cada librería tiene su propia convención para "unwrapping" de moléculas, lo que puede generar pequeñas discrepancias si no se configura igual en las tres.
- **Convenciones de selección**: una misma selección puede no ser trivialmente equivalente entre sintaxis (por ejemplo, `backbone` en MDAnalysis vs. una máscara Amber en pytraj), por lo que conviene validar explícitamente que los átomos seleccionados coincidan antes de comparar métricas.
- **Rendimiento en trayectorias muy grandes**: MDTraj suele ser más rápido en operaciones puramente vectorizadas sobre NumPy; pytraj es competitivo cuando el flujo de trabajo es nativamente Amber; MDAnalysis ofrece el balance más flexible cuando se trabaja con múltiples formatos de origen dentro del mismo proyecto.

Para una comparación cuantitativa propia, un patrón recomendado es correr el mismo análisis (por ejemplo RMSD de backbone) sobre el mismo par topología/trayectoria en las tres librerías, y comparar los arreglos resultantes con `numpy.allclose` para detectar discrepancias por encima de una tolerancia razonable.

## Capacidades principales

Permite hacer tareas típicas de análisis de dinámica molecular como cálculo de RMSD, RMSF, radios de giro, funciones de distribución radial, distancias, ángulos, dihedros, densidades y contactos intermoleculares, entre muchas otras.

## Uso en biofísica y simulación

En biofísica computacional se utiliza para estudiar la flexibilidad de proteínas, cambios conformacionales, estabilidad de complejos, interacción proteína-ligando y organización de membranas, todo a partir de trayectorias de dinámica molecular. Su diseño orientado a objetos hace que sea sencillo escribir scripts personalizados de análisis, lo que permite a grupos de investigación implementar métricas y protocolos específicos sin reescribir manejadores de formatos o utilidades básicas.

## Ecosistema y extensibilidad

La librería sirve como base para herramientas adicionales ("MDAKits") de análisis más especializado, que aprovechan su infraestructura para implementar métodos avanzados sin duplicar el manejo de datos. Gracias a su naturaleza en Python y su integración con el ecosistema científico (NumPy, SciPy, matplotlib), se adapta bien a flujos de trabajo interactivos en cuadernos Jupyter y al procesamiento de sistemas con millones de partículas.

## Referencias

- [1] https://www.mdanalysis.org
- [2] https://www.mdanalysis.org/pages/mdakits/
- [3] https://pmc.ncbi.nlm.nih.gov/articles/PMC3144279/
- [4] https://bio.tools/mdanalysis
- [5] https://academic.oup.com/bioinformatics/article/33/17/2768/3859177
- [6] https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.21787
- [7] https://www.youtube.com/watch?v=p3OUUnHXQjU
- Guía de instalación oficial: https://userguide.mdanalysis.org/stable/installation.html
