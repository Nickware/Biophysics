# Pymol

PyMOL es un software de visualización molecular ampliamente utilizado en biología estructural y bioquímica. Fue creado por Warren Lyford Delano y es desarrollado por DeLano Scientific LLC. Este programa permite generar imágenes 3D de alta calidad de moléculas pequeñas y macromoléculas biológicas, como proteínas y ácidos nucleicos[^1].

## Características principales

- **Código abierto**: PyMOL es extensible gracias a su integración con el lenguaje de programación Python, lo que permite realizar análisis avanzados utilizando bibliotecas como NumPy o pylab[^1].
- **Visualización profesional**: Es capaz de representar sitios activos, interacciones moleculares y estructuras cristalinas, facilitando el estudio de interacciones biomoleculares a nivel atómico[^2].
- **Creación de películas**: PyMOL permite generar animaciones a partir de escenas moleculares, útiles para presentaciones científicas[^3].
- **Exportación de figuras**: Ofrece herramientas para crear imágenes con fondo transparente y ajustar representaciones moleculares para trabajos científicos[^2][^5].


## Uso y funciones destacadas

- **Manipulación molecular**: Permite limpiar estructuras cristalinas, resaltar sitios de unión, representar contactos polares y exportar imágenes[^2].
- **Selección avanzada**: Los usuarios pueden seleccionar aminoácidos específicos, cadenas o regiones de una proteína mediante comandos personalizados[^4].
- **Análisis estructural**: Es útil para identificar propiedades físicas y químicas en estructuras moleculares[^8].


PyMOL es una herramienta esencial para investigadores en bioinformática, biología estructural y química computacional debido a su versatilidad y capacidad para producir visualizaciones detalladas y científicamente precisas.

## Integración de PyMOL con Otros Paquetes y Softwares

**PyMOL** es el **"visualizador universal"** de química y biofísica, con **APIs ricas** (Python, plugins, scripts) que lo convierten en el hub central de workflows computacionales.

## **1. API Python Nativa (Pymol Python Module)**

```python
import pymol  # pymol.cmd
from pymol import cmd

# Cargar estructuras
cmd.load("protein.pdb", "prot")
cmd.load("ligand.mol2", "lig")

# Análisis automático
cmd.select("active_site", "prot and resi 100-150")
cmd.show("cartoon", "prot")
cmd.show("sticks", "lig")
cmd.ray(1)  # Renderizado
```

## **2. Integración con Paquetes de Simulación**

| Software | Método | Uso típico |
|----------|--------|------------|
| **GROMACS** | `gmx trjconv -pbc mol -o prot.pdb` → PyMOL | Trayectorias MD |
| **LAMMPS** | Dump → VMD → PyMOL | Coloides/membranas |
| **NAMD** | DCD → VMD → PyMOL | Biomoléculas grandes |
| **ESPResSo** | VMD export → PyMOL | Materia blanda CG |
| **pDynamo** | XYZ/PDB directo | QC/MM optimizaciones |

## **3. MDAnalysis ↔ PyMOL (Análisis + Visualización)**

```python
import MDAnalysis as mda
import pymol
from pymol import cmd

# MDAnalysis: RMSD + clustering
u = mda.Universe("traj.dcd", "prot.pdb")
rmsd = mda.analysis.rms.RMSD(u).run()
clusters = mda.analysis.cluster.ClusterAnalysis(u).run()

# PyMOL: Visualizar cluster representativo
cmd.load("cluster_0.pdb", "cluster_rep")
cmd.show("cartoon", "cluster_rep")
cmd.zoom()
cmd.png("cluster_representative.png")
```

## **4. Plugins Oficiales PyMOL**

| Plugin | Función | Integración |
|--------|---------|-------------|
| **APBS** | Cargas electrostáticas | `Plugin → Electrostatics` |
| **CEalign** | Superposición estructural | `Wizard → CEalign` |
| **ColorBlindness** | Paletas accesibles | CBM + VMD colors |
| **Caver** | Túneles en proteínas | Análisis de cavidades |
| **FABP** | Clustering conformacional | Análisis de trayectorias |

## **5. Scripts PyMOL + Análisis Automatizado**

```python
# Script maestro: MD → Clustering → PyMOL
import subprocess, pymol, MDAnalysis as mda

# 1. Clustering con MDAnalysis
u = mda.Universe("traj.xtc", "prot.gro")
clusters = mda.analysis.cluster.ClusterAnalysis(u).run()

# 2. PyMOL para cada cluster
for i, cluster in enumerate(clusters.results.clustercenters):
    cmd.load(f"cluster_{i}.pdb")
    cmd.show("cartoon")
    cmd.color("spectrum", "cluster_%d" % i)
    cmd.group("cluster_%d" % i)

cmd.zoom()
cmd.png("all_clusters.png")
```

## **6. Pipeline Biofísica Completa**

```
Simulación (LAMMPS/GROMACS) 
    ↓ PDB/DCD
MDAnalysis (RMSD/RMSF/clustering)
    ↓ Cluster reps PDB
PyMOL (Visualización + películas)
    ↓ PNG/MP4 para paper
```

## **7. Integración con Jupyter Notebooks**

```python
# Jupyter + PyMOL (interactivo)
from pymol import cmd
import ipywidgets as widgets

def show_structure(pdb_file):
    cmd.delete("all")
    cmd.load(pdb_file)
    cmd.show("cartoon")
    cmd.zoom()

widgets.interact(show_structure, 
                pdb_file=["1abc.pdb", "2xyz.pdb"])
```

## **8. Plugins y Extensiones Avanzados**

### **a) PyMOL + Machine Learning**
```python
# AlphaFold + PyMOL
cmd.fetch("AF-P12345-F1-model_v4", "alphafold")
cmd.show("cartoon", "alphafold")
cmd.color("blue", "alphafold")
```

### **b) PyMOL + Quantum Chemistry**
```python
# ORCA/Gaussian → PyMOL (densidades)
cmd.load("density.cube", "electron_density")
cmd.isosurface("density_surf", "electron_density", 0.02)
cmd.color("yellow", "density_surf")
```

## **9. Workflows Especializados**

| Aplicación | Pipeline |
|------------|----------|
| **Docking** | AutoDock → PyMOL (poses) + APBS (energías) |
| **MD Analysis** | GROMACS → MDAnalysis → PyMOL clusters |
| **QC/MM** | pDynamo → PyMOL (reactivos/productos) |
| **Cryo-EM** | Relion → PyMOL (modelos atómicos) |

## **10. Automatización de Películas**

```python
# Morphing entre estados
cmd.create("state1", "resi 1-100")
cmd.create("state2", "resi 101-200")
cmd.morph("state1", "state2", "morph_traj")
cmd.mplay("morph_traj")
cmd.mpng("morph_", 50)  # 50 frames PNG
```

## **Ventajas Únicas PyMOL**

| Característica | Beneficio |
|----------------|-----------|
| **API Python completa** | Scripts reproducibles |
| **Plugins ecosystem** | 100+ extensiones oficiales |
| **Formato universal** | PDB/DCD/XYZ/VMD/Trajectory |
| **Raytracing nativo** | Figuras publication-ready |
| **Multi-state** | Conformers + morphing |

## **Workflow "Paper-Ready" (5 líneas)**

```python
import pymol, MDAnalysis as mda
u = mda.Universe("traj.dcd", "prot.pdb")
rmsd = mda.analysis.rms.RMSD(u).run()  # Análisis
cmd.load("prot.pdb"); cmd.cartoon()    # Visualizar
cmd.png("paper_figure.png", dpi=300)   # ¡Listo!
```

## **Ecosistema Completo Recomendado**

```
PyMOL (Visualización) ←→ 
MDAnalysis (Análisis) ←→ 
GROMACS/LAMMPS (Simulación) ←→ 
pDynamo (QC/MM) ←→ 
Jupyter (Reporte)
```

**PyMOL es el "Photoshop de la biología estructural"** - conecta todos los software y genera figuras profesionales con un solo comando. ¡Indispensable para publicaciones! 

[^1]: https://es.wikipedia.org/wiki/PyMOL

[^2]: https://www.youtube.com/watch?v=3QFz5g2iwqo

[^3]: https://www.youtube.com/watch?v=OoQy7a_o4Gc

[^4]: https://www.youtube.com/watch?v=EedCvWBVCJ4

[^5]: https://www.youtube.com/watch?v=hVouZg5OKvc

[^6]: https://www.youtube.com/watch?v=EjmNiDwt5Us

[^7]: https://www.youtube.com/watch?v=QgCFjFgMzlM

[^8]: https://www.youtube.com/watch?v=NaSP9s-OIj4

[^9]: https://legadoweb.minciencias.gov.co/nexoglobal/blog/primer-video-utilizando-pymol

[^10]: https://www.instagram.com/recompsci/reel/DFTP0QGolWY/-descubre-cómo-hacer-películas-moleculares-con-pymol-3-en-este-reel-te-muestro-p/

