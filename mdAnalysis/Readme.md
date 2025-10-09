# Análisis de RMSD con MDAnalysis: Tutorial Paso a Paso

Este proyecto demuestra cómo usar **MDAnalysis** en Python para:

1. Crear una topología molecular dummy (archivo `.gro`).
2. Generar una trayectoria simulada (archivo `.xtc`).
3. Calcular el **RMSD (Root Mean Square Deviation)**.
4. Graficar los resultados.

Es ideal para aprender los fundamentos de MDAnalysis sin necesidad de simulaciones reales.

## Estructura del Proyecto

proyecto-md/
├── topology.gro        # Topología estructural (3 átomos)
├── trajectory.xtc      # Trayectoria simulada (10 frames)
├── analyze_rmsd.py     # Script principal de análisis
├── generate_trajectory.py  # Genera la trayectoria dummy
└── rmsd_plot.png       # Salida del análisis

## Paso 1: Crear una topología dummy (`topology.gro`)

Crear un archivo `.gro` de GROMACS con una molécula ficticia (tipo agua simple):

```gro
Simple Dummy System
3
    1SOL    OW1    1   0.000   0.000   0.000
    1SOL    HW2    2   0.100   0.000   0.000
    1SOL    HW3    3   0.000   0.100   0.000
   1.00000   1.00000   1.00000
```

- **Residuo**: `SOL` (no `protein`)
- **Átomos**: 3 (OW1, HW2, HW3)
- **Unidad**: nanómetros
- **Caja**: 1×1×1 nm

> Importante: Como no hay residuos llamados `protein`, no puedes usar `select_atoms('protein')`.

---

## Paso 2: Generar una trayectoria simulada (`trajectory.xtc`)

Ejecutar un script para crear una trayectoria con 10 frames, donde las posiciones varían ligeramente (simulando dinámica).

### Script: `generate_trajectory.py`

```python
import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCWriter
import numpy as np

# Cargar topología
u = mda.Universe('topology.gro')
n_frames = 10

# Escribir archivo XTC
with XTCWriter('trajectory.xtc', n_atoms=len(u.atoms)) as W:
    for i in range(n_frames):
        # Añadir pequeño desplazamiento
        u.atoms.positions += 0.02 * np.sin(i * 0.5) * np.array([[1, 0, 0], [0, 1, 0], [-1, -1, 0]])
        
        # Actualizar timestep
        ts = u.trajectory.ts
        ts._pos = u.atoms.positions.astype(np.float32)
        ts.time = i * 100.0  # tiempo en ps
        ts.frame = i
        
        # Escribir frame
        W.write(u)  # escribe el Universe, no el timestep
```

Este script:
- Modifica las posiciones en cada paso.
- Asigna un tiempo real (0, 100, 200, ..., 900 ps).
- Guarda todo en formato `.xtc`.

---

## Paso 3: Calcular el RMSD

Usar el módulo `rms` de MDAnalysis para calcular cuánto cambia la estructura respecto al primer frame.

### Script: `analyze_rmsd.py`

```python
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import matplotlib.pyplot as plt

# Cargar simulación
u = mda.Universe('topology.gro', 'trajectory.xtc')

# Selección correcta de átomos
# Como el residuo es SOL, no se puede usar 'protein'
atom_group = u.select_atoms('resname SOL')  # o u.atoms

print(f"Átomos seleccionados: {len(atom_group)}")
print(f"Residuos: {[r.resname for r in u.residues]}")

# Calcular RMSD
R = rms.RMSD(atom_group, ref=atom_group, ref_frame=0).run()

# Mostrar resultados
print("Datos RMSD (Frame, Time, RMSD):")
print(R.results.rmsd)
```

---

## Paso 4: Graficar el RMSD

Finalmente, graficar el RMSD a lo largo del tiempo.

```python
import numpy as np

# Verificar que hay datos
if len(R.results.rmsd) > 0:
    plt.figure(figsize=(9, 5))
    plt.plot(R.results.rmsd[:, 0], R.results.rmsd[:, 2], 'o-', color='darkblue')
    plt.xlabel('Frame', fontsize=12)
    plt.ylabel('RMSD (Å)', fontsize=12)
    plt.title('RMSD respecto al frame inicial', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('rmsd_plot.png', dpi=150, bbox_inches='tight')
    print(" Gráfica guardada: rmsd_plot.png")
else:
    print(" No se generaron datos de RMSD.")
```

---

## Requisitos

- Python 3.7+
- MDAnalysis ≥ 2.0
- NumPy
- Matplotlib

### Instalación:

```bash
pip install mdanalysis numpy matplotlib
```

---

## Cómo ejecutar

```bash
# 1. Generar la trayectoria
python generate_trajectory.py

# 2. Analizar y graficar
python analyze_rmsd.py
```

---

## Notas importantes

- El archivo `.gro` solo contiene estructura, no topología química completa.
- El RMSD será pequeño porque el sistema es rígido y las fluctuaciones son suaves.
- Si se usa `select_atoms('protein')` en este sistema, obtendrás 0 átomos → RMSD vacío.
- Si no se observan gráficas, **guardar la figura con `plt.savefig()`** en lugar de confiar en `plt.show()`.

---

## Perspectivas

- Usar una proteína real (como ubiquitina) con `pdb` y `xtc`.
- Analizar RMSF, contacto entre residuos o estructura secundaria.
- Integrar con Jupyter Notebook para visualización interactiva.
