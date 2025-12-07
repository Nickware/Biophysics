# generate_dummy_trajectory.py
import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCWriter
import numpy as np

# Cargar universo
u = mda.Universe('topology.gro')

n_frames = 10
n_atoms = len(u.atoms)

# Crear el writer
with XTCWriter('trajectory.xtc', n_atoms) as W:
    for i in range(n_frames):
        # Modificar posiciones (ejemplo: pequeño movimiento)
        displacement = 0.02 * np.sin(i * 0.5) * np.array([[1, 0, 0], [0, 1, 0], [-1, -1, 0]])
        u.atoms.positions = u.atoms.positions + displacement

        # Acceder al timestep actual
        ts = u.trajectory.ts

        # Asegurarse de que las posiciones estén en float32
        ts._pos = u.atoms.positions.astype(np.float32)

        # Modificar atributos del timestep
        ts.frame = i
        ts.time = i * 100.0  # asignar tiempo en ps

        # Escribir el universo (o AtomGroup), NO el timestep
        W.write(u)  # <-- Esto SÍ funciona
