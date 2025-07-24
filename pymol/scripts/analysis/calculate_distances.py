import MDAnalysis as mda
import numpy as np

# Cargar estructura
u = mda.Universe("1YKT.pdb")

# Seleccionar átomos (ej: residuos 10 y 20, átomos CA)
atom1 = u.select_atoms("resid 10 and name CA")[0]
atom2 = u.select_atoms("resid 20 and name CA")[0]

# Calcular distancia
distancia = np.linalg.norm(atom1.position - atom2.position)

# Guardar en .dat
with open("distancias.dat", "w") as f:
    f.write(f"# Distancia entre residuo 10 y 20\n")
    f.write(f"0.0 {distancia:.2f}\n")  # Tiempo 0 ps
