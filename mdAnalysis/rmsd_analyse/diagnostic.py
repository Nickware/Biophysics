import MDAnalysis as mda
from MDAnalysis.analysis import rms
import matplotlib.pyplot as plt

# Cargar simulación
u = mda.Universe('topology.gro', 'trajectory.xtc')

#  Diagnóstico: ¿qué hay en el sistema?
print("Residuos:", [res.resname for res in u.residues])
print("Nombres de átomos:", [atom.name for atom in u.atoms])
print("Número total de átomos:", len(u.atoms))

#  Seleccionar correctamente
# Como no hay 'protein', seleccionamos por resname o todos los átomos
atom_group = u.select_atoms('resname SOL')  # o 'all'

print("Átomos seleccionados:", len(atom_group))
if len(atom_group) == 0:
    print(" ERROR: No se seleccionaron átomos. Revisa los nombres.")
else:
    print("Nombres de átomos seleccionados:", [a.name for a in atom_group])
