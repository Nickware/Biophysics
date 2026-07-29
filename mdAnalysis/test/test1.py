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
