from pymol import cmd

# Cargar una estructura desde el PDB (ej: ubiquitina)
cmd.fetch("1UBQ", async_=0)

# Mostrar diferentes estilos
cmd.show("cartoon", "1UBQ")
cmd.color("cyan", "1UBQ")
cmd.show("sticks", "resn GLY+LYS")  # Muestra bastones en glicinas y lisinas

# Seleccionar dos átomos específicos para medir la distancia
cmd.select("atom1", "1UBQ and resi 1 and name CA")  # carbono alfa del residuo 1
cmd.select("atom2", "1UBQ and resi 50 and name CA") # carbono alfa del residuo 50

# Medir distancia entre esos dos átomos
cmd.distance("dist_CA", "atom1", "atom2")

# Ajustar la vista
cmd.zoom()

# Exportar imagen
cmd.png("proteina_1ubq.png", width=800, height=600, dpi=300, ray=1)
