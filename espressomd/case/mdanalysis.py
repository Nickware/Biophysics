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
