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
