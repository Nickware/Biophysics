import MDAnalysis as mda
from MDAnalysis.analysis import rms
import matplotlib.pyplot as plt

# Cargar
u = mda.Universe('topology.gro', 'trajectory.xtc')

# Diagnóstico
print("Residuos:", [r.resname for r in u.residues])
atom_group = u.select_atoms('resname SOL')  # correcto
print("Átomos seleccionados:", len(atom_group))

# RMSD
R = rms.RMSD(atom_group, ref=atom_group, ref_frame=0)
R.run()

print("Datos RMSD:\n", R.results.rmsd)

# Gráfica
if len(R.results.rmsd) > 0:
    plt.plot(R.results.rmsd[:, 0], R.results.rmsd[:, 2], 's-')
    plt.xlabel("Frame")
    plt.ylabel("RMSD (Å)")
    plt.title("RMSD - Sistema Dummy")
    plt.grid(True)
    plt.savefig("rmsd_plot.png", dpi=150)
    print(" Gráfica guardada: rmsd_plot.png")
else:
    print("Sin datos de RMSD.")
