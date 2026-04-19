import MDAnalysis as mda
from MDAnalysis.analysis import rdf
from MDAnalysis.analysis.msd import EinsteinMSD  # ← CAMBIO AQUÍ
import matplotlib.pyplot as plt

# Cargar LAMMPS DATA + DUMP
u = mda.Universe("argon_liquid.data", "argon_liquid.dump",
                 format="LAMMPS", atom_style="id type x y z")

argon = u.select_atoms("type 1")  # Átomos de argón

# 1. RDF (distribución radial) - SIN CAMBIOS
rdf_analyzer = rdf.InterRDF(argon, argon, nbins=50, range=(0.0, 10.0))
rdf_analyzer.run()

# 2. MSD + Difusión (MÉTODO CORREGIDO)
msd_analyzer = EinsteinMSD(argon, msd_type='xyz', fft=True)
msd_analyzer.run(start=0, stop=-1000, step=50, n_frames=20)

# Coeficiente D = MSD/(6t) en 3D
times = msd_analyzer.results.times
msd = msd_analyzer.results.msd[:, 0]  # Componente x (o promedio)
D = np.polyfit(times[5:], msd[5:], 1)[0] / 6  # Ajuste lineal

print(f"Coeficiente de difusión D: {D:.4f} Å²/ps")

# Visualización
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# RDF
ax1.plot(rdf_analyzer.bins, rdf_analyzer.rdf, 'r-', lw=2)
ax1.axhline(1.0, color='k', ls='--', alpha=0.7, label='Gas ideal')
ax1.set_xlabel('Distancia (Å)')
ax1.set_ylabel('g(r)')
ax1.set_title('Función RDF Argón')
ax1.legend()
ax1.grid(True, alpha=0.3)

# MSD
ax2.plot(times, msd, 'bo-', label='MSD simulada')
ax2.plot(times, 6*D*times, 'r--', lw=2, label=f'Ajuste D={D:.3f}')
ax2.set_xlabel('Tiempo (ps)')
ax2.set_ylabel('MSD (Å²)')
ax2.set_title('Mean Square Displacement')
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('argon_rdf_msd.png', dpi=300)
plt.show()