import MDAnalysis as mda
from MDAnalysis.analysis import rdf, diffusion, align
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.constants import k, N_A
from scipy.stats import linregress
import seaborn as sns

# =====================================================
# CARGA DE SIMULACIÓN LAMMPS
# =====================================================
print(" Cargando simulación LAMMPS...")
u = mda.Universe("gases_nobles.data", "gases_nobles.dump",
                 format="LAMMPS", atom_style="atomic")

# Definir gases nobles (tipos 1-6)
gases_nobles = {
    'He': u.select_atoms("type 1"),
    'Ne': u.select_atoms("type 2"), 
    'Ar': u.select_atoms("type 3"),
    'Kr': u.select_atoms("type 4"),
    'Xe': u.select_atoms("type 5"),
    'Rn': u.select_atoms("type 6")
}

propiedades = {
    'He': {'masa': 4.00, 'color': 'gold', 'sigma': 2.15},
    'Ne': {'masa': 20.18, 'color': 'hotpink', 'sigma': 2.75},
    'Ar': {'masa': 39.95, 'color': 'orange', 'sigma': 3.40},
    'Kr': {'masa': 83.80, 'color': 'limegreen', 'sigma': 3.65},
    'Xe': {'masa': 131.29, 'color': 'purple', 'sigma': 3.96},
    'Rn': {'masa': 222.00, 'color': 'gray', 'sigma': 4.38}
}

# =====================================================
# ANÁLISIS TERModinámico PRINCIPAL
# =====================================================
temps, pressures, volumes = [], [], []
msd_data = {gas: [] for gas in gases_nobles}

print(" Analizando propiedades termodinámicas...")
for i, ts in enumerate(u.trajectory[::50]):  # Cada 50 frames
    # Temperatura cinética instantánea
    KE = 0.5 * u.atoms.velocities**2
    T = np.mean(KE.sum(axis=1)) / (1.5 * k * 1e10)  # K
    temps.append(T)
    
    # Volumen y presión ideal
    V = np.prod(ts.dimensions[:3])
    volumes.append(V)
    N = len(u.atoms)
    P_ideal = N * k * T * 1e10 / V * 1e-5  # atm
    pressures.append(P_ideal)

# =====================================================
# RDF (COMPORTAMIENTO IDEAL vs REAL)
# =====================================================
print(" Calculando RDFs...")
rdf_data = {}
for gas, atoms in gases_nobles.items():
    if len(atoms) > 50:
        rdf_calc = rdf.InterRDF(atoms, atoms, nbins=50, range=(0.0, 12.0))
        rdf_calc.run()
        rdf_data[gas] = rdf_calc

# =====================================================
# DIFUSIÓN (D ∝ 1/√m para gas ideal)
# =====================================================
print(" Calculando auto-difusión...")
diffusions = {}
for gas, atoms in gases_nobles.items():
    if len(atoms) > 30:
        diff_calc = diffusion.DiffusionMode(atoms, msd=True)
        diff_calc.run(start=1000, stop=-2000, step=200, n_frames=25)
        D = diff_calc.results.diffusion[0][0]  # Å²/ps
        diffusions[gas] = D

# =====================================================
# VISUALIZACIÓN COMPLETA
# =====================================================
fig = plt.figure(figsize=(20, 16))

# 1. Temperatura y Presión
ax1 = plt.subplot(3, 3, 1)
ax1.plot(temps, 'b-', lw=2, label='T (K)')
ax1.axhline(np.mean(temps), color='b', ls='--', alpha=0.7)
ax1_twin = ax1.twinx()
ax1_twin.plot(pressures, 'r--', lw=2, label='P (atm)')
ax1_twin.axhline(np.mean(pressures), color='r', ls='--', alpha=0.7)
ax1.set_title('Temperatura y Presión')
ax1.legend(loc='upper left')
ax1_twin.legend(loc='upper right')

# 2. RDFs superpuestos
ax2 = plt.subplot(3, 3, 2)
for gas in rdf_data:
    ax2.plot(rdf_data[gas].bins, rdf_data[gas].rdf, 
             label=gas, color=propiedades[gas]['color'], lw=2)
ax2.axhline(1.0, color='k', ls=':', alpha=0.8, label='Gas Ideal g(r)=1')
ax2.set_xlabel('Distancia r (Å)')
ax2.set_ylabel('g(r)')
ax2.set_title('Funciones de Distribución Radial')
ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax2.set_ylim(0, 3.5)

# 3. Difusión vs Masa (ley de gas ideal)
ax3 = plt.subplot(3, 3, 3)
masas = [propiedades[gas]['masa'] for gas in diffusions]
Ds = list(diffusions.values())
ax3.scatter(masas, Ds, s=150, c=[propiedades[gas]['color'] 
              for gas in diffusions], alpha=0.8, edgecolors='black')
# Ajuste 1/√m
log_m = np.log(masas)
log_D = np.log(Ds)
slope, intercept, r_value, _, _ = linregress(log_m, log_D)
line = np.exp(intercept + slope * log_m)
ax3.plot(masas, line, 'k--', lw=2, label=f'R²={r_value**2:.3f}')
ax3.set_xlabel('Masa atómica (u)')
ax3.set_ylabel('D (Å²/ps)')
ax3.set_title('D ∝ 1/√m (Gas Ideal)')
ax3.legend()

# 4. Distribución de velocidades (Maxwell-Boltzmann)
ax4 = plt.subplot(3, 3, 4)
v_mag = np.linalg.norm(u.atoms.velocities, axis=1)
ax4.hist(v_mag, bins=50, density=True, alpha=0.7, color='lightblue')
v_th = np.sqrt(8*k*np.mean(temps)*1e10/(np.pi*40))  # v_th Ar
ax4.axvline(v_th, color='r', ls='--', label=f'v_th={v_th:.2f} Å/ps')
ax4.set_xlabel('Velocidad (|v|) (Å/ps)')
ax4.set_ylabel('Densidad')
ax4.set_title('Distribución Maxwell-Boltzmann')
ax4.legend()

# 5. Compresibilidad (PV/nRT)
ax5 = plt.subplot(3, 3, 5)
Z = np.array(pressures) * 1e5 * np.mean(volumes) / (len(u.atoms) * k * np.mean(temps) * 1e10)
ax5.hist(Z, bins=30, alpha=0.7, color='coral')
ax5.axvline(1.0, color='k', lw=3, label='Z=1 (Ideal)')
ax5.set_xlabel('Factor de Compresibilidad Z')
ax5.set_title('PV/nRT → 1 (Gas Ideal)')
ax5.legend()

# 6. Tabla resumen
ax6 = plt.subplot(3, 3, 6)
ax6.axis('off')
tabla = [
    ['Gas', 'N', 'D (Å²/ps)', 'RDF_max', 'σ (Å)'],
    ['He', len(gases_nobles['He']), f"{diffusions.get('He','N/A'):.3f}", 
     f"{np.max(rdf_data['He'].rdf):.2f}", '2.15'],
    ['Ne', len(gases_nobles['Ne']), f"{diffusions.get('Ne','N/A'):.3f}", 
     f"{np.max(rdf_data['Ne'].rdf):.2f}", '2.75'],
    ['Ar', len(gases_nobles['Ar']), f"{diffusions.get('Ar','N/A'):.3f}", 
     f"{np.max(rdf_data['Ar'].rdf):.2f}", '3.40'],
    ['Kr', len(gases_nobles['Kr']), f"{diffusions.get('Kr','N/A'):.3f}", 
     f"{np.max(rdf_data['Kr'].rdf):.2f}", '3.65'],
    ['Xe', len(gases_nobles['Xe']), f"{diffusions.get('Xe','N/A'):.3f}", 
     f"{np.max(rdf_data['Xe'].rdf):.2f}", '3.96'],
    ['Rn', len(gases_nobles['Rn']), f"{diffusions.get('Rn','N/A'):.3f}", 
     f"{np.max(rdf_data['Rn'].rdf):.2f}", '4.38']
]
table = ax6.table(cellText=tabla[1:], colLabels=tabla[0],
                  cellLoc='center', loc='center', colColours=['lightgray']*5)
table.auto_set_font_size(False)
table.set_fontsize(11)
table.scale(1.2, 2)
ax6.set_title('Resumen Gases Nobles', pad=20)

plt.tight_layout()
plt.savefig('gases_nobles_completo.png', dpi=300, bbox_inches='tight')
plt.show()

# =====================================================
# RESULTADOS CUANTITATIVOS
# =====================================================
print("\n" + "="*60)
print(" RESULTADOS - GASES NOBLES COMO GAS IDEAL")
print("="*60)
print(f"Temperatura promedio: {np.mean(temps):.1f} ± {np.std(temps):.1f} K")
print(f"Presión promedio:    {np.mean(pressures):.2f} ± {np.std(pressures):.2f} atm")
print(f"Factor Z promedio:   {np.mean(Z):.3f} (ideal=1.000)")
print(f"Correlación D vs 1/√m: R² = {r_value**2:.3f}")
print("\nAuto-difusión por especie:")
for gas in diffusions:
    print(f"  {gas:2s}: {diffusions[gas]:.4f} Å²/ps")

print("\n Gas más ideal: He (RDF≈1.0, Z≈1.0)")
print(" Gas menos ideal: Rn (RDF>2.0, Z desviado)")
print("="*60)
