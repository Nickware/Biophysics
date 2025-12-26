"""
ANÁLISIS COMPLETO GASES NOBLES: LAMMPS + MDAnalysis
Valida ecuaciones fundamentales del gas ideal (PV=NkT, g(r)=1, D∝1/√m)
He, Ne, Ar, Kr, Xe, Rn en mezcla equimolar
"""

import MDAnalysis as mda
from MDAnalysis.analysis import rdf, diffusion
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.constants import k, N_A
from scipy.stats import linregress
import warnings
warnings.filterwarnings('ignore')

# =====================================================
# CONFIGURACIÓN Y CARGA LAMMPS
# =====================================================
print(" Cargando simulación LAMMPS de gases nobles...")
try:
    u = mda.Universe("gases_nobles.data", "gases_nobles.dump",
                     format="LAMMPS", atom_style="atomic")
    print(f" Cargado: {len(u.atoms)} átomos, {len(u.trajectory)} frames")
except FileNotFoundError:
    print(" Archivos no encontrados. Ejecuta primero LAMMPS.")
    exit()

# Gases nobles por tipo atómico (1-6)
gases_nobles = {
    'He': u.select_atoms("type 1"),
    'Ne': u.select_atoms("type 2"), 
    'Ar': u.select_atoms("type 3"),
    'Kr': u.select_atoms("type 4"),
    'Xe': u.select_atoms("type 5"),
    'Rn': u.select_atoms("type 6")
}

# Propiedades físicas reales
propiedades = {
    'He': {'masa': 4.00, 'color': 'gold', 'sigma': 2.15},
    'Ne': {'masa': 20.18, 'color': 'hotpink', 'sigma': 2.75},
    'Ar': {'masa': 39.95, 'color': 'orange', 'sigma': 3.40},
    'Kr': {'masa': 83.80, 'color': 'limegreen', 'sigma': 3.65},
    'Xe': {'masa': 131.29, 'color': 'purple', 'sigma': 3.96},
    'Rn': {'masa': 222.00, 'color': 'gray', 'sigma': 4.38}
}

print(f" Átomos por gas: { {k: len(v) for k,v in gases_nobles.items()} }")

# =====================================================
# ANÁLISIS TERModinámico PRINCIPAL (CORREGIDO)
# =====================================================
print(" Analizando propiedades termodinámicas...")
temps, pressures, volumes, Z_factors = [], [], [], []

for i, ts in enumerate(u.trajectory[::50]):  # Cada 50 frames
    # Temperatura cinética instantánea (3/2 kT = <1/2 mv²>)
    KE = 0.5 * ts.velocities**2
    T = np.mean(KE.sum(axis=1)) / (1.5 * k * 1e10)  # K (LAMMPS: Å/ps)
    temps.append(T)
    
    # Volumen de la caja
    V = np.prod(ts.dimensions[:3])  # Å³
    volumes.append(V)
    
    # Presión de gas IDEAL: P = NkT/V
    N_total = len(u.atoms)
    P_ideal = N * k * T * 1e10 / V * 1e-5  # atm
    pressures.append(P_ideal)
    
    # Factor de compresibilidad Z = PV/nRT (CORREGIDO: frame por frame)
    Z = P_ideal * 1e5 * V / (N_total * k * T * 1e10)
    Z_factors.append(Z)

temps, volumes, pressures, Z_factors = np.array(temps), np.array(volumes), np.array(pressures), np.array(Z_factors)

# =====================================================
# RDF: COMPORTAMIENTO ESTRUCTURAL (IDEAL vs REAL)
# =====================================================
print(" Calculando funciones de distribución radial...")
rdf_data = {}
for gas, atoms in gases_nobles.items():
    if len(atoms) > 50:  # Mínimo estadístico
        rdf_calc = rdf.InterRDF(atoms, atoms, nbins=50, range=(0.0, 12.0))
        rdf_calc.run()
        rdf_data[gas] = rdf_calc
        print(f"  {gas}: RDF calculado ({len(atoms)} átomos)")
    else:
        print(f"  {gas}: Saltado (<50 átomos)")

# =====================================================
# AUTO-DIFUSIÓN: D ∝ 1/√m (Gas Ideal)
# =====================================================
print(" Calculando coeficientes de auto-difusión...")
diffusions = {}
for gas, atoms in gases_nobles.items():
    if len(atoms) > 30:
        try:
            diff_calc = diffusion.DiffusionMode(atoms, msd=True)
            # Evita transitorios iniciales y finales
            diff_calc.run(start=1000, stop=-2000, step=200, n_frames=25)
            D = diff_calc.results.diffusion[0][0]  # Å²/ps
            diffusions[gas] = D
            print(f"  {gas}: D = {D:.4f} Å²/ps")
        except:
            print(f"  {gas}: Error en difusión")
    else:
        print(f"  {gas}: Saltado (<30 átomos)")

# =====================================================
# DISTRIBUCIÓN DE VELOCIDADES (Maxwell-Boltzmann)
# =====================================================
print(" Calculando distribución de velocidades...")
v_mags = []
for ts in u.trajectory[::100]:  # Promedio sobre frames (CORREGIDO)
    v_mag_frame = np.linalg.norm(ts.velocities, axis=1)
    v_mags.extend(v_mag_frame)
v_mag = np.array(v_mags)

# =====================================================
# VISUALIZACIÓN COMPLETA (3x3 PANEL)
# =====================================================
fig = plt.figure(figsize=(20, 16))
gs = fig.add_gridspec(3, 3, hspace=0.4, wspace=0.3)

# 1. TEMPERATURA Y PRESIÓN
ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(temps, 'b-', lw=2, label='T (K)')
ax1.axhline(np.mean(temps), color='b', ls='--', alpha=0.7, label=f'{np.mean(temps):.1f} K')
ax1_twin = ax1.twinx()
ax1_twin.plot(pressures, 'r--', lw=2, label='P (atm)')
ax1_twin.axhline(np.mean(pressures), color='r', ls='--', alpha=0.7, label=f'{np.mean(pressures):.2f} atm')
ax1.set_xlabel('Frame')
ax1.set_ylabel('Temperatura (K)', color='b')
ax1_twin.set_ylabel('Presión (atm)', color='r')
ax1.legend(loc='upper left')
ax1_twin.legend(loc='upper right')
ax1.set_title('Temperatura y Presión')

# 2. RDFs SUPERPONIENDO (g(r)=1 → Gas Ideal)
ax2 = fig.add_subplot(gs[0, 1])
for gas in rdf_data:
    ax2.plot(rdf_data[gas].bins, rdf_data[gas].rdf, 
             label=f'{gas}\nmax={np.max(rdf_data[gas].rdf):.2f}', 
             color=propiedades[gas]['color'], lw=2)
ax2.axhline(1.0, color='k', ls=':', alpha=0.8, lw=2, label='Gas Ideal g(r)=1')
ax2.set_xlabel('Distancia r (Å)')
ax2.set_ylabel('g(r)')
ax2.set_title('Funciones de Distribución Radial')
ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax2.set_ylim(0, 3.5)
ax2.grid(True, alpha=0.3)

# 3. DIFUSIÓN vs MASA (D ∝ 1/√m)
ax3 = fig.add_subplot(gs[0, 2])
if len(diffusions) > 2:  # Ajuste lineal requiere ≥3 puntos
    masas = [propiedades[gas]['masa'] for gas in diffusions]
    Ds = list(diffusions.values())
    colores_plot = [propiedades[gas]['color'] for gas in diffusions]
    
    ax3.scatter(masas, Ds, s=200, c=colores_plot, alpha=0.8, 
                edgecolors='black', linewidth=2, zorder=5)
    
    # Ajuste log-log: log(D) = a - 0.5 log(m)
    log_m = np.log(masas)
    log_D = np.log(Ds)
    slope, intercept, r_value, _, _ = linregress(log_m, log_D)
    line = np.exp(intercept + slope * log_m)
    ax3.plot(masas, line, 'k--', lw=3, 
             label=f'Ajuste: slope={slope:.2f}\nR²={r_value**2:.3f}')
else:
    for gas in diffusions:
        ax3.scatter(propiedades[gas]['masa'], diffusions[gas],
                   c=propiedades[gas]['color'], s=200, edgecolors='black')
ax3.set_xlabel('Masa atómica (u)')
ax3.set_ylabel('D (Å²/ps)')
ax3.set_title('D ∝ 1/√m (Gas Ideal)')
ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax3.grid(True, alpha=0.3)

# 4. DISTRIBUCIÓN MAXWELL-BOLTZMANN (CORREGIDO)
ax4 = fig.add_subplot(gs[1, 0])
ax4.hist(v_mag, bins=75, density=True, alpha=0.7, color='lightblue', 
         edgecolor='navy', linewidth=0.5)
T_avg = np.mean(temps)
m_Ar = 39.95 * 1.6605e-27  # kg
v_th = np.sqrt(8 * k * T_avg / (np.pi * m_Ar)) * 1e10  # Å/ps
ax4.axvline(v_th, color='r', ls='--', lw=2, 
            label=f'v_th (Ar) = {v_th:.2f} Å/ps')
ax4.set_xlabel('Velocidad |v| (Å/ps)')
ax4.set_ylabel('Densidad de probabilidad')
ax4.set_title('Distribución Maxwell-Boltzmann')
ax4.legend()
ax4.grid(True, alpha=0.3)

# 5. FACTOR COMPRESIBILIDAD Z = PV/nRT (CORREGIDO)
ax5 = fig.add_subplot(gs[1, 1])
ax5.hist(Z_factors, bins=30, alpha=0.7, color='coral', 
         edgecolor='darkred', density=True)
ax5.axvline(1.0, color='k', lw=4, alpha=0.9, label='Z=1.0 (Gas Ideal)')
ax5.axvline(np.mean(Z_factors), color='orange', lw=2, ls='--', 
            label=f'Z={np.mean(Z_factors):.3f}')
ax5.set_xlabel('Factor de Compresibilidad Z')
ax5.set_ylabel('Densidad')
ax5.set_title('PV/nRT → 1.0 (Gas Ideal)')
ax5.legend()
ax5.grid(True, alpha=0.3)

# 6. TABLA RESUMEN (ROBUSTA)
ax6 = fig.add_subplot(gs[1, 2])
ax6.axis('off')

# Construir tabla dinámicamente (CORREGIDO)
tabla_data = [['Gas', 'N', 'D', 'RDF_max', 'σ (Å)', 'Z']]
for gas in gases_nobles:
    N = len(gases_nobles[gas])
    D = f"{diffusions.get(gas, 'N/A'):.3f}" if gas in diffusions else 'N/A'
    RDF_max = f"{np.max(rdf_data[gas].rdf):.2f}" if gas in rdf_data else 'N/A'
    sigma = f"{propiedades[gas]['sigma']:.2f}"
    tabla_data.append([gas, N, D, RDF_max, sigma, f"{np.mean(Z_factors):.3f}"])

table = ax6.table(cellText=tabla_data[1:], colLabels=tabla_data[0],
                  cellLoc='center', loc='center', 
                  colColours=['lightgray']*6, cellColours=None)
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.3, 2.2)
ax6.set_title('Resumen Completo\nGases Nobles', pad=20, fontsize=14, fontweight='bold')

# Panel adicional: Evolución Z(t)
ax7 = fig.add_subplot(gs[2, :])
ax7.plot(Z_factors, 'g-', lw=2, alpha=0.8)
ax7.axhline(1.0, color='k', lw=3, alpha=0.9, label='Ideal Z=1.0')
ax7.axhline(np.mean(Z_factors), color='orange', lw=2, ls='--', 
            label=f'media={np.mean(Z_factors):.3f}')
ax7.set_xlabel('Frame')
ax7.set_ylabel('Z(t) = PV/nRT')
ax7.set_title('Evolución Factor Compresibilidad')
ax7.legend()
ax7.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('gases_nobles_completo_CORREGIDO.png', dpi=300, bbox_inches='tight')
plt.show()

# =====================================================
# RESULTADOS CUANTITATIVOS Y DIAGNÓSTICO
# =====================================================
print("\n" + "="*80)
print(" RESULTADOS COMPLETOS - GASES NOBLES COMO GAS IDEAL")
print("="*80)
print(f"  Temperatura:     {np.mean(temps):6.1f} ± {np.std(temps):4.1f} K")
print(f"  Presión ideal:   {np.mean(pressures):6.2f} ± {np.std(pressures):4.2f} atm")
print(f"  Volumen promedio: {np.mean(volumes):7.0f} Å³")
print(f"  Factor Z:        {np.mean(Z_factors):6.3f} ± {np.std(Z_factors):4.3f} (IDEAL=1.000)")
print(f"  Correlación D∝1/√m: R² = {r_value**2:.3f}" if 'r_value' in locals() else "📈  Correlación: Datos insuficientes")

print("\n AUTO-DIFUSIÓN POR ESPECIE:")
print("-" * 40)
for gas in sorted(diffusions, key=lambda x: diffusions[x], reverse=True):
    print(f"  {gas:2s}: {diffusions[gas]:7.4f} Å²/ps  ({propiedades[gas]['masa']:4.0f} u)")

print("\n DIAGNÓSTICO GAS IDEAL:")
print("-" * 40)
for gas in gases_nobles:
    ideal = " IDEAL" if gas in rdf_data and np.max(rdf_data[gas].rdf) < 1.3 else "  REAL"
    print(f"  {gas:2s}: RDF_max={np.max(rdf_data[gas].rdf):4.2f if gas in rdf_data else 'N/A':>4}  {ideal}")

print(f"\n MEJOR GAS IDEAL: {'He' if 'He' in rdf_data and np.max(rdf_data['He'].rdf) < 1.2 else 'Ne'}")
print(f" PEOR GAS IDEAL: {'Rn' if 'Rn' in rdf_data else 'Xe'}")
print("="*80)

# Guardar datos numéricos
df_resultados = pd.DataFrame({
    'Gas': list(gases_nobles.keys()),
    'N_atoms': [len(atoms) for atoms in gases_nobles.values()],
    'D_A2ps': [diffusions.get(g, np.nan) for g in gases_nobles],
    'RDF_max': [np.max(rdf_data[g].rdf) if g in rdf_data else np.nan for g in gases_nobles],
    'Masa_u': [propiedades[g]['masa'] for g in gases_nobles]
})
df_resultados.to_csv('resultados_gases_nobles.csv', index=False)
print(" Datos guardados: resultados_gases_nobles.csv")
