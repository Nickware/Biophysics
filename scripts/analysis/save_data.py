import MDAnalysis as mda
import numpy as np

# Cargar estructura
u = mda.Universe("1YKT.pdb")

# Definir los pares de residuos/átomos a analizar
pares_residuos = [
    (10, 20),  # Residuo 10 vs Residuo 20 (CA)
    (15, 30),  # Residuo 15 vs Residuo 30
    (5, 25),   # Residuo 5 vs Residuo 25
    # ¡Añade más pares aquí!
]

# Abrir archivo para guardar datos
with open("distancias.dat", "w") as f:
    # Escribir cabecera (nombres de columnas)
    f.write("# Tiempo(ps) " + " ".join([f"Distancia_{r1}-{r2}(A)" for r1, r2 in pares_residuos]) + "\n")
    
    # Calcular distancias para cada par
    distancias = []
    for r1, r2 in pares_residuos:
        atom1 = u.select_atoms(f"resid {r1} and name CA")[0]  # Átomo CA del residuo r1
        atom2 = u.select_atoms(f"resid {r2} and name CA")[0]  # Átomo CA del residuo r2
        distancia = np.linalg.norm(atom1.position - atom2.position)
        distancias.append(distancia)
    
    # Guardar datos (tiempo=0 para estructura estática)
    f.write(f"0.0 " + " ".join([f"{d:.2f}" for d in distancias]) + "\n")
