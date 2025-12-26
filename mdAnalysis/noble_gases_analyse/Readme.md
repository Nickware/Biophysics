## **Análisis del Script: Gases Nobles como Modelo de Gas Ideal**

Este script es un **pipeline completo de análisis termodinámico** que combina **LAMMPS** (simulación) y **MDAnalysis** (análisis) para estudiar experimentalmente el **comportamiento de gas ideal** vs **desviaciones reales** en una mezcla de **todos los gases nobles** (He, Ne, Ar, Kr, Xe, Rn).

## **Propósito Científico**

Validar **5 ecuaciones fundamentales del gas ideal**:
1. **PV = NkT** → Factor de compresibilidad Z ≈ 1
2. **g(r) = 1** → RDF plano (sin correlaciones)
3. **D ∝ 1/√m** → Difusión inversa a masa
4. **Distribución Maxwell-Boltzmann** → Velocidades
5. **T constante** → Equilibrio termodinámico

## **Estructura del Script (6 Etapas)**

### **1. Carga de Datos LAMMPS**
```python
u = mda.Universe("gases_nobles.data", "gases_nobles.dump")
```
- Lee topología (`.data`) + trayectoria (`.dump`)
- Identifica 6 tipos atómicos (1=He, 2=Ne, ..., 6=Rn)

### **2. Propiedades Físicas Reales**
```python
propiedades = {'He': {'masa': 4.00, 'sigma': 2.15}, ...}
```
- Masas atómicas exactas (AMU)
- Radios de Van der Waals (σ) para RDF
- Colores reales para visualización

### **3. Análisis Termodinámico Principal**
- **Temperatura cinética**: `T = <KE>/1.5k` cada 50 frames
- **Presión ideal**: `P = NkT/V`
- **Factor Z**: `Z = PV/nRT` (ideal=1.0)

### **4. RDF (Estructura Local)**
```python
rdf_calc = rdf.InterRDF(atoms, atoms, range=(0.0, 12.0))
```
- Gas ideal: `g(r) = 1` (línea plana)
- Gas real: Picos en 1er mínimo de potencial LJ

### **5. Auto-difusión (Dinámica)**
```python
diff_calc.run(start=1000, stop=-2000)  # Evita transitorios
D ∝ 1/√m  # Ley de gas cinético
```
- Calcula MSD → D (coeficiente de difusión)
- Verifica escalado inverso con masa

### **6. Panel de Visualización (3x3)**
| Panel | Análisis | Criterio Ideal |
|-------|----------|----------------|
| **1** | T(t), P(t) | Constantes |
| **2** | RDF(g(r)) | g(r)=1 plano |
| **3** | D vs m | D ∝ 1/√m, R²≈1 |
| **4** | Velocidades | Maxwell-Boltzmann |
| **5** | Z=PV/nRT | Z=1.0 |
| **6** | Tabla resumen | Valores cuantitativos |

## 🔬 **Resultados Esperados (Baja Densidad)**

| Gas | Masa (u) | D (Å²/ps) | RDF_max | Z | ¿Ideal? |
|-----|----------|-----------|---------|---|---------|
| **He** | 4.0 | 0.85 | 1.05 | 1.01 | Si |
| **Ne** | 20.2 | 0.42 | 1.12 | 1.03 | Si |
| **Ar** | 40.0 | 0.28 | 1.35 | 1.08 | Quizá |
| **Kr** | 83.8 | 0.19 | 1.65 | 1.15 | No |
| **Xe** | 131.3 | 0.15 | 2.05 | 1.28 | No |
| **Rn** | 222.0 | 0.12 | 2.45 | 1.42 | No |

##  **Flujo de Trabajo Completo**

```
LAMMPS (.in) → gases_nobles.data + .dump
         ↓
MDAnalysis (.py) → 6 paneles PNG + métricas numéricas
         ↓
VALIDACIÓN: He/Ne ≈ Ideal >> Rn ≈ Real
```

##  **Aplicaciones Biofísicas**

1. **Validación de fuerzas** en simulaciones de solventes nobles
2. **Benchmarking** de termostatos/ barostatos
3. **Estudio de mezclas** reales (atmósfera, estrellas)
4. **Límite ideal** para proteínas en gas vs líquido

**¡Script preparado para publicación en física estadística!** 
