# Simulación de un polímero de polietileno (PE) en pDynamo usando mecánica molecular 

## Objetivos de la simulación

- **Estudiar la estructura y conformación a nivel molecular del polietileno (PE)** bajo condiciones ambientales típicas (temperatura y presión ambientales), observando su flexibilidad y comportamiento dinámico.
- **Explorar la energía potencial y optimización de la geometría** para entender las conformaciones preferidas y posibles transiciones conformacionales entre diferentes estados.
- **Investigar las interacciones intramoleculares y no enlazantes** (van der Waals, electrostáticas) que definen la estabilidad del polímero.
- **Obtener parámetros moleculares atomísticos** como energías, radios de giro y ángulos entre enlaces, que puedan ser usados posteriormente como insumos para simulaciones mesoscópicas.
- **Analizar la dinámica molecular clásica** para describir movimientos, difusiones y cambios conformacionales en escalas temporales cortas.

A continuación se describe la metodología para realizar una simulación de un polímero de polietileno (PE) en pDynamo usando mecánica molecular clásica, siguiendo pasos similares a los de la glucosa pero adaptados a polímeros:

## Metodología para simulación de polímero PE en pDynamo

1. **Preparación de la estructura del polímero**

   - Generar o adquirir la estructura tridimensional de un fragmento representativo del polietileno (una cadena de unidades repetidas de etileno).
   - Puede estar en formato .pdb o .mol. Asegúrate de que la cadena tenga la longitud adecuada para estudio (ej. decámero).

2. **Importar la estructura en pDynamo**

   ```
   pythonfrom pDynamo import *
   polymer = ImportSystem('polyethylene.pdb')
   ```

3. **Asignar campo de fuerza**

   - Escoger un campo de fuerza clásico adecuado para polímeros hidrocarbonados, como OPLS-AA.

   ```
   pythonmmModel = MMModelOPLS.WithParameterSet('bookSmallExamples')
   polymer.DefineMMModel(mmModel)
   ```

4. **Definir modelo de interacciones no enlazantes**

   ```
   pythonnbModel = NBModelFull.WithDefaults()
   polymer.DefineNBModel(nbModel)
   ```

5. **Resumen e inspección inicial**

   ```
   pythonpolymer.Summary()
   initialEnergy = polymer.Energy()
   print(f"Energía inicial PE: {initialEnergy}")
   ```

6. **Optimización de geometría**

   - Minimizar la energía para relajar la estructura inicial.

   ```
   pythonoptimizer = ConjugateGradientMinimizer()
   optimizer.Iterate(polymer, terminationEnergyDelta=1.0e-4, maxIterations=1000)
   optimizedEnergy = polymer.Energy()
   print(f"Energía optimizada PE: {optimizedEnergy}")
   ```

7. **Configuración y ejecución de dinámica molecular**

   - Puede hacer simulaciones para estudiar la flexibilidad y dinámica de las cadenas de PE.

   ```
   pythonfrom pDynamo.Models import IntegratorLeapFrog
   
   timestep = 1.0  # fs
   temperature = 300  # K
   integrator = IntegratorLeapFrog(timestep)
   
   for step in range(1000):
       integrator.Integrate(polymer)
       energy = polymer.Energy()
       print(f"Paso {step}, energía: {energy}")
   ```

8. **Análisis**

   - Estudiar propiedades como el radio de giro, conformaciones, energía, ángulos entre enlaces C-C.
   - Guardar resultados para visualización y análisis detallado.

------

## Notas adicionales

- Para polímeros como PE, es recomendable usar fragmentos medianamente largos para capturar comportamiento real de cadena.
- Explorar diferentes condiciones de temperatura y presión para simular propiedades termodinámicas.
- También es posible agregar solventes (como agua) para estudiar interacción polímero-solvente.

Esta metodología básica permitirá analizar el polietileno desde la mecánica molecular clásica en pDynamo, brindando información útil sobre estructura y dinámica del polímero.[pdynamo+1](https://www.pdynamo.org/)

1. [https://www.pdynamo.org](https://www.pdynamo.org/)
2. https://www.pdynamo.org/tutorials/molecular-dynamics