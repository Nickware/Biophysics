# ESPResSo Installation Suite (Debian/Ubuntu/Deepin)

Este repositorio contiene un script de automatización avanzado para la compilación e instalación de **ESPResSo** (Extensible Simulation Package for Research on Soft matter). El script está diseñado para gestionar todo el ciclo de vida de la instalación, desde la resolución de dependencias del sistema hasta la configuración de entornos virtuales de Python.

##  Características Principales

* **Soporte Multi-Distribución:** Compatible con Debian, Ubuntu, Deepin, Zorin OS y Linux Mint.
* **Aislamiento de Entorno:** Utiliza `venv` para evitar conflictos con las librerías de Python del sistema.
* **Gestión Inteligente de Permisos:** Aunque el script requiere `sudo` para instalar paquetes de sistema, las carpetas de usuario y el entorno virtual mantienen la propiedad del usuario real para evitar problemas de permisos.
* **Aceleración por Hardware:** Soporte integrado para **NVIDIA CUDA** mediante el flag `--with-cuda`.
* **Compilación en Paralelo:** Detecta automáticamente el número de núcleos disponibles (`nproc`) para minimizar el tiempo de compilación.
* **Validación Post-Instalación:** Realiza un test de importación automático para asegurar que el motor de simulación está listo.

##  Requisitos Previos

* Sistema operativo basado en Debian (64 bits recomendado).
* Conexión a internet estable.
* Controladores de NVIDIA instalados (solo si se planea usar `--with-cuda`).

##  Instrucciones de Uso

### 1. Preparación

Clona este repositorio o guarda el contenido del script en un archivo:

```bash
chmod +x install_espresso_debian.sh

```

### 2. Instalación Estándar (CPU)

Ideal para simulaciones generales o desarrollo de scripts de análisis:

```bash
sudo ./install_espresso_debian.sh

```

### 3. Instalación con Soporte GPU (CUDA)

Para simulaciones de dinámica molecular a gran escala que requieran aceleración por GPU:

```bash
sudo ./install_espresso_debian.sh --with-cuda

```

---

##  Estructura del Entorno Creado

Tras la ejecución, el script organiza los archivos en la carpeta `~/espresso` de la siguiente manera:

* `/ESPResSo`: Código fuente original clonado desde el repositorio oficial.
* `/build`: Binarios compilados y la interfaz de Python.
* `/espresso_env`: Entorno virtual de Python con todas las dependencias científicas (`numpy`, `scipy`, `cython`, etc.).
* `activate_espresso.sh`: Script de entrada rápida.

---

##  Activación del Entorno

Para comenzar a trabajar en tus simulaciones, no es necesario configurar variables de entorno manualmente. Simplemente ejecuta:

```bash
source ~/espresso/activate_espresso.sh

```

Esto activará el entorno virtual y configurará el `PYTHONPATH` necesario para que Python reconozca el módulo `espressomd`.

### Verificación manual

Una vez activado, puedes comprobar el estado con:

```bash
python3 -c "import espressomd; print(espressomd.features())"

```

---

## Dependencias de Python incluidas

El script instala automáticamente un stack científico optimizado:

* **Core:** `cython`, `numpy`, `scipy`.
* **Visualización & Análisis:** `matplotlib`, `PyOpenGL`, `tqdm`, `pint`.
* **Notebooks:** `jupyterlab`, `nbconvert`.
* **I/O:** `h5py`, `lxml`.

---

## Solución de Problemas (FAQ)

**1. ¿Por qué mis archivos pertenecen a root?**
El script actual ha sido mejorado para usar `SUDO_USER`. Si usaste una versión antigua, puedes arreglarlo con: `sudo chown -R $USER:$USER ~/espresso`.

**2. Error de compilación en Deepin/Debian antiguo:**
ESPResSo requiere un compilador compatible con C++17. Asegúrate de que `gcc --version` sea superior a 7.0.

**3. CUDA no detectado:**
Asegúrate de que los drivers de NVIDIA estén cargados (`nvidia-smi`). Si instalaste CUDA en una ruta no estándar, podrías necesitar ajustar el flag `-DCUDA_TOOLKIT_ROOT_DIR` en el script.

---

**Autores:** N.Torres | **Versión:** 0.0.2 | **Licencia:** GPLv3
