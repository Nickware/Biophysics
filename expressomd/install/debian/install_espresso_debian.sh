# #!/bin/bash
# Script de instalación de ESPResSo para distribuciones Debian/Ubuntu
# Este script instala dependencias y configura ESPResSo desde el código fuente
# Autores: N.Torres | Versión: 0.0.1 | Fecha: 18-09-2025

set -e  # Detener el script en caso de error

# Colores para output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Función para imprimir mensajes
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Verificar si es Ubuntu/Debian
if ! command -v lsb_release &> /dev/null || ! lsb_release -i | grep -E "Ubuntu|Debian" &> /dev/null; then
    print_error "Este script está diseñado para distribuciones Debian/Ubuntu"
    exit 1
fi

# Verificar versión de Ubuntu
if lsb_release -i | grep "Ubuntu" &> /dev/null; then
    UBUNTU_VERSION=$(lsb_release -r | awk '{print $2}')
    print_status "Detectado Ubuntu $UBUNTU_VERSION"
fi

# Verificar si es sudo
if [ "$EUID" -ne 0 ]; then
    print_error "Por favor ejecuta este script con sudo:"
    echo "sudo $0"
    exit 1
fi

# Configuración
ESPRESSO_DIR="${HOME}/espresso"
BUILD_DIR="${ESPRESSO_DIR}/build"
VENV_DIR="${ESPRESSO_DIR}/espresso_env"
REQUIREMENTS_FILE="${ESPRESSO_DIR}/requirements.txt"

# Crear archivo de requisitos temporal
cat > /tmp/requirements.txt << 'EOF'
cmake>=3.18
cython>=0.29
numpy>=1.18
scipy>=1.5
packaging
setuptools>=60.0
h5py>=3.0
matplotlib>=3.3
pint
tqdm
PyOpenGL
jupyterlab>=4.3
nbformat
nbconvert
lxml[html_clean]
jupyter_console
sphinx
sphinxcontrib-bibtex
sphinx-toggleprompt
EOF

print_status "Actualizando lista de paquetes..."
apt update

print_status "Instalando dependencias del sistema..."
apt install -y build-essential cmake cmake-curses-gui python3-dev \
    openmpi-bin libboost-all-dev libfftw3-dev libfftw3-mpi-dev \
    libhdf5-dev libhdf5-openmpi-dev python3-pip libgsl-dev \
    freeglut3-dev ffmpeg doxygen graphviz

# Opcional: Instalar CUDA toolkit si se solicita
if [ "$1" == "--with-cuda" ]; then
    print_status "Instalando CUDA toolkit..."
    apt install -y nvidia-cuda-toolkit
fi

print_status "Creando directorio de trabajo..."
mkdir -p "${ESPRESSO_DIR}"
mkdir -p "${BUILD_DIR}"
cp /tmp/requirements.txt "${REQUIREMENTS_FILE}"

print_status "Creando entorno virtual Python..."
python3 -m venv "${VENV_DIR}"
source "${VENV_DIR}/bin/activate"

print_status "Actualizando pip y setuptools..."
pip install --upgrade pip setuptools

print_status "Instalando dependencias Python..."
pip install -r "${REQUIREMENTS_FILE}"

print_status "Descargando ESPResSo..."
cd "${ESPRESSO_DIR}"
if [ ! -d "ESPResSo" ]; then
    git clone https://github.com/espressomd/espresso.git ESPResSo
else
    print_warning "Directorio ESPResSo ya existe, omitiendo clonación"
fi

cd ESPResSo
git submodule update --init --recursive

print_status "Configurando build de ESPResSo..."
cd "${BUILD_DIR}"

# Configurar flags de compilación
CMAKE_FLAGS=".. -DCMAKE_BUILD_TYPE=Release -DPYTHON_EXECUTABLE=$(which python3)"

if [ "$1" == "--with-cuda" ]; then
    CMAKE_FLAGS="${CMAKE_FLAGS} -DESPRESSO_BUILD_WITH_CUDA=ON"
    print_status "Configurando con soporte CUDA..."
fi

cmake ${CMAKE_FLAGS}

print_status "Compilando ESPResSo (esto puede tomar varios minutos)..."
make -j$(nproc)

print_status "Instalando ESPResSo..."
make install

print_status "Verificando instalación..."
python3 -c "import espressomd; print('ESPResSo version:', espressomd.version())"

# Crear script de activación
cat > "${ESPRESSO_DIR}/activate_espresso.sh" << 'EOF'
#!/bin/bash
source ${VENV_DIR}/bin/activate
export PYTHONPATH="${BUILD_DIR}/python:${PYTHONPATH}"
echo "Entorno ESPResSo activado. Versión:"
python -c "import espressomd; print(espressomd.version())"
EOF

chmod +x "${ESPRESSO_DIR}/activate_espresso.sh"

print_status "Instalación completada exitosamente!"
echo ""
echo "Para usar ESPResSo:"
echo "  source ${ESPRESSO_DIR}/activate_espresso.sh"
echo ""
echo "O para activar manualmente:"
echo "  source ${VENV_DIR}/bin/activate"
echo "  export PYTHONPATH=\"${BUILD_DIR}/python:\${PYTHONPATH}\""
echo ""
echo "Para comprobar la instalación:"
echo "  python -c \"import espressomd; print(espressomd.version())\""

# Limpiar archivo temporal
rm -f /tmp/requirements.txt
