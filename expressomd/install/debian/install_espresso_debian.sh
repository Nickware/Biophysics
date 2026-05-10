#!/bin/bash
# Script de instalación de ESPResSo para Debian, Ubuntu, Deepin y derivados
# Autores: N.Torres | Versión: 0.0.2 | Fecha: 09-05-2026

set -e  # Detener en caso de error

# Colores para output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Funciones de logs
print_status() { echo -e "${GREEN}[INFO]${NC} $1"; }
print_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
print_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# 1. Verificación de Distribución (Incluye Deepin)
if [ -f /etc/os-release ]; then
    . /etc/os-release
    OS_NAME=$ID
    OS_PRETTY=$PRETTY_NAME
else
    print_error "No se pudo determinar la distribución."
    exit 1
fi

# Validar si es una familia compatible (Debian, Ubuntu, Deepin, Zorin, etc.)
case "$OS_NAME" in
    ubuntu|debian|deepin|zorin|linuxmint)
        print_status "Distribución detectada: $OS_PRETTY"
        ;;
    *)
        print_error "Distribución '$OS_NAME' no soportada oficialmente por este script."
        exit 1
        ;;
esac

# 2. Verificación de privilegios
if [ "$EUID" -ne 0 ]; then
    print_error "Por favor ejecuta este script con sudo o como root."
    exit 1
fi

# 3. Configuración de rutas (Asegurando que no se usen rutas relativas de root)
# Si se ejecuta con sudo, queremos el HOME del usuario real, no de /root
REAL_USER=${SUDO_USER:-$USER}
REAL_HOME=$(getent passwd "$REAL_USER" | cut -d: -f6)

ESPRESSO_DIR="${REAL_HOME}/espresso"
BUILD_DIR="${ESPRESSO_DIR}/build"
VENV_DIR="${ESPRESSO_DIR}/espresso_env"
REQUIREMENTS_FILE="${ESPRESSO_DIR}/requirements.txt"

# 4. Dependencias del Sistema
print_status "Actualizando repositorios..."
apt update

print_status "Instalando dependencias base y científicas..."
# Se agrega python3-venv y git explícitamente
apt install -y build-essential cmake cmake-curses-gui python3-dev python3-venv git \
    openmpi-bin libboost-all-dev libfftw3-dev libfftw3-mpi-dev \
    libhdf5-dev libhdf5-openmpi-dev python3-pip libgsl-dev \
    freeglut3-dev ffmpeg doxygen graphviz libhdf5-mpi-dev

# Soporte para CUDA
if [[ "$*" == *"--with-cuda"* ]]; then
    print_status "Instalando CUDA toolkit..."
    apt install -y nvidia-cuda-toolkit
fi

# 5. Preparación del entorno Python
print_status "Configurando entorno en $ESPRESSO_DIR"
mkdir -p "$ESPRESSO_DIR" "$BUILD_DIR"

cat > "$REQUIREMENTS_FILE" << 'EOF'
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
jupyterlab
EOF

# Crear VENV con el usuario real para evitar problemas de permisos después
sudo -u "$REAL_USER" python3 -m venv "$VENV_DIR"
# Usamos un subshell para las operaciones de Python del usuario
sudo -u "$REAL_USER" bash << EOF
    source "$VENV_DIR/bin/activate"
    pip install --upgrade pip setuptools wheel
    pip install -r "$REQUIREMENTS_FILE"
EOF

# 6. Descarga y compilación
print_status "Obteniendo código fuente de ESPResSo..."
if [ ! -d "${ESPRESSO_DIR}/ESPResSo" ]; then
    sudo -u "$REAL_USER" git clone https://github.com/espressomd/espresso.git "${ESPRESSO_DIR}/ESPResSo"
fi

cd "${ESPRESSO_DIR}/ESPResSo"
sudo -u "$REAL_USER" git submodule update --init --recursive

print_status "Iniciando configuración de CMake..."
cd "$BUILD_DIR"

CMAKE_OPTS="-DCMAKE_BUILD_TYPE=Release -DPYTHON_EXECUTABLE=$VENV_DIR/bin/python3"

if [[ "$*" == *"--with-cuda"* ]]; then
    CMAKE_OPTS="$CMAKE_OPTS -DESPRESSO_BUILD_WITH_CUDA=ON"
fi

sudo -u "$REAL_USER" cmake $CMAKE_OPTS ../ESPResSo

print_status "Compilando con $(nproc) núcleos..."
sudo -u "$REAL_USER" make -j$(nproc)

# 7. Finalización y Scripts de activación
print_status "Creando script de acceso rápido..."
cat > "${ESPRESSO_DIR}/activate_espresso.sh" << EOF
#!/bin/bash
source "$VENV_DIR/bin/activate"
export PYTHONPATH="$BUILD_DIR/python:\$PYTHONPATH"
echo "Entorno ESPResSo (Deepin/Debian Ready) activado."
EOF

chown "$REAL_USER:$REAL_USER" -R "$ESPRESSO_DIR"
chmod +x "${ESPRESSO_DIR}/activate_espresso.sh"

print_status "¡Instalación finalizada!"
echo -e "Para comenzar, ejecuta: ${YELLOW}source ${ESPRESSO_DIR}/activate_espresso.sh${NC}"
# Limpiar archivo temporal
rm -f /tmp/requirements.txt
