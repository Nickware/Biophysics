# #!/bin/bash
# Script de instalación de ESPResSo para distribuciones Fedora/RHEL
# Este script instala dependencias y configura ESPResSo desde el código fuente
# Autores: N.Torres | Versión: 0.0.1 | Fecha: 18-09-2025

set -e  # Detener el script en caso de error

# Colores para output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
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

print_step() {
    echo -e "${BLUE}[STEP]${NC} $1"
}

# Detectar distribución
detect_distro() {
    if [ -f /etc/fedora-release ]; then
        DISTRO="fedora"
        VERSION=$(grep -oE '[0-9]+' /etc/fedora-release)
    elif [ -f /etc/redhat-release ]; then
        if grep -q "CentOS" /etc/redhat-release; then
            DISTRO="centos"
        elif grep -q "Rocky" /etc/redhat-release; then
            DISTRO="rocky"
        elif grep -q "AlmaLinux" /etc/redhat-release; then
            DISTRO="almalinux"
        else
            DISTRO="rhel"
        fi
        VERSION=$(grep -oE '[0-9]+' /etc/redhat-release | head -1)
    else
        print_error "Distribución no compatible. Este script es para Fedora/RHEL/CentOS/Rocky/AlmaLinux"
        exit 1
    fi
}

# Verificar si es sudo
check_root() {
    if [ "$EUID" -ne 0 ]; then
        print_error "Por favor ejecuta este script con sudo:"
        echo "sudo $0"
        exit 1
    fi
}

# Instalar dependencias según la distribución
install_dependencies() {
    print_step "Instalando dependencias del sistema..."
    
    case $DISTRO in
        fedora)
            dnf install -y epel-release
            dnf groupinstall -y "Development Tools" "Development Libraries"
            dnf install -y cmake cmake-gui python3-devel openmpi-devel \
                boost-devel fftw-devel fftw-static fftw-openmpi-devel \
                hdf5-devel hdf5-openmpi-devel python3-pip gsl-devel \
                freeglut-devel ffmpeg doxygen graphviz \
                environment-modules redhat-rpm-config \
                libicu-devel openssl-devel libxml2-devel \
                openblas-devel lapack-devel blas-devel
            ;;
        centos|rocky|almalinux|rhel)
            if [ $VERSION -eq 7 ]; then
                yum install -y epel-release
                yum groupinstall -y "Development Tools"
                yum install -y cmake3 cmake3-gui python3-devel openmpi-devel \
                    boost-devel fftw-devel fftw-static fftw-openmpi-devel \
                    hdf5-devel hdf5-openmpi-devel python3-pip gsl-devel \
                    freeglut-devel ffmpeg doxygen graphviz \
                    environment-modules redhat-rpm-config \
                    libicu-devel openssl-devel libxml2-devel
                # Crear symlink para cmake3 -> cmake
                ln -sf /usr/bin/cmake3 /usr/local/bin/cmake
            else
                dnf install -y epel-release
                dnf groupinstall -y "Development Tools" "Development Libraries"
                dnf install -y cmake cmake-gui python3-devel openmpi-devel \
                    boost-devel fftw-devel fftw-static fftw-openmpi-devel \
                    hdf5-devel hdf5-openmpi-devel python3-pip gsl-devel \
                    freeglut-devel ffmpeg doxygen graphviz \
                    environment-modules redhat-rpm-config \
                    libicu-devel openssl-devel libxml2-devel \
                    openblas-devel lapack-devel blas-devel
            fi
            ;;
    esac
}

# Instalar CUDA toolkit si se solicita
install_cuda() {
    if [ "$WITH_CUDA" = true ]; then
        print_step "Instalando CUDA toolkit..."
        
        case $DISTRO in
            fedora)
                dnf install -y kernel-devel kernel-headers
                dnf config-manager --add-repo https://developer.download.nvidia.com/compute/cuda/repos/fedora${VERSION}/x86_64/cuda-fedora${VERSION}.repo
                dnf install -y cuda-toolkit
                ;;
            centos|rocky|almalinux|rhel)
                if [ $VERSION -eq 7 ]; then
                    yum install -y kernel-devel kernel-headers
                    yum config-manager --add-repo https://developer.download.nvidia.com/compute/cuda/repos/rhel7/x86_64/cuda-rhel7.repo
                    yum install -y cuda-toolkit
                else
                    dnf install -y kernel-devel kernel-headers
                    dnf config-manager --add-repo https://developer.download.nvidia.com/compute/cuda/repos/rhel${VERSION}/x86_64/cuda-rhel${VERSION}.repo
                    dnf install -y cuda-toolkit
                fi
                ;;
        esac
        
        # Configurar variables de entorno para CUDA
        echo 'export PATH=/usr/local/cuda/bin:$PATH' >> /etc/profile.d/cuda.sh
        echo 'export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH' >> /etc/profile.d/cuda.sh
        source /etc/profile.d/cuda.sh
    fi
}

# Configurar entorno Python
setup_python_env() {
    print_step "Configurando entorno Python..."
    
    # Instalar pip y herramientas básicas
    pip3 install --upgrade pip setuptools wheel
    
    # Crear archivo de requisitos
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
lxml
jupyter_console
sphinx
sphinxcontrib-bibtex
sphinx-toggleprompt
scikit-build
pybind11
EOF

    # Instalar dependencias Python
    pip3 install -r /tmp/requirements.txt
}

# Compilar e instalar ESPResSo
compile_espresso() {
    print_step "Compilando ESPResSo..."
    
    cd "$ESPRESSO_DIR"
    
    if [ ! -d "ESPResSo" ]; then
        git clone https://github.com/espressomd/espresso.git ESPResSo
        cd ESPResSo
        git submodule update --init --recursive
    else
        print_warning "Directorio ESPResSo ya existe, actualizando..."
        cd ESPResSo
        git pull
        git submodule update --init --recursive
    fi
    
    mkdir -p build
    cd build
    
    # Configurar flags de compilación
    local CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=Release -DPYTHON_EXECUTABLE=$(which python3)"
    
    if [ "$WITH_CUDA" = true ]; then
        CMAKE_FLAGS="$CMAKE_FLAGS -DESPRESSO_BUILD_WITH_CUDA=ON"
        print_status "Configurando con soporte CUDA..."
        
        # Configurar compilador para CUDA en Fedora
        if [ "$DISTRO" = "fedora" ] && [ $VERSION -ge 36 ]; then
            export CC=gcc
            export CXX=g++
            CMAKE_FLAGS="$CMAKE_FLAGS -DCMAKE_CUDA_FLAGS='--compiler-bindir=/usr/bin/g++'"
        fi
    fi
    
    # Configurar cmake
    if [ "$DISTRO" = "centos" ] && [ $VERSION -eq 7 ]; then
        cmake3 $CMAKE_FLAGS ..
    else
        cmake $CMAKE_FLAGS ..
    fi
    
    # Compilar
    make -j$(nproc)
    
    # Instalar
    make install
    
    # Verificar instalación
    print_step "Verificando instalación..."
    python3 -c "import espressomd; print('ESPResSo version:', espressomd.version())"
}

# Crear entorno de usuario
setup_user_environment() {
    print_step "Configurando entorno de usuario..."
    
    # Crear script de activación
    cat > /usr/local/bin/activate_espresso << 'EOF'
#!/bin/bash
export PYTHONPATH="/usr/local/lib/python3.$(python3 -c "import sys; print(sys.version_info[1])")/site-packages:${PYTHONPATH}"
export LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"
echo "Entorno ESPResSo configurado"
python3 -c "import espressomd; print('ESPResSo version:', espressomd.version())"
EOF
    
    chmod +x /usr/local/bin/activate_espresso
    
    # Crear módulo environment (opcional)
    mkdir -p /etc/modulefiles/espresso
    cat > /etc/modulefiles/espresso/4.3.0 << 'EOF'
#%Module1.0
proc ModulesHelp { } {
    puts stderr "ESPResSo Molecular Dynamics software"
}

set version 4.3.0
set prefix /usr/local

prepend-path PATH $prefix/bin
prepend-path LD_LIBRARY_PATH $prefix/lib
prepend-path PYTHONPATH $prefix/lib/python3.$(python3 -c "import sys; print(sys.version_info[1])")/site-packages
EOF
}

# Main execution
main() {
    print_status "Iniciando instalación de ESPResSo para Fedora/RHEL"
    
    # Configuración
    WITH_CUDA=false
    if [ "$1" == "--with-cuda" ]; then
        WITH_CUDA=true
    fi
    
    ESPRESSO_DIR="/opt/espresso"
    
    # Detectar distribución
    detect_distro
    print_status "Detectada distribución: $DISTRO $VERSION"
    
    # Verificar privilegios
    check_root
    
    # Instalar dependencias
    install_dependencies
    
    # Instalar CUDA si se solicita
    install_cuda
    
    # Configurar Python
    setup_python_env
    
    # Compilar ESPResSo
    compile_espresso
    
    # Configurar entorno
    setup_user_environment
    
    print_status "Instalación completada exitosamente!"
    echo ""
    echo "Para usar ESPResSo:"
    echo "  activate_espresso"
    echo ""
    echo "O manualmente:"
    echo "  export PYTHONPATH=\"/usr/local/lib/python3.\$(python3 -c \"import sys; print(sys.version_info[1])\")/site-packages:\${PYTHONPATH}\""
    echo "  export LD_LIBRARY_PATH=\"/usr/local/lib:\${LD_LIBRARY_PATH}\""
    echo ""
    echo "Para verificar:"
    echo "  python3 -c \"import espressomd; print(espressomd.version())\""
}

# Ejecutar main
main "$@"
