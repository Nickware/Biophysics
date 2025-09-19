## Script de instalación de ESPResSo para Distribuciones Fedora/RHEL
## Instrucciones de uso:

1. **Guardar el script** como `install_espresso_fedora.sh`

2. **Hacerlo ejecutable:**
   ```bash
   chmod +x install_espresso_fedora.sh
   ```

3. **Ejecutar con sudo:**
   ```bash
   sudo ./install_espresso_fedora.sh
   ```

4. **Para instalación con soporte CUDA:**
   ```bash
   sudo ./install_espresso_fedora.sh --with-cuda
   ```

## Características específicas para Fedora/RHEL:

- **Detección automática** de distribución (Fedora, RHEL, CentOS, Rocky, AlmaLinux)
- **Soporte para múltiples versiones** (Fedora 36+, RHEL/CentOS 7+)
- **Manejo de repositorios EPEL** automático
- **Soporte para cmake3** en CentOS/RHEL 7
- **Configuración de CUDA** con repositorios oficiales NVIDIA
- **Módulos environment** para gestión de entornos
- **Dependencias específicas** de Fedora/RHEL

## Script de activación rápido:

```bash
#!/bin/bash
# espresso_env.sh
export PYTHONPATH="/usr/local/lib/python3.$(python3 -c 'import sys; print(sys.version_info[1])')/site-packages:${PYTHONPATH}"
export LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"
echo "Entorno ESPResSo configurado para Fedora/RHEL"
python3 -c "import espressomd; print('Versión:', espressomd.version())"
```

Este script maneja las particularidades de las distribuciones basadas en Fedora/RHEL y proporciona una instalación completa y optimizada de ESPResSo.