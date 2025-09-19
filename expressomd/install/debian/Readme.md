

## Script de instalación de ESPResSo para Distribuciones Debian/Ubuntu 
## Instrucciones de uso:

1. **Guardar el script** como `install_espresso_debian.sh`

2. **Hacerlo ejecutable:**
   ```bash
   chmod +x install_espresso_debian.sh
   ```

3. **Ejecutar con sudo:**
   ```bash
   sudo ./install_espresso_debian.sh
   ```

4. **Para instalación con soporte CUDA:**
   ```bash
   sudo ./install_espresso_debian.sh --with-cuda
   ```

## Características del script:

- ✅ **Verificación del sistema**: Confirma que es Debian/Ubuntu
- ✅ **Instalación automática**: Todas las dependencias del sistema
- ✅ **Entorno virtual**: Aísla las dependencias Python
- ✅ **Soporte CUDA opcional**: Parámetro `--with-cuda`
- ✅ **Compilación optimizada**: Usa todos los cores disponibles (`nproc`)
- ✅ **Verificación**: Confirma que la instalación fue exitosa
- ✅ **Script de activación**: Crea un script fácil para usar el entorno

## Script de activación simplificado:

También se puede crear un script más simple para usuarios que solo quieran usar ESPResSo:

```bash
#!/bin/bash
# activate_espresso.sh
source ~/espresso/espresso_env/bin/activate
export PYTHONPATH=~/espresso/build/python:$PYTHONPATH
echo "ESPResSo environment activated"
python -c "import espressomd; print('Version:', espressomd.version())"
```

Este script automatiza completamente el proceso de instalación siguiendo las instrucciones oficiales de ESPResSo para Ubuntu/Debian.