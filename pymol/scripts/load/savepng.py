import pymol
from pymol import cmd # import PyMOL commands


#iniciar PyMol en modo sin GUI (headless)
pymol.finish_launching(['pymol', '-qc'])

# Comandos de PyMol
cmd.load('1aa7.cif')
cmd.show('cartoon')
cmd.ray(800, 600)
cmd.png('1aa7.png')

# Salir de PyMol
cmd.quit()