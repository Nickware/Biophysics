from pymol import cmd
import numpy as np
distances = np.loadtxt("distancias.dat")
cmd.distance("medida", "residue 10", "residue 20")
