# test_mdanalysis_basico.py
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD

def test_carga_universo():
    """Verifica que un Universe se construye con el número esperado de átomos."""
    u = mda.Universe(PSF, DCD)
    assert len(u.atoms) > 0
    assert u.trajectory.n_frames > 0

def test_seleccion_basica():
    """Verifica que una selección de átomos de backbone no esté vacía."""
    u = mda.Universe(PSF, DCD)
    backbone = u.select_atoms("backbone")
    assert len(backbone) > 0
