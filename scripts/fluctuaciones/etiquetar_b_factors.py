from pymol import cmd

def etiquetar_b_factors(objeto):
    """
    Este script elimina moleculas de agua, ajusta el fondo,
    colorea segun los B-factors y etiqueta los carbonos alfa
    """
    # Eliminar moléculas de agua
    cmd.remove(f"{objeto} and resn HOH")
    
    # Ajustar el fondo
    cmd.bg_color("white")

    # Obtener la lista de átomos del objeto
    atomos = cmd.get_model(objeto).atom

    # Filtrar los carbonos alfa y obtener sus B-factors
    b_factors = [atom.b for atom in atomos if atom.name == "CA"]

    # Calcular el B-factor mínimo y máximo
    b_min = min(b_factors)
    b_max = max(b_factors)

    # Crear rampa de color
    cmd.ramp_new('color_bar', objeto, [b_min, b_max], 'rainbow')

    # Colorear el objeto según los B-factors
    cmd.spectrum("b", selection=objeto, quiet=0)

    # Añadir etiquetas con los B-factors
    for atom in atomos:
        if atom.name == "CA":
            etiqueta = f"{atom.b:.2f}"
            cmd.label(f"{objeto} and chain {atom.chain} and resi {atom.resi} and name CA", f'"{etiqueta}"')
    
    print(f"Rampa de color creada con mínimo B-factor: {b_min} y máximo B-factor: {b_max}")

# Registrar la función en PyMOL
cmd.extend("etiquetar_b_factors", etiquetar_b_factors)