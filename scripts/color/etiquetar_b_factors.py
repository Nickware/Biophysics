from pymol import cmd

def etiquetar_b_factors(objeto):
    """
    Este script elimina moléculas de agua, ajusta el fondo,
    colorea según los B-factors y etiqueta los carbonos alfa.
    """
    # Eliminar moléculas de agua
    cmd.remove(f"{objeto} and resn HOH")

    # Cambiar el fondo a blanco
    cmd.bg_color("white")

    # Obtener la lista de átomos del objeto
    atomos = cmd.get_model(objeto).atom

    # Filtrar los carbonos alfa y obtener sus B-factors
    b_factors = [atomo.b for atomo in atomos if atomo.name == 'CA']

    # Calcular el B-factor mínimo y máximo
    min_b_factor = min(b_factors)
    max_b_factor = max(b_factors)

    # Crear la rampa de color
    cmd.ramp_new('color_bar', objeto, [min_b_factor, max_b_factor], 'rainbow')

    # Colorear el objeto según el B-factor
    cmd.spectrum("b", selection=objeto, quiet=0)

    # Añadir etiquetas con los B-factors
    for atomo in atomos:
        if atomo.name == 'CA':
            etiqueta = f"{atomo.b:.2f}"
            cmd.label(f"{objeto} and chain {atomo.chain} and resi {atomo.resi} and name CA", f'"{etiqueta}"')

    print(f"Rampa de color creada con min B-factor: {min_b_factor} y max B-factor: {max_b_factor}")

# Registrar la función en PyMOL
cmd.extend("etiquetar_b_factors", etiquetar_b_factors)
