from pymol import cmd

# Descargar un complejo real desde la RCSB PDB: 3PTB = tripsina bovina
# con el inhibidor benzamidina (código de residuo BEN) ya presente en el
# mismo archivo — no hay que inventar un .mol2 aparte para el ligando.
cmd.fetch("3ptb", "complejo", type="pdb")

# Separar receptor y ligando en objetos distintos a partir del complejo descargado
cmd.create("prot", "complejo and polymer.protein")
cmd.create("lig", "complejo and resn BEN")

# Selección y representación
cmd.select("active_site", "prot within 5 of lig")
cmd.show("cartoon", "prot")
cmd.show("sticks", "lig")
cmd.color("cyan", "prot")

# Renderizado y exportación
cmd.bg_color("white")
cmd.ray(1200, 900)
cmd.png("figura.png", dpi=300)