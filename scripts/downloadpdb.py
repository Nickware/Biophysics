#from pymol import cmd
#for i in range(1, 10):
#    cmd.create(f"copy_{i}",'/tmp/1aa7.pdb')

from pymol import cmd
cmd.reinitialize()
for pdb_id in ["1aa7", "3MD2"]:  
    cmd.fetch(pdb_id)  
    cmd.spectrum("count", selection=pdb_id)  
    cmd.png(f"{pdb_id}.png")