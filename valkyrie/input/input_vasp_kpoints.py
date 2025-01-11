import numpy as np
from ase.io import read

def kpoints(ka, kb, kc):
    kpoints = """Automatic generation 
0 
Gamma 
{} {} {}
0 0 0
""".format(ka, kb, kc)
    with open ("KPOINTS", "w") as file:
        print(kpoints, file = file)
    return 0


def kpoints_byhand(poscar, ceiling):
    atom_size = [ max(poscar.get_positions()[:,i]) - min(poscar.get_positions()[:,i]) for i in [0, 1, 2] ]
    cell_size = [ np.linalg.norm(poscar.cell[i]) for i in [0, 1, 2] ]
    k = []
    for i in [0, 1, 2]:
        if cell_size[i] - atom_size[i] >= 8: # If vacuum size >= 8A, then set the K-mesh as 1.
            k.append(1)
        else:
            k.append( int( max(1, np.ceil(ceiling / cell_size[i]))) )
    kpoints(k[0], k[1], k[2])
    return 0


if __name__ == "__main__":
    poscar = read("POSCAR")
    kpoints_byhand(poscar)
