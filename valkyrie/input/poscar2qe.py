from ase.io import read, write
import os

def read_poscar(poscar):
    f = read(poscar)
    f.write("pw.in", format = "espresso-in")
    ntyp = 0
    nat = 0
    with open("pw.in", "r") as file:
        lines = file.readlines()
    # Read nat and ntyp
    for i in lines:
        if "ntyp" in i:
            ntyp = int(i.split()[2])
        if "nat" in i:
            nat = int(i.split()[2])
    if ntyp == 0 or nat == 0:
        raise Exception("Error to read POSCAR, check pw.in.")

    # Read CELL and POSITIONS
    for i, line in enumerate(lines):
        if "CELL_PARAMETERS" in line:
            cell_start = i
            cell_end = i + 3
        if "ATOMIC_POSITIONS" in line:
            pos_start = i
            pos_end = i + nat
    cell = ""
    pos = ""
    for i in range(cell_start, cell_end + 1):
        cell = cell + lines[i]
    for i in range(pos_start, pos_end + 1):
        pos = pos + lines[i]

    os.remove("pw.in")
    return ntyp, nat, cell, pos
