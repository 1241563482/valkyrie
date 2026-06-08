# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 10:28:31 2023

@author: zyj
"""

from ase.io import read
import os
import numpy as np


pos = read("POSCAR")
atoms = set(pos.get_chemical_symbols())
if os.path.exists("DOS.dat"):
    os.remove("DOS.dat")
stacked_data = []

#Energy s  py pz px dxy dyz dz2 dxz dx2 tot
#  0    1  2  3  4  5   6   7   8   9   10

# Header
header = "Energy".ljust(8) + "Total".ljust(8)
for i in atoms:
    for j in ["s", "p", "d", "tot"]:
        header = header + f"{i}_{j}".ljust(8)
    #header = header + f"{i}_s".ljust(28)+f"{i}_p".ljust(28)+f"{i}_d".ljust(28)+f"{i}_tot".ljust(28)

# TDOS
for i in [0, 1]:
    data = np.genfromtxt("TDOS.dat", usecols = i)
    stacked_data.append(data)

# PDOS
for i in atoms:
    s = np.genfromtxt("PDOS_{}.dat".format(i), usecols = 1)
    py = np.genfromtxt("PDOS_{}.dat".format(i), usecols = 2)
    pz = np.genfromtxt("PDOS_{}.dat".format(i), usecols = 3)
    px = np.genfromtxt("PDOS_{}.dat".format(i), usecols = 4)
    dxy = np.genfromtxt("PDOS_{}.dat".format(i), usecols = 5)
    dyz = np.genfromtxt("PDOS_{}.dat".format(i), usecols = 6)
    dz2 = np.genfromtxt("PDOS_{}.dat".format(i), usecols = 7)
    dxz = np.genfromtxt("PDOS_{}.dat".format(i), usecols = 8)
    dx2 = np.genfromtxt("PDOS_{}.dat".format(i), usecols = 9)
    tot = np.genfromtxt("PDOS_{}.dat".format(i), usecols = 10)
    
    p = py + px + pz
    d = dxy + dyz + dz2 + dxz + dx2
    
    stacked_data.append(s)
    stacked_data.append(p)
    stacked_data.append(d)
    stacked_data.append(tot)
    
# Output
stacked_data = np.transpose(stacked_data)
stacked_data = np.vstack(stacked_data)
#stacked_data = np.round(stacked_data, decimals=4)
np.savetxt("DOS.dat", stacked_data, delimiter = "  ", fmt='%.4f', header = header)

os.system(f"sed -i 's/# /#/g' DOS.dat")
