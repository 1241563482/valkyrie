# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 20:45:48 2023

@author: zyj
"""

import shutil
from Input import get_info


def ph(encut, dim, spin, fd, symmetry, fun, u_atom, u_value, lmaxmix):
    incar = """# INCAR for ph
ISTART = 0        
#LREAL = Auto    
ENCUT = {}   
PREC = Accurate     
LWAVE = .FALSE.  
LCHARG = .FALSE.    
ADDGRID = .TRUE. 
ISMEAR = 0          
SIGMA = 0.06   
NELM = 150         
EDIFF = 1E-07
EDIFFG = -1E-03
ALGO = N
NPAR = 4
SYMPREC = 1E-12""".format(encut)
    with open("INCAR", "w") as file:
        file.write(incar)
    # Spin
    if spin == True:
        with open("INCAR", "a") as file:
            print("\nISPIN = 2", file = file)
    # Fermi-dirac smearing
    if fd != None:
        with open('INCAR', 'r') as file1, open('INCAR-tmp', 'w') as file2:
            for line in file1:
                if 'ISMEAR' not in line and 'SIGMA' not in line:
                    file2.write(line)
            kB = 8.617330337217213e-05
            print("\nISMEAR = -1\nSIGMA = {:.5f}".format(kB * fd), file = file2)
            print("<=> Valkyrie: Considering electron enthalpy at {} K".format(fd))
        shutil.move("INCAR-tmp", "INCAR")
    # GGA+U
    if fun == "ggau" and u_atom is not None:
        u_value = float(u_value)
        symbols = get_info.get_symbols("POSCAR")
        ldaul = ""
        ldauu = ""
        ldauj = ""
        for i in range(len(symbols)):
            if symbols[i] == u_atom:
                ldaul = ldaul + "2 "
                ldauu = ldauu + str(u_value + 1) + " "
                ldauj = ldauj + "1 "
            else:
                ldaul = ldaul + "-1 "
                ldauu = ldauu + "0 "
                ldauj = ldauj + "0 " 
        ggau_part = """
# GGA+U part
LASPH = .T.
LDAU = .T.
LDAUTYPE = 2
LDAUL = {}  # 2 for +U, -1 for not +U
LDAUU = {}  # Ueff = U-J
LDAUJ = {}
LMAXMIX = {}  # If f electron, use 6
    """.format(ldaul, ldauu, ldauj, lmaxmix)
        with open("INCAR", "a") as file:
            file.write(ggau_part)
    # HSE06
    elif fun == "hse":
        hse_part = """
# HSE part
LHFCALC= .TRUE.
HFSCREEN = 0.2
ALGO = A
TIME = 0.4
AEXX = 0.25
"""
        with open("INCAR", "a") as file:
            file.write(hse_part)


    # band.conf
    band_conf = """
ATOM_NAME = ZYJ
DIM = {} {} {}
PRIMITIVE_AXES = Auto
BAND = Auto
FORCE_SETS = READ
FORCE_CONSTANTS= WRITE""".format(*dim)
    with open("band.conf", "w") as file:
        file.write(band_conf)

    # mesh.conf
    mesh_conf = """
ATOM_NAME = Fe Si
DIM = {} {} {}
#FORCE_CONSTANTS = READ  # For DFPT
MP = 15 15 15
TPROP = T
TMIN = 0
TMAX = 10000
TSTEP = 10""".format(*dim)
    with open("mesh.conf", "w") as file:
        file.write(mesh_conf)

    # loop-out.sh
    loop_out = """
#! /bin/bash
mkdir -p output
for i in `ls POSCAR-* | sed 's/POSCAR-//g'`
do
cp $i/vasprun.xml output/"vasprun.xml-"$i
done
phonopy -f output/*
phonopy -s -p band.conf --tolerance={}
phonopy-bandplot --gnuplot > ph.out
phonopy -t mesh.conf --tolerance={}
python get_free_energy.py
#phonopy --dim="2 2 2" --irreps="0 0 0 1e-3"
""".format(symmetry, symmetry)
    with open("loop-out.sh", "w") as file:
        file.write(loop_out)
        
        
