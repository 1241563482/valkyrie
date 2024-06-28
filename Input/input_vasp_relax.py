# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 16:42:36 2023

@author: zyj
"""

import shutil
from Input import get_info

def relax(encut, pressure, spin, fd, u_atom, u_value, lmaxmix, optcell):
    pressure = pressure * 10

    if optcell == True:
        incar = """# INCAR for relax-optcell
PREC = Accurate
EDIFF = 1e-6
EDIFFG = -1e-3
ENCUT = {}
KSPACING = 0.157
SIGMA = 0.030
POTIM = 0.020
IBRION = 2
ISIF = 3
NSW = 150
ISMEAR = 0
LWAVE = F
LCHARG = F
NPAR = 2
KPAR = 2
IVDW = 11
""".format(encut)
        
        OPTCELL = """100
110
000"""
        with open("OPTCELL", "w") as file:
            print(OPTCELL, file = file)
        
    else:
        incar = """# INCAR for relax
PREC = Accurate
EDIFF = 1e-6
EDIFFG = -1e-3
ENCUT = {}
KSPACING = 0.157
SIGMA = 0.030
POTIM = 0.020
IBRION = 2
ISIF = 3
NSW = 150
ISMEAR = 0
LWAVE = F
LCHARG = F
PSTRESS = {}
NPAR = 2
KPAR = 2""".format(encut, pressure)
    with open("INCAR", "w")as file:
        file.write(incar)


    # Spin
    if spin == True:
        with open("INCAR", "a") as file:
            print("\nISPIN = 2", file = file)
            print("\nLORBIT = 11", file = file)
            print("<=> Valkyrie: Default is FM, add MAGMOM in INCAR for other mag-configurations.")
    
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
    
    # GGA+U part
    if u_atom is not None:
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

        
        
