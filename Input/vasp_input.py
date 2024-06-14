# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 13:27:57 2023

@author: zyj
"""

import subprocess, platform, shutil
from Input import get_info


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


    
    
def dos_ggau(encut, spin):
    incar = """# INCAR for dos-gga+U
EDIFF = 1E-7
EDIFFG = -0.005
ISTART = 1
ICHARG = 11
ISMEAR = -5
SIGMA = 0.05
IBRION = -1
PREC = Accurate
ENCUT = {}
ISIF = 2
NCORE = 6
LORBIT = 11
KSPACING = 0.157
NEDOS =  2001
LWAVE = .FALSE.
LCHARG = .FALSE.
# GGA+U part
LASPH = .T.
LDAU = .T.
LDAUTYPE = 2
LDAUL = 2 -1 -1 -1 -1  # 2 for +U, -1 for not +U
LDAUU = 3 0 0 0 0  # Ueff = U-J
LDAUJ = 1 0 0 0 0
LMAXMIX = 4  # If f electron, use 6""".format(encut)
    with open("INCAR", "w") as file:
        file.write(incar)
    if spin == True:
        with open("INCAR", "a") as file:
            print("\nISPIN = 2", file = file)

def band_hse(encut, spin):
    incar = """# INCAR for band-hse
EDIFF = 1E-5
EDIFFG = -0.01
ISTART = 1
ICHARG = 11
ISMEAR = 0
SIGMA = 0.03
IBRION = -1
PREC = Accurate
ENCUT = {}
ISIF = 2
NCORE = 6
LORBIT = 11
LWAVE = .FALSE.
LCHARG = .FALSE.
# HSE06
LVHAR = .T.
LHFCALC = .TRUE.
HFSCREEN = 0.2
TIME = 0.4
AEXX = 0.25
ALGO = A
# 2D
#IVDW = 11
#VDW_S8 = 0.722
#VDW_SR = 1.217""".format(encut)
    with open("INCAR", "w") as file:
        file.write(incar)
    if spin == True:
        with open("INCAR", "a") as file:
            print("\nISPIN = 2", file = file)            

def dos_hse(encut, spin):
    incar = """# INCAR for dos-hse
EDIFF = 1E-5
EDIFFG = -0.01
ISTART = 1
ICHARG = 11
ISMEAR = 0
SIGMA = 0.03
IBRION = -1
PREC = Accurate
ENCUT = {}
ISIF = 2
NCORE = 6
LORBIT = 11
LWAVE = .FALSE.
LCHARG = .FALSE.
# HSE06
LVHAR = .T.
LHFCALC = .TRUE.
HFSCREEN = 0.2
TIME = 0.4
AEXX = 0.25
ALGO = A
# 2D
#IVDW = 11
#VDW_S8 = 0.722
#VDW_SR = 1.217""".format(encut)
    with open("INCAR", "w") as file:
        file.write(incar)
    if spin == True:
        with open("INCAR", "a") as file:
            print("\nISPIN = 2", file = file)


def elf(encut):
    incar = """# INCAR for ELF
system = ELF
ISTART = 0
EDIFF = 1E-7
EDIFFG = -1E-3
ICHARG = 2
ISMEAR = 0
ENCUT = {}
IBRION = -1
ISIF = 2
PREC = Accurate
NSW = 0
NPAR = 4
LORBIT = 11
LELF = .True.
KSPACING = 0.157
LCHARG = .False.
LWAVE = .False.""".format(encut)
    with open("INCAR", "w") as file:
        file.write(incar)




