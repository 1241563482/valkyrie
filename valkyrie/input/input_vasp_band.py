import shutil
from Input import get_info


def band(encut, pressure, spin, fd, fun, u_atom, u_value, lmaxmix):
    incar = """# INCAR for band
EDIFF = 1E-6
EDIFFG = -0.005
ISTART = 1
ICHARG = 11
ISMEAR = 0
SIGMA = 0.03
IBRION = -1
PREC = Accurate
ENCUT = {}
ISIF = 2
LORBIT = 11
LWAVE  = .FALSE.
LCHARG = .FALSE.
NPAR = 2
KPAR = 2""".format(encut)
    with open("INCAR", "w") as file:
        file.write(incar)
    if fd != None:
        with open('INCAR', 'r') as file1, open('INCAR-tmp', 'w') as file2:
            for line in file1:
                if 'ISMEAR' not in line and 'SIGMA' not in line:
                    file2.write(line)
            kB = 8.617330337217213e-05
            print("\nISMEAR = -1\nSIGMA = {:.5f}".format(kB * fd), file = file2)
            print("<=> Valkyrie: Considering electron enthalpy at {} K".format(fd))
        shutil.move("INCAR-tmp", "INCAR")
    # Spin
    if spin == True:
        with open("INCAR", "a") as file:
            print("\nISPIN = 2", file = file)
    # Fermi-dirac
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
