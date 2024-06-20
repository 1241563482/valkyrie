import shutil
from Input import get_info

def elf(encut, pressure, spin, fun, u_atom, u_value, lmaxmix):
    pressure = pressure * 10
    incar = """# INCAR for elf
ISMEAR = 0
EDIFF = 1E-7
EDIFFG = -0.001
SIGMA = 0.05
IBRION = -1
NSW = 0
ISIF = 2
ENCUT = {}
PSTRESS = {}
PREC = Accurate
KSPACING = 0.157
#NPAR = 2
KPAR = 2
NELM = 120
LELF = .True.
LCHARG = .False.
LWAVE = .False.
""".format(encut, pressure)
    with open("INCAR", "w") as file:
        file.write(incar)
    # Spin
    if spin == True:
        with open("INCAR", "a") as file:
            print("\nISPIN = 2", file = file)
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
    
        
