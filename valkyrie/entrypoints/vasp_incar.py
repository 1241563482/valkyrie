import shutil
from ..input import get_info


def addSpin(f, spin):
    with open(f, "a") as file:
        if spin:
            file.write("\nISPIN = 2\nLORBIT = 11")
            print(f"<=> Valkyrie: Add spin to {f}.")
    return 0


def addFermiDirac(f, fermiDirac):
    if fermiDirac > 0:
        with open(f, 'r') as file1, open('INCAR-tmp', 'w') as file2:
            for line in file1:
                if 'ISMEAR' not in line and 'SIGMA' not in line:
                    file2.write(line)
            kB = 8.617330337217213e-05
            print("\nISMEAR = -1\nSIGMA = {:.5f}".format(kB * fermiDirac), file=file2)
            print(f"<=> Valkyrie: Considering electron enthalpy at {fermiDirac} K ({f}).")
        shutil.move("INCAR-tmp", f)
    return 0


def ggau(f, u, fElectron):
    try:
        u_atom = u[0]
        u_value = float(u[1])
    except:
        raise Exception("\n\n<=> ERROR: Input u wrong, e.g. -u Li 3.")
    
    if float(u_value) <= 0:
        raise Exception("Ueff must be positive.")
    
    if fElectron:
        print("<=> Valkyrie: f electron is considered, LMAXMIX = 6 ({f}).")
        lmaxmix = 6
    else:
        lmaxmix = 4

    symbols = get_info.get_symbols("POSCAR")
    ldaul = ""
    ldauu = ""
    ldauj = ""
    for i in range(len(symbols)):
        if symbols[i] == u_atom:
            ldaul = ldaul + "2 "
            ldauu = ldauu + str(u_value ) + " "
            ldauj = ldauj + "0 "
        else:
            ldaul = ldaul + "-1 "
            ldauu = ldauu + "0 "
            ldauj = ldauj + "0 " 
    ggau_part = """
# GGA+U part
LASPH = .T.
LDAU = .T.
LDAUTYPE = 2
LDAUL = {} # 2 for +U, -1 for not +U
LDAUU = {} # Ueff = U-J
LDAUJ = {}
LMAXMIX = {} # If f electron, use 6
""".format(ldaul, ldauu, ldauj, lmaxmix)
    with open(f, "a") as file:
        file.write(ggau_part)
    print(f"<=> Valkyrie: GGA+U for {u_atom}, Ueff={u_value} ({f}).")
    return 0


def hse(f):
    hse_part = """
# HSE part
LHFCALC= .TRUE.
HFSCREEN = 0.2
ALGO = A
TIME = 0.4
AEXX = 0.25
"""
    with open(f, "a") as file:
        file.write(hse_part)
    print(f"<=> Valkyrie: HSE06 calculation ({f}).")
    return 0


def modify_INCAR(*args, **kwargs):
    files = kwargs.get("file", ["INCAR"])
    fun = kwargs.get("fun", "gga")
    u = kwargs.get("u", None)
    spin = kwargs.get("spin", False)
    fermiDirac = kwargs.get("fermiDirac", -1)
    fElectron = kwargs.get("fElectron", False)

    for f in files:
        addSpin(f, spin)
        addFermiDirac(f, fermiDirac)

        if fun!= "ggau" and u is not None:
            fun = "ggau"
            print("<=> WARNING: Get input u is not None, fun is changed to ggau.")
        elif fun == "ggau":
            ggau(f, u, fElectron)
        elif fun == "hse":
            hse(f)
        
    return 0
