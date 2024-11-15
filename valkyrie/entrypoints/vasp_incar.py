import shutil
from ..input import get_info

def modify_INCAR(file, spin, fermiDirac, fun, u, fElectron, **kwargs):
    for f in file:
        with open(f, "a") as file:
            if spin:
                file.write("\nISPIN = 2\nLORBIT = 11")
    
        if fermiDirac > 0:
            with open(f, 'r') as file1, open('INCAR-tmp', 'w') as file2:
                for line in file1:
                    if 'ISMEAR' not in line and 'SIGMA' not in line:
                        file2.write(line)
                kB = 8.617330337217213e-05
                print("\nISMEAR = -1\nSIGMA = {:.5f}".format(kB * fermiDirac), file=file2)
                print(f"<=> Valkyrie: Considering electron enthalpy at {fermiDirac} K")
            shutil.move("INCAR-tmp", f)



        if fun != "ggau" and u is not None:
            raise Exception("u only works for fun=ggau")
        
        if fun == "ggau" and u is not None:
            u_atom = u[0]
            u_value = float(u[1])
            if float(u_value) <= 0:
                raise Exception("Ueff must be positive.")
            lmaxmix = 6 if fElectron else 4

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
            with open(f, "a") as file:
                file.write(ggau_part)
            print(f"<=> Valkyrie: GGA+U for {u_atom}, Ueff = {u_value}.")

        elif fun == "hse":
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
            print(f"<=> Valkyrie: HSE06 calculation.")
        
    return 0