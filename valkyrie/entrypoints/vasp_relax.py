import os, shutil
from ..input import input_vasp_potcar, get_info, input_vasp_kpoints
from ..shell_scripts import job
from .. import __shell__, run_vasp, run_vasp_opt
from ase.io import read

def gen_INCAR(encut, pressure, spin, fd, u_atom, u_value, lmaxmix, optcell):
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
    if fd != False:
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
    return 0



def main(*args, pressure = 0, pot = "auto", spin = False, not_sub = False, fermiDirac = False,
         queue = "9242opa!", nodes = 24, comment = "relax", symmetry = False, optcell = False,
         fun = "gga", u = "None", encut = 600, fElectron = False,
         **kwargs):

    # Fun and gga+U part
    if fun == "ggau":
        if args.u is not None:
            u_atom = args.u[0]
            u_value = args.u[1]
            if args.f_electron == True:
                lmaxmix = 6 
            else:
                lmaxmix = 4
        else:
            raise Exception("Missing -u tag from input")
    else:
        u_atom = None
        u_value = None
        lmaxmix = None
         
    ##### End Parameters #####


    # POSCAR
    poscar = read("POSCAR")
    chemical_formula = poscar.get_chemical_formula()
    if not os.path.exists("POSCAR"):
        raise Exception("No POSCAR file found at current path.")
    #if symmetry != None:
    #    os.system("phonopy --symmetry --tolerance={} | grep space | head -1".format(symmetry))
    #    shutil.copy2("PPOSCAR", "POSCAR")
    

    input_vasp_potcar.potcar(pot) # POTCAR
    encut = get_info.get_encut(encut) # ENCUT
    gen_INCAR(encut, pressure, spin, fermiDirac, u_atom, u_value, lmaxmix, optcell) # INCAR

    # KPOINTS
    if optcell == True:
        input_vasp_kpoints.kpoints_byhand(poscar)
        os.system("sed -i 's/KSPACING/#KSPACING/g' INCAR")

    # Job and Sub
    if optcell == False:
        job.gen_job(job = "job", job_file = "{}/job_vasp_relax".format(__shell__), task = run_vasp)
    else:
        job.gen_job(job = "job", job_file = "{}/job_vasp_relax_optcell".format(__shell__), task = run_vasp_opt)
    job.control_job("job", queue, nodes, comment)
    print("<=> Valkyrie: Only generate input file.") if not_sub == True else job.sub("job")
        

    # Output
    if optcell == True:
        print("<=> Valkyrie: Relax for {} with vasprelax version, ENCUT = {}, POTCAR = {}."\
            .format(chemical_formula, encut, pot))
        if pressure != 0.0:
            print("\033[0;31;40m", "\b<=> WARNING: Pressure and optcell tags can not use together!", "\033[0m")
            print("\033[0;31;40m", "\b<=> WARNING: The Pressure tag is ignored!", "\033[0m")
    else:
        print("<=> Valkyrie: Relax for {} under the pressure of {} GPa, ENCUT = {}, POTCAR = {}."\
            .format(chemical_formula, pressure, encut, pot))

    return 0



if __name__ == "__main__":
    main()