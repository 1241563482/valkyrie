import os
from ..input import input_vasp_potcar, get_info, input_vasp_kpoints
from ..shell_scripts import job
from .. import __shell__, run_vasp, run_vasp_opt
from ase.io import read
from . import vasp_incar

def gen_input(encut, pressure, spin, fermiDirac, fun, u, fElectron, optcell, **kwargs):
    pressure = pressure * 10
    incar = f"""# INCAR for relax
PREC = Accurate
EDIFF = 1e-6
EDIFFG = -1e-3
ENCUT = {encut}
PSTRESS = {pressure}
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
"""
    with open("INCAR", "w") as file:
        file.write(incar)

    if optcell == True:
        OPTCELL = """100
110
000"""
        with open("OPTCELL", "w") as file:
            print(OPTCELL, file = file)
        with open("INCAR", "a") as file:
            print("IVDW = 11", file = file)
        
    vasp_incar.modify_INCAR(file       =   ["INCAR"],
                            encut      =   encut, 
                            pressure   =   pressure, 
                            spin       =   spin, 
                            fermiDirac =   fermiDirac,
                            fun        =   fun, 
                            u          =   u,
                            fElectron  =   fElectron
                            )
    return 0


def main(*args, pressure = 0, pot = "auto", spin = False, notSub = False, fermiDirac = -1,
         queue = "9242opa!", nodes = 24, comment = "relax", symmetry = False, optcell = False,
         fun = "gga", u = "None", encut = 0, fElectron = False,
         **kwargs):

    try:
        poscar = read("POSCAR").get_chemical_formula()
    except:
        raise Exception("No POSCAR file found at current path.")
    
    input_vasp_potcar.potcar(pot) # POTCAR
    encut = get_info.get_encut(encut) # ENCUT
    gen_input(**locals()) # INCAR


    # KPOINTS and job
    if optcell:
        input_vasp_kpoints.kpoints_byhand(read("POSCAR"), ceiling = 50)
        os.system("sed -i 's/KSPACING/#KSPACING/g' INCAR")
        job.gen_job(job = "job", job_file = "{}/job_vasp_relax_optcell".format(__shell__), task = run_vasp_opt)
        job.control_job("job", queue, nodes, comment)
        print("<=> Valkyrie: Only generate input file.") if notSub else job.sub("job")
        print("<=> Valkyrie: Relax for {} with vasprelax version, ENCUT = {}, POTCAR = {}."\
            .format(poscar, encut, pot))
        if pressure != 0.0:
            print("\033[0;31;40m", "\b<=> WARNING: Pressure and optcell tags can not use together!", "\033[0m")
            print("\033[0;31;40m", "\b<=> WARNING: The Pressure tag is ignored!", "\033[0m")
    else:
        job.gen_job(job = "job", job_file = "{}/job_vasp_relax".format(__shell__), task = run_vasp)
        job.control_job("job", queue, nodes, comment)
        print("<=> Valkyrie: Only generate input file.") if notSub else job.sub("job")
        print("<=> Valkyrie: Relax for {} under the pressure of {} GPa, ENCUT = {}, POTCAR = {}."\
            .format(poscar, pressure, encut, pot))
 
    return 0



if __name__ == "__main__":
    main()