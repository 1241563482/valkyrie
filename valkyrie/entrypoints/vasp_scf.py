import os
from ..input import input_vasp_potcar, get_info
from ..shell_scripts import job
from .. import __shell__, run_vasp
from ase.io import read
from . import vasp_incar


def gen_INCAR(encut, pressure, spin, fermiDirac, fun, u, fElectron, **kwargs):
    pressure = pressure * 10
    incar = f"""# INCAR for scf
ISMEAR = 0
EDIFF = 1E-6
EDIFFG = -0.001
SIGMA = 0.05
IBRION = -1
NSW = 0
ISIF = 2
ENCUT = {encut}
PSTRESS = {pressure}
PREC = Accurate
KSPACING = 0.157
NPAR = 2
KPAR = 2
NELM = 120
"""
    with open("INCAR", "w") as file:
        file.write(incar)
    
    vasp_incar.modify_INCAR(file       =   ["INCAR"],
                            spin       =   spin, 
                            fermiDirac =   fermiDirac,
                            fun        =   fun, 
                            u          =   u,
                            fElectron  =   fElectron
                            )
    return 0



def main(*args, queue = "9242opa!", nodes = 24, comment = "band", notSub = False,
         pressure = 0, pot = "auto", encut = 0,
         spin = False,
         fermiDirac = -1,
         fun = "gga", u = "None", fElectron = False,
         **kwargs):
    try:
        poscar = read("POSCAR").get_chemical_formula()
    except:
        raise Exception("No POSCAR file found at current path.")
    input_vasp_potcar.potcar(pot) # POTCAR
    encut = get_info.get_encut(encut) # ENCUT
    gen_INCAR(**locals()) # INCAR


    # Job and sub
    job.gen_job(job = "job", job_file = "{}/job_scf".format(__shell__), task = run_vasp)
    job.control_job("job", queue, nodes, comment)
    print("<=> Valkyrie: Only generate input file.") if notSub else job.sub("job")
    print("<=> Valkyrie: Scf for {} under the pressure of {} GPa, ENCUT = {}, POTCAR = {}."\
        .format(poscar, pressure, encut, pot))

    return 0


if __name__ == "__main__":
    main()