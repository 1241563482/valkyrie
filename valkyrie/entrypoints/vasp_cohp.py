import os, shutil
from ..input import input_vasp_potcar, get_info, input_vasp_kpoints
from ..shell_scripts import job
from ase.io import read
from .. import __shell__, __python__, run_vasp
from . import vasp_incar


def gen_INCAR(encut, spin, fun, u, fElectron, **kwargs):
    INCAR1 = f"""# COHP step1
EDIFF = 1E-6
EDIFFG = -0.01
SIGMA = 0.05
IBRION = -1
NSW = 0
ISIF = 2
ENCUT = {encut}
PREC = Accurate
ISMEAR = -5
NEDOS = 2001
ISYM = -1
LWAVE = .F.
LCHARG = .F.
NELM = 120"""
    
    INCAR2 = f"""# COHP step2
EDIFF = 1E-6
EDIFFG = -0.01
SIGMA = 0.05
IBRION = -1
NSW = 0
ISIF = 2
ENCUT = {encut}
PREC = Accurate
ISMEAR = -5
NEDOS = 2001
ISYM = -1
LCHARG = .F.
NELM = 120"""
    
    with open("INCAR1", "w") as file1, open("INCAR2", "w") as file2:
        file1.write(INCAR1)
        file2.write(INCAR2)
    
    vasp_incar.modify_INCAR(file       =   ["INCAR1", "INCAR2"],
                            spin       =   spin, 
                            fun        =   fun, 
                            u          =   u,
                            fElectron  =   fElectron
                            )

    return 0

def main(*args, queue = "9242opa!", nodes = 24, comment = "cohp", notSub = False,
         pot = "auto", encut = 0,
         spin = False,
         fun = "gga", u = "None", fElectron = False,
         **kwargs):
    
    # POSCAR
    poscar = read("POSCAR")
    chemical_formula = poscar.get_chemical_formula()
    if not os.path.exists("POSCAR"):
        raise Exception("No POSCAR file found at current path.")
    input_vasp_potcar.potcar(pot) # POTCAR
    encut = get_info.get_encut(encut) # ENCUT
    gen_INCAR(encut, spin, fun, u, fElectron) # INCAR
    input_vasp_kpoints.kpoints_byhand(poscar, ceiling = 50) # KPOINTS


    # Job and Sub
    job.gen_job(job = "job", job_file = "{}/job_cohp".format(__shell__), task = run_vasp)
    job.control_job("job", queue, nodes, comment)
    print(f"<=> Valkyrie: COHP for {chemical_formula} , ENCUT = {encut}, POTCAR = {pot}.")
    print("<=> Valkyrie: Only generate input file.") if notSub else job.sub("job")


    # Scripts
    shutil.copy(f"{__python__}/vasp_cohp.py", "COHP.py")
    shutil.copy(f"{__python__}/vasp_cohp_plot.py", "plot.py")

    return 0




if __name__ == "__main__":
    main()
