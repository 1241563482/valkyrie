from ..input import input_vasp_potcar, get_info, input_vasp_kpoints
from ..shell_scripts import job
from .. import __shell__, __python__, run_vasp_gam, run_vasp
from ase.io import read
from . import vasp_incar
import os, shutil


def gen_INCAR(encut, pressure, spin, fermiDirac, fun, u, fElectron, tBegin, tEnd, **kwargs):
    pressure = pressure * 10
    incar = f"""SYSTEM = AIMD
LWAVE = .FALSE.
LCHARG = .FALSE.
ISTART = 0
ICHARG = 2
ISYM = 0
PREC = Normal
ALGO = N
ENCUT = {encut}
ISMEAR = 0
SIGMA = 0.05
EDIFF = 1e-5
NELMIN = 8
MDALGO = 2
ISIF = 2
IBRION = 0
#LREAL = A

NSW = 10000
POTIM = 1
SMASS = 2
NELM = 150
PSTRESS = {pressure}
TEBEG = {tBegin}
TEEND = {tEnd}

NBLOCK = 5
KBLOCK = 4
NPAR = 8
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



def main(*args, queue = "9242opa!", nodes = 48, comment = "vasp_MD", notSub = False,
         pressure = 0, pot = "auto", encut = 0,
         temperature = [300, 300],
         spin = False,
         fermiDirac = -1,
         fun = "gga", u = "None", fElectron = False,
         kpoints = [1, 1, 1], dim = [2, 2, 2],
         **kwargs):
    

    try:
        os.system("phonopy -d --dim {} {} {} > /dev/null".format(*dim))
        os.system("mv POSCAR POSCAR0")
        os.system("mv SPOSCAR POSCAR")
        poscar = read("POSCAR")
        formula = poscar.get_chemical_formula()
    except:
        raise Exception("No POSCAR file found at current path.")


    input_vasp_potcar.potcar(pot) # POTCAR
    encut = get_info.get_encut(encut) # ENCUT
    [tBegin, tEnd] = temperature
    gen_INCAR(**locals()) # INCAR
    input_vasp_kpoints.kpoints(*kpoints)
    

    # Other input files
    shutil.copy2("{}/vasp_md_MDdata.py".format(__python__), "MDdata.py")
    shutil.copy2("{}/vasp_md_reprocess.py".format(__python__), "reprocess.py")

    # Job and sub
    task = run_vasp_gam if kpoints == [1, 1, 1] else run_vasp
    job.gen_job(job = "job", job_file = "{}/job_vasp_md".format(__shell__), task = task)
    job.control_job("job", queue, nodes, comment)

    print("<=> Valkyrie: Only generate input file.") if notSub else job.sub("job")
    print(f"<=> Valkyrie: MD for {formula} under the pressure of {pressure} GPa, ENCUT = {encut}, POTCAR = {pot}")
    
    print("<=> Valkyrie: Super cell is {} {} {}. KPOINTS are {} {} {}.".format(*dim, *kpoints))
    print(f"<=> Valkyrie: There are {len(poscar)} atoms in the cell.")
    print(f"<=> Valkyrie: TEBEG = {tBegin}, TEEND = {tEnd}.")

    return 0


if __name__ == "__main__":
    main()