from ..input import input_vasp_potcar, get_info, input_vasp_kpoints
from ..shell_scripts import job
from .. import run_vasp_gam, run_vasp
from ..io_utils import run_command
from ..resources import copy_python_script, shell_template
from ase.io import read
from . import vasp_incar
import subprocess
from pathlib import Path


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

#NBLOCK = 5
#KBLOCK = 4
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
    

    if not Path("POSCAR").exists():
        raise Exception("No POSCAR file found at current path.")
    run_command(["phonopy", "-d", "--dim", *map(str, dim)], stdout=subprocess.DEVNULL)
    Path("POSCAR").rename("POSCAR0")
    Path("SPOSCAR").rename("POSCAR")
    poscar = read("POSCAR")
    formula = poscar.get_chemical_formula()


    input_vasp_potcar.potcar(pot) # POTCAR
    encut = get_info.get_encut(encut) # ENCUT
    [tBegin, tEnd] = temperature
    gen_INCAR(**locals()) # INCAR
    input_vasp_kpoints.kpoints(*kpoints)
    

    # Other input files
    copy_python_script("vasp_md_MDdata.py", "MDdata.py")
    copy_python_script("vasp_md_reprocess.py", "reprocess.py")

    # Job and sub
    task = run_vasp_gam if kpoints == [1, 1, 1] else run_vasp
    job.gen_job(job = "job", job_file = shell_template("job_vasp_md"), task = task)
    job.control_job("job", queue, nodes, comment)

    print("<=> Valkyrie: Only generate input file.") if notSub else job.sub("job")
    print(f"<=> Valkyrie: MD for {formula} under the pressure of {pressure} GPa, ENCUT = {encut}, POTCAR = {pot}")
    
    print("<=> Valkyrie: Super cell is {} {} {}. KPOINTS are {} {} {}.".format(*dim, *kpoints))
    print(f"<=> Valkyrie: There are {len(poscar)} atoms in the cell.")
    print(f"<=> Valkyrie: TEBEG = {tBegin}, TEEND = {tEnd}.")

    return 0


if __name__ == "__main__":
    main()
