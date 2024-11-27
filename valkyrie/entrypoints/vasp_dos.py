import os, shutil
from ..input import input_vasp_potcar, get_info, input_vasp_kpoints
from ..shell_scripts import job
from ase.io import read
from .. import __shell__, __python__, run_vasp
from . import vasp_incar

def gen_input(encut, spin, fun, fermiDirac, u, fElectron, **kwargs):
    incar_scf=f"""# INCAR for scf
ISMEAR = 0
EDIFF = 1E-6
EDIFFG = -0.001
SIGMA = 0.05
IBRION = -1
NSW = 0
ISIF = 2
ENCUT = {encut}
PREC = Accurate
KSPACING = 0.157
NPAR = 2
KPAR = 2
NELM = 120
"""
    with open("INCAR_scf", "w") as file:
        file.write(incar_scf)

    incar_dos = """# INCAR for dos
EDIFF = 1E-7
EDIFFG = -0.001
ISTART = 1
ICHARG = 11
ISMEAR = -5
SIGMA = 0.05
IBRION = -1
PREC = Accurate
ENCUT = {}
ISIF = 2
LORBIT = 11

KSPACING = 0.157
NEDOS  =  2001
LWAVE  = .FALSE.
LCHARG = .FALSE.
NPAR = 2
KPAR = 2""".format(encut)
    with open("INCAR_dos", "w") as file:
        file.write(incar_dos)
    
    if fun == "hse":
        fun1="gga"
        fun2="hse"
    else:
        fun1=fun
        fun2=fun
    vasp_incar.modify_INCAR(file        =  ["INCAR_scf"],
                            encut       =  encut, 
                            spin        =  spin, 
                            fermiDirac  =  fermiDirac,
                            fun         =  fun1, 
                            u           =  u,
                            fElectron   =  fElectron
                            ) 
    vasp_incar.modify_INCAR(file        =  ["INCAR_dos"],
                            encut       =  encut, 
                            spin        =  spin, 
                            fermiDirac  =  fermiDirac,
                            fun         =  fun2, 
                            u           =  u,
                            fElectron   =  fElectron
                            ) 
    return 0

def main(*args, queue="9242opa!", nodes=24, comment="band", notSub=False,
         pot="auto", encut=0,
         spin=False, fermiDirac=-1,
         fun="gga", u="None", fElectron=False,
         **kwargs):

    try:
        poscar=read("POSCAR").get_chemical_formula()
    except:
        raise Exception("No POSCAR file found at current path.")
    
    input_vasp_potcar.potcar(pot) # POTCAR
    encut=get_info.get_encut(encut) # ENCUT
    gen_input(**locals()) # INCAR


    # job
    job.gen_job(job="job", job_file="{}/job_dos".format(__shell__), task=run_vasp)
    shutil.copy2("{}/vasp_dos_sum_spd.py".format(__python__), "vasp_dos_sum_spd.py")
    job.control_job("job", queue, nodes, comment)
    print("<=> Valkyrie: Only generate input file.") if notSub else job.sub("job")

    print("<=> Valkyrie: Running dos for {}, ENCUT = {}, POTCAR = {}, fun = {}.".format(poscar, encut, pot, fun))
    os.system("rm -f Brillouin_Zone_3D.jpg INCAR KPOINTS PLOT.in SYMMETRY PRIMCELL.vasp")
    return 0

if __name__ == "__main__":
    main()