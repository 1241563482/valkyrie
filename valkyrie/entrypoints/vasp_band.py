import os, shutil
from ..input import input_vasp_potcar, get_info, input_vasp_kpoints
from ..shell_scripts import job
from .. import __shell__, run_vasp, run_vasp_opt
from ase.io import read
from . import vasp_incar


def gen_input(encut, pressure, spin, fermiDirac, fun, u, fElectron, optcell, **kwargs):
    pressure = pressure * 10

    incar_scf = f"""# INCAR for scf
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
    with open("INCAR_scf", "w") as file:
        file.write(incar_scf)

    incar_band = """# INCAR for band
EDIFF = 1E-6
EDIFFG = -0.005
ISTART = 1
ICHARG = 11
ISMEAR = 0
SIGMA = 0.03
IBRION = -1
PREC = Accurate
ENCUT = {}
ISIF = 2
LORBIT = 11
LWAVE  = .FALSE.
LCHARG = .FALSE.
NPAR = 2
KPAR = 2""".format(encut)
    with open("INCAR_band", "w") as file:
        file.write(incar_band)

    vasp_incar.modify_INCAR(file       =   ["INCAR_scf", "INCAR_band"],
                            encut      =   encut, 
                            pressure   =   pressure, 
                            spin       =   spin, 
                            fermiDirac =   fermiDirac,
                            fun        =   fun, 
                            u          =   u,
                            fElectron  =   fElectron
                            ) 
    return 0



def main(*args, queue = "9242opa!", nodes = 24, comment = "band", notSub = False,
         pot = "auto", encut = 0,
         spin = False,
         fun = "gga", u = "None", fElectron = False,
         **kwargs):


    try:
        poscar = read("POSCAR").get_chemical_formula()
    except:
        raise Exception("No POSCAR file found at current path.")
    
    input_vasp_potcar.potcar(pot) # POTCAR
    encut = get_info.get_encut(encut) # ENCUT
    gen_input(**locals()) # INCAR


    
    # KPOINTS and job
    if fun == "gga" or fun == "ggau":
        os.system("vaspkit -task 303 > /dev/null 2>&1")
        job.gen_job(job = "job", job_file = "{}/job_band".format(__shell__), task = run_vasp)
    
    elif fun == "hse":
        if not os.path.isfile("KPATH.in"):
            os.system("vaspkit -task 303 > /dev/null 2>&1") # KPOINTS for gga band
        os.system("echo -e '251\n2\n0.04\n0.06' | vaspkit > /dev/null") # KPOINTS for hse band
        os.rename("KPOINTS", "KPOINTS-band")
        with open("KPOINTS-band", "r") as file:
            numbers = file.readline().strip().split()

        input_vasp_kpoints.kpoints(numbers[1], numbers[2], numbers[3]) # KPOINTS for scf
        os.rename("KPOINTS", "KPOINTS-scf")
        os.system("sed -i 's/KSPACING/#KSPACING/g' INCAR-scf")
        job.gen_job(job = "job", job_file = "{}/job_band_hse".format(__shell__), task = run_vasp)
    

    # Job and Sub
    job.control_job("job", queue, nodes, comment)
    print("<=> Valkyrie: Only generate input file.") if notSub else job.sub("job")
    print("<=> Valkyrie: Running band for {}, ENCUT = {}, POTCAR = {}, fun = {}."\
          .format(poscar, encut, pot, fun))
    os.system("rm -f Brillouin_Zone_3D.jpg INCAR PLOT.in SYMMETRY PRIMCELL.vasp")
    return 0


if __name__ == "__main__":
    main()