import os, shutil
from ..input import input_vasp_potcar, get_info, input_vasp_kpoints
from ..shell_scripts import job
from .. import __shell__, run_vasp, __python__
from ase.io import read
from . import vasp_incar

def gen_input(encut, dim, spin, fermiDirac, fun, u, fElectron, symmetry, **kwargs):
    incar = f"""# INCAR for ph
ISTART = 0        
#LREAL = Auto    
ENCUT = {encut}   
PREC = Accurate     
LWAVE = .FALSE.  
LCHARG = .FALSE.    
ADDGRID = .TRUE. 
ISMEAR = 0          
SIGMA = 0.06   
NELM = 150         
EDIFF = 1E-07
EDIFFG = -1E-03
ALGO = N
NPAR = 4
SYMPREC = 1E-12"""
    with open("INCAR", "w") as file:
        file.write(incar)
    
    vasp_incar.modify_INCAR(file        =   ["INCAR"],
                            spin        =   spin, 
                            fermiDirac  =   fermiDirac,
                            fun         =   fun, 
                            u           =   u,
                            fElectron   =   fElectron
                            )
    
    # band.conf
    band_conf = f"""
ATOM_NAME = Bilibili@yijiezhu
DIM = {dim[0]} {dim[1]} {dim[2]}
PRIMITIVE_AXES = Auto
BAND = Auto
FORCE_SETS = READ
FORCE_CONSTANTS= WRITE"""
    with open("band.conf", "w") as file:
        file.write(band_conf)

    # mesh.conf
    mesh_conf = f"""
ATOM_NAME = Bilibili@yijiezhu
DIM = {dim[0]} {dim[1]} {dim[2]}
#FORCE_CONSTANTS = READ  # For DFPT
MP = 15 15 15
TPROP = T
TMIN = 0
TMAX = 10000
TSTEP = 10"""
    with open("mesh.conf", "w") as file:
        file.write(mesh_conf)

    # loop-out.sh
    loop_out = f"""
#! /bin/bash
mkdir -p output
for i in `ls POSCAR-* | sed 's/POSCAR-//g'`
do
cp $i/vasprun.xml output/"vasprun.xml-"$i
done
phonopy -f output/*
phonopy -s -p band.conf --tolerance={symmetry}
phonopy-bandplot --gnuplot > ph.out
#phonopy -t mesh.conf --tolerance={symmetry}
#python get_free_energy.py
#phonopy --dim="2 2 2" --irreps="0 0 0 1e-3"
"""
    with open("loop-out.sh", "w") as file:
        file.write(loop_out)

    return 0



def main(*args, queue = "9242opa!", nodes = 24, comment = "ph", notSub = False,
         pressure = 0, pot = "auto", encut = 0,
         spin = False,
         fermiDirac = -1,
         fun = "gga", u = "None",  fElectron = False, 
         kpoints = [3, 3, 3], dim = [2, 2, 2], symmetry = 0.1,
         **kwargs):
    

    # POSCAR
    try:
        poscar = read("POSCAR").get_chemical_formula()
    except:
        raise Exception("No POSCAR file found at current path.")
    
    # Super cell
    os.system("phonopy -d --dim {} {} {} --tolerance {} > /dev/null".format(*dim, symmetry))

    input_vasp_potcar.potcar(pot) # POTCAR
    encut = get_info.get_encut(encut) # ENCUT
    gen_input(**locals()) # input

    input_vasp_kpoints.kpoints(*kpoints)


    # Other input files
    shutil.copy2("{}/get_free_energy.py".format(__python__), "get_free_energy.py")
    shutil.copy2("{}/vasp_ph_loop_in.py".format(__python__), "loop-in.py")
    

    # Job and Sub
    job.gen_job(job = "job", job_file = "{}/job_ph".format(__shell__), task = run_vasp)
    job.control_job("job", queue, nodes, comment)
    print("<=> Valkyrie: Only generate input file.") if notSub else os.system("python loop-in.py")

    print("<=> Valkyrie: Ph for {}, ENCUT = {}, POTCAR = {}.".format(poscar, encut, pot))  
    print("<=> Valkyrie: Super cell is {} {} {}. KPOINTS are {} {} {}.".format(*dim, *kpoints))


    return 0


if __name__ == "__main__":
    main()