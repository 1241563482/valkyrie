import os
import numpy as np
from ..input import get_info
from ..shell_scripts import job
from .. import __shell__, run_pw, qe_pot_path
from ase.io import read
from ..input import poscar2qe, input_qe_pot

def gen_INPUT(encut, pressure, spin, k, nat, ntyp, cell, pos, pot, optcell, **kwargs):
    if optcell == False:
        CELL = f"""&CELL
    press = {pressure}  !kbar
    press_conv_thr = 0.02
    cell_dynamics = 'bfgs'
/
"""
    else:
        CELL = f"""&CELL
    cell_dofree = 2Dxy
    cell_dynamics = 'bfgs'
/
"""
    
    relax_in = f"""&CONTROL
    calculation = 'vc-relax' , prefix = 'pwscf'
    pseudo_dir = '/fsa/home/js_zhuyj/mypps/QE/SSSP_1.3.0_PBE_efficiency'
    outdir = './tmp'
    forc_conv_thr = 1.0d-7
    etot_conv_thr = 1.0d-7
/
&SYSTEM
    ibrav= 0, nat= {nat}, ntyp= {ntyp},  
    occupations = 'smearing', smearing = 'mp', degauss = 0.02
    ecutwfc = {encut}
/
&ELECTRONS
    electron_maxstep = 100
    conv_thr = 1.0d-11
!    diagonalization = 'cg'
/
&IONS
    ion_dynamics='bfgs'
/
{CELL}
ATOMIC_SPECIES
{pot}
K_POINTS automatic
{k[0]} {k[1]} {k[2]} 0 0 0

{cell}
{pos}
"""

    with open("relax.in", "w") as file:
        print(relax_in, file = file)
    print("<=> Valkyrie: Kmesh for relax.in: {} {} {}".format(k[0], k[1], k[2]))
    return 0


def main(*args, queue = "9242opa!", nodes = 24, comment = "qe_relax", notSub = False,
         pressure = 0, pot = "auto", encut = 100,
         spin = False,
         optcell = False,
         fermiDirac = -1,
         fun = "gga", u = "None", fElectron = False,
         **kwargs):
    
    try:
        ntyp, nat, cell, pos = poscar2qe.read_poscar("POSCAR")
        poscar = read("POSCAR")
        formula = poscar.get_chemical_formula()
    except:
        raise Exception("No POSCAR file found at current path.")
    

    # K points
    atom_size = [ max(poscar.get_positions()[:,i]) - min(poscar.get_positions()[:,i]) for i in [0, 1, 2] ]
    cell_size = [ np.linalg.norm(poscar.cell[i]) for i in [0, 1, 2] ]
    k = []
    for i in [0, 1, 2]:
        if cell_size[i] - atom_size[i] >= 8: # If vacuum size >= 8A, then set the K-mesh as 1.
            k.append(1)
        else:
            k.append( int( max(1, np.ceil(40 / cell_size[i]))) )

    # POTCAR
    pot = input_qe_pot.get_pot("POSCAR")
    # INCAR
    gen_INPUT(**locals())
    

    # Job and Sub
    job.gen_job(job="job", job_file="{}/job_qe_relax".format(__shell__), task=run_pw)
    job.control_job("job", queue, nodes, comment)
    print("<=> Valkyrie: Only generate input file.") if notSub else job.sub("job")
    print(f"<=> Valkyrie: QE relax for {formula} under the pressure of {pressure} GPa, ENCUT={encut}Ry.")
    print(f"<=> Valkyrie: Potential: {pot.split()[-1]}")

    # Output
    if optcell == True:
        print("<=> Valkyrie: QE relax for {} with cell_dofree, ENCUT={}Ry, Potential={}.".format(poscar, encut, pot))
        if pressure != 0.0:
            print("\033[0;31;40m", "\b<=> WARNING: Pressure and cell_dofree tags can not use together!", "\033[0m")
            print("\033[0;31;40m", "\b<=> WARNING: The Pressure tag is ignored!", "\033[0m")
        
    return 0

if __name__ == "__main__":
    main()