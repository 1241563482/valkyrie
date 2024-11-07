from ase.dft import kpoints
from ase.io import read
import shutil
from Calculators import job

def relax(pressure: int,
          encut: int,
          k: list,
          nat: int,
          ntyp: int,
          cell: str,
          pos: str,
          pot: str,
          optcell: bool):
    

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
