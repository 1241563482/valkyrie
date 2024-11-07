from ase.dft import kpoints
from ase.io import read
import shutil
from Calculators import job

def scf(pressure: int,
          encut: int,
          k: list,
          nat: int,
          ntyp: int,
          cell: str,
          pos: str,
          pot: str
          ):
    
    scf_in = f"""&CONTROL
    calculation = 'scf' , prefix = 'pwscf'
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
&CELL
    press = {pressure}  !kbar
    press_conv_thr = 0.02
    cell_dynamics = 'bfgs'

ATOMIC_SPECIES
{pot}
K_POINTS automatic
{k[0]} {k[1]} {k[2]} 0 0 0

{cell}
{pos}
"""

    with open("scf.in", "w") as file:
        print(scf_in, file = file)
    print("<=> Valkyrie: Kmesh for scf.in: {} {} {}".format(k[0], k[1], k[2]))
