# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 13:08:07 2023

@author: zyj
"""

import os, shutil, sys, math
from Input import input_qe_pot, input_vasp_relax
from Input import get_info, poscar2qe, input_qe_relax
from Calculators import job
from ase.io import read
import numpy as np


def relax(args, __shell__, __python__, __work__):
    ##### Start Parameters #####
    pressure = args.pressure
    spin = args.spin
    not_sub = args.not_sub
    fd = args.fermi_dirac
    q = args.q
    n = 24 if args.n is None else args.n
    comment = "relax" if args.comment is None else args.comment
    symmetry = args.symmetry
    optcell = args.optcell
    ##### End Parameters #####


    # POSCAR
    poscar = read("POSCAR")
    formula = poscar.get_chemical_formula()
    if not os.path.exists("POSCAR"):
        raise Exception("No POSCAR file found at current path.")
    if symmetry != None:
        os.system("phonopy --symmetry --tolerance={} | grep space | head -1".format(symmetry))
        shutil.copy2("PPOSCAR", "POSCAR")
    ntyp, nat, cell, pos = poscar2qe.read_poscar("POSCAR")

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
    # ENCUT
    encut = args.encut
    # INCAR
    input_qe_relax.relax(pressure * 10, encut, k, nat, ntyp, cell, pos, pot, optcell)
    

    # Job and Sub
    job.gen_job("job", "{}/job_qe_relax".format(__shell__))
    job.control_job("job", q, n, comment)
    if not_sub == True:
        print("<=> Valkyrie: Only generate input file.")
    else:
        job.sub("job")

    # Output
    if optcell == True:
        print("<=> Valkyrie: Relax for {} with cell_dofree, ENCUT = {}, POTCAR = {}."\
            .format(formula, encut, pot.split()[-1]))
        if pressure != 0.0:
            print("\033[0;31;40m", "\b<=> WARNING: Pressure and cell_dofree tags can not use together!", "\033[0m")
            print("\033[0;31;40m", "\b<=> WARNING: The Pressure tag is ignored!", "\033[0m")
    else:
        print("<=> Valkyrie: Relax for {} under the pressure of {} GPa, ENCUT = {}, POTCAR = {}."\
                .format(formula, pressure, encut, pot.split()[-1]))
        
    return 0

