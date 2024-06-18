# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 13:08:07 2023

@author: zyj
"""

import os, shutil, sys
from Input import input_vasp_potcar, input_vasp_relax
from Input import get_info
from Calculators import job
from ase.io import read


def relax(args, __shell__, __python__, __work__):
    ##### Start Parameters #####
    pressure = args.pressure
    pot = args.pot
    spin = args.spin
    not_sub = args.not_sub
    fd = args.fermi_dirac
    q = args.q
    n = args.n
    comment = args.comment
    symmetry = args.symmetry

    # Fun and gga+U part
    fun = args.fun
    if fun == "ggau":
        if args.u is not None:
            u_atom = args.u[0]
            u_value = args.u[1]
            if args.f_electron == True:
                lmaxmix = 6 
            else:
                lmaxmix = 4
        else:
            raise Exception("Missing -u tag from input")
    else:
        u_atom = None
        u_value = None
        lmaxmix = None
         
    ##### End Parameters #####


    # POSCAR
    poscar = read("POSCAR").get_chemical_formula()
    if not os.path.exists("POSCAR"):
        raise Exception("No POSCAR file found at current path.")
    if symmetry != None:
        os.system("phonopy --symmetry --tolerance={} | grep space | head -1".format(symmetry))
        shutil.copy2("PPOSCAR", "POSCAR")
    # POTCAR
    input_vasp_potcar.potcar(pot)
    # ENCUT
    encut = get_info.get_encut(args.encut)
    # INCAR
    input_vasp_relax.relax(encut, pressure, spin, fd, u_atom, u_value, lmaxmix)

    # Job and Sub
    job.gen_job("job", "{}/job_relax".format(__shell__))
    job.control_job("job", q, n, comment)
    if not_sub == True:
        print("<=> Valkyrie: Only generate input file.")
    else:
        job.sub("job")

    
    print("<=> Valkyrie: Relax for {} under the pressure of {} GPa, ENCUT = {}, POTCAR = {}."\
        .format(poscar, pressure, encut, pot))

    return 0

