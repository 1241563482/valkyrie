# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 20:42:21 2023

@author: zyj
"""

import os, platform, shutil, sys
from Input import input_vasp_potcar, input_vasp_ph
from Input import input_vasp_kpoints
from Input import get_info
from Calculators import job
from ase.io import read



def ph(args, __shell__, __python__, __work__):
    ##### Start Parameters #####
    pot = args.pot
    spin = args.spin
    not_sub = args.not_sub
    fd = args.fermi_dirac
    kpoints = args.kpoints
    encut = args.encut
    dim = args.dim
    symmetry = args.symmetry
    q = args.q
    n = args.n
    comment = args.comment
    # Fun and gga+U part
    fun = args.fun
    if fun == "ggau":
        if args.u is not None:
            u_atom = args.u[0]
            u_value = args.u[1]
            if float(u_value) <= 0:
                raise Exception("Ueff must be positive.")
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
        os.system("phonopy -d --dim {} {} {} --tolerance {} > /dev/null".format(*dim, symmetry))
    else:
        symmetry = 1e-5
        os.system("phonopy -d --dim {} {} {} --tolerance {} > /dev/null".format(*dim, symmetry))
    # POTCAR
    input_vasp_potcar.potcar(pot)
    # ENCUT
    encut = get_info.get_encut(args.encut)
    # INCAR
    input_vasp_ph.ph(encut, dim, spin, fd, symmetry, fun, u_atom, u_value, lmaxmix)
    # KPOINTS
    input_vasp_kpoints.kpoints(*kpoints)
    # Other input files
    shutil.copy2("{}/get_free_energy.py".format(__python__), "get_free_energy.py")
    shutil.copy2("{}/vasp_ph_loop_in.py".format(__python__), "loop-in.py")
    shutil.copy2("{}/../Calculators/job.py".format(__python__), "job.py")
    shutil.copy2("{}/../set_up.py".format(__python__), "set_up.py")
    
    # Job and Sub
    job.gen_job("job", "{}/job_ph".format(__shell__))
    job.control_job("job", q, n, comment)
    if not_sub == True:
        print("<=> Valkyrie: Only generate input file.")
    else:
        os.system("python loop-in.py")
    
    print("<=> Valkyrie: Ph for {}, ENCUT = {}, POTCAR = {}.".format(poscar, encut, pot))  
    print("<=> Valkyrie: Super cell is {} {} {}. KPOINTS are {} {} {}.".format(*dim, *kpoints))


    return 0


