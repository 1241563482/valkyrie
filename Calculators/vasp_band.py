# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 10:49:30 2023

@author: zyj
"""

import os, platform, shutil, sys
from Input import input_vasp_band, input_vasp_potcar
from Input import input_vasp_kpoints, input_vasp_scf
from Input import get_info
from Calculators import job
from ase.io import read



def band(args, __shell__, __python__, __work__):
    ##### Start Parameters #####
    pressure = args.pressure
    pot = args.pot
    spin = args.spin
    not_sub = args.not_sub
    fd = args.fermi_dirac
    q = args.q
    n = 24 if args.n is None else args.n
    comment = "band" if args.comment is None else args.comment
    symmetry = args.symmetry
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
    # POTCAR
    input_vasp_potcar.potcar(pot)
    # ENCUT
    encut = get_info.get_encut(args.encut)
    # INCAR
    input_vasp_scf.scf(encut, pressure, spin, fd, "gga", u_atom, u_value, lmaxmix)
    os.rename("INCAR", "INCAR-scf")
    input_vasp_band.band(encut, pressure, spin, fd, fun, u_atom, u_value, lmaxmix)
    os.rename("INCAR", "INCAR-band")
    
    
    # KPOINTS and job
    if fun == "gga" or fun == "ggau":
        os.system("vaspkit -task 303 > /dev/null 2>&1")
        job.gen_job("job", "{}/job_band".format(__shell__))
    
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
        job.gen_job("job", "{}/job_band_hse".format(__shell__))
    
    # Job and Sub
    job.control_job("job", q, n, comment)
    if not_sub == True:
        print("<=> Valkyrie: Only generate input file.")
    else:
        job.sub("job")
    
    print("<=> Valkyrie: Running band for {}, ENCUT = {}, POTCAR = {}, fun = {}.".format(poscar, encut, pot, fun))
    os.system("rm -f Brillouin_Zone_3D.jpg INCAR PLOT.in SYMMETRY PRIMCELL.vasp")
    return 0


