from ase.io import read, write
import os, subprocess

def dmft(args, __shell__, __python__, __work__):
    ##### Start Parameters #####
    not_sub = args.not_sub
    fd = args.fermi_dirac
    q = args.q
    n = 24 if args.n is None else args.n
    comment = "dmft" if args.comment is None else args.comment
    jobname = os.path.basename(os.getcwd())
    ##### End Parameters #####


    poscar = read("POSCAR")
    write("POSCAR.cif", poscar, format="cif")
    os.system("cif2struct.py POSCAR.cif")
    os.system("mv POSCAR.struct {}.struct".format(jobname))

    output = subprocess.check_output("init_lapw -b")
    if "ERROR" in output:
        os.system("mv {}.struct_sgroup {}.struct".format(jobname, jobname))
        output = subprocess.check_output("init_lapw -b")
    if "init_lapw finished ok" in output:
        print("<=> Valkyrie: Successfully generate input file for wien2k.")
    

