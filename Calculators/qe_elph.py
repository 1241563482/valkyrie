from ase.io import read, write
from Input import poscar2qe, input_qe_elph
from Input import input_qe_pot
from Calculators import sub_job


def elph(args, __shell__, __python__, __work__):
    # Parameters:
    poscar = args.f
    encut = args.encut
    qmesh = args.qmesh
    k = args.k
    not_sub = args.not_sub
    q = args.q
    n = args.n
    comment = args.comment
    if q == None:
        q = "9242opa!"
    if n == None:
        n = "48"

    # Read POSCAR
    ntyp, nat, cell, pos = poscar2qe.read_poscar(poscar)
    pot = input_qe_pot.get_pot(poscar)

    # Write input
    print("<=> Valkyrie: Running elph for {}, ecutwfc = {}.".format(poscar, encut))
    print("<=> Valkyrie: Pot for elph are \n{}".format(pot), end = "")
    input_qe_elph.scf1(encut, k, nat, ntyp, cell, pos, pot)
    input_qe_elph.scf2(encut, k, nat, ntyp, cell, pos, pot)
    input_qe_elph.loop_elph(qmesh, n, q)
    input_qe_elph.matdyn_in(poscar, pot, k)
    input_qe_elph.other_input(__shell__, __python__, __work__)

    # Sub job
    if not_sub == True:
        sub_job.control_job(q, n, comment)
        print("<=> Valkyrie: Only generate input file.")
    else:
        sub_job.control_job(q, n, comment)
        sub_job.sub("job")
        print("<=> Valkyrie: Running ELPH ... ...")

