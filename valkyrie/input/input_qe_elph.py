from ase.dft import kpoints
from ase.io import read
import shutil
from Calculators import job
from set_up import sub_command

def scf1(encut, k, nat, ntyp, cell, pos, pot):
    if k[2] == 1:
        k2 = 1
    else:
        k2 = 2 * k[2]

    scf1_in = """&CONTROL
    calculation = 'scf' , prefix = 'pwscf'
    pseudo_dir = '/fsa/home/js_zhuyj/mypps/QE/SSSP_1.3.0_PBE_efficiency'
    outdir = './tmp'
    forc_conv_thr = 1.0d-7
    etot_conv_thr = 1.0d-7
/
&SYSTEM
    ibrav= 0, nat= {}, ntyp= {},  
    occupations = 'smearing', smearing = 'mp', degauss = 0.02
    ecutwfc = {}, la2F = .true.
/
&ELECTRONS
    electron_maxstep = 100
    conv_thr = 1.0d-11
/
&IONS
    ion_dynamics='bfgs'
/
ATOMIC_SPECIES
{}
K_POINTS automatic
{} {} {} 0 0 0

{}
{}
""".format(nat, ntyp, encut, pot, 2 * k[0], 2 * k[1], k2, cell, pos)
    with open("scf1.in", "w") as file:
        print(scf1_in, file = file)
    print("<=> Valkyrie: Kmesh for scf1.in: {} {} {}".format(2 * k[0], 2 * k[1], k2))

def scf2(encut, k, nat, ntyp, cell, pos, pot):
    scf2 = """&CONTROL
    calculation = 'scf', prefix = 'pwscf'
    pseudo_dir = '/fsa/home/js_zhuyj/mypps/QE/SSSP_1.3.0_PBE_efficiency' 
    outdir = './tmp'
    forc_conv_thr = 1.0d-7
    etot_conv_thr = 1.0d-7
/
&SYSTEM
    ibrav= 0, nat= {}, ntyp= {},  
    occupations = 'smearing', smearing = 'mp', degauss = 0.02
    ecutwfc = {}
/
&ELECTRONS
    electron_maxstep = 100
    conv_thr = 1.0d-11
/
&IONS
    ion_dynamics='bfgs'
/
ATOMIC_SPECIES
{}
K_POINTS automatic
{} {} {} 0 0 0

{}
{}
""".format(nat, ntyp, encut, pot, k[0], k[1], k[2], cell, pos)
    with open("scf2.in", "w") as file:
        print(scf2, file = file)
    print("<=> Valkyrie: Kmesh for scf2.in: {} {} {}".format(k[0], k[1], k[2]))


def loop_elph(qmesh, n, q):
    loop_elph = """#!/bin/bash
PREFIX="pwscf"

# q-vectors parallel params.
nbeg=1
nqs=8

for ((iq=$nbeg; iq<=$nqs; ++iq))  
do
    # create temp dir for each q-vector.
    phq_tmp=tmp/phq_$iq
    mkdir -p $phq_tmp
    cd $phq_tmp
    ln -sf  ../$PREFIX.save $PREFIX.save
    ln -sf  ../$PREFIX.xml $PREFIX.xml
    ln -sf  ../$PREFIX.a2Fsave $PREFIX.a2Fsave
    cd ../../
    
    cat > elph.in_$iq << EOF
&INPUTPH
    tr2_ph = 1.0e-16,
    prefix = '$PREFIX',
    outdir = '$phq_tmp',
    fildvscf = 'dvscf',
    fildyn = '$PREFIX.dyn',
    ldisp = .true.,
    trans = .true.,   ! Not recalculate the phonon.
    asr = .true.
!    alpha_mix = 0.3
!    recover = .true.,
    nq1={}, nq2={}, nq3={}
    electron_phonon = 'interpolated'
    start_q = $iq,
    last_q = $iq,
    el_ph_sigma = 0.02, 
    el_ph_nsigma = 30
/
EOF

    cp job-elph job-tmp
    sed -i "s/NUMBER/$iq/g" job-tmp
    cat job-tmp
    {} job-tmp
    rm -f job-tmp
done
""".format(qmesh[0], qmesh[1], qmesh[2], sub_command)
    with open("loop-elph.sh", "w") as file:
        print(loop_elph, file = file)
    print("<=> Valkyrie: Qmesh for elph: {} {} {}".format(qmesh[0], qmesh[1], qmesh[2]))

def matdyn_in(poscar, pot, k):
    f = read(poscar)
    # Kpoints part
    k_points = kpoints.get_special_points(f.cell)
    k_points_matdyn = ""
    for i, point in enumerate(k_points):
        if i == len(k_points) - 1:
            n = "1"
        else:
            n = "30"
        for j in k_points[point]:
            j = "{:.12f}".format(j)
            k_points_matdyn += str(j)
            k_points_matdyn += " "
        k_points_matdyn += n
        k_points_matdyn += "\n"
    
    # Mass part
    nat = len(pot.splitlines())
    mass_matdyn = ""
    for i in range(nat):
        mass = pot.splitlines()[i].split()[1]
        mass_matdyn += "amass({}) = {}".format(str(i+1), str(mass))
        if i == nat - 1:
            continue
        mass_matdyn += "\n"
    
    matdyn_in_freq = """&input
asr = 'simple'
{}
flfrc = 'ph.fc'
flfrq = 'ph.freq'
q_in_band_form = .true.
q_in_cryst_coord = .true.
la2F = .true.
dos = .false.
/
{}
{}
""".format(mass_matdyn, str(len(k_points)), k_points_matdyn)
    with open("matdyn.in.freq", "w") as file:
        print(matdyn_in_freq, file = file)

    matdyn_in_dos = """&input
asr = 'simple'
{}
flfrc = 'ph.fc'
flfrq = 'ph.freq'
q_in_band_form = .true.
la2F = .true.
dos = .true.
fldos = 'phonon.dos'
nk1={}, nk2={}, nk3={}
ndos = 500
""".format(mass_matdyn, k[0], k[1], k[2])
    with open("matdyn.in.dos", "w") as file:
        print(matdyn_in_dos, file = file)


def other_input(__shell__, __python__, __work__):
    q2r_in = """&input
fildyn = 'pwscf.dyn'
zasr = 'simple'
flfrc = 'ph.fc'
la2f = .true.
/
"""
    with open("q2r.in", "w") as file:
        print(q2r_in, file = file)

    plotband_in = """ph.freq
-200 1200
freq.plot
freq.ps
0.0
100.0 0.0"""
    with open("plotband.in", "w") as file:
        print(plotband_in, file = file)

    shutil.copy("{}/qe_elph_alpha2F.py".format(__python__), "qe_alpha2F.py")
    shutil.copy("{}/qe_elph_sumlambda.py".format(__python__), "sumlambda.py")
    shutil.copy("{}/qe_elph_linewidth.sh".format(__shell__), "qe_elph_linewidth.sh")
    job.gen_job("job", "{}/job_qe_elph_1".format(__shell__))
    job.gen_job("job-elph", "{}/job_qe_elph_2".format(__shell__))
    job.gen_job("job-re", "{}/job_qe_elph_reprocess".format(__shell__))


