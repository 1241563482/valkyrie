from periodictable import elements
from ase.io import read, write
from ase.data import chemical_symbols, atomic_numbers

potcar = ["",
    "H.pbe-rrkjus_psl.1.0.0.UPF",
    "He_ONCV_PBE-1.0.oncvpsp.upf",
    "li_pbe_v1.4.uspp.F.UPF",
    "be_pbe_v1.4.uspp.F.UPF",
    "b_pbe_v1.4.uspp.F.UPF",
    "C.pbe-n-kjpaw_psl.1.0.0.UPF",
    "N.pbe-n-radius_5.UPF",
    "O.pbe-n-kjpaw_psl.0.1.UPF",
    "f_pbe_v1.4.uspp.F.UPF",
    "Ne_ONCV_PBE-1.0.oncvpsp.upf",
    "na_pbe_v1.5.uspp.F.UPF",
    "Mg.pbe-n-kjpaw_psl.0.3.0.UPF",
    "Al.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Si.pbe-n-rrkjus_psl.1.0.0.UPF",
    "P.pbe-n-rrkjus_psl.1.0.0.UPF",
    "s_pbe_v1.4.uspp.F.UPF",
    "cl_pbe_v1.4.uspp.F.UPF",
    "Ar_ONCV_PBE-1.1.oncvpsp.upf",
    "K.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Ca_pbe_v1.uspp.F.UPF",
    "Sc_ONCV_PBE-1.0.oncvpsp.upf",
    "ti_pbe_v1.4.uspp.F.UPF",
    "v_pbe_v1.4.uspp.F.UPF",
    "cr_pbe_v1.5.uspp.F.UPF",
    "mn_pbe_v1.5.uspp.F.UPF",
    "Fe.pbe-spn-kjpaw_psl.0.2.1.UPF",
    "Co_pbe_v1.2.uspp.F.UPF",
    "ni_pbe_v1.4.uspp.F.UPF",
    "Cu.paw.z_11.ld1.psl.v1.0.0-low.upf",
    "Zn_pbe_v1.uspp.F.UPF",
    "Ga.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "ge_pbe_v1.4.uspp.F.UPF",
    "As.pbe-n-rrkjus_psl.0.2.UPF",
    "Se_pbe_v1.uspp.F.UPF",
    "br_pbe_v1.4.uspp.F.UPF",
    "Kr_ONCV_PBE-1.0.oncvpsp.upf",
    "Rb_ONCV_PBE-1.0.oncvpsp.upf",
    "Sr_pbe_v1.uspp.F.UPF",
    "Y_pbe_v1.uspp.F.UPF",
    "Zr_pbe_v1.uspp.F.UPF",
    "Nb.pbe-spn-kjpaw_psl.0.3.0.UPF",
    "Mo_ONCV_PBE-1.0.oncvpsp.upf",
    "Tc_ONCV_PBE-1.0.oncvpsp.upf",
    "Ru_ONCV_PBE-1.0.oncvpsp.upf",
    "Rh_ONCV_PBE-1.0.oncvpsp.upf",
    "Pd_ONCV_PBE-1.0.oncvpsp.upf",
    "Ag_ONCV_PBE-1.0.oncvpsp.upf",
    "Cd.pbe-dn-rrkjus_psl.0.3.1.UPF",
    "In.pbe-dn-rrkjus_psl.0.2.2.UPF",
    "Sn_pbe_v1.uspp.F.UPF",
    "sb_pbe_v1.4.uspp.F.UPF",
    "Te_pbe_v1.uspp.F.UPF",
    "I.pbe-n-kjpaw_psl.0.2.UPF",
    "Xe_ONCV_PBE-1.1.oncvpsp.upf",
    "Cs_pbe_v1.uspp.F.UPF",
    "Ba.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "La.paw.z_11.atompaw.wentzcovitch.v1.2.upf",
    "Ce.paw.z_12.atompaw.wentzcovitch.v1.2.upf",
    "Pr.paw.z_13.atompaw.wentzcovitch.v1.2.upf",
    "Nd.paw.z_14.atompaw.wentzcovitch.v1.2.upf",
    "Pm.paw.z_15.atompaw.wentzcovitch.v1.2.upf",
    "Sm.paw.z_16.atompaw.wentzcovitch.v1.2.upf",
    "Eu.paw.z_17.atompaw.wentzcovitch.v1.2.upf",
    "Gd.paw.z_18.atompaw.wentzcovitch.v1.2.upf",
    "Tb.paw.z_19.atompaw.wentzcovitch.v1.2.upf",
    "Dy.paw.z_20.atompaw.wentzcovitch.v1.2.upf",
    "Ho.paw.z_21.atompaw.wentzcovitch.v1.2.upf",
    "Er.paw.z_22.atompaw.wentzcovitch.v1.2.upf",
    "Tm.paw.z_23.atompaw.wentzcovitch.v1.2.upf",
    "Yb.paw.z_24.atompaw.wentzcovitch.v1.2.upf",
    "Lu.paw.z_25.atompaw.wentzcovitch.v1.2.upf",
    "Hf-sp.oncvpsp.upf",
    "Ta_pbe_v1.uspp.F.UPF",
    "W_pbe_v1.2.uspp.F.UPF",
    "Re_pbe_v1.2.uspp.F.UPF",
    "Os_pbe_v1.2.uspp.F.UPF",
    "Ir_pbe_v1.2.uspp.F.UPF",
    "pt_pbe_v1.4.uspp.F.UPF",
    "Au_ONCV_PBE-1.0.oncvpsp.upf",
    "Hg_ONCV_PBE-1.0.oncvpsp.upf",
    "Tl_pbe_v1.2.uspp.F.UPF",
    "Pb.pbe-dn-kjpaw_psl.0.2.2.UPF",
    "Bi_pbe_v1.uspp.F.UPF",
    "Po.pbe-dn-rrkjus_psl.1.0.0.UPF",
    "At.us.z_17.ld1.psl.v1.0.0-high.upf",
    "Rn.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "Fr.paw.z_19.ld1.psl.v1.0.0-high.upf",
    "Ra.paw.z_20.ld1.psl.v1.0.0-high.upf",
    "Ac.us.z_11.ld1.psl.v1.0.0-high.upf",
    "Th.paw.z_12.ld1.uni-marburg.v0.upf",
    "Pa.paw.z_13.ld1.uni-marburg.v0.upf",
    "U.paw.z_14.ld1.uni-marburg.v0.upf",
    "Np.paw.z_15.ld1.uni-marburg.v0.upf",
    "Pu.paw.z_16.ld1.uni-marburg.v0.upf",
    "Am.paw.z_17.ld1.uni-marburg.v0.upf",
    "Cm.paw.z_18.ld1.uni-marburg.v0.upf",
    "Bk.paw.z_19.ld1.uni-marburg.v0.upf",
    "Cf.paw.z_20.ld1.uni-marburg.v0.upf",
    "Es.paw.z_21.ld1.uni-marburg.v0.upf",
    "Fm.paw.z_22.ld1.uni-marburg.v0.upf",
    "Md.paw.z_23.ld1.uni-marburg.v0.upf",
    "No.paw.z_24.ld1.uni-marburg.v0.upf",
    "Lr.paw.z_25.ld1.uni-marburg.v0.upf"
]

def get_pot(poscar):
    poscar = read(poscar)
    atoms = poscar.get_chemical_symbols()
    uniq_atoms = []
    for i in atoms:
        if i not in uniq_atoms:
            uniq_atoms.append(i)
    pot = ""
    for atom in uniq_atoms:
        pot += atom
        pot += " "
        pot += str(elements.symbol(atom).mass)
        pot += " "
        pot += potcar[atomic_numbers[atom]]
        pot += "\n"

    return pot