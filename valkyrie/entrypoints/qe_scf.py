import shutil
import subprocess
from pathlib import Path

import numpy as np
from ase.io import read

from .. import run_pw
from ..input import input_qe_pot, input_qe_scf, poscar2qe
from ..resources import shell_template
from ..shell_scripts import job


def _kmesh_from_poscar(poscar):
    atom_size = [max(poscar.get_positions()[:, i]) - min(poscar.get_positions()[:, i]) for i in [0, 1, 2]]
    cell_size = [np.linalg.norm(poscar.cell[i]) for i in [0, 1, 2]]
    k = []
    for i in [0, 1, 2]:
        if cell_size[i] - atom_size[i] >= 8:
            k.append(1)
        else:
            k.append(int(max(1, np.ceil(40 / cell_size[i]))))
    return k


def main(*args, queue="9242opa!", nodes=24, comment="qe_scf", notSub=False,
         pressure=0, encut=100, spin=False, fermiDirac=-1, symmetry=None,
         **kwargs):
    if not Path("POSCAR").exists():
        raise Exception("No POSCAR file found at current path.")

    poscar = read("POSCAR")
    formula = poscar.get_chemical_formula()
    if symmetry is not None:
        result = subprocess.run(
            ["phonopy", "--symmetry", f"--tolerance={symmetry}"],
            check=True,
            capture_output=True,
            text=True,
        )
        for line in result.stdout.splitlines():
            if "space" in line.lower():
                print(line)
                break
        shutil.copy2("PPOSCAR", "POSCAR")
        poscar = read("POSCAR")

    ntyp, nat, cell, pos = poscar2qe.read_poscar("POSCAR")
    k = _kmesh_from_poscar(poscar)
    pot = input_qe_pot.get_pot("POSCAR")
    input_qe_scf.scf(pressure * 10, encut, k, nat, ntyp, cell, pos, pot)

    job.gen_job(job="job", job_file=shell_template("job_qe_scf"), task=run_pw, inputFile="< scf.in")
    job.control_job("job", queue, nodes, comment)
    print("<=> Valkyrie: Only generate input file.") if notSub else job.sub("job")

    print("<=> Valkyrie: Scf for {} under the pressure of {} GPa, ENCUT = {} Ry, POTCAR = {}."\
        .format(formula, pressure, encut, pot.split()[-1]))
    return 0


if __name__ == "__main__":
    main()
