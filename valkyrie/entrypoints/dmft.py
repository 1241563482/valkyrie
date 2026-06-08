import os
import shutil
import subprocess
from pathlib import Path

from ase.io import read, write


def main(*args, queue="9242opa!", nodes=24, comment="dmft", notSub=False,
         fermiDirac=-1, fun="gga", **kwargs):
    if not Path("POSCAR").exists():
        raise Exception("No POSCAR file found at current path.")

    jobname = os.path.basename(os.getcwd())
    poscar = read("POSCAR")
    write("POSCAR.cif", poscar, format="cif")

    subprocess.run(["cif2struct.py", "POSCAR.cif"], check=True)
    shutil.move("POSCAR.struct", f"{jobname}.struct")

    output = subprocess.run(["init_lapw", "-b"], check=True, capture_output=True, text=True).stdout
    if "ERROR" in output:
        shutil.move(f"{jobname}.struct_sgroup", f"{jobname}.struct")
        output = subprocess.run(["init_lapw", "-b"], check=True, capture_output=True, text=True).stdout
    if "init_lapw finished ok" in output:
        print("<=> Valkyrie: Successfully generate input file for wien2k.")
    return 0


if __name__ == "__main__":
    main()
