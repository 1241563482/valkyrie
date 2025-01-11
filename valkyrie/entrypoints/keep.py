import os, fnmatch

def removeFiles(keepFiles):
    files = os.listdir()
    for file in files:
        if not any(fnmatch.fnmatch(file, pattern) for pattern in keepFiles):
            os.remove(file)
            print(f"Delete: {file}")

def main(*args, task = None, keepFiles = [], **kwargs):
    vaspStan = ["INCAR", "KPOINTS", "POSCAR", "POTCAR", "OUTCAR", "job*", "*py", "*png"]
    if task is None and keepFiles == []:
        print("<=> Valkyrie: Enter the files to keep!")
        return 0
    elif task == "vasp_relax":
        keepFiles = vaspStan + ["CONTCAR", "XDATCAR"]
    elif task == "vasp_scf":
        keepFiles = vaspStan
    elif task == "vasp_elf":
        keepFiles = vaspStan + ["ELFCAR"]
    elif task == "vasp_md":
        keepFiles = vaspStan + ["XDATCAR"]
    
    removeFiles(keepFiles)
    return 0


if __name__ == "__main__":
    main()