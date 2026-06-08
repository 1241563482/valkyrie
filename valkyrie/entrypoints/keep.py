import fnmatch
from pathlib import Path


def plan_removals(paths, keepFiles):
    """Return files that are not matched by any keep pattern.

    Directories are intentionally skipped; this command only removes files and
    symlinks in the current directory.
    """
    removals = []
    for path in paths:
        path = Path(path)
        if not (path.is_file() or path.is_symlink()):
            continue
        if not any(fnmatch.fnmatch(path.name, pattern) for pattern in keepFiles):
            removals.append(path)
    return removals


def removeFiles(keepFiles, dry_run=False):
    for path in plan_removals(Path(".").iterdir(), keepFiles):
        if dry_run:
            print(f"Would delete: {path.name}")
        else:
            path.unlink()
            print(f"Delete: {path.name}")

    for path in Path(".").iterdir():
        if path.is_dir() and not any(fnmatch.fnmatch(path.name, pattern) for pattern in keepFiles):
            print(f"Skip directory: {path.name}")


def main(*args, task=None, keepFiles=[], dry_run=False, **kwargs):
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

    removeFiles(keepFiles, dry_run=dry_run)
    return 0


if __name__ == "__main__":
    main()
