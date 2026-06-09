"""Environment and input checks for Valkyrie."""

import importlib
import json
import shlex
import shutil
import sys
from pathlib import Path

from .. import config
from ..resources import python_script, shell_template


REQUIRED_MODULES = ["numpy", "ase", "matplotlib", "phonopy", "periodictable"]
SHELL_TEMPLATES = [
    "job_vasp_relax",
    "job_scf",
    "job_band",
    "job_dos",
    "job_elf",
    "job_ph",
    "job_vasp_md",
    "job_cohp",
    "job_qe_relax",
    "job_qe_scf",
    "job_dmft",
]
PYTHON_SCRIPTS = [
    "vasp_dos_sum_spd.py",
    "vasp_elf_electride.py",
    "vasp_md_MDdata.py",
    "vasp_ph_loop_in.py",
    "vasp_cohp.py",
    "get_free_energy.py",
]
EXTERNAL_TOOLS = ["bsub", "vaspkit", "phonopy", "lobster", "cif2struct.py", "init_lapw"]


def check_item(name, status, message, detail=None):
    return {"name": name, "status": status, "message": message, "detail": detail}


def command_executable(command):
    """Return the executable token from a shell-like command string.

    MPI launchers are skipped so the actual scientific executable is checked.
    """
    parts = shlex.split(command)
    skip_next = False
    for part in parts:
        if skip_next:
            skip_next = False
            continue
        if part in {"mpiexec", "mpiexec.hydra", "mpirun", "srun"}:
            continue
        if part in {"-np", "-n", "-npool", "--np", "--ntasks"}:
            skip_next = True
            continue
        if part.startswith("-"):
            continue
        return part
    return parts[0] if parts else None


def check_executable(name, executable, required=False):
    if not executable:
        return check_item(name, "warn", "no executable configured")
    path = Path(executable)
    if path.is_absolute():
        if path.exists():
            return check_item(name, "ok", "configured executable exists", str(path))
        return check_item(name, "error" if required else "warn", "configured executable not found", str(path))
    found = shutil.which(executable)
    if found:
        return check_item(name, "ok", "executable found on PATH", found)
    return check_item(name, "error" if required else "warn", "executable not found on PATH", executable)


def check_python_version():
    version = sys.version_info
    if version >= (3, 8):
        return check_item("python", "ok", f"Python {version.major}.{version.minor}.{version.micro}")
    return check_item("python", "error", "Valkyrie requires Python >= 3.8", sys.version.split()[0])


def check_imports():
    results = []
    for module_name in REQUIRED_MODULES:
        try:
            module = importlib.import_module(module_name)
            results.append(check_item(f"python module {module_name}", "ok", "importable", getattr(module, "__version__", None)))
        except Exception as exc:
            results.append(check_item(f"python module {module_name}", "error", "not importable", str(exc)))
    return results


def check_resources():
    results = []
    for name in SHELL_TEMPLATES:
        path = Path(shell_template(name))
        status = "ok" if path.exists() else "error"
        message = "resource exists" if path.exists() else "resource missing"
        results.append(check_item(f"shell resource {name}", status, message, str(path)))
    for name in PYTHON_SCRIPTS:
        path = Path(python_script(name))
        status = "ok" if path.exists() else "error"
        message = "resource exists" if path.exists() else "resource missing"
        results.append(check_item(f"python resource {name}", status, message, str(path)))
    return results


def check_configured_paths():
    results = []
    for key in ["run_vasp", "run_vasp_opt", "run_vasp_gam", "run_pw"]:
        executable = command_executable(getattr(config, key))
        results.append(check_executable(key, executable))
    qe_path = Path(config.qe_pot_path)
    if qe_path.exists():
        results.append(check_item("qe_pot_path", "ok", "path exists", str(qe_path)))
    else:
        results.append(check_item("qe_pot_path", "warn", "path not found", str(qe_path)))
    return results


def check_external_tools():
    return [check_executable(f"tool {tool}", tool) for tool in EXTERNAL_TOOLS]


def check_local_inputs(path="."):
    directory = Path(path)
    results = []
    poscar = directory / "POSCAR"
    if poscar.exists():
        try:
            from ase.io import read

            atoms = read(str(poscar))
            results.append(check_item("POSCAR", "ok", "readable by ASE", atoms.get_chemical_formula()))
        except Exception as exc:
            results.append(check_item("POSCAR", "error", "exists but is not readable by ASE", str(exc)))
    else:
        results.append(check_item("POSCAR", "warn", "not found in current directory", str(poscar)))

    potcar = directory / "POTCAR"
    if potcar.exists():
        try:
            text = potcar.read_text(errors="ignore")
            if "ENMAX" in text:
                results.append(check_item("POTCAR", "ok", "contains ENMAX entries", str(potcar)))
            else:
                results.append(check_item("POTCAR", "warn", "exists but no ENMAX entries were found", str(potcar)))
        except OSError as exc:
            results.append(check_item("POTCAR", "error", "could not read POTCAR", str(exc)))
    else:
        results.append(check_item("POTCAR", "warn", "not found in current directory", str(potcar)))
    return results


def run_checks(path="."):
    checks = [check_python_version()]
    checks.extend(check_imports())
    checks.extend(check_resources())
    checks.extend(check_configured_paths())
    checks.extend(check_external_tools())
    checks.extend(check_local_inputs(path))
    return checks


def summarize(checks, strict=False):
    errors = sum(1 for item in checks if item["status"] == "error")
    warnings = sum(1 for item in checks if item["status"] == "warn")
    ok = sum(1 for item in checks if item["status"] == "ok")
    failed = errors > 0 or (strict and warnings > 0)
    return {"ok": ok, "warnings": warnings, "errors": errors, "failed": failed}


def print_text(checks, summary):
    print("Valkyrie check")
    for item in checks:
        label = item["status"].upper()
        print(f"[{label}] {item['name']}: {item['message']}")
        if item.get("detail"):
            print(f"      {item['detail']}")
    print(f"\nSummary: {summary['ok']} ok, {summary['warnings']} warnings, {summary['errors']} errors")


def main(*args, json_output=False, strict=False, **kwargs):
    checks = run_checks()
    summary = summarize(checks, strict=strict)
    data = {"checks": checks, "summary": summary}
    if json_output:
        print(json.dumps(data, indent=2, sort_keys=True))
    else:
        print_text(checks, summary)
    return 1 if summary["failed"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
