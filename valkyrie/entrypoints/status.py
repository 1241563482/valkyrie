"""Local calculation status reporting."""

import json
import re
from pathlib import Path

from ..run_record import load_history, load_latest_run


FLOAT_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
VASP_OUTCAR_NAMES = ["OUTCAR", "OUTCAR_scf", "OUTCAR_dos", "OUTCAR_band"]
QE_OUTPUT_NAMES = ["scf.out", "relax.out", "pw.out"]


def _read_text(path):
    return Path(path).read_text(errors="ignore")


def parse_vasp_outcar(path):
    """Parse a VASP OUTCAR-like file."""
    path = Path(path)
    text = _read_text(path)
    energies = [float(match) for match in re.findall(r"free\s+energy\s+TOTEN\s+=\s+(%s)\s+eV" % FLOAT_RE, text)]
    sigma0 = [
        float(match)
        for match in re.findall(r"energy\(sigma->0\)\s*=\s*(%s)" % FLOAT_RE, text)
    ]
    pressures = [
        float(match)
        for match in re.findall(r"(?:external|total)\s+pressure\s*=\s*(%s)" % FLOAT_RE, text, flags=re.IGNORECASE)
    ]
    elapsed = re.findall(r"Elapsed time \(sec\):\s*(%s)" % FLOAT_RE, text)
    completed = "General timing and accounting" in text or "Voluntary context switches" in text
    converged = "reached required accuracy" in text or "aborting loop because EDIFF is reached" in text
    result = {
        "file": path.name,
        "kind": "vasp_outcar",
        "completed": completed,
        "converged": converged,
        "steps": len(energies),
    }
    if energies:
        result["final_energy"] = {"value": energies[-1], "unit": "eV"}
    if sigma0:
        result["sigma_zero_energy"] = {"value": sigma0[-1], "unit": "eV"}
    if pressures:
        result["pressure"] = {"value": pressures[-1], "unit": "kB"}
    if elapsed:
        result["elapsed_time_sec"] = float(elapsed[-1])
    return result


def parse_qe_output(path):
    """Parse a Quantum ESPRESSO output file."""
    path = Path(path)
    text = _read_text(path)
    energies = [
        float(match)
        for match in re.findall(r"!\s+total energy\s+=\s+(%s)\s+Ry" % FLOAT_RE, text)
    ]
    forces = [
        float(match)
        for match in re.findall(r"Total force\s+=\s+(%s)" % FLOAT_RE, text, flags=re.IGNORECASE)
    ]
    completed = "JOB DONE." in text
    lower_text = text.lower()
    converged = "convergence has been achieved" in lower_text or "bfgs converged" in lower_text
    result = {
        "file": path.name,
        "kind": "qe_output",
        "completed": completed,
        "converged": converged,
        "steps": len(re.findall(r"iteration\s+#", text, flags=re.IGNORECASE)),
    }
    if energies:
        result["final_energy"] = {"value": energies[-1], "unit": "Ry"}
    if forces:
        result["total_force"] = {"value": forces[-1], "unit": "Ry/Bohr"}
    return result


def read_structure_summary(path):
    """Read POSCAR/CONTCAR formula if ASE can parse it."""
    directory = Path(path)
    for name in ["CONTCAR", "POSCAR"]:
        structure = directory / name
        if not structure.exists() or structure.stat().st_size == 0:
            continue
        try:
            from ase.io import read

            atoms = read(str(structure))
            return {"file": name, "formula": atoms.get_chemical_formula(), "natoms": len(atoms)}
        except Exception as exc:
            return {"file": name, "error": str(exc)}
    return None


def find_vasp_outputs(path):
    directory = Path(path)
    outputs = []
    for name in VASP_OUTCAR_NAMES:
        candidate = directory / name
        if candidate.exists():
            outputs.append(candidate)
    for candidate in sorted(directory.glob("OUTCAR*")):
        if candidate not in outputs and candidate.is_file():
            outputs.append(candidate)
    return outputs


def find_qe_outputs(path):
    directory = Path(path)
    outputs = []
    for name in QE_OUTPUT_NAMES:
        candidate = directory / name
        if candidate.exists():
            outputs.append(candidate)
    for candidate in sorted(directory.glob("*.out")):
        if candidate not in outputs and candidate.is_file():
            outputs.append(candidate)
    return outputs


def detect_engine(path):
    if find_vasp_outputs(path):
        return "vasp"
    if find_qe_outputs(path):
        return "qe"
    latest = load_latest_run(path)
    if latest:
        command = latest.get("command", "")
        if command.startswith("vasp"):
            return "vasp"
        if command.startswith("qe"):
            return "qe"
    return "unknown"


def collect_status(path=".", engine="auto", history=0):
    directory = Path(path).resolve()
    latest = load_latest_run(directory)
    status = {
        "path": str(directory),
        "run": latest,
        "history": load_history(directory, history) if history else [],
        "engine": engine,
        "state": "unknown",
        "structure": read_structure_summary(directory),
        "outputs": [],
        "notes": [],
    }
    if not directory.exists():
        status["state"] = "missing_path"
        status["notes"].append("directory does not exist")
        return status

    detected = detect_engine(directory) if engine == "auto" else engine
    status["engine"] = detected

    if detected == "vasp":
        outputs = find_vasp_outputs(directory)
        if not outputs:
            status["state"] = "no_output"
            status["notes"].append("no VASP OUTCAR files found")
            return status
        status["outputs"] = [parse_vasp_outcar(output) for output in outputs]
    elif detected == "qe":
        outputs = find_qe_outputs(directory)
        if not outputs:
            status["state"] = "no_output"
            status["notes"].append("no QE .out files found")
            return status
        status["outputs"] = [parse_qe_output(output) for output in outputs]
    else:
        status["state"] = "no_output"
        status["notes"].append("no VASP or QE output files found")
        return status

    primary = status["outputs"][0] if status["outputs"] else {}
    if primary.get("completed"):
        status["state"] = "completed"
    else:
        status["state"] = "running_or_incomplete"
    return status


def print_text(status):
    print(f"Directory: {status['path']}")
    run = status.get("run")
    if run:
        print(f"Last val command: {run.get('command')} at {run.get('recorded_at')}")
    else:
        print("Last val command: none recorded")
    print(f"Engine: {status['engine']}")
    print(f"State: {status['state']}")
    structure = status.get("structure")
    if structure:
        if "formula" in structure:
            print(f"Structure: {structure['formula']} ({structure['natoms']} atoms, {structure['file']})")
        else:
            print(f"Structure: could not read {structure['file']} ({structure['error']})")
    for output in status.get("outputs", []):
        print(f"\nOutput: {output['file']}")
        print(f"  completed: {output['completed']}")
        print(f"  converged: {output['converged']}")
        print(f"  steps: {output['steps']}")
        if "final_energy" in output:
            energy = output["final_energy"]
            print(f"  final energy: {energy['value']} {energy['unit']}")
        if "pressure" in output:
            pressure = output["pressure"]
            print(f"  pressure: {pressure['value']} {pressure['unit']}")
        if "total_force" in output:
            force = output["total_force"]
            print(f"  total force: {force['value']} {force['unit']}")
        if "elapsed_time_sec" in output:
            print(f"  elapsed time: {output['elapsed_time_sec']} sec")
    for note in status.get("notes", []):
        print(f"Note: {note}")
    if status.get("history"):
        print("\nHistory:")
        for record in status["history"]:
            print(f"  {record.get('recorded_at')}  {record.get('command')}")


def main(*args, path=".", json_output=False, engine="auto", history=0, **kwargs):
    status = collect_status(path=path, engine=engine, history=history)
    if json_output:
        print(json.dumps(status, indent=2, sort_keys=True))
    else:
        print_text(status)
    return 0


if __name__ == "__main__":
    main()
