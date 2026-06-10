from pathlib import Path

from .. import config, run_abacus
from ..abacus_ase_io import poscar_to_stru, read_poscar, write_input, write_kpt
from ..resources import shell_template
from ..shell_scripts import job


FALLBACK_INPUT = {
    "calculation": "scf",
    "ecutwfc": 120,
    "scf_thr": "1e-06",
    "scf_nmax": 400,
    "out_chg": 0,
    "cal_force": 1,
    "smearing_method": "gauss",
    "smearing_sigma": 0.01,
    "mixing_type": "broyden",
    "mixing_ndim": 10,
    "ks_solver": "genelpa",
    "basis_type": "lcao",
    "symmetry": -1,
    "dft_functional": "pbe",
}

SOC_DEFAULTS = {
    "noncolin": 1,
    "lspinorb": 1,
    "nspin": 4,
    "out_mul": True,
    "sc_mag_switch": 1,
    "sc_thr": "1e-05",
    "nsc": 400,
    "nsc_min": 2,
    "alpha_trial": 0.01,
    "sccut": 3,
}


def parse_value(value):
    if value in {"True", "False"}:
        return value == "True"
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def read_abacus_input(path):
    path = Path(path)
    if not path.exists():
        return {}
    params = {}
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or stripped == "INPUT_PARAMETERS":
            continue
        parts = stripped.split(maxsplit=1)
        if len(parts) == 1:
            params[parts[0]] = ""
        else:
            params[parts[0]] = parse_value(parts[1])
    return params


def read_abacus_kpt(path):
    path = Path(path)
    if not path.exists():
        return {}
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if len(lines) < 4 or lines[0] != "K_POINTS":
        return {}
    values = [int(float(value)) for value in lines[3].split()]
    return {
        "kmode": lines[2],
        "kpts": values[:3],
        "koffset": values[3:6] if len(values) >= 6 else [0, 0, 0],
    }


def _is_section_header(line):
    stripped = line.strip()
    return bool(stripped) and stripped.replace("_", "").isupper()


def read_stru_pp_basis(path):
    path = Path(path)
    if not path.exists():
        return {}, {}

    lines = path.read_text().splitlines()
    species = []
    pp = {}
    basis_files = []
    section = None
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if _is_section_header(stripped):
            section = stripped
            continue
        if section == "ATOMIC_SPECIES":
            parts = stripped.split()
            if len(parts) >= 3:
                species.append(parts[0])
                pp[parts[0]] = parts[2]
        elif section == "NUMERICAL_ORBITAL":
            basis_files.append(stripped.split()[0])

    basis = {symbol: basis_files[i] for i, symbol in enumerate(species) if i < len(basis_files)}
    return pp, basis


def parse_key_value_map(items):
    result = {}
    for item in items or []:
        if "=" not in item:
            raise ValueError(f"Expected ELEMENT=FILE mapping, got {item!r}")
        key, value = item.split("=", 1)
        key = key.strip()
        value = value.strip()
        if not key or not value:
            raise ValueError(f"Expected ELEMENT=FILE mapping, got {item!r}")
        result[key] = value
    return result


def _missing_elements(atoms, mapping):
    return [symbol for symbol in dict.fromkeys(atoms.get_chemical_symbols()) if symbol not in mapping]


def _apply_no_soc(params):
    params["noncolin"] = 0
    params["lspinorb"] = 0
    params["nspin"] = 1
    return params


def build_input_parameters(reference_dir, pseudo_dir=None, orbital_dir=None, ecutwfc=None,
                           basis_type="lcao", no_soc=False):
    params = dict(FALLBACK_INPUT)
    params.update(SOC_DEFAULTS)
    params.update(read_abacus_input(Path(reference_dir) / "INPUT"))
    params.setdefault("calculation", "scf")
    params.setdefault("pseudo_dir", config.abacus_pseudo_dir)
    params.setdefault("orbital_dir", config.abacus_orbital_dir)
    if pseudo_dir is not None:
        params["pseudo_dir"] = pseudo_dir
    if orbital_dir is not None:
        params["orbital_dir"] = orbital_dir
    if ecutwfc is not None:
        params["ecutwfc"] = ecutwfc
    if basis_type is not None:
        params["basis_type"] = basis_type
    if no_soc:
        _apply_no_soc(params)
    return params


def build_kpt_parameters(reference_dir, kpoints=None):
    params = {"kmode": "Gamma", "kpts": [12, 15, 1], "koffset": [0, 0, 0]}
    params.update(read_abacus_kpt(Path(reference_dir) / "KPT"))
    if kpoints is not None:
        params["kpts"] = kpoints
    return params


def main(*args, poscar="POSCAR", reference_dir="1", pseudo_dir=None, orbital_dir=None,
         pp=None, basis=None, kpoints=None, ecutwfc=None, basis_type="lcao",
         no_soc=False, queue="9242opa!", nodes=24, comment="abacus_scf",
         notSub=False, command=None, no_record=False, **kwargs):
    poscar_path = Path(poscar)
    if not poscar_path.exists():
        raise Exception("No POSCAR file found at current path.")

    atoms = read_poscar(poscar_path)
    reference_path = Path(reference_dir)

    ref_pp, ref_basis = read_stru_pp_basis(reference_path / "STRU")
    pp_map = {**ref_pp, **parse_key_value_map(pp)}
    basis_map = {**ref_basis, **parse_key_value_map(basis)}

    missing_pp = _missing_elements(atoms, pp_map)
    if missing_pp:
        missing = ", ".join(missing_pp)
        raise ValueError(f"No ABACUS pseudopotential mapping for element(s): {missing}. Use --pp ELEMENT=FILE.")

    input_params = build_input_parameters(
        reference_path,
        pseudo_dir=pseudo_dir,
        orbital_dir=orbital_dir,
        ecutwfc=ecutwfc,
        basis_type=basis_type,
        no_soc=no_soc,
    )
    if input_params.get("basis_type") == "lcao":
        missing_basis = _missing_elements(atoms, basis_map)
        if missing_basis:
            missing = ", ".join(missing_basis)
            raise ValueError(f"No ABACUS orbital mapping for element(s): {missing}. Use --basis ELEMENT=FILE.")
    else:
        basis_map = None

    kpt_params = build_kpt_parameters(reference_path, kpoints=kpoints)

    poscar_to_stru(poscar_path, "STRU", pp=pp_map, basis=basis_map, scaled=True)
    write_input("INPUT", parameters=input_params)
    write_kpt("KPT", parameters=kpt_params, atoms=atoms)

    job.gen_job(job="job", job_file=shell_template("job_abacus_scf"), task=run_abacus)
    job.control_job("job", queue, nodes, comment)
    print("<=> Valkyrie: Only generate input file.") if notSub else job.sub("job")

    print(f"<=> Valkyrie: ABACUS scf for {atoms.get_chemical_formula()}.")
    print(f"<=> Valkyrie: KPOINTS = {' '.join(map(str, kpt_params['kpts']))}, ecutwfc = {input_params.get('ecutwfc')}.")
    print(f"<=> Valkyrie: Pseudopotentials = {pp_map}.")
    if basis_map:
        print(f"<=> Valkyrie: Orbitals = {basis_map}.")
    return 0


if __name__ == "__main__":
    main()
