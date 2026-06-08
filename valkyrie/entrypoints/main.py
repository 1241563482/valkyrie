import argparse
import importlib

from .. import __picture__, __version__


COMMAND_MODULES = {
    "keep": "keep",
    "vasp_relax": "vasp_relax",
    "vasp_scf": "vasp_scf",
    "vasp_band": "vasp_band",
    "vasp_dos": "vasp_dos",
    "vasp_elf": "vasp_elf",
    "vasp_ph": "vasp_ph",
    "vasp_md": "vasp_md",
    "vasp_cohp": "vasp_cohp",
    "qe_relax": "qe_relax",
    "qescf": "qe_scf",
    "qe_scf": "qe_scf",
    "dmft": "dmft",
}


def add_common_job_args(parser, comment):
    parser.add_argument(
        "-ns",
        "--notSub",
        action="store_true",
        help="Only generate the input file, not sub.",
    )
    parser.add_argument(
        "-q",
        "--queue",
        type=str,
        default="9242opa!",
        help="Node name.",
    )
    parser.add_argument(
        "-n",
        "--nodes",
        type=int,
        default=24,
        help="Number of cores.",
    )
    parser.add_argument(
        "-c",
        "--comment",
        type=str,
        default=comment,
        help="Comment.",
    )


def add_vasp_args(parser):
    parser.add_argument(
        "--encut",
        type=int,
        default=0,
        help="ENCUT",
    )
    parser.add_argument(
        "--pot",
        type=str,
        default="auto",
        nargs="+",
        help="POTCAR type, eg Li_sv or V",
    )
    parser.add_argument(
        "--spin",
        action="store_true",
        help="Add spin (default FM).",
    )
    parser.add_argument(
        "--fun",
        type=str,
        default="gga",
        choices=["gga", "ggau", "hse"],
        help="Functional.",
    )
    parser.add_argument(
        "-u",
        type=str,
        nargs=2,
        default=None,
        help="Atom for GGA+U and Ueff.",
    )
    parser.add_argument(
        "--fElectron",
        action="store_true",
        default=None,
        help="f electron for GGA+U.",
    )
    parser.add_argument(
        "-fd",
        "--fermiDirac",
        type=int,
        default=-1,
        help="F-D smearing for electron enthalpy.",
    )


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Python based first-principles calculations tools",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-v",
        "--version",
        help="print version",
        action="version",
        version=__version__,
    )

    subparser = parser.add_subparsers(title="Valid subcommands", dest="command")
    vasp_parsers = []
    common_parsers = []

    keep_parser = subparser.add_parser(
        "keep",
        help="Keep some files, delete others.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    keep_parser.add_argument(
        "--keepFiles",
        type=str,
        default=[],
        nargs="+",
        help="Files to keep.",
    )
    keep_parser.add_argument(
        "--task",
        type=str,
        default=None,
        choices=["vasp_relax", "vasp_scf", "vasp_elf", "vasp_md"],
        help="Task to handle.",
    )
    keep_parser.add_argument(
        "--dry-run",
        dest="dry_run",
        action="store_true",
        help="Print files that would be deleted without deleting them.",
    )

    vasp_relax_parser = subparser.add_parser(
        "vasp_relax",
        help="Relax by vasp.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    common_parsers.append((vasp_relax_parser, "vasp_relax"))
    vasp_parsers.append(vasp_relax_parser)
    vasp_relax_parser.add_argument("-p", "--pressure", type=float, default=0.0, help="Pressure")
    vasp_relax_parser.add_argument("--optcell", action="store_true", help="Constrained relax.")

    vasp_scf_parser = subparser.add_parser(
        "vasp_scf",
        help="Scf by vasp.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    common_parsers.append((vasp_scf_parser, "vasp_scf"))
    vasp_parsers.append(vasp_scf_parser)
    vasp_scf_parser.add_argument("-p", "--pressure", type=float, default=0.0, help="Pressure")

    vasp_band_parser = subparser.add_parser(
        "vasp_band",
        help="Band by vasp, default is gga.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    common_parsers.append((vasp_band_parser, "vasp_band"))
    vasp_parsers.append(vasp_band_parser)

    vasp_dos_parser = subparser.add_parser(
        "vasp_dos",
        help="Dos by vasp, default is gga.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    common_parsers.append((vasp_dos_parser, "vasp_dos"))
    vasp_parsers.append(vasp_dos_parser)

    vasp_elf_parser = subparser.add_parser(
        "vasp_elf",
        help="Elf by vasp.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    common_parsers.append((vasp_elf_parser, "vasp_elf"))
    vasp_parsers.append(vasp_elf_parser)

    vasp_ph_parser = subparser.add_parser(
        "vasp_ph",
        help="Ph by vasp + phonopy.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    common_parsers.append((vasp_ph_parser, "vasp_ph"))
    vasp_parsers.append(vasp_ph_parser)
    vasp_ph_parser.add_argument("-k", "--kpoints", type=int, nargs=3, default=[3, 3, 3], help="KPOINTS mesh.")
    vasp_ph_parser.add_argument("-d", "--dim", type=int, nargs=3, default=[2, 2, 2], help="Super cell size.")

    vasp_md_parser = subparser.add_parser(
        "vasp_md",
        help="AIMD by vasp.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    common_parsers.append((vasp_md_parser, "vasp_md"))
    vasp_parsers.append(vasp_md_parser)
    vasp_md_parser.add_argument("-t", "--temperature", type=float, nargs=2, default=[300, 300], help="TEBEG and TEEND for vasp.")
    vasp_md_parser.add_argument("-k", "--kpoints", type=int, nargs=3, default=[1, 1, 1], help="KPOINTS mesh.")
    vasp_md_parser.add_argument("-d", "--dim", type=int, nargs=3, default=[3, 3, 3], help="Super cell size.")

    vasp_cohp_parser = subparser.add_parser(
        "vasp_cohp",
        help="COHP by vasp and lobster.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    common_parsers.append((vasp_cohp_parser, "vasp_cohp"))
    vasp_parsers.append(vasp_cohp_parser)

    qe_relax_parser = subparser.add_parser(
        "qe_relax",
        help="Relax by qe.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    common_parsers.append((qe_relax_parser, "qe_relax"))
    qe_relax_parser.add_argument("-p", "--pressure", type=float, default=0.0, help="Pressure")
    qe_relax_parser.add_argument("--encut", type=int, default=100, help="ENCUT (Ry)")
    qe_relax_parser.add_argument("--spin", action="store_true", help="Add spin (default FM).")
    qe_relax_parser.add_argument("-fd", "--fermiDirac", type=int, default=-1, help="F-D smearing for electron enthalpy.")
    qe_relax_parser.add_argument("-opt", "--optcell", action="store_true", help="Constrained relax.")

    for command_name in ["qescf", "qe_scf"]:
        qescf_parser = subparser.add_parser(
            command_name,
            help="Scf by QE.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        common_parsers.append((qescf_parser, command_name))
        qescf_parser.add_argument("--encut", type=int, default=100, help="ENCUT for QE")
        qescf_parser.add_argument("--spin", action="store_true", help="Add spin (default FM).")
        qescf_parser.add_argument("-p", "--pressure", type=float, default=0.0, help="Pressure")
        qescf_parser.add_argument("-fd", "--fermiDirac", type=int, default=-1, help="F-D smearing for electron enthalpy.")
        qescf_parser.add_argument("--symmetry", type=float, default=None, help="Symmetry tolerance for phonopy preprocessing.")

    dmft_parser = subparser.add_parser(
        "dmft",
        help="Dmft by EDMFT",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    common_parsers.append((dmft_parser, "dmft"))
    dmft_parser.add_argument("-fd", "--fermiDirac", type=int, default=-1, help="F-D smearing for electron enthalpy.")
    dmft_parser.add_argument("--fun", type=str, default="gga", choices=["gga", "lda"], help="Functional.")

    for parser_item, comment in common_parsers:
        add_common_job_args(parser_item, comment)
    for parser_item in vasp_parsers:
        add_vasp_args(parser_item)

    parsed_args = parser.parse_args(argv)
    if parsed_args.command is None:
        parser.print_help()
    return parsed_args


def load_command(command):
    module_name = COMMAND_MODULES.get(command)
    if module_name is None:
        raise RuntimeError(f"unknown command {command}")
    try:
        module = importlib.import_module(f"valkyrie.entrypoints.{module_name}")
    except ModuleNotFoundError as exc:
        if exc.name == f"valkyrie.entrypoints.{module_name}":
            raise RuntimeError(f"unknown command {command}") from exc
        raise
    try:
        return getattr(module, "main")
    except AttributeError as exc:
        raise RuntimeError(f"command {command} has no main()") from exc


def main(argv=None):
    args = parse_args(argv)
    if args.command:
        print(__picture__)
        command = load_command(args.command)
        command(**vars(args))


if __name__ == "__main__":
    main()
