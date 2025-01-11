import argparse, importlib, logging, os, subprocess
from .. import __picture__, __version__
from ..logger import set_logger

def parse_args():
    parser = argparse.ArgumentParser(
        description="Python based first-principles calculations tools",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-v",
        "--version",
        help = "print version",
        action = 'version',
        version = __version__
    )

    parser_log = argparse.ArgumentParser(
        add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser_log.add_argument(
        "-ll",
        "--log-level",
        choices = ["DEBUG", "INFO", "WARNING", "ERROR"],
        default = "INFO",
        help = "set verbosity level by strings: ERROR, WARNING, INFO and DEBUG",
    )
    parser_log.add_argument(
        "-lp",
        "--log-path",
        type = str,
        default = "log.txt",
        help = "set log file to log messages to disk",
    )



    subparser = parser.add_subparsers(title="Valid subcommands", dest="command")
    subparserList = []

    keep_parser = subparser.add_parser(
        "keep",
        help = "Keep some files, delete others.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    keep_parser.add_argument(
        "--keepFiles",
        type = str,
        default = [],
        nargs = "+",
        help = "Files to keep."
    ) 
    keep_parser.add_argument(
        "--task",
        type = str,
        default = None,
        choices = ["vasp_relax", "vasp_scf", "vasp_elf", "vasp_md"],
        help = "Task to handle."
    )



    # Relax by vasp
    vasp_relax_parser = subparser.add_parser(
        "vasp_relax",
        help = "Relax by vasp.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparserList.append(vasp_relax_parser)
    vasp_relax_parser.add_argument(
        "-p",
        "--pressure",
        type = float,
        default = 0.0,
        help = "Pressure"
    )
    vasp_relax_parser.add_argument(
        "--optcell",
        action = "store_true",
        help = "Constrained relax."
    )


    
    # Scf by vasp
    vasp_scf_parser = subparser.add_parser(
        "vasp_scf",
        help = "Scf by vasp.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparserList.append(vasp_scf_parser)
    vasp_scf_parser.add_argument(
        "-p",
        "--pressure",
        type = float,
        default = 0.0,
        help = "Pressure"
    )
    

    # Band by vasp
    vasp_band_parser = subparser.add_parser(
        "vasp_band",
        help = "Band by vasp, default is gga.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparserList.append(vasp_band_parser)
    

    # DOS by vasp
    vasp_dos_parser = subparser.add_parser(
        "vasp_dos",
        help = "Dos by vasp, default is gga.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparserList.append(vasp_dos_parser)
    
    
    # ELF by vasp
    vasp_elf_parser = subparser.add_parser(
        "vasp_elf",
        help = "Elf by vasp.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparserList.append(vasp_elf_parser)


    # Ph by vasp + phonopy
    vasp_ph_parser = subparser.add_parser(
        "vasp_ph",
        help = "Ph by vasp + phonopy.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparserList.append(vasp_ph_parser)
    vasp_ph_parser.add_argument(
        "-k",
        "--kpoints",
        type = int,
        nargs = 3,
        default = [3, 3, 3],
        help = "KPOINTS mesh."
    )
    vasp_ph_parser.add_argument(
        "-d",
        "--dim",
        type = int,
        nargs = 3,
        default = [2, 2, 2],
        help = "Super cell size."
    )


    # MD by vasp
    vasp_md_parser = subparser.add_parser(
        "vasp_md",
        help = "AIMD by vasp.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparserList.append(vasp_md_parser)
    vasp_md_parser.add_argument(
        "-t",
        "--temperature",
        type = float,
        nargs = 2,
        default = [300, 300],
        help = "TEBEG and TEEND for vasp."
        )
    vasp_md_parser.add_argument(
        "-k",
        "--kpoints",
        type = int,
        nargs = 3,
        default = [1, 1, 1],
        help = "KPOINTS mesh."
    )
    vasp_md_parser.add_argument(
        "-d",
        "--dim",
        type = int,
        nargs = 3,
        default = [3, 3, 3],
        help = "Super cell size."
    )


    # COHP by vasp and lobster
    vasp_cohp_parser = subparser.add_parser(
        "vasp_cohp",
        help = "COHP by vasp and lobster.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparserList.append(vasp_cohp_parser)

    
    ## ELPH by QE
    #elph_parser = subparser.add_parser(
    #    "elph",
    #    help = "Elph by QE.",
    #    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    #)
    #subparserList.append(elph_parser)

    #elph_parser.add_argument(
    #    "-f",
    #    type = str,
    #    default = "POSCAR",
    #    help="POSCAR file."
    #)
    #elph_parser.add_argument(
    #    "--qmesh",
    #    type = int,
    #    nargs = 3,
    #    default = [4, 4, 4],
    #    help = "Q mesh."
    #)
    #elph_parser.add_argument(
    #    "-k",
    #    type = int,
    #    nargs = 3,
    #    default = [16, 16, 16],
    #    help = "K mesh."
    #)
    #elph_parser.add_argument(
    #    "--encut",
    #    type = float,
    #    default = 100.0,
    #    help = "ENCUT (Ry)."
    #)
    

    # Relax by QE
    qe_relax_parser = subparser.add_parser(
        "qe_relax",
        help = "Relax by qe.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparserList.append(qe_relax_parser)
    qe_relax_parser.add_argument(
        "-p",
        "--pressure",
        type = float,
        default = 0.0,
        help = "Pressure"
    )
    qe_relax_parser.add_argument(
        "--encut",
        type = int,
        default = 100,
        help = "ENCUT (Ry)"
    )
    qe_relax_parser.add_argument(
        "--spin",
        action = "store_true",
        help = "Add spin (default FM)."
    )
    qe_relax_parser.add_argument(
        "-fd",
        "--fermiDirac",
        type = int,
        default = -1,
        help = "F-D smearing for electron enthalpy."
    )
    qe_relax_parser.add_argument(
        "-opt",
        "--optcell",
        action = "store_true",
        help = "Constrained relax."
    )


    # Scf by QE
    qescf_parser = subparser.add_parser(
        "qescf",
        help = "Scf by QE.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    qescf_parser.add_argument(
        "--encut",
        type = int,
        default = 100,
        help = "ENCUT for QE"
    )
    qescf_parser.add_argument(
        "--spin",
        action = "store_true",
        help = "Add spin (default FM)."
    )
    qescf_parser.add_argument(
        "-p",
        "--pressure",
        type = float,
        default = 0.0,
        help = "Pressure"
    )
    qescf_parser.add_argument(
        "-fd",
        "--fermiDirac",
        type = int,
        default = -1,
        help = "F-D smearing for electron enthalpy."
    )


    # DMFT by EDMFT
    dmft_parser = subparser.add_parser(
        "dmft",
        help = "Dmft by EDMFT",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    dmft_parser.add_argument(
        "-fd",
        "--fermiDirac",
        type = int,
        default = -1,
        help = "F-D smearing for electron enthalpy."
    )
    dmft_parser.add_argument(
        "--fun",
        type = str,
        default = "gga",
        choices = ["gga", "lda"],
        help = "Functional."
    )



    for item in subparserList:
        comment = item.prog.split()[-1]
        item.add_argument(
            "-ns",
            "--notSub",
            action = "store_true",
            help = "Only generate the input file, not sub."
        )
        item.add_argument(
            "-q",
            "--queue",
            type = str,
            default = "9242opa!",
            help = "Node name."
        )
        item.add_argument(
            "-n",
            "--nodes",
            type = int,
            default = "24",
            help = "Number of cores."
        )
        item.add_argument(
            "-c",
            "--comment",
            type = str,
            default = comment,
            help = "Comment."
        )
        if comment.startswith("vasp"):
            item.add_argument(
                "--encut",
                type = int,
                default = 0,
                help = "ENCUT"
            )
            item.add_argument(
                "--pot",
                type = str,
                default = "auto",
                nargs = "+",
                help = "POTCAR type, eg Li_sv or V"
            )
            item.add_argument(
                "--spin",
                action = "store_true",
                help = "Add spin (default FM)."
            )
            item.add_argument(
                "--fun",
                type = str,
                default = "gga",
                choices = ["gga", "ggau", "hse"],
                help = "Functional."
            )
            item.add_argument(
                "-u",
                type = str,
                nargs = 2,
                default = None,
                help = "Atom for GGA+U and Ueff."
            )
            item.add_argument(
                "--fElectron",
                action = "store_true",
                default = None,
                help = "f electron for GGA+U."
            )
            item.add_argument(
                "-fd",
                "--fermiDirac",
                type = int,
                default = -1,
                help = "F-D smearing for electron enthalpy."
            )

    parsed_args = parser.parse_args()
    if parsed_args.command is None:
        parser.print_help()
    return parsed_args


def main():
    args = parse_args()
    dict_args = vars(args)
    print(dict_args)
    if args.command:
        try:
            print(__picture__)
            f = getattr(importlib.import_module('valkyrie.entrypoints.{}'.format(args.command)), "main")
        except:
            raise RuntimeError(f"unknown command {args.command}")
        f(**dict_args)

if __name__ == "__main__":
    main()