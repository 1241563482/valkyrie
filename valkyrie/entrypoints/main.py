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
    # relax by vasp
    relax_parser = subparser.add_parser(
        "vasp_relax",
        help = "Relax by vasp.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    relax_parser.add_argument(
        "-p",
        "--pressure",
        type = float,
        default = 0.0,
        help = "Pressure"
    )
    relax_parser.add_argument(
        "--encut",
        type = int,
        default = 0,
        help = "ENCUT"
    )
    relax_parser.add_argument(
        "--pot",
        type = str,
        default = "auto",
        nargs = "+",
        help = "POTCAR type, eg Li_sv or V"
    )
    relax_parser.add_argument(
        "--spin",
        action = "store_true",
        help = "Add spin (default FM)."
    )
    relax_parser.add_argument(
        "-ns",
        "--notSub",
        action = "store_true",
        help = "Only generate the input file, not sub."
    )
    relax_parser.add_argument(
        "-fd",
        "--fermiDirac",
        type = int,
        default = -1,
        help = "F-D smearing for electron enthalpy."
    )
    relax_parser.add_argument(
        "--fun",
        type = str,
        default = "gga",
        choices = ["gga", "ggau"],
        help = "Functional."
    )
    relax_parser.add_argument(
        "-u",
        type = str,
        nargs = 2,
        default = None,
        help = "Atom for GGA+U and Ueff."
    )
    relax_parser.add_argument(
        "--fElectron",
        action = "store_true",
        default = None,
        help = "f electron for GGA+U."
    )
    relax_parser.add_argument(
        "--optcell",
        action = "store_true",
        help = "Constrained relax."
    )
    relax_parser.add_argument(
        "-q",
        "--queue",
        type = str,
        default = "9242opa!",
        help = "Node name."
    )
    relax_parser.add_argument(
        "-n",
        "--nodes",
        type = int,
        default = "24",
        help = "Number of cores."
    )
    relax_parser.add_argument(
        "-c",
        "--comment",
        type = str,
        default = "relax",
        help = "Comment."
    )

    
    # scf
    scf_parser = subparser.add_parser(
        "vasp_scf",
        help = "Scf by vasp.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    scf_parser.add_argument(
        "--encut",
        type = int,
        default = 0,
        help = "ENCUT"
    )
    scf_parser.add_argument(
        "--pot",
        type = str,
        default = "auto",
        nargs = "+",
        help = "POTCAR type, eg Li_sv or V"
    )
    scf_parser.add_argument(
        "--spin",
        action = "store_true",
        help = "Add spin (default FM)."
    )
    scf_parser.add_argument(
        "-ns",
        "--notSub",
        action = "store_true",
        help = "Only generate the input file, not sub."
    )
    scf_parser.add_argument(
        "-p",
        "--pressure",
        type = float,
        default = 0.0,
        help = "Pressure"
    )
    scf_parser.add_argument(
        "-fd",
        "--fermiDirac",
        type = float,
        default = -1,
        help = "F-D smearing for electron enthalpy."
    )
    scf_parser.add_argument(
        "--fun",
        type = str,
        default = "gga",
        choices = ["gga", "ggau", "hse"],
        help = "Functional."
    )
    scf_parser.add_argument(
        "-u",
        type = str,
        nargs = 2,
        default = None,
        help = "Atom for GGA+U and Ueff."
    )
    scf_parser.add_argument(
        "--fElectron",
        action = "store_true",
        default = None,
        help = "f electron for GGA+U."
    )
    scf_parser.add_argument(
        "-q",
        "--queue",
        type = str,
        default = "9242opa!",
        help = "Node name."
    )
    scf_parser.add_argument(
        "-n",
        "--nodes",
        type = int,
        default = "24",
        help = "Number of cores."
    )
    scf_parser.add_argument(
        "-c",
        "--comment",
        type = str,
        default = "scf",
        help = "Comment."
    )

    
    # band
    band_parser = subparser.add_parser(
        "band",
        help = "Band by vasp, default is gga.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    band_parser.add_argument(
        "--encut",
        type = int,
        default = 0,
        help = "ENCUT"
    )
    band_parser.add_argument(
        "--pot",
        type = str,
        default = "auto",
        nargs = "+",
        help = "POTCAR type, eg Li_sv or V"
    )
    band_parser.add_argument(
        "-p",
        "--pressure",
        type = float,
        default = 0.0,
        help = "Pressure"
    ) 
    band_parser.add_argument(
        "-ns",
        "--notSub",
        action = "store_true",
        help = "Only generate the input file, not sub."
    )  
    band_parser.add_argument(
        "--spin",
        action = "store_true",
        help = "Add spin (only FM)."
    )
    band_parser.add_argument(
        "--fun",
        type = str,
        default = "gga",
        choices = ["gga", "ggau", "hse"],
        help = "Functional for band."
    )
    band_parser.add_argument(
        "--fElectron",
        action = "store_true",
        help = "f electron for GGA+U."
    )
    band_parser.add_argument(
        "-u",
        type = str,
        nargs = 2,
        help = "Atom for GGA+U and Ueff."
    )
    band_parser.add_argument(
        "-fd",
        "--fermiDirac",
        type = int,
        default = -1,
        help = "F-D smearing for electron enthalpy."
    )
    

    # dos
    dos_parser = subparser.add_parser(
        "dos",
        help = "Dos by vasp, default is gga.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    dos_parser.add_argument(
        "--encut",
        type = int,
        default = 0,
        help = "ENCUT"
    )
    dos_parser.add_argument(
        "--pot",
        type = str,
        default = "auto",
        nargs = "+",
        help = "POTCAR type, eg 'Li_sv' or 'V Se'"
    )
    dos_parser.add_argument(
        "-p",
        "--pressure",
        type = float,
        default = 0.0,
        help = "Pressure"
    ) 
    dos_parser.add_argument(
        "-ns",
        "--notSub",
        action = "store_true",
        help = "Only generate the input file, not sub."
    )  
    dos_parser.add_argument(
        "--spin",
        action = "store_true",
        help = "Add spin (only FM)."
    )
    dos_parser.add_argument(
        "--fun",
        type = str,
        default = "gga",
        choices = ["gga", "ggau", "hse"],
        help = "Functional for dos."
    )
    dos_parser.add_argument(
        "--fElectron",
        action = "store_true",
        help = "f electron for GGA+U."
    )
    dos_parser.add_argument(
        "-u",
        type = str,
        nargs = 2,
        help = "Atom for GGA+U and Ueff."
    )
    dos_parser.add_argument(
        "-fd",
        "--fermiDirac",
        type = int,
        default = -1,
        help = "F-D smearing for electron enthalpy."
    )
    
    
    # elf
    elf_parser = subparser.add_parser(
        "elf",
        help = "Elf by vasp.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    elf_parser.add_argument(
        "--encut",
        type = int,
        default = 0,
        help = "ENCUT"
    )
    elf_parser.add_argument(
        "--pot",
        type = str,
        default = "auto",
        nargs = "+",
        help = "POTCAR type, eg Li_sv or V"
    )
    elf_parser.add_argument(
        "--spin",
        action = "store_true",
        help = "Add spin (default FM)."
    )
    elf_parser.add_argument(
        "-ns",
        "--notSub",
        action = "store_true",
        help = "Only generate the input file, not sub."
    )
    elf_parser.add_argument(
        "-p",
        "--pressure",
        type = float,
        default = 0.0,
        help = "Pressure"
    )
    elf_parser.add_argument(
        "--fun",
        type = str,
        default = "gga",
        choices = ["gga", "ggau"],
        help = "Functional."
    )
    elf_parser.add_argument(
        "-u",
        type = str,
        nargs = 2,
        default = None,
        help = "Atom for GGA+U and Ueff."
    )
    elf_parser.add_argument(
        "--fElectron",
        action = "store_true",
        default = None,
        help = "f electron for GGA+U."
    )
    
    # ph
    ph_parser = subparser.add_parser(
        "vasp_ph",
        help = "Ph by vasp + phonopy.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    ph_parser.add_argument(
        "--encut",
        type = int,
        default = 0,
        help = "ENCUT"
    )
    ph_parser.add_argument(
        "--pot",
        type = str,
        default = "auto",
        nargs = "+",
        help = "POTCAR type, eg Li_sv or V"
    )
    ph_parser.add_argument(
        "--spin",
        action = "store_true",
        help = "Add spin (default FM)."
    )
    ph_parser.add_argument(
        "-ns",
        "--notSub",
        action = "store_true",
        help = "Only generate the input file, not sub."
    )
    ph_parser.add_argument(
        "-fd",
        "--fermiDirac",
        type = int,
        default = -1,
        help = "F-D smearing for electron enthalpy."
    )
    ph_parser.add_argument(
        "-k",
        "--kpoints",
        type = int,
        nargs = 3,
        default = [3, 3, 3],
        help = "KPOINTS mesh."
    )
    ph_parser.add_argument(
        "-d",
        "--dim",
        type = int,
        nargs = 3,
        default = [2, 2, 2],
        help = "Super cell size."
    )
    ph_parser.add_argument(
        "--fun",
        type = str,
        default = "gga",
        choices = ["gga", "ggau", "hse"],
        help = "Functional."
    )
    ph_parser.add_argument(
        "-u",
        type = str,
        nargs = 2,
        default = None,
        help = "Atom for GGA+U and Ueff."
    )
    ph_parser.add_argument(
        "--fElectron",
        action = "store_true",
        default = None,
        help = "f electron for GGA+U."
    )
    ph_parser.add_argument(
        "-q",
        "--queue",
        type = str,
        default = "9242opa!",
        help = "Node name."
    )
    ph_parser.add_argument(
        "-n",
        "--nodes",
        type = int,
        default = "24",
        help = "Number of cores."
    )
    ph_parser.add_argument(
        "-c",
        "--comment",
        type = str,
        default = "ph",
        help = "Comment."
    )

    # md
    # COHP
    cohp_parser = subparser.add_parser(
        "vasp_cohp",
        help = "COHP by vasp and lobster.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    cohp_parser.add_argument(
        "--spin",
        action = "store_true",
        help = "Add spin (default FM)."
    )
    cohp_parser.add_argument(
        "-ns",
        "--notSub",
        action = "store_true",
        help = "Only generate the input file, not sub."
    )
    cohp_parser.add_argument(
        "--pot",
        type = str,
        default = "auto",
        nargs = "+",
        help = "POTCAR type, eg Li_sv / Li_sv In_d"
    )
    cohp_parser.add_argument(
        "--encut",
        type = int,
        default = 0,
        help = "ENCUT"
    )
    cohp_parser.add_argument(
        "-q",
        "--queue",
        type = str,
        default = "9242opa!",
        help = "Node name."
    )
    cohp_parser.add_argument(
        "-n",
        "--nodes",
        type = int,
        default = "24",
        help = "Number of cores."
    )
    cohp_parser.add_argument(
        "-c",
        "--comment",
        type = str,
        default = "cohp",
        help = "Comment."
    )
    cohp_parser.add_argument(
        "--fun",
        type = str,
        default = "gga",
        choices = ["gga", "ggau", "hse"],
        help = "Functional."
    )
    cohp_parser.add_argument(
        "-u",
        type = str,
        nargs = 2,
        default = None,
        help = "Atom for GGA+U and Ueff."
    )
    cohp_parser.add_argument(
        "--fElectron",
        action = "store_true",
        default = None,
        help = "f electron for GGA+U."
    )

    # fermi
    # summary (relax, zpe, electride, free energy)
    
    # ELPH 
    elph_parser = subparser.add_parser(
        "elph",
        help = "Elph by QE.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    elph_parser.add_argument(
        "-f",
        type = str,
        default = "POSCAR",
        help="POSCAR file."
    )
    elph_parser.add_argument(
        "--qmesh",
        type = int,
        nargs = 3,
        default = [4, 4, 4],
        help = "Q mesh."
    )
    elph_parser.add_argument(
        "-k",
        type = int,
        nargs = 3,
        default = [16, 16, 16],
        help = "K mesh."
    )
    elph_parser.add_argument(
        "--encut",
        type = float,
        default = 100.0,
        help = "ENCUT (Ry)."
    )
    elph_parser.add_argument(
        "-ns",
        "--notSub",
        action = "store_true",
        help = "Only generate the input file, not sub."
    )
    

    # qe_relax
    qerelax_parser = subparser.add_parser(
        "qerelax",
        help = "Relax by qe.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    qerelax_parser.add_argument(
        "-p",
        "--pressure",
        type = float,
        default = 0.0,
        help = "Pressure"
    )
    qerelax_parser.add_argument(
        "--encut",
        type = int,
        default = 100,
        help = "ENCUT (Ry)"
    )
    qerelax_parser.add_argument(
        "--spin",
        action = "store_true",
        help = "Add spin (default FM)."
    )
    qerelax_parser.add_argument(
        "-ns",
        "--notSub",
        action = "store_true",
        help = "Only generate the input file, not sub."
    )
    qerelax_parser.add_argument(
        "-fd",
        "--fermiDirac",
        type = int,
        default = -1,
        help = "F-D smearing for electron enthalpy."
    )
    qerelax_parser.add_argument(
        "-opt",
        "--optcell",
        action = "store_true",
        help = "Constrained relax."
    )


    # QE scf
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
        "-ns",
        "--notSub",
        action = "store_true",
        help = "Only generate the input file, not sub."
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


    # dmft
    dmft_parser = subparser.add_parser(
        "dmft",
        help = "Dmft by EDMFT",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    dmft_parser.add_argument(
        "-ns",
        "--notSub",
        action = "store_true",
        help = "Only generate the input file, not sub."
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


    parsed_args = parser.parse_args()
    if parsed_args.command is None:
        parser.print_help()
    return parsed_args


def main():
    args = parse_args()
    dict_args = vars(args)

    if args.command:
        try:
            print(__picture__)
            f = getattr(importlib.import_module('valkyrie.entrypoints.{}'.format(args.command)), "main")
        except:
            raise RuntimeError(f"unknown command {args.command}")
        f(**dict_args)

if __name__ == "__main__":
    main()