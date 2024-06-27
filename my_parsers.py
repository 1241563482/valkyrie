# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 18:48:08 2023

@author: zyj
"""

import argparse

class MyParsers():
    def __init__(self):
        # main
        self.parser = argparse.ArgumentParser(
            description="Python based first-principles calculations tools",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.subparser = self.parser.add_subparsers(
            title="Valid subcommands",
            dest="command"
        )
        self.parser.add_argument(
            "-q",
            type = str,
            default = None,
            help = "Node name."
        )
        self.parser.add_argument(
            "-n",
            type = int,
            default = None,
            help = "Number of cores."
        )
        self.parser.add_argument(
            "--comment",
            type = str,
            default = None,
            help = "Comment."
        )
        self.parser.add_argument(
            "--symmetry",
            type = float,
            default = None,
            help = "Find the symmetry before calculations."
        )
        self.parser.add_argument(
            "-i",
            "--input",
            type = str,
            default = None,
            help = "Input folder name for other input fils."
        )



        # relax
        self.relax_parser = self.subparser.add_parser(
            "relax",
            help = "Relax by vasp.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.relax_parser.add_argument(
            "-p",
            "--pressure",
            type = float,
            default = 0.0,
            help = "Pressure"
        )
        self.relax_parser.add_argument(
            "--encut",
            type = int,
            default = 0,
            help = "ENCUT"
        )
        self.relax_parser.add_argument(
            "--pot",
            type = str,
            default = "auto",
            nargs = "+",
            help = "POTCAR type, eg Li_sv or V"
        )
        self.relax_parser.add_argument(
            "--spin",
            action = "store_true",
            help = "Add spin (default FM)."
        )
        self.relax_parser.add_argument(
            "-ns",
            "--not-sub",
            action = "store_true",
            help = "Only generate the input file, not sub."
        )
        self.relax_parser.add_argument(
            "-fd",
            "--fermi-dirac",
            type = int,
            default = None,
            help = "F-D smearing for electron enthalpy."
        )
        self.relax_parser.add_argument(
            "--fun",
            type = str,
            default = "gga",
            choices = ["gga", "ggau"],
            help = "Functional."
        )
        self.relax_parser.add_argument(
            "-u",
            type = str,
            nargs = 2,
            default = None,
            help = "Atom for GGA+U and Ueff."
        )
        self.relax_parser.add_argument(
            "--f-electron",
            action = "store_true",
            default = None,
            help = "f electron for GGA+U."
        )
        
            
        # scf
        self.scf_parser = self.subparser.add_parser(
            "scf",
            help = "Scf by vasp.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.scf_parser.add_argument(
            "--encut",
            type = int,
            default = 0,
            help = "ENCUT"
        )
        self.scf_parser.add_argument(
            "--pot",
            type = str,
            default = "auto",
            nargs = "+",
            help = "POTCAR type, eg Li_sv or V"
        )
        self.scf_parser.add_argument(
            "--spin",
            action = "store_true",
            help = "Add spin (default FM)."
        )
        self.scf_parser.add_argument(
            "-ns",
            "--not-sub",
            action = "store_true",
            help = "Only generate the input file, not sub."
        )
        self.scf_parser.add_argument(
            "-p",
            "--pressure",
            type = float,
            default = 0.0,
            help = "Pressure"
        )
        self.scf_parser.add_argument(
            "-fd",
            "--fermi-dirac",
            type = int,
            default = None,
            help = "F-D smearing for electron enthalpy."
        )
        self.scf_parser.add_argument(
            "--fun",
            type = str,
            default = "gga",
            choices = ["gga", "ggau", "hse"],
            help = "Functional."
        )
        self.scf_parser.add_argument(
            "-u",
            type = str,
            nargs = 2,
            default = None,
            help = "Atom for GGA+U and Ueff."
        )
        self.scf_parser.add_argument(
            "--f-electron",
            action = "store_true",
            default = None,
            help = "f electron for GGA+U."
        )

        
        # band
        self.band_parser = self.subparser.add_parser(
            "band",
            help = "Band by vasp, default is gga.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.band_parser.add_argument(
            "--encut",
            type = int,
            default = 0,
            help = "ENCUT"
        )
        self.band_parser.add_argument(
            "--pot",
            type = str,
            default = "auto",
            nargs = "+",
            help = "POTCAR type, eg Li_sv or V"
        )
        self.band_parser.add_argument(
            "-p",
            "--pressure",
            type = float,
            default = 0.0,
            help = "Pressure"
        ) 
        self.band_parser.add_argument(
            "-ns",
            "--not-sub",
            action = "store_true",
            help = "Only generate the input file, not sub."
        )  
        self.band_parser.add_argument(
            "--spin",
            action = "store_true",
            help = "Add spin (only FM)."
        )
        self.band_parser.add_argument(
            "--fun",
            type = str,
            default = "gga",
            choices = ["gga", "ggau", "hse"],
            help = "Functional for band."
        )
        self.band_parser.add_argument(
            "--f-electron",
            action = "store_true",
            help = "f electron for GGA+U."
        )
        self.band_parser.add_argument(
            "-u",
            type = str,
            nargs = 2,
            help = "Atom for GGA+U and Ueff."
        )
        self.band_parser.add_argument(
            "-fd",
            "--fermi-dirac",
            type = int,
            default = None,
            help = "F-D smearing for electron enthalpy."
        )
        

        # dos
        self.dos_parser = self.subparser.add_parser(
            "dos",
            help = "Dos by vasp, default is gga.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.dos_parser.add_argument(
            "--encut",
            type = int,
            default = 0,
            help = "ENCUT"
        )
        self.dos_parser.add_argument(
            "--pot",
            type = str,
            default = "auto",
            nargs = "+",
            help = "POTCAR type, eg 'Li_sv' or 'V Se'"
        )
        self.dos_parser.add_argument(
            "-p",
            "--pressure",
            type = float,
            default = 0.0,
            help = "Pressure"
        ) 
        self.dos_parser.add_argument(
            "-ns",
            "--not-sub",
            action = "store_true",
            help = "Only generate the input file, not sub."
        )  
        self.dos_parser.add_argument(
            "--spin",
            action = "store_true",
            help = "Add spin (only FM)."
        )
        self.dos_parser.add_argument(
            "--fun",
            type = str,
            default = "gga",
            choices = ["gga", "ggau", "hse"],
            help = "Functional for dos."
        )
        self.dos_parser.add_argument(
            "--f-electron",
            action = "store_true",
            help = "f electron for GGA+U."
        )
        self.dos_parser.add_argument(
            "-u",
            type = str,
            nargs = 2,
            help = "Atom for GGA+U and Ueff."
        )
        self.dos_parser.add_argument(
            "-fd",
            "--fermi-dirac",
            type = int,
            default = None,
            help = "F-D smearing for electron enthalpy."
        )
        
        
        # elf
        self.elf_parser = self.subparser.add_parser(
            "elf",
            help = "Elf by vasp.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.elf_parser.add_argument(
            "--encut",
            type = int,
            default = 0,
            help = "ENCUT"
        )
        self.elf_parser.add_argument(
            "--pot",
            type = str,
            default = "auto",
            nargs = "+",
            help = "POTCAR type, eg Li_sv or V"
        )
        self.elf_parser.add_argument(
            "--spin",
            action = "store_true",
            help = "Add spin (default FM)."
        )
        self.elf_parser.add_argument(
            "-ns",
            "--not-sub",
            action = "store_true",
            help = "Only generate the input file, not sub."
        )
        self.elf_parser.add_argument(
            "-p",
            "--pressure",
            type = float,
            default = 0.0,
            help = "Pressure"
        )
        self.elf_parser.add_argument(
            "--fun",
            type = str,
            default = "gga",
            choices = ["gga", "ggau"],
            help = "Functional."
        )
        self.elf_parser.add_argument(
            "-u",
            type = str,
            nargs = 2,
            default = None,
            help = "Atom for GGA+U and Ueff."
        )
        self.elf_parser.add_argument(
            "--f-electron",
            action = "store_true",
            default = None,
            help = "f electron for GGA+U."
        )
        
        # ph
        self.ph_parser = self.subparser.add_parser(
            "ph",
            help = "Ph by vasp + phonopy.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.ph_parser.add_argument(
            "-f",
            "--folder",
            type = str,
            nargs = "+",
            help="A relaxed folder for ph."
        )
        self.ph_parser.add_argument(
            "--encut",
            type = int,
            default = 0,
            help = "ENCUT"
        )
        self.ph_parser.add_argument(
            "--pot",
            type = str,
            default = "auto",
            nargs = "+",
            help = "POTCAR type, eg Li_sv or V"
        )
        self.ph_parser.add_argument(
            "--spin",
            action = "store_true",
            help = "Add spin (default FM)."
        )
        self.ph_parser.add_argument(
            "-ns",
            "--not-sub",
            action = "store_true",
            help = "Only generate the input file, not sub."
        )
        self.ph_parser.add_argument(
            "-fd",
            "--fermi-dirac",
            type = int,
            default = None,
            help = "F-D smearing for electron enthalpy."
        )
        self.ph_parser.add_argument(
            "-k",
            "--kpoints",
            type = int,
            nargs = 3,
            default = [3, 3, 3],
            help = "KPOINTS mesh."
        )
        self.ph_parser.add_argument(
            "-d",
            "--dim",
            type = int,
            nargs = 3,
            default = [2, 2, 2],
            help = "Super cell size."
        )
        self.ph_parser.add_argument(
            "--fun",
            type = str,
            default = "gga",
            choices = ["gga", "ggau", "hse"],
            help = "Functional."
        )
        self.ph_parser.add_argument(
            "-u",
            type = str,
            nargs = 2,
            default = None,
            help = "Atom for GGA+U and Ueff."
        )
        self.ph_parser.add_argument(
            "--f-electron",
            action = "store_true",
            default = None,
            help = "f electron for GGA+U."
        )

        # md
        # cohp
        # fermi
        # summary (relax, zpe, electride, free energy)
        
        # ELPH 
        self.elph_parser = self.subparser.add_parser(
            "elph",
            help = "Elph by QE.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.elph_parser.add_argument(
            "-f",
            type = str,
            default = "POSCAR",
            help="POSCAR file."
        )
        self.elph_parser.add_argument(
            "--qmesh",
            type = int,
            nargs = 3,
            default = [4, 4, 4],
            help = "Q mesh."
        )
        self.elph_parser.add_argument(
            "-k",
            type = int,
            nargs = 3,
            default = [16, 16, 16],
            help = "K mesh."
        )
        self.elph_parser.add_argument(
            "--encut",
            type = float,
            default = 100.0,
            help = "ENCUT (Ry)."
        )
        self.elph_parser.add_argument(
            "-ns",
            "--not-sub",
            action = "store_true",
            help = "Only generate the input file, not sub."
        )
        

        # qe_relax
        self.qerelax_parser = self.subparser.add_parser(
            "qerelax",
            help = "Relax by qe.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.qerelax_parser.add_argument(
            "-p",
            "--pressure",
            type = float,
            default = 0.0,
            help = "Pressure"
        )
        self.qerelax_parser.add_argument(
            "--encut",
            type = int,
            default = 100,
            help = "ENCUT (Ry)"
        )
        self.qerelax_parser.add_argument(
            "--spin",
            action = "store_true",
            help = "Add spin (default FM)."
        )
        self.qerelax_parser.add_argument(
            "-ns",
            "--not-sub",
            action = "store_true",
            help = "Only generate the input file, not sub."
        )
        self.qerelax_parser.add_argument(
            "-fd",
            "--fermi-dirac",
            type = int,
            default = None,
            help = "F-D smearing for electron enthalpy."
        )
