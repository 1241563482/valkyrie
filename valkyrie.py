#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 14:18:04 2023

@author: zyj
"""


from my_parsers import MyParsers
from Calculators import vasp_relax, vasp_scf, vasp_band, vasp_ph, vasp_dos, vasp_elf
from Calculators import qe_elph, qe_relax
from set_up import __shell__, __python__, __work__
import os

valkyrie = """
#--------------------------------------------------#
|                 _  _                  _          |
|  /\   /\  __ _ | || | __ _   _  _ __ (_)  ___    |
|  \ \ / / / _` || || |/ /| | | || '__|| | / _ \   |
|   \ V / | (_| || ||   < | |_| || |   | ||  __/   |
|    \_/   \__,_||_||_|\_\ \__, ||_|   |_| \___|   |
|                           |___/                  |
|                                     by Yijie Zhu |
#--------------------------------------------------#
"""
print(valkyrie)

# main parser
args = MyParsers().parser.parse_args()

if args.command == "relax":
    vasp_relax.relax(args, __shell__, __python__, __work__)
elif args.command == "scf":
    vasp_scf.scf(args, __shell__, __python__, __work__)
elif args.command == "band":
    vasp_band.band(args, __shell__, __python__, __work__)
elif args.command == "dos":
    vasp_dos.dos(args, __shell__, __python__, __work__)
elif args.command == "elf":
    vasp_elf.elf(args, __shell__, __python__, __work__)
elif args.command == "ph":
    vasp_ph.ph(args, __shell__, __python__, __work__)    
elif args.command == "elph":
    qe_elph.elph(args, __shell__, __python__, __work__)
elif args.command == "qerelax":
    qe_relax.relax(args, __shell__, __python__, __work__)


if args.input != None:
    out = []
    for file in os.listdir(args.input):
        os.system("cp {}/{} ./".format(args.input, file))
        out.append(file)
    print("<=> Valkyrie: {} are copied from {}.".format(out, args.input))
