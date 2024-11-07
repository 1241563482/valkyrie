# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 21:36:50 2023

@author: Yijie Zhu
"""

import subprocess, os, platform
import shutil
import ase.io


def get_symbols(poscar):
    symbol_list = []
    poscar = poscar
    atoms = ase.io.read(poscar)
    for ii in atoms.symbols:
        if ii not in symbol_list:
            symbol_list.append(ii)
    return symbol_list


def check_vaspkit():
    paths = os.environ["PATH"].split(os.pathsep)
    for path in paths:
        executable_path = os.path.join(path, "vaspkit")
        if os.path.isfile(executable_path) and os.access(executable_path, os.X_OK):
            return 0
    raise FileNotFoundError("Vaspkit not found.")
        

def get_wavecar():
    if os.path.exists("../scf/CHG") and os.path.exists("../scf/CHGCAR")\
        and os.path.exists("../scf/WAVECAR") and os.path.isdir("../scf"):
        shutil.copy2("../scf/CHG", "CHG")
        shutil.copy2("../scf/CHGCAR", "CHGCAR")
        shutil.copy2("../scf/WAVECAR", "WAVECAR")
    else:
        raise FileNotFoundError("Output of scf not found.")        
        
        
def get_encut(encut_input):
    enmax_list = []
    if encut_input == 0:
        if os.path.exists("POTCAR"):
            with open('POTCAR', 'r') as file:
                for line in file:
                    if 'ENMAX' in line:
                        enmax = float(line.split('=')[1].split(';')[0].strip())
                        enmax_list.append(enmax)
            if enmax_list:
                encut_output = max(enmax_list) * 1.5  # ENCUT = 1.5 * max(ENMAX)
            else:
                raise FileNotFoundError("No ENMAX found in POTCAR file.")
        else:
            raise FileNotFoundError("No POTCAR file found.")
        return encut_output
    elif type(encut_input) == float or type(encut_input) == int:
        return encut_input
    else:
        raise ValueError("ENCUT input error.")


def get_input_manually(__work__, cal_path):
    files = os.listdir(__work__)
    for file in files:
        if file.startwith("INCAR") or file.startwith("KPOINTS") or file.startwith("POTCAR")\
        or file.startwith("job") or file.endwith(".sh") or file.endwith(".py"):                                                           
            source_path = os.path.join(__work__, file)
            destination_path = os.path.join(cal_path, file)
            shutil.copy2(source_path, destination_path)
            print("Copy input file from {} to {}".format(source_path, destination_path))

            
            
        




        
        