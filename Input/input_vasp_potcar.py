# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 16:40:52 2023

@author: zyj
"""

import subprocess, platform
from Input import get_info

def potcar(pot):
    # Generate POTCAR.
    if pot == "auto":
        if platform.system() == "Linux":
            get_info.check_vaspkit()
            subprocess.run("vaspkit -task 103 > /dev/null", shell = True)
            print("<=> Valkyrie: Auto generate POTCAR by vaspkit.")
        else:
            print("<=> Valkyrie: Check under windows, POTCAR is generated autoly")
    else:
        tag = ""
        for i in pot:
            tag = tag + i + "\n"
        if platform.system() == "Linux":
            get_info.check_vaspkit()
            subprocess.run("echo -e '{}'|vaspkit -task 104 > /dev/null".format(tag), shell = True)
            print("<=> Valkyrie: Generate POTCAR {} by vaspkit.".format(pot))
        else:
            print("<=> Valkyrie: Check under windows, POTCAR is generated by {}".format(pot))
    return 0
