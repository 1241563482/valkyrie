# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 16:49:47 2023

@author: zyj
"""

def kpoints(ka, kb, kc):
    kpoints = """Automatic generation 
0 
Gamma 
{} {} {}
0 0 0
""".format(ka, kb, kc)
    with open ("KPOINTS", "w") as file:
        print(kpoints, file = file)
    return 0