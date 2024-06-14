# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 21:16:20 2023

@author: zyj
"""

import numpy as np
import yaml


def get_info(ff):
    with open(ff) as f:
        y = yaml.safe_load(f)
        thermal=y['thermal_properties']
    a=np.zeros([1000,4])
    for i in range(1000):
        t= thermal[i]['temperature']
        s= thermal[i]['entropy'] * 0.0103636 / y["natom"] /1000 
        f= thermal[i]['free_energy'] * 0.0103636 / y["natom"] 
        ts= t*s
        a[i,0]=t
        a[i,1]=s
        a[i,2]=f
        a[i,3]=ts
    return a

data=get_info("thermal_properties.yaml")
np.savetxt('out.dat', data, header = "T\tS\tF\tTS")
