import platform, os
import subprocess

comment_startswith = "#BSUB -J" 
core_number_startswith = "#BSUB -n" 
node_name_startswith = "#BSUB -q" 
sub_command = "bsub < "
job_head = """#!/bin/sh
#BSUB -J ph
#BSUB -n 12 
#BSUB -q 9242opa!
"""
environment = """
source activate magus-new
module load ips/2018u4
export I_MPI_ADJUST_REDUCE=3
export I_MPI_FABRICS=shm
"""
run_vasp = r"mpirun -np ${LSB_DJOB_NUMPROC} vasp_std"




if platform.system() == "Linux":
    host_name = subprocess.run("echo $USER", shell=True, capture_output=True)
    host_name = host_name.stdout.strip().decode('utf-8')
    if host_name == "jiansun":
        __shell__ = "/share/home/jiansun/zyj/val/shell_scripts"
        __python__ = "/share/home/jiansun/zyj/val/python_scripts"
    elif host_name == "js_zhuyj":
        __shell__ = "/fsa/home/js_zhuyj/valkyrie/shell_scripts" # To change
        __python__ = "/fsa/home/js_zhuyj/valkyrie/python_scripts" # To change
    elif host_name == "yijiezhu":
        __shell__ = "/home/yijiezhu/valkyrie/shell_scripts"
        __python__ = "/home/yijiezhu/valkyrie/python_scripts"
    __work__ = os.getcwd()
else:
    __shell__ = os.getcwd() + "\\shell_scripts"
    __python__ = os.getcwd() + "\\python_scripts"
    __work__ = os.getcwd()
