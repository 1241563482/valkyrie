import platform, os
import subprocess

comment_startswith = "#BSUB -J" 
core_number_startswith = "#BSUB -n" 
node_name_startswith = "#BSUB -q" 
sub_command = "bsub <"


job_head = """#!/bin/sh
#BSUB -J Bilibili@yijiezhu
#BSUB -n 114514
#BSUB -q APEX legends never die
"""


environment = """
source activate magus-new
module load ips/2018u4
export I_MPI_ADJUST_REDUCE=3
export I_MPI_FABRICS=shm
"""


run_vasp = "mpiexec.hydra vasp_std"
run_vasp_opt = "mpiexec.hydra vasprelax"
run_vasp_gam = "mpiexec.hydra vasp_gam"
run_pw = "mpiexec.hydra /fsa/home/js_zhuyj/software/QE/qe-7.3.1/bin/pw.x -npool 4 "

qe_pot_path = '/fsa/home/js_zhuyj/mypps/QE/SSSP_1.3.0_PBE_efficiency'
