import platform, os
import subprocess

comment_startswith = "#BSUB -J" 
core_number_startswith = "#BSUB -n" 
node_name_startswith = "#BSUB -q" 
sub_command = "bsub <"


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

run_vasp = "mpiexec.hydra vasp_std"
#__run_vasp__ = "mpirun -np \$LSB_DJOB_NUMPROC vasp_std"   # A \ needs to set before the $