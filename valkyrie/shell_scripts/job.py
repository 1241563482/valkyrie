import sys, os
sys.path.append('../')
from set_up import *

def gen_job(*args, job = "job", job_file = None, task = run_vasp, **kwargs):
    with open(job_file, "r") as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if "__TASK__" in line:
                lines[i] = lines[i].replace("__TASK__", task)
        body = ''.join(lines)
    with open(job, "w") as file:
        print(job_head, file = file)
        print(environment, file = file)
        print(body, file = file)
    
        
def sub(job):
    os.system(sub_command + job)
    

def control_job(job, q, n, comment):
    with open(job, "r") as file:
        lines = file.readlines()
    line_number = 0
    
    for line in lines:
        if comment is not None:
            if line.startswith(comment_startswith) or line.startswith(comment_startswith.lower()):
                lines[line_number] = comment_startswith + " " + str(comment) + "\n"
        if n is not None:
            if line.startswith(core_number_startswith) or line.startswith(core_number_startswith.lower()):
                lines[line_number] = core_number_startswith + " " + str(n) + "\n"
        if q is not None:
            if line.startswith(node_name_startswith) or line.startswith(node_name_startswith.lower()):
                lines[line_number] = node_name_startswith + " " + str(q) + "\n"
        line_number = line_number + 1 
    
    with open(job, 'w') as file:
        file.writelines(lines)

    
if __name__ == "__main__":
    gen_job("job", "/home/yijiezhu/valkyrie/shell_scripts/job_relax")
