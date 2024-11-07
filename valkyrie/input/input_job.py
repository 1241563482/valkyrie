import sys
sys.path.append('../')
from set_up import comment_startswith
from set_up import core_number_startswith
from set_up import node_name_startswith


def job_head(q, n, comment):
    with open("job", "r") as file:
        lines = file.readlines()
    line_number = 0
    
    for line in lines:
        if comment is not None:
            if line.startswith(comment_startswith) or line.startswith(comment_startswith.lower()):
                lines[line_number] = comment_startswith + " " + comment + "\n"
        if n is not None:
            if line.startswith(core_number_startswith) or line.startswith(core_number_startswith.lower()):
                lines[line_number] = core_number_startswith + " " + n + "\n"
        if q is not None:
            if line.startswith(node_name_startswith) or line.startswith(node_name_startswith.lower()):
                lines[line_number] = node_name_startswith + " " + q + "\n"
        line_number = line_number + 1 
    
    with open("job", 'w') as file:
        file.writelines(lines)
        


