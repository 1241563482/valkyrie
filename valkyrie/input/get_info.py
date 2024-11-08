import subprocess, os, platform
import shutil
import ase.io


def check_vaspkit():
    paths = os.environ["PATH"].split(os.pathsep)
    for path in paths:
        executable_path = os.path.join(path, "vaspkit")
        if os.path.isfile(executable_path) and os.access(executable_path, os.X_OK):
            return 0
    raise FileNotFoundError("Vaspkit not found.")
        

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
        return round(encut_output, 2)
    elif type(encut_input) == float or type(encut_input) == int:
        return round(encut_input, 2)


def get_input_manually(__work__, cal_path):
    files = os.listdir(__work__)
    for file in files:
        if file.startwith("INCAR") or file.startwith("KPOINTS") or file.startwith("POTCAR")\
        or file.startwith("job") or file.endwith(".sh") or file.endwith(".py"):                                                           
            source_path = os.path.join(__work__, file)
            destination_path = os.path.join(cal_path, file)
            shutil.copy2(source_path, destination_path)
            print("Copy input file from {} to {}".format(source_path, destination_path))

            
            
        




        
        