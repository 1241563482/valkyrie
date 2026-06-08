import os
import shutil

import valkyrie.shell_scripts.job as job


def main():
    dir_list = []
    for file in os.listdir():
        if file.startswith("POSCAR-"):
            file = file[7:]
            dir_list.append(file)
    dir_list = sorted(dir_list)
    for i in dir_list:
        if os.path.exists(i):
            os.chdir(i)
        else:
            os.mkdir(i)
            os.chdir(i)
        shutil.copy2("../INCAR", "INCAR")
        shutil.copy2("../KPOINTS", "KPOINTS")
        shutil.copy2("../POSCAR-{}".format(i), "POSCAR")
        shutil.copy2("../POTCAR", "POTCAR")
        shutil.copy2("../job", "job")
        job.control_job("job", None, None, i)
        job.sub("job")
        os.chdir("../")


if __name__ == "__main__":
    main()
