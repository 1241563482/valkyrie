from ase.io import read
import glob


def main():
    file = glob.glob("*struct")[0]
    atoms = read(file)
    arr = atoms.get_chemical_symbols()
    number = [index for index, element in enumerate(arr) if element == "Fe"]

    command = "init_dmft.py -ca "
    for i in number:
        command = command + str(i + 1)
        if i != number[-1]:
            command = command + ","
    command = command + " -ot "
    for i in number:
        command = command + "d"
        if i != number[-1]:
            command = command + ","
    command = command + " -qs "
    for i in number:
        command = command + "2"
        if i != number[-1]:
            command = command + ","
    command = command + " -p 5 -us "
    for i in number:
        command = command + str(i + 1)
        if i != number[-1]:
            command = command + ","
    print(command)


if __name__ == "__main__":
    main()
