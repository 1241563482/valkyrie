import ase.io
import numpy as np


def main():
    sposcar = ase.io.read('POSCAR')
    #sposcar = make_supercell(sposcar, [[2, 0, 0], [0, 2, 0], [0, 0, 2]])

    elements = sposcar.get_chemical_symbols()
    symbols = []
    for element in elements:
        if element not in symbols:
            symbols.append(element)
    atom_pos = sposcar.get_positions()

    atom_number = {}
    for symbol in symbols:
        index = []
        for j, sym in enumerate(elements):
            if sym == symbol:
                index.append(j)
        atom_number.update({symbol: index})

    distance = []
    atom_pair = []
    atom_type = []
    for i, symbol1 in enumerate(symbols):
        for j, symbol2 in enumerate(symbols):
            if i <= j: # generate posible bond pairs bwteeen different atom types
            # if i < j: # generate all posible bond pairs
                continue
            print(symbol1, symbol2, i, j)
            min_dis = 999
            for aa in atom_number[symbol1]:
                for bb in atom_number[symbol2]:
                    dis = np.linalg.norm(atom_pos[aa] - atom_pos[bb])
                    if dis < min_dis and dis > 0:
                        min_a = aa
                        min_b = bb
                        min_dis = dis
            print(min_dis, min_a, min_b, symbol1, symbol2)

            distance.append(min_dis)
            atom_pair.append([min_a, min_b])
            atom_type.append([symbol1, symbol2])

    bondpairs = atom_pair

    """bondpairs_raw = get_bondpairs(sposcar)
    bondpairs = []
    for bp in bondpairs_raw:
        tmp = list(bp[:2])
        tmp.sort()
        if tmp not in bondpairs and tmp[0] != tmp[1]:
            bondpairs.append(tmp)"""
    with open('labels', 'w') as f:
        for bp in bondpairs:
            f.write(f'{sposcar.get_chemical_symbols()[bp[0]]}[{bp[0]}]-{sposcar.get_chemical_symbols()[bp[1]]}[{bp[1]}]\n')
    with open('lobsterin', 'w') as f:
        f.write('COHPstartEnergy  -10\n')
        f.write('COHPendEnergy  10\n')
        f.write('usebasisset pbeVaspFit2015\n')
        f.write('useRecommendedBasisFunctions\n')
        f.write('gaussianSmearingWidth 0.05\n')
        for bp in bondpairs:
            f.write(f'cohpbetween atom {bp[0]+1} and atom {bp[1]+1}\n')


if __name__ == "__main__":
    main()
