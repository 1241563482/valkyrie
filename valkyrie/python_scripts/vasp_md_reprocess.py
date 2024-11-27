from MDdata import MDdata
import numpy as np
import os
import matplotlib.pyplot as plt

# read data
msd = MDdata(file='XDATCAR', file_route='./', s=0, e=-1, temperature=900, timestep=1, filetype='vasp', symbol2number=None, smooth_window_size=None)
data = msd.MSD()
lengths = [len(v) for k, v in data.items()]
n = min(lengths)
result = np.column_stack([v[:n] for k, v in data.items()])
shape = result.shape
Atom=[]

for key, values in data.items():
    Atom.append(key)
with open('data.dat', 'w') as f:
    for element in Atom:
        f.write(element + '\t\t\t\t\t')
    f.write('\n')
    np.savetxt(f, result)

# read tag
with open('INCAR', 'r') as f:
    lines = f.readlines()
for line in lines:
    if 'NBLOCK' in line:
        nblock = int(line.split('=')[1].strip())
        print(f"NBLOCK = {nblock}")
    if 'NSW' in line:
        nsw= int(line.split('=')[1].strip())
        print(f"NSW = {nsw}")
    if 'POTIM' in line:
        potim = int(line.split('=')[1].strip())
        print(f"POTIM = {potim}")

# plot
time = np.linspace(0,nsw*potim/1000,int(nsw/nblock-1))
time = time[:len(result)]

for i in range(shape[1]):
    plt.plot(time, result[:, i], label=Atom[i])
plt.ylim(0, np.max(result))
plt.title(os.getcwd()[24:])
plt.legend()

plt.savefig('fig.png',dpi=300)

plt.yscale('log')
plt.ylim(0.01,300)
plt.legend()
plt.savefig('fig-log.png',dpi=300)
