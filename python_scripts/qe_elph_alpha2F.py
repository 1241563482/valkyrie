import numpy as np
import os

inputfile="alpha2F.dat-0.10"

with open(inputfile, 'r') as infile, open('output.txt', 'w') as outfile:
    lines = infile.readlines()
    lines = lines[3:]
    for i in range(len(lines)):
        outfile.write(lines[i].rstrip())
        if i%3 == 2:
            outfile.write("\n")    
infile.close()
outfile.close()
with open('output.txt', 'r') as outfile:
    lines = outfile.readlines()
header = np.arange(0, 0.62, 0.02)
column_files = {column: open(f'alpha2F_sigma={column:.2f}.dat', 'w') for column in header}

for line in lines:
    data = line.strip().split('  ')
    for i, column in enumerate(header):
        column_files[column].write(data[i] + '\n')

for file in column_files.values():
    file.close()
os.system("mv alpha2F_sigma=0.00.dat freq.dat")

sigma=np.arange(0.02, 0.6, 0.02)
for i in sigma:
    with open('freq.dat', 'r') as freq_file:
        freq_data = [float(line.strip()) for line in freq_file]
    with open('alpha2F_sigma={:.2f}.dat'.format(i), 'r') as sigma_file:
        sigma_data = [float(line.strip()) for line in sigma_file]

    Lambda_list = []
    Lambda = 0
    dw = freq_data[1]
    for j in range(len(freq_data)):
        if freq_data[j]==0:
            continue
        else:
            Lambda = Lambda + 2*sigma_data[j]/freq_data[j]*dw
            Lambda_list.append(Lambda)
    with open('lambda_sigma={:.2f}.dat'.format(i), 'w') as f:
        f.write(str(0.0)+'\n')
        for k in range(len(Lambda_list)):
            f.write(str(Lambda_list[k])+'\n')
    freq_file.close()
    sigma_file.close()
    f.close()

os.system("rm -f output.txt")
os.system("mkdir -p alpha2F-out")
os.system("mv alpha2F_sigma* lambda_sigma* freq.dat alpha2F-out")
