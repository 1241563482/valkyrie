# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 20:23:08 2023

@author: zyj
"""

import numpy as np
import math
import sys, os


def get_atom_information():
    with open('ELFCAR', 'r') as file:
        lines = file.readlines()    
    print("Reading ELFCAR...")

    if len(lines) >= 7:
     
        # Get the atom numbers and grid for ELFCAR.
        atom_numbers = lines[6].strip()
        numbers = [int(num) for num in atom_numbers.split() if num.isdigit()]
        total_atom_numbers = sum(numbers)
        print("There are {} atoms in this ELFCAR.".format(total_atom_numbers))
        
        grid = lines[total_atom_numbers + 9]
        [Nx, Ny, Nz] = [int(num) for num in grid.split() if num.isdigit()]
        print("The grids for this ELFCAR is {} * {} * {}.\n".format(Nx, Ny, Nz))
        
        # Get the atom positions.
        print("Reading the atom positions...")
        atom_position_list = np.empty((0, 3), dtype=float)
        for i in range(8, total_atom_numbers+8):
            line = lines[i].strip().split()
            if len(line) == 3:
                x, y, z = map(float, line)
                atom = np.array([[x, y, z]])
                atom_position_list = np.vstack((atom_position_list, atom))
        for i, coord in enumerate(atom_position_list, start=1):
            x, y, z = coord
            print(f'Atom {i}: ({x:.6f}, {y:.6f}, {z:.6f})')    
    else:
        print("Input error, check your ELFCAR.")
        
    print("")
    
    return total_atom_numbers, atom_position_list, Nx, Ny, Nz

def pbc(site):
    if site[0]<0:
        site[0]=site[0]+Nx
    if site[1]<0:
        site[1]=site[1]+Ny
    if site[2]<0:
        site[2]=site[2]+Nz
    if site[0]>=Nx:
        site[0]=site[0]-Nx
    if site[1]>=Ny:
        site[1]=site[1]-Ny
    if site[2]>=Nz:
        site[2]=site[2]-Nz
    return site


def get_elf():
    with open('ELFCAR', 'r') as file:
        for _ in range(total_atom_numbers + 10):
            next(file)
        data_lines = file.readlines()

    # Reshape the ELF data.
    data = []
    for line in data_lines:
        values = line.split()
        for value in values:
            data.append(float(value)) 
    data=np.array(data)
    data=data.reshape((Nz, Ny, Nx))
    data=data.swapaxes(0, 2)
    Max = np.max(data)
    data = data/Max
    
    return data
    

def distance(positions, site):
    limit = 0.2
    PBC = [(-1, 0, 0), (0, -1, 0), (0, 0, -1), (-1, -1, 0), (-1, 0, -1), 
           (0, -1, -1), (-1, -1, -1), (1, 0, 0), (0, 1, 0), (0, 0, 1), 
           (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1), (0, 0, 0)]
    for i in positions:
        point_a = np.array(i)
        neighbors = [point_a + np.array(offset) for offset in PBC]
        for neighbor in neighbors:
            
            dist = np.linalg.norm(neighbor - np.array(site))
            if dist < limit:
                return 0    
    return 1


def correction_x(data,site):
    Rx = round(Nx/10)
    site = pbc([round(site[0]*Nx), round(site[1]*Ny), round(site[2]*Nz)])
    max = data[site[0]][site[1]][site[2]]
    max_x = site[0]
    for i in range(-Rx,Rx):
        x=site[0]+i 
        if x < 0:
            x = x + Rx
        elif x >= Nx:
            x = x - Rx
        if max < data[x][site[1]][site[2]]:
            max = data[x][site[1]][site[2]]
            max_x = x
    site[0] = max_x/Nx
    site[1] = site[1]/Ny
    site[2] = site[2]/Nz
    return site


def correction_y(data,site):
    Ry = round(Ny/10)
    site = pbc([round(site[0]*Nx), round(site[1]*Ny), round(site[2]*Nz)])
    max = data[site[0]][site[1]][site[2]]
    max_y = site[1]
    for i in range(-Ry,Ry):
        y=site[1]+i 
        if y < 0:
            y = y + Ry
        elif y >= Ny:
            y = y - Ry
        if max < data[site[0]][y][site[2]]:
            max = data[site[0]][y][site[2]]
            max_y = y
    site[0] = site[0]/Nx
    site[1] = max_y/Ny
    site[2] = site[2]/Nz
    return site


def correction_z(data,site):
    Rz = round(Nz/10)
    site = pbc([round(site[0]*Nx), round(site[1]*Ny), round(site[2]*Nz)])
    max = data[site[0]][site[1]][site[2]]
    max_z = site[2]
    for i in range(-Rz,Rz):
        z=site[2]+i 
        if z < 0:
            z = z + Rz
        elif z >= Nz:
            z = z - Rz
        if max < data[site[0]][site[1]][z]:
            max = data[site[0]][site[1]][z]
            max_z = z
    site[0] = site[0]/Nx
    site[1] = site[1]/Ny
    site[2] = max_z/Nz
    return site


def get_axis_length():
    
    with open('ELFCAR', 'r') as file:
        lines = file.readlines()  
    [a1,b1,c1] = lines[2].split()
    a1 = round(float(a1),8)
    b1 = round(float(b1),8)
    c1 = round(float(c1),8)
    Lx = math.sqrt(a1**2 + b1**2 + c1**2)
    
    [a2,b2,c2] = lines[3].split()
    a2 = round(float(a2),8)
    b2 = round(float(b2),8)
    c2 = round(float(c2),8)
    Ly = math.sqrt(a2**2 + b2**2 + c2**2)
    
    [a3,b3,c3] = lines[4].split()
    a3 = round(float(a3),8)
    b3 = round(float(b3),8)
    c3 = round(float(c3),8)
    Lz = math.sqrt(a3**2 + b3**2 + c3**2)
    
    A = np.array([a1,b1,c1])
    B = np.array([a2,b2,c2])
    C = np.array([a3,b3,c3])
    cross_product = np.cross(B, C)
    dot_product = np.dot(A, cross_product)
    volume = abs(dot_product)
    
    return Lx, Ly, Lz, volume


def get_radius_x(data,site,min_elf):
    x0 = round(site[0]*Nx)
    [x,y,z] = pbc([round(site[0]*Nx), round(site[1]*Ny), round(site[2]*Nz)])
    while data[x][y][z] > min_elf:
        x = x + 1
        [x,y,z]=pbc([x,y,z])
        if x == x0:
            return 0, Nx-1, length_x
    if x > x0:
        x_right = x
    else:
        x_right = x+Nx
        
    x0 = round(site[0]*Nx)
    [x,y,z] = pbc([round(site[0]*Nx), round(site[1]*Ny), round(site[2]*Nz)])
    while data[x][y][z] > min_elf:
        x = x-1
        [x,y,z]=pbc([x,y,z])
        if x == x0:
            return 0, Nx-1, length_x
    if x < x0:
        x_left = x
    else:
        x_left = x-Nx
        
    return x_left, x_right, (x_right-x_left-2)/Nx*length_x


def get_radius_y(data,site,min_elf):
    y0 = round(site[1]*Ny)
    [x,y,z] = pbc([round(site[0]*Nx), round(site[1]*Ny), round(site[2]*Nz)])
    while data[x,y,z] > min_elf:
        y = y+1
        [x,y,z]=pbc([x,y,z])
        if y == y0:
            return 0, Ny-1, length_y
    if y > y0:
        y_right = y
    else:
        y_right = y+Ny
        
    y0 = round(site[1]*Ny)
    [x,y,z] = pbc([round(site[0]*Nx), round(site[1]*Ny), round(site[2]*Nz)])
    while data[x][y][z] > min_elf:
        y = y-1
        [x,y,z]=pbc([x,y,z])
        if y == y0:
            return 0, Ny-1, length_y
    if y < y0:
        y_left = y
    else:
        y_left = y-Ny
        
    return y_left, y_right, (y_right-y_left-2)/Ny*length_y


def get_radius_z(data,site,min_elf):
    z0 = round(site[2]*Nz)
    [x,y,z] = pbc([round(site[0]*Nx), round(site[1]*Ny), round(site[2]*Nz)])
    while data[x][y][z] > min_elf:
        z = z+1
        [x,y,z]=pbc([x,y,z])
        if z == z0:
            return 0, Nz-1, length_z
    if z > z0:
        z_right = z
    else:
        z_right = z+Nz
        
    z0 = round(site[2]*Nz)
    [x,y,z] = pbc([round(site[0]*Nx), round(site[1]*Ny), round(site[2]*Nz)])
    while data[x][y][z] > min_elf:
        z = z-1
        [x,y,z]=pbc([x,y,z])
        if z == z0:
            return 0, Nz-1, length_z
    if z < z0:
        z_left = z
    else:
        z_left = z-Nz
        
    return z_left, z_right, (z_right-z_left-2)/Nz*length_z


def get_volume(data,min_elf):
    count = 0
    density = V/Nx/Ny/Nz
    for ii in range(x_left,x_right+1):
        for jj in range(y_left, y_right+1):
            for kk in range(z_left, z_right+1):
                # For PBC
                [a, b, c] = pbc([ii, jj, kk])
                # Count for ELF > min_elf
                if data[a, b, c] > max(min_elf) and [a, b, c] not in already_counted:
                    count = count + 1
                    already_counted.append([a, b, c])
    return count*density


def if_electride_candidate(data,site):
    if data[site[0],site[1],site[2]] < criteria_for_position:
        return 0
    offset=[[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1],\
            [1,1,0],[-1,-1,0],[1,0,1],[-1,0,-1],[0,1,1],[0,-1,-1],\
            [1,1,1],[-1,-1,-1]]
    data_compare=[]
    for i in range(14):
        site_new=pbc([site[j] - offset[i][j] for j in range(3)])
        data_compare.append(data[site_new[0],site_new[1],site_new[2]])
    if data[site[0],site[1],site[2]] >= np.max(data_compare):
        # print(site,site[0]/Nx,site[1]/Ny,site[2]/Nz)
        # print(data[site[0],site[1],site[2]],data_compare)
        # print()
        return 1
    else:
        return 0
        
def get_min_elf(data,site):
    x = round(site[0]*Nx)
    y = round(site[1]*Ny)
    z = round(site[2]*Nz)
    max_value = data[x,y,z]
    min_value = []
    min_value.append(min(data[:,y,z]))
    min_value.append(min(data[x,:,z]))
    min_value.append(min(data[x,y,:]))
    min_elf = []
    min_elf.append((max_value - min_value[0])/2)
    min_elf.append((max_value - min_value[1])/2)
    min_elf.append((max_value - min_value[2])/2)
    return min_elf                 
    
    


##################
###### Main ######
##################

if os.path.exists("summary.dat"):
    os.remove("summary.dat")
if os.path.exists("Electride.dat"):
    os.remove("Electride.dat")

# Get the POSCAR information and ELF data.
total_atom_numbers, atom_position_list, Nx, Ny, Nz = get_atom_information()
data = get_elf()
criteria_for_position = 0.8

# Find the sites with ELF > criteria_for_volume.
positions_candidate = [] 
for i in range(Nx):
    for j in range(Ny):
        for k in range(Nz):
            if if_electride_candidate(data, [i,j,k]) == 1:
                positions_candidate.append([i/(Nx-1),j/(Ny-1),k/(Nz-1)])


length_x, length_y, length_z, V = get_axis_length()

# Find the position of the electride.
electride=[]
for i in positions_candidate:
    if distance(atom_position_list, i) == 1:
        i = correction_x(data,i)
        i = correction_y(data,i)
        i = correction_z(data,i)
        atom_position_list = np.vstack((atom_position_list, np.array(i)))
        electride.append(i)
if len(electride)==0:
    print("No electride found in this system!")
    with open("Electride.dat","a") as file:
        print("Summary:\t0 0 0 0", file=file)
    with open("summary.dat", "a") as g:
        print(0, V, 0, file = g)
    sys.exit()


# Calculate the radius and volume of each electride.
print("\nThe x,y,z axis lengths are: (Ansgtrom)\n {}, {}, {}\n".format(length_x,length_y,length_z))
volume=[]
electride_V_checked=[]
already_counted = []
for i in electride:
    criteria_for_volume = get_min_elf(data,i)
    x_left, x_right, radius_x = get_radius_x(data, i, criteria_for_volume[0])
    y_left, y_right, radius_y = get_radius_y(data, i, criteria_for_volume[1])
    z_left, z_right, radius_z = get_radius_z(data, i, criteria_for_volume[2])
    volume_i = get_volume(data, criteria_for_volume)
    print(Nx*Ny*Nz, len(already_counted))
    if volume_i > V/3:  # Exclude the free electron gas.
        continue
    volume.append(volume_i)
    electride_V_checked.append(i)
    print("For electride at {}, Rx, Ry, Rz are:".format(i), "(Ansgtrom)")
    print(radius_x, radius_y, radius_z)
    print("The volume is",volume[-1],"(Ansgtrom^3)\n")

if len(electride_V_checked)==0:
    print("No electride found in this system!")
    with open("Electride.dat","a") as file:
        print("Summary:\t0 0 0 0", file=file)
    with open("summary.dat", "a") as g:
        print(0, V, 0, file = g)
    sys.exit()


print("!!!NOTE!!!\nThe criteria for judging the electride volume is ELF >", max(criteria_for_volume))


if os.path.exists("Electride.dat"):
    os.remove("Electride.dat")

with open("Electride.dat","a") as file:
    print("\t\t\tNumber\t\tTotal_volume\t\tAverage_volume", file = file)
    print("Summary:\t", len(electride_V_checked), '\t\t', sum(volume), '\t\t',np.mean(volume), file = file)
    print("The volumes for each electrides are:", end=" ", file=file)
    for i in volume:
        print(i, end=' ', file=file)
    print("",file=file)
    print("The volume of the cell is {} (Angstrom^3)".format(V), file = file)
    
with open("summary.dat", "a") as g:
    print(sum(volume), V, sum(volume)/V, file = g)    

