# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 16:47:33 2020

@author: PSN
"""


import re
from collections import Counter

import ase.io
import ase.neighborlist
import numpy as np
from scipy import optimize, constants
from scipy.special import erfc
from scipy.integrate import simps


class MDdata():
    def __init__(self,file,file_route,s=0,e=-1, temperature=1, timestep=1, filetype="vasp", symbol2number=None, smooth_window_size=None):
        self.file_route = file_route
        self.filetype = filetype
        
        if filetype == "vasp":
            traj = ase.io.read(file, index='%d:%d'%(s,e))
            self.traj=traj
        if filetype == "lammps":
            traj = ase.io.read(file,
                               format="lammps-dump-text", index='%d:%d'%(s,e))
            self.traj=traj
            # self.set_element(symbol2number)
            
            
            
        if filetype == "xyz":
            self.traj = ase.io.read(file, index='%d:%d'%(s,e))
            # c = np.loadtxt(file_route+".cell")
            # x = c.T[2]
            # y = c.T[6]
            # z = c.T[10]
            # for ii in range(len(self.traj)):
            #     self.traj[ii].set_cell([x[ii], y[ii], z[ii]])
                
                
        [ii.set_pbc([0,0,0]) for ii in self.traj]
        self.step = len(self.traj)
        self.atom_number = len(self.traj[0].get_chemical_symbols())
        self.atom_name = self.traj[0].get_chemical_symbols()
        self.atom_mass = self.traj[0].get_masses() * constants.atomic_mass
        self.cell_const = np.zeros((self.step, 3, 3))
        # direct pos
        self.position_d = np.zeros((self.step, self.atom_number, 3))
        self.volume = np.zeros(self.step)
        for time, ii in enumerate(self.traj):
            self.position_d[time] = ii.get_scaled_positions()
            self.cell_const[time] = ii.get_cell()#[:]
            self.volume[time] = ii.get_volume()
        self.move_position_pbc()
        # fix mass center
        self.fix_center_of_mass()
        # Cartesian pos
        self.position = self.position_d @ self.cell_const
        # velocity calculation unit : a / fs
        self.index = self.get_index()
        self.temperature = temperature
        self.time_step = timestep
        # read sth from INCAR
        # with open(self.file_route+r'\INCAR', 'r') as file:
        #     text = file.read()
        # self.temperature = float(re.findall(r'TEBEG =\s*(\d+)\s*', text)[0])
        # self.time_step = float(re.findall(r'POTIM =\s*(\d+)\s*', text)[0])
        self.beta = 1/(self.temperature * constants.Boltzmann)
        self.velocity = self.position[2:, :, :] - self.position[:-2, :, :]
        self.velocity /= 2*self.time_step

    def set_element(self, type2number):
        '''
        set proper chemical symbol for lammpstrj
        '''
        
        wrong_list = self.traj[0].get_atomic_numbers()
        symbols = self.traj[0].get_chemical_symbols()
        for key, value in type2number.items():
            for index, jj in enumerate(wrong_list):
                if jj == value+1:
                    symbols[index] = key
        for ii in self.traj:
            ii.set_chemical_symbols(symbols)

    def move_position_pbc(self):
        '''
        If some atoms move to another side of the cell (pbc), move them back.
        '''
        for time in range(self.step-1):
            self.position_d[time+1] -= \
            np.rint(self.position_d[time+1] - self.position_d[time])
    
    def fix_center_of_mass(self):
        total_mass = self.atom_mass.sum()
        center_of_mass = (self.position_d * self.atom_mass[None, :, None]).sum(axis=1)/total_mass
        self.position_d -= (center_of_mass-center_of_mass[0])[:, None, :]

    def get_index(self):
        '''
        Return a dict, key is element name, value is its list of index.
        For ndarray slicing.
        '''
        ret = {}
        for ii in Counter(self.atom_name):
            ret[ii] = [index for index, name in enumerate(self.atom_name) \
                       if name == ii]
        return ret

    def get_force(self):
        traj = ase.io.read(self.file_route+r'\OUTCAR', index=':')
        force = np.zeros((self.step, self.atom_number, 3))
        if len(traj) > self.step:
            traj =  traj[0:self.step]
        for time, ii in enumerate(traj):
            force[time] = ii.get_forces()
        return force

    def VDOS(self, multi=False):
        nframe = self.velocity.shape[0]
        spectrum = np.fft.rfft(self.velocity, axis=0)
        spectrum = self.atom_mass[None, :, None] * abs(spectrum)**2
        x = np.fft.rfftfreq(nframe, d=self.time_step*1e-15)
        if multi == False:
            spectrum = spectrum.sum(axis=2).sum(axis=1)
            spectrum -= min(spectrum)
            spectrum *= 2e-5*self.beta*self.time_step / (self.atom_number*nframe)
            return [x, spectrum]
        else:
            spectrum = spectrum.sum(axis=2)
            ret = {}
            for ii in self.index:
                s = spectrum[:, self.index[ii]].sum(axis=1)
                s -= min(s)
                s *= 2e-5*self.beta*self.time_step / (len(self.index[ii])*nframe)
                ret[ii] = [x, s]
            return ret

    def VACF(self):
        x, spectrum = self.VDOS()
        vacf = np.fft.irfft(spectrum)/(2*self.beta)
        t = np.arange(0, vacf.shape[0], 1)
        t = t/x[-1]
        return [t, vacf]

    def twoPT(self, multi=False, freq=None, spectrum=None, 
              density=None, mass=None, number=None):
        '''

        Use 2PT-MF method
        Must eliminate the high freq part of the spectrum

        '''        
        n = density
        if multi == False:
            n = self.atom_number/self.volume.mean()
            mass = self.atom_mass[0]
            number = self.atom_number
            freq, spectrum = self.VDOS()
        for ii in np.arange(0, spectrum.shape[0], 1):
            if freq[ii] >= 1e14:
                spectrum[ii] = 0
        df = freq[1] - freq[0]
        # normalized diffusivity
        delta = 1e10 * 2*spectrum[0]/9 * (6/np.pi)**(2/3) *\
                (np.pi/(mass*self.beta))**0.5 * n**(1/3)
        gamma = optimize.fsolve(lambda x: 2*(1-x)**3/(2-x) - x**0.4*delta**0.6,
                                0)[0]
        f = gamma**0.4*delta**0.6
        alpha = 12*f/spectrum[0]
        # twoPT_S_HS = 12*f*alpha/(alpha**2 + 4*(np.pi*freq)**2)
        # twoPT_S_solid = spectrum - twoPT_S_HS
        def W_entropy_solid(nu):
            a = self.beta * constants.Planck * nu
            return a/(np.exp(a)-1) - np.log(1-np.exp(-a))
        # IG_component = 2.5+np.log((2*np.pi*self.atom_mass[0] / \
        #                           (constants.Planck**2*self.beta))**1.5 * \
        #                           1e-30/(f*n))
        # W_entropy_HS = IG_component
        # W_entropy_HS += np.log((1+gamma+gamma**2-gamma**3)/(1-gamma)**3)
        # W_entropy_HS += gamma*(3*gamma-4)/(1-gamma)**2
        # W_entropy_HS /= 3
        # entropy_HS = constants.Boltzmann*twoPT_S_HS.sum()*df*W_entropy_HS
        # the first term in freq is 0, can not be included, or divide by zero
        # entropy_solid = constants.Boltzmann*\
        # (twoPT_S_solid[1:]*W_entropy_solid(freq[1:])).sum()*df

        # 2PT-MF
        M2 = simps((2*np.pi*freq)**2 * spectrum, dx=df)/3
        M4 = simps((2*np.pi*freq)**4 * spectrum, dx=df)/3

        def func_for_A_B_f(B):
            A = 4*B/(2+(np.pi*(1+4*B/alpha**2))**0.5)
            f = A*spectrum[0]/(12*(4*B/np.pi)**0.5)
            ret = B - (f*(2*M2*A-M4-A**2)+M4-M2**2)/(2*f*(1-f)*A)
            return ret
        B = optimize.fsolve(func_for_A_B_f, 1e28)[0]
        A = 4*B/(2+(np.pi*(1+4*B/alpha**2))**0.5)
        MF_f = A*spectrum[0]/(12*(4*B/np.pi)**0.5)

        def K(s):
            ret = A * (np.pi/(4*B))**0.5
            ret *= np.exp(s**2/(4*B))
            ret *= erfc(s/(2*B**0.5))
            return ret
        
        IG_component = 2.5+np.log((2*np.pi*mass / \
                                  (constants.Planck**2*self.beta))**1.5 * \
                                  1e-30/(MF_f*n))
        W_entropy_HS = IG_component
        W_entropy_HS += np.log((1+gamma+gamma**2-gamma**3)/(1-gamma)**3)
        W_entropy_HS += gamma*(3*gamma-4)/(1-gamma)**2
        W_entropy_HS /= 3

        temp = 1j*2*np.pi*freq
        MF_S_gas = np.abs(MF_f*6*(1/(K(temp) + temp) + 1/(K(-temp) - temp)))
        MF_S_solid = spectrum - MF_S_gas
        MF_entropy_gas = constants.Boltzmann*MF_S_gas.sum()*df*W_entropy_HS
        MF_entropy_solid = constants.Boltzmann*\
                           (MF_S_solid[1:]*W_entropy_solid(freq[1:])).sum()*df
        MF_entropy = MF_entropy_gas + MF_entropy_solid
        
        if self.filetype == "vasp":
            E = self.getPTE('E')*constants.eV/number
        elif self.filetype == "lammps":
            E = self.getPTE_lmp('E')*constants.eV/number
        else:
            E = 0
        # PV = x.getPTE('P')*1e5*x.volume.mean()*1e-30/x.atom_number
        TS = MF_entropy*self.temperature
        F = (E - TS)/constants.eV
        return {"MF_S_gas":MF_S_gas,
                "MF_S_solid":MF_S_solid,
                "MF_entropy_gas":MF_entropy_gas,
                "MF_entropy_solid": MF_entropy_solid,
                "MF_entropy":MF_entropy,
                "MF_f": MF_f,
                "free_energy": F,
                "temperature": self.temperature}
    
    def multi_twoPT(self):
        l = {}
        d = {}
        VDOS = self.VDOS(multi=True)
        r = {'H': 53, 
             'O': 48,}
        total = 0
        for ii in self.index:
            total += len(self.index[ii])*r[ii]**3            
        for ii in self.index:
            d[ii] = total/(self.volume.mean()*r[ii]**3)
            f, s = VDOS[ii]
            a = (self.twoPT(multi=True, freq=f, spectrum=s, 
                            density=d[ii], mass=self.atom_mass[self.index[ii][0]], 
                            number=len(self.index[ii])))
            l[ii] = a
        return l

    def getPTE(self, target='pte'):
        '''
        get free energy, P and T from OUTCAR
        '''
        ret = []
        with open(self.file_route+r'\OUTCAR', 'r') as file:
            text = file.read()
        for ii in target:
            if ii in ['p', 'P']:
                a = re.findall(r'total pressure\s+=\s+(\d+\.\d+)', text)
            if ii in ['t', 'T']:
                a = re.findall(r'temperature\s+(\d+\.\d+)', text)
            if ii in ['e', 'E']:
                a = re.findall(r'free  energy   TOTEN\s+=\s+(\d+\.\d+)', text)
            ret.append(np.array(a, dtype='float').mean())
        if len(target) == 1:
            ret = ret[0]
        return ret
    
    def getPTE_lmp(self, target='pte'):
        ret = []
        with open(self.file_route+r'\log.lammps', 'r') as file:
            text = file.read()
        for ii in target:
            if ii in ['p', 'P']:
                a = re.findall(r'Press\s+=\s+(\d+\.\d+)', text)
            if ii in ['t', 'T']:
                a = re.findall(r'Temp\s+=\s+(\d+\.\d+)', text)
            if ii in ['e', 'E']:
                a = re.findall(r'PotEng\s+=\s+(\d+\.\d+)', text)
            ret.append(np.array(a, dtype='float').mean())
        if len(target) == 1:
            ret = ret[0]
        return ret

    # get T from velocity
    # def temperature(self):
    #     E = 0.5*self.velocity**2
    #     E = E.sum(axis=2)
    #     E *= self.atom_mass
    #     E = E.mean(axis=1)
    #     T = E*1e10 / (1.5*constants.Boltzmann)
    #     return T

    def MSD(self, start=0, end=None):
        '''
        Mean-square displacement of different elements.
        The return is a dict, key is chemical symbal.
        Three axes are seperated, add sum(axis=1) to get total msd
        '''
        d = (self.position[start:end, :, :] - self.position[start, :, :])**2
        ret = {}
        for name, l in self.index.items():
            #ret[name] = d[:, l, :].mean(axis=1)
            ret[name] = d[:, l, :].mean(axis=1).sum(axis=1)
        return ret
    
    def van_Hove_coorelation_self(self, ele_index, start, end, dt_max=1000):
        ret = np.zeros([dt_max, 200])
        for dt in range(dt_max):
            dr_all = []
            for t0 in range(start, end, 20):
                dr = self.position[dt+t0, ele_index, :] - self.position[t0, ele_index, :] 
                dr = np.sqrt((dr**2).sum())
                dr_all.append(dr)
            vh_y, bin_edge = np.histogram(dr_all, bins=200, range=(0, 4))
            bin_width = bin_edge[1] - bin_edge[0]
            vh_x = bin_edge[1:] - bin_width*0.5
            # vh_y = vh_y/(4*np.pi*vh_x**2*bin_width)
            vh_y = vh_y/bin_width
            vh_y = vh_y/len(dr_all)
            ret[dt] = vh_y
        return ret

    def CN(self, e1, e2, cutoff, frame=-1):
        n = ase.neighborlist.neighbor_list('i', self.traj[frame], {(e1, e2): cutoff})
        cn = np.bincount(np.bincount(n, minlength=self.atom_number)[self.index[e1]], minlength=12)
        return cn/cn.sum()
    
    def RDF(self, e1, e2, r_max=5, frame=-1):
        a = self.traj[frame]
        d = ase.neighborlist.neighbor_list('d', a, {(e1, e2): r_max})
        h, bin_edges = np.histogram(d, bins=500, range=(0, r_max))
        rdf = h/(4*np.pi/3*(bin_edges[1:]**3 - bin_edges[:-1]**3)) * a.get_volume()/(len(self.index[e1])*len(self.index[e2]))
        if e1 != e2:
            rdf /= 2
        return bin_edges[0:-1], rdf
    
    def avg_CN(self, e1, e2, cutoff, dt=100):
        ret = 0
        sample_range = range(0, self.step, dt)
        for ii in sample_range:
            ret += self.CN(e1, e2, cutoff, frame=ii)
        return ret/len(sample_range)
    
    def avg_RDF(self, e1, e2, r_max=5, dt=100):
        retx, rety = 0, 0
        sample_range = range(0, self.step, dt)
        for ii in sample_range:
            rdf = self.RDF(e1, e2, r_max, frame=ii)
            retx = rdf[0]
            rety += rdf[1]
        return retx, rety/len(sample_range)
    
    def adf(self, e1, e2, cutoff, sample_range=[-1, ]):
        x = np.arange(0, 180,0.5)
        rety = 0
        for t in sample_range:
            degs = list(self.get_angles(self.traj[t], {(e1, e2): cutoff}))
            y, _ = np.histogram(degs, bins=360, range=(0, 180), density=True)
            rety += y
        return x, rety/len(sample_range)
        
    def angle(self, e1, e2, dt=20):
        rdf = self.position[::dt, : , None, :]-self.position[::dt, None, :, :]
        index = np.argsort(np.sqrt((rdf**2).sum(axis=3)), axis=2)
        vector = np.take_along_axis(rdf, index[..., None], axis=2)[..., 1:3, :]
        cos = (vector[..., 0, :]*vector[..., 1, :]).sum(axis=2)
        cos = cos/(np.linalg.norm(vector[..., 0, :], axis=2)*\
                   np.linalg.norm(vector[..., 1, :], axis=2))
        angle = np.degrees(np.arccos(cos))
        histy, bin_edge = np.histogram(angle, bins=180, range=(0, 180))
        bin_width = bin_edge[1] - bin_edge[0]
        retx = bin_edge[1:] - bin_width*0.5
        return retx, histy/(angle.size*bin_width)
    
    def smooth_traj(self, window_size=500):
        window = np.full([window_size], 1/window_size)[:, None, None]
        return convolve(self.position, window, mode="nearest")
    
    def distance_distribution(self, e1, e2, cutoff=5, frame=-1):
        i, d = ase.neighborlist.neighbor_list('id', self.traj[frame], {(e1, e2): cutoff})
        ds, ret = [], []
        for ii in self.index[e1]:
            ds.append(np.sort(d[i==ii])[0:12])
        ds = np.array(ds)
        # print(ds.shape)
        for ii in range(ds.shape[1]):
            y, x = np.histogram(ds[:, ii], bins=100, range=(0.5, 1.5))
            ret.append(y)
        retx = x[1:] - (x[1]-x[0])*0.5
        return retx, ret
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

