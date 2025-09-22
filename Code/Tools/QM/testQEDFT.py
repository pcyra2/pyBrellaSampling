import numpy as np
# from fit_functional import *
from pyscf import gto, scf, dft, ci, fci
# from pyscf.pbc.tools.pyscf_ase import ase_atoms_to_pyscf,pyscf_to_ase_atoms
from pyscf.tools import cubegen
# from ase.io import read
import density_functional_approximation_dm21 as dm21
import json
import tensorflow as tf
import os
from pprint import pprint
import pyBrellaSampling.Code.Tools.QM.pyscf_tools as pyscf_tools
import matplotlib.pyplot as plt

gpus = tf.config.list_physical_devices('GPU')
if gpus:
    try:    
        # Currently, memory growth needs to be the same across GPUs
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
    except RuntimeError as e:
        # Memory growth must be set before GPUs have been initialized
        print(e)

DataFile = "Water_Energies.json"

def json_write(data:dict):
    try:
        with open(DataFile, "w") as f:
            json.dump(data, f,  indent="\t",sort_keys=True, )
    except TypeError:
        with open(DataFile, "w") as f:
            json.dump(data, f, indent="\t",)
def json_read():
    def parse_float_keys(dct):
        rval = dict()
        for key, val in dct.items():
            try:
                # Convert the key to an integer
                int_key = float(key)
                # Assign value to the integer key in the new dict
                rval[int_key] = val
            except ValueError:
                # Couldn't convert key to an integer; Use original key
                rval[key] = val
        return rval

    data = {}
    if os.path.isfile(DataFile):
        with open(DataFile, "r") as f:
            data = json.load(f, object_hook=parse_float_keys)
            f.close()
    return data



try:
    data=json_read()
except FileNotFoundError:
    data = {}

pprint(data)
x=np.arange(-0.4,3.0,0.1)
for Delta in x:
    Delta = round(Delta,2)
    if Delta not in data.keys():
        data[Delta] = {}

        #create vector alongwhcih to move H
        vector = np.array([0.0,0.0,0.174])-np.array([0.00,0.770,-0.505])
        print('vector',vector)
        vector=vector/np.linalg.norm(vector)
        print('vector',vector)
        H_pos=np.array([0.00,0.770,-0.505])-vector*Delta
        print('H_pos',H_pos)
    
        mol = gto.Mole()
        mol.verbose = 5
        mol.output = 'out_H2O_'+str(Delta)
        mol.atom =  'O 0.000 0.000 0.174; H 0.000 -0.770 -0.505; H '+str(list(H_pos))
        mol.basis = 'sto3g' # this is a small basis set we would wan tot explore the use of larger basis sets
        mol.spin= 0
        #mol.symmetry = True
        mol.build()
        
        #do some testing stuff
        print('Full CI')
        mf = scf.UHF(mol)
        e = mf.kernel()
        mf.conv_tol = 1e-11
        cisolver = fci.FCI(mol, mf.mo_coeff)
        e, fcivec = cisolver.kernel()

        dm0, dm1 = cisolver.make_rdm12(fcivec=fcivec, norb=mf.mo_coeff.shape[1], nelec=mol.nelectron)
        data[Delta]["FCI"]=e

        # print('LDA, VWN')
        # mf = pyscf_tools.DFT(mol, "LDA, VWN", "None",False, 5, False, None, None,dm0)
        # data[Delta]["PBE"]=mf.e_tot


        # print('PBE, VWN3')
        # mf = pyscf_tools.DFT(mol, "PBE, VWN3", "None",False, 5, False, None, None,dm0)
        # data[Delta]["PBE,VWN3"]=mf.e_tot

        print('PBE, PBE')
        mf = pyscf_tools.DFT(mol, "PBE, PBE", "None",True, 5, False, None, None,dm0)
        data[Delta]["PBE"]=mf.e_tot

        print('B3LYP')
        mf_1 = pyscf_tools.DFT(mol, "B3LYP", "None",True, 5, False, None, None,dm0)
        data[Delta]["B3LYP"]=mf_1.e_tot
        # dm0 = mf.make_rdm1()
        
        # print('DM21')
        # mf = pyscf_tools.NN_MF(mol, "None", 3, None, None, dm0=dm0)
        # data[Delta]["DM21"]=mf.e_tot

        print("QE-DFT")
        mf_2 = pyscf_tools.DFT(mol, "gga_5p", "None",True, 5, False, None, None,dm0)
        data[Delta]["QE-DFT"]=mf_2.e_tot
        json_write(data)
        # pprint(data)



# plt.plot(x,e_pbelda,':',label='PBE|VWN3')
# plt.plot(x,e_pbe,':',label='PBE|PBE')
# plt.plot(x,e_lda,':',color='red',label='LDA|VWN3')
# plt.plot(x,e_fci,'--',color='purple',label='DM21')
# plt.plot(x,e_b3lyp,'--',color='purple',label='B3LYP')

# plt.xlabel("Distance from Equilibrium / ${\AA}$" )
# plt.ylabel("Energy / $E_{h}$")
# #plt.xlim(0.75,2.5)
# plt.ylim(-75.8,-74)

# plt.legend(loc='upper right')

# plt.savefig('fun_eval.png')
# plt.show()
def zero_data(dat:list):
    min_val = min(dat)
    return [i - min_val for i in dat]

def zero_dict(dat:dict)->dict:
    x = list(dat.keys())
    keys = list(dat[x[0]].keys())
    for key in keys:
        minima =  min([dat[i][key] for i in x])

        for i in x:
            dat[i][key] = dat[i][key]-minima
    return dat

data = zero_dict(data)
x = list(data.keys())
print(data)
DataFile = [None]*(len(x)+1)
DataFile[0] = "Distance\tFCI\tPBE\tB3LYP\QEDFT"
for j,i in enumerate(x):#range(len(x)):
    DataFile[j+1] = f"{round(i,2)}\t{round(data[i]['FCI'],6)}\t{round(data[i]['PBE'],6)}\t{round(data[i]['B3LYP'],6)}\t{round(data[i]['QE-DFT'],6)}"

with open("water_dissoc.dat", "w") as f:
    for line in DataFile:
        print(line, file=f)


plt.plot(x,zero_data([data[i]['FCI'] for i in x]),':',label="FCI")
plt.plot(x,zero_data([data[i]['PBE'] for i in x]),':',label="PBE")
plt.plot(x,zero_data([data[i]['B3LYP'] for i in x]),':',label="B3LYP")
plt.plot(x,zero_data([data[i]['QE-DFT'] for i in x]),':',label="QE-DFT")

plt.xlabel("Distance from Equilibrium / ${\AA}$" )
plt.ylabel("Energy / $E_{h}$")
#plt.xlim(0.75,2.5)
# plt.ylim(-75.8,-74)

plt.legend(loc='upper right')

plt.savefig('fun_eval.png')
plt.show()

# fig = plotly.GraphObjects.Figure()
# fig.add_scatter(x=x, y=zero_data([data[i]['FCI'] for i in x]), mode='lines+markers', name='FCI')
# fig.add_scatter(x=x, y=zero_data([data[i]['PBE'] for i in x]), mode='lines+markers', name='PBE')
# fig.add_scatter(x=x, y=zero_data([data[i]['B3LYP'] for i in x]), mode='lines+markers', name='B3LYP')
# fig.add_scatter(x=x, y=zero_data([data[i]['QE-DFT'] for i in x]), mode='lines+markers', name='QE-DFT')
# fig.show()