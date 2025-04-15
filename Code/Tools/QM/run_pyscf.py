#!/usr/bin/python3
import pyBrellaSampling.Code.Tools.io as io
import pyBrellaSampling.Code.Tools.classes as classes
import pyBrellaSampling.Code.Tools.QM.pyscf_tools as pyscf_tools

import sys
import os
from pprint import pprint
from pyscf import grad as pygrad
from pyscf import qmmm

def main():
    # inputFilename = sys.argv[1]
    inputFilename = "qmmm_0.input"
    inpFile = io.textRead(inputFilename)
    nat = int(inpFile[0].split()[0])
    ncharges = int(inpFile[0].split()[1])
    atoms = [str]*nat
    charge = int(os.getenv("charge"))
    spin = int(os.getenv("spin"))
    method = os.getenv("method")
    basis = os.getenv("basis")
    for index, line in enumerate(inpFile[1:nat+1]):
        words = line.split()
        atoms[index] = f"{words[3]} {words[0]} {words[1]} {words[2]}"
    charges = [float]*ncharges
    charge_loc = [tuple]*ncharges
    for i, line in enumerate(inpFile[nat+1:]):
        words = line.split()
        charges[i] = float(words[0]) #### WARNING, this should be 3. changed for testing!
        charge_loc[i] = (float(words[1]), float(words[2]), float(words[3]))
    io.textDump(charges, "charges")
    mol = pyscf_tools.genMol(atoms, charge, spin, basis,True)
    if method.casefold() == "hf":
        mf, grad = pyscf_tools.UHF(mol, charges, charge_loc)
        qm_grad = pygrad.uhf(mf)
        # pprint(vars(mf.mm_mol))
    elif method.casefold() == "dm21":
        mf, grad = pyscf_tools.NN_MF(mol, "None", 8, charges, charge_loc)
    else:
        mf, grad = pyscf_tools.DFT(mol, method, "None", True, 8, False, charges, charge_loc)
        qm_grad = pygrad.UKS(mf)
    muliken, dipole = mf.analyze()
    pc_grad = qmmm.QMMMGrad(qm_grad).grad_nuc_mm()

    result = [str]*(nat+ncharges+1)
    result[0] = f"{mf.e_tot} {ncharges}"
    for i in range(nat):
        result[i+1] = f"{grad[i][0]} {grad[i][1]} {grad[i][2]} {muliken[1][i]}"
    for i in range(ncharges):
        result[i+nat+1] = f"{pc_grad[i][0]} {pc_grad[i][1]} {pc_grad[i][2]}"
    io.textDump(result, f"{inputFilename}.result" )
    exit(0)

if __name__ == "__main__":
    main()
