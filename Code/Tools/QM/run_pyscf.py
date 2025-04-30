#!/usr/bin/python3
# import pyBrellaSampling.Code.Tools.io as io
# import pyBrellaSampling.Code.Tools.classes as classes
# import pyBrellaSampling.Code.Tools.QM.pyscf_tools as pyscf_tools

# import sys
# import os
# from pprint import pprint
# from pyscf import grad as pygrad
# # from gpu4pyscf import qmmm
# from pyscf.qmmm.itrf import _QMMMGrad
# from pyscf import lib
import numpy
import time
import os

# def grad_nuc_mm(qmmm, mol=None): # Credit to pySCF. This is pulled from v.2.8.0 © Copyright 2025, The PySCF Developers. Pulled purely for compatibility between DM21 and qmmm
#         '''Nuclear gradients of the QM-MM nuclear energy
#         (in the form of point charge Coulomb interactions)
#         with respect to MM atoms.
#         '''
#         if mol is None:
#             mol = mol
#         mm_mol = qmmm.base.mm_mol
#         coords = mm_mol.atom_coords()
#         charges = mm_mol.atom_charges()
#         g_mm = numpy.zeros_like(coords)
#         for i in range(mol.natm):
#             q1 = mol.atom_charge(i)
#             r1 = mol.atom_coord(i)
#             r = lib.norm(r1-coords, axis=1)
#             g_mm += q1 * numpy.einsum('i,ix,i->ix', charges, r1-coords, 1/r**3)
#             print
#         return g_mm

def main():
    # # inputFilename = sys.argv[1]
    inputFilename = "qmmm_0.input"

    #### TEST
    # inpFile = io.textRead(inputFilename)
    # nat = int(inpFile[0].split()[0])
    # ncharges = int(inpFile[0].split()[1])
    # atoms = [str]*nat
    # charge = int(os.getenv("charge"))
    # spin = int(os.getenv("spin"))
    # method = os.getenv("method")
    # basis = os.getenv("basis")
    # for index, line in enumerate(inpFile[1:nat+1]):
    #     words = line.split()
    #     atoms[index] = f"{words[3]} {round(float(words[0]),6)} {round(float(words[1]),6)} {round(float(words[2]),6)}"
    # charges = [float]*ncharges
    # charge_loc = [tuple]*ncharges
    # pntchrg_line = [str]*(ncharges+1)
    # pntchrg_line[0] = ncharges
    # for i, line in enumerate(inpFile[nat+1:]):
    #     words = line.split()
    #     charges[i] = float(words[3]) #### WARNING, this should be 3. changed for testing!
    #     charge_loc[i] = (float(words[0]), float(words[1]), float(words[2]))
    #     pntchrg_line[i+1] = f"{words[3]} {words[0]} {words[1]} {words[2]}"
    # io.textDump(pntchrg_line, f"{inputFilename}.pntchrg")
    # mol = pyscf_tools.genMol(atoms, charge, spin, basis,False)
    # pprint(vars(mol))
    # print(mol.atom_coords())
    # print(mol.atom_charges())
    # print(mol.energy_nuc())
    # mf, grad= pyscf_tools.DFT(mol, "PBE", "None", True, 9, False,charges, charge_loc )
    # # if method.casefold() == "hf":
    # #     mf, grad = pyscf_tools.UHF(mol, charges, charge_loc)
    # #     qm_grad = pygrad.uhf(mf)
    # #     # pprint(vars(mf.mm_mol))
    # # elif method.casefold() == "dm21":
    # #     mf, grad = pyscf_tools.NN_MF(mol, "None", 8, charges, charge_loc)
    # # else:
    # #     mf, grad = pyscf_tools.DFT(mol, method, "None",True, 8, False, charges, charge_loc)
    # qm_grad = pygrad.UKS(mf)
    # muliken, dipole = mf.analyze()
    # pc_grad = grad_nuc_mm(qm_grad, mol)

    # result = [str]*(nat+ncharges+1)
    # result[0] = f"{mf.e_tot*627.5} {ncharges}"
    # for i in range(nat):
    #     result[i+1] = f"{grad[i][0]*-1185.82151} {grad[i][1]*-1185.82151} {grad[i][2]*-1185.82151} {muliken[1][i]}"
    # for i in range(ncharges):
    #     result[i+nat+1] = f"{pc_grad[i][0]} {pc_grad[i][1]} {pc_grad[i][2]}"
    # io.textDump(result, f"{inputFilename}.result" )
    while True:
        time.sleep(2)
        if os.path.isfile(inputFilename+".result"):
            exit(0)
    

if __name__ == "__main__":
    main()
