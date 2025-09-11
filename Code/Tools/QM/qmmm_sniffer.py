import os
import time
from pyscf import grad as pygrad
import numpy

import pyBrellaSampling.Code.Tools.io as io
import pyBrellaSampling.Code.Tools.classes as classes
import pyBrellaSampling.Code.Tools.QM.pyscf_tools as pyscf_tools
from pyscf import lib, gto, df
inputFilename = "qmmm_0.input"

eh2kcal = 627.509469
grad_fix = -1185.82151

def grad_hcore_mm(qmmm, qm_mol, dm):# Credit to pySCF. This is pulled from v.2.8.0 © Copyright 2025, The PySCF Developers. Pulled purely for compatibility between DM21 and qmmm
    mol = qm_mol
    mm_mol = qmmm.mm_mol

    coords = mm_mol.atom_coords()
    charges = mm_mol.atom_charges()

    intor = 'int3c2e_ip2'
    nao = mol.nao
    max_memory = qmmm.max_memory - lib.current_memory()[0]
    blksize = int(min(max_memory*1e6/8/nao**2/3, 200))
    blksize = max(blksize, 1)
    cintopt = gto.moleintor.make_cintopt(mol._atm, mol._bas,
                                            mol._env, intor)
    g = numpy.empty_like(coords)
    for i0, i1 in lib.prange(0, charges.size, blksize):
        fakemol = gto.fakemol_for_charges(coords[i0:i1])
        j3c = df.incore.aux_e2(mol, fakemol, intor, aosym='s1',
                                comp=3, cintopt=cintopt)
        g[i0:i1] = numpy.einsum('ipqk,qp->ik', j3c * charges[i0:i1], dm).T
    return g

def grad_nuc_mm(qmmm, mol,dm): # Credit to pySCF. This is pulled from v.2.8.0 © Copyright 2025, The PySCF Developers. Pulled purely for compatibility between DM21 and qmmm
        '''Nuclear gradients of the QM-MM nuclear energy
        (in the form of point charge Coulomb interactions)
        with respect to MM atoms.
        '''
        mm_mol = qmmm.mm_mol
        coords = mm_mol.atom_coords()
        charges = mm_mol.atom_charges()
        # g_mm = numpy.zeros_like(coords)
        g_mm = grad_hcore_mm(qmmm, mol, dm)
        for i in range(mol.natm):
            # q1 = qm_charges[i]
            q1 = mol.atom_charge(i)
            r1 = mol.atom_coord(i)
            r = lib.norm(coords -r1, axis=1)
            g_mm -= q1 * numpy.einsum('i,ix,i->ix', charges, coords-r1, 1/r**3)
        return g_mm



def run_qmmm(dm=None):
    start = time.perf_counter()
    inpFile = io.textRead(inputFilename)
    os.remove(inputFilename)
    io.textDump(inpFile, inputFilename.replace("input", "RossTemp"))
    if os.path.isfile(inputFilename+".result"):
        os.remove(inputFilename+".result")
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
        try:
            words = line.split()
            charges[i] = float(words[3]) #### WARNING, this should be 3. changed for testing!
            charge_loc[i] = (float(words[0]), float(words[1]), float(words[2]))
        except IndexError:
            print(f"ERROR: Problem reading charge {i} from input file. Check that there are {ncharges} charges in the input file."  )
            print(f"Line is : {line}")
            charges[i] = 0.0
            charge_loc[i] = (0.0, 0.0, 0.0)
        # print(f"{charge_loc[i]=}")
        # print(f"{charges[i]=}")
    
    mol = pyscf_tools.genMol(atoms, charge, spin, basis,False) # Symmetry breaks this! 
    if method.casefold() == "hf":
        mf, grad = pyscf_tools.UHF(mol, charges, charge_loc)
        # qm_grad = pygrad.uhf(mf)
        # pprint(vars(mf.mm_mol))
    elif method.casefold() == "dm21":
        mf = pyscf_tools.NN_MF(mol, "None", 8, charges, charge_loc, dm0=dm)
        dm0 = mf.make_rdm1()
        grad = pyscf_tools.fdiff_forces(mol, "None", 3, charges, charge_loc, 0.05, dm0)
    else:
        mf, grad = pyscf_tools.DFT(mol, method, "None",False, 3, False, charges, charge_loc)
        # qm_grad = pygrad.UKS(mf)
    # muliken, dipole = mf.analyze()
    muliken = mf.mulliken_pop(verbose=1)
    dm = mf.make_rdm1()
    pc_grad = grad_nuc_mm(mf, mol, dm)

    result = [str]*(nat+ncharges+1)
    result[0] = f"{mf.e_tot*eh2kcal} {ncharges}"
    for i in range(nat):
        result[i+1] = f"{grad[i][0]*grad_fix} {grad[i][1]*grad_fix} {grad[i][2]*grad_fix} {muliken[1][i]}"
    for i in range(ncharges):
        result[i+nat+1] = f"{pc_grad[i][0]*grad_fix} {pc_grad[i][1]*grad_fix} {pc_grad[i][2]*grad_fix}"
    io.textDump(result, f"{inputFilename}.result" )
    io.textDump(result, f"{inputFilename}.result.OLD")
    end = time.perf_counter()
    print(f"INFO: QM calculation took {end - start} s")
    return dm

def main():
    dm = None
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    while True:
        if os.path.isfile(inputFilename):
            dm = run_qmmm()
        if os.path.isfile("./kill") == True:
            exit(0)
        time.sleep(2)

if __name__ == "__main__":
    main()