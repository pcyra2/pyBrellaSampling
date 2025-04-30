import pyscf
from pyscf import qmmm as qmmm
import pyscf.fci as fci
import pyscf.grad
import pyscf.tools.cubegen as cubegen
# from pyscf.symm import msym
from pyscf import cc
from pyscf import lib
from pyscf.geomopt.berny_solver import optimize
from pyscf.hessian import thermo
try:
    from pyscf.qsdopt.qsd_optimizer import QSD
except:
    print("WARNING: QSD optimizer not found... Dont try to perform a TS search")

import libmsym 
import density_functional_approximation_dm21 as dm21

import pyBrellaSampling.Code.Tools.classes as classes

from pprint import pprint

try: 
    from gpu4pyscf.dft.rks import RKS as RKSG
    from gpu4pyscf.dft.uks import UKS as UKSG
    print("INFO: gpu4pyscf is imported successfully")
except ModuleNotFoundError: 
    print("ERROR: gpu4pyscf is not installed... Using standard pyscf instead. GPU Toggle will not work. ")
    from pyscf.dft import RKS as RKSG
    from pyscf.dft import UKS as UKSG

#from pyscf.dft import RKS 
#from pyscf.dft import UKS 

import os
import numpy



pyscf.symm.geom.TOLERANCE = 5e-3 ## Needs to be called outside of the function to be initialised first. 

XC_ALIAS = {
    # Conventional name : name in XC_CODES
    'BLYP'              : 'B88,LYP',
    'BP86'              : 'B88,P86',
    'PW91'              : 'PW91,PW91',
    'PBE'               : 'PBE,PBE',
    'REVPBE'            : 'PBE_R,PBE',
    'PBESOL'            : 'PBE_SOL,PBE_SOL',
    'PKZB'              : 'PKZB,PKZB',
    'TPSS'              : 'TPSS,TPSS',
    'REVTPSS'           : 'REVTPSS,REVTPSS',
    'SCAN'              : 'SCAN,SCAN',
    'RSCAN'             : 'RSCAN,RSCAN',
    'R2SCAN'            : 'R2SCAN,R2SCAN',
    'SCANL'             : 'SCANL,SCANL',
    'R2SCANL'           : 'R2SCANL,R2SCANL',
    'SOGGA'             : 'SOGGA,PBE',
    'BLOC'              : 'BLOC,TPSSLOC',
    'OLYP'              : 'OPTX,LYP',
    'OPBE'              : 'OPTX,PBE',
    'RPBE'              : 'RPBE,PBE',
    'BPBE'              : 'B88,PBE',
    'MPW91'             : 'MPW91,PW91',
    'HFLYP'             : 'HF,LYP',
    'HFPW92'            : 'HF,PW_MOD',
    'SPW92'             : 'SLATER,PW_MOD',
    'SVWN'              : 'SLATER,VWN',
    'MS0'               : 'MS0,REGTPSS',
    'MS1'               : 'MS1,REGTPSS',
    'MS2'               : 'MS2,REGTPSS',
    'MS2H'              : 'MS2H,REGTPSS',
    'MVS'               : 'MVS,REGTPSS',
    'MVSH'              : 'MVSH,REGTPSS',
    'SOGGA11'           : 'SOGGA11,SOGGA11',
    'SOGGA11_X'         : 'SOGGA11_X,SOGGA11_X',
    'KT1'               : 'KT1,VWN',
    'KT2'               : 'GGA_XC_KT2',
    'KT3'               : 'GGA_XC_KT3',
    'DLDF'              : 'DLDF,DLDF',
    'GAM'               : 'GAM,GAM',
    'M06_L'             : 'M06_L,M06_L',
    'M06_SX'            : 'M06_SX,M06_SX',
    'M11_L'             : 'M11_L,M11_L',
    'MN12_L'            : 'MN12_L,MN12_L',
    'MN15_L'            : 'MN15_L,MN15_L',
    'N12'               : 'N12,N12',
    'N12_SX'            : 'N12_SX,N12_SX',
    'MN12_SX'           : 'MN12_SX,MN12_SX',
    'MN15'              : 'MN15,MN15',
    'MBEEF'             : 'MBEEF,PBE_SOL',
    'SCAN0'             : 'SCAN0,SCAN',
    'PBEOP'             : 'PBE,OP_PBE',
    'BOP'               : 'B88,OP_B88',
    # new in libxc-4.2.3
    'REVSCAN'           : 'MGGA_X_REVSCAN,MGGA_C_REVSCAN',
    'REVSCAN_VV10'      : 'MGGA_X_REVSCAN,MGGA_C_REVSCAN_VV10',
    'SCAN_VV10'         : 'MGGA_X_SCAN,MGGA_C_SCAN_VV10',
    'SCAN_RVV10'        : 'MGGA_X_SCAN,MGGA_C_SCAN_RVV10',
    'M05'               : 'HYB_MGGA_X_M05,MGGA_C_M05',
    'M06'               : 'HYB_MGGA_X_M06,MGGA_C_M06',
    'M05_2X'            : 'HYB_MGGA_X_M05_2X,MGGA_C_M05_2X',
    'M06_2X'            : 'HYB_MGGA_X_M06_2X,MGGA_C_M06_2X',
    # extra aliases
    'SOGGA11X'          : 'SOGGA11_X',
    'M06L'              : 'M06_L',
    'M11L'              : 'M11_L',
    'MN12L'             : 'MN12_L',
    'MN15L'             : 'MN15_L',
    'N12SX'             : 'N12_SX',
    'MN12SX'            : 'MN12_SX',
    'M052X'             : 'M05_2X',
    'M062X'             : 'M06_2X',
    'GGA_C_BCGP'    : 'GGA_C_BCGP',
    'LDA'           : 'LDA' ,
    'SLATER'        : 'SLATER' ,
    'VWN3'          : 'VWN3'     ,
    'VWNRPA'        : 'VWNRPA'   ,
    'VWN5'          : 'VWN5'     ,
    'B88'           : 'B88'      ,
    'PBE0'          : 'PBE0'     ,
    'PBE1PBE'       : 'PBE1PBE'  ,
    'OPTXCORR'      : 'OPTXCORR' ,
    'B3LYP'         : 'B3LYP'    ,
    'B3LYPG'        : 'B3LYPG'   ,  # used by Gaussian
    'B3LYP5'        : 'B3LYP5'   , # VWN5 version
    'B3P86'         : 'B3P86'    ,
    'B3P86G'        : 'B3P86G'   ,  # used by Gaussian
    'B3P86V5'       : 'B3P86V5'  , # VWN5 version
    #'O3LYP5'       : 'O3LYP5'  #'.1161*HF + .9262*SLATER + .8133*OPTXCORR, .81*LYP + .19*VWN5',
    #'O3LYPG'       : 'O3LYPG'  #'.1161*HF + .9262*SLATER + .8133*OPTXCORR, .81*LYP + .19*VWNRPA',
    'O3LYP'         : 'O3LYP'    , # in libxc == '.1161*HF + 0.071006917*SLATER + .8133*OPTX, .81*LYP + .19*VWN5', may be erroreous
    'MPW3PW'        : 'MPW3PW'   ,  # VWN5 version
    'MPW3PW5'       : 'MPW3PW5'  ,
    'MPW3PWG'       : 'MPW3PWG'  ,  # used by Gaussian
    'MPW3LYP'       : 'MPW3LYP'  ,  # VWN5 version
    'MPW3LYP5'      : 'MPW3LYP5' ,
    'MPW3LYPG'      : 'MPW3LYPG' ,  # used by Gaussian
    'REVB3LYP'      : 'REVB3LYP' ,  # VWN5 version
    'REVB3LYP5'     : 'REVB3LYP5',
    'REVB3LYPG'     : 'REVB3LYPG',  # used by Gaussian
    'X3LYP'         : 'X3LYP'    ,
    'X3LYPG'        : 'X3LYPG'   ,  # used by Gaussian
    'X3LYP5'        : 'X3LYP5'   ,
    'CAMB3LYP'      : 'CAMB3LYP' ,
    'CAMYBLYP'      : 'CAMYBLYP' ,
    'CAMYB3LYP'     : 'CAMYB3LYP',
    'B5050LYP'      : 'B5050LYP' ,
    'MPW1LYP'       : 'MPW1LYP'  ,
    'MPW1PBE'       : 'MPW1PBE'  ,
    'PBE50'         : 'PBE50'    ,
    'REVPBE0'       : 'REVPBE0'  ,
    'B1B95'         : 'B1B95'    ,
    'TPSS0'         : 'TPSS0'    ,
}

def genMol(atoms, charge: int, spin:int, basis: str, symmetry:bool)->pyscf.M:
    """Generates a pySCF molecule from a .xyz file

    Args:
        path (str|list): Path to file, including extension or list of atoms. This list can either be "atom x y z" or [atom, x, y, z]
        charge (int): Net charge
        spin (int): Net spin in pySCF format (2S)
        basis (str): Basis set to generate the molecule in 
        symmetry (bool): Whether to use symmetry

    Returns:
        mol (pyscf.M): pySCF molecule object. 
    """
    # print(type(atoms))
    if type(atoms) == str:
        if ".xyz" in atoms:
            assert os.path.isfile(atoms), "Coordinate file does not exist."
            mol = pyscf.gto.Mole(atom=atoms, unit="Ang")
        elif ";" in atoms:
            mol = pyscf.gto.Mole(atom=atoms, unit="Ang")
        else:
            raise Exception("Unknown atoms")
    elif type(atoms) == list:
        text = ""
        for atom in atoms:
            if len(atom) == 4:
                text = text + "; "+ str(atom[0]) + " "  + str(atom[1]) + " " + str(atom[2]) + " " + str(atom[3]) 
            else:
                text = text + "; " + atom
        mol = pyscf.gto.Mole(atom=text, unit="Ang")
    elif type(atoms) == classes.molecule:
        text = ""
        for atom in atoms.atoms:
            text += f"{atom.element} {atom.x} {atom.y} {atom.z} ;"
        mol = pyscf.gto.Mole(atom=text, unit="Ang")
    mol.basis = basis
    mol.charge = charge
    mol.spin = spin
    mol.symmetry = symmetry
    
    # mol.unit="Ang"
    mol.build()
    # sym = pyscf.symm.geom.detect_symm(mol._atom)
    # print("INFO Symetry is currently at " +mol.topgroup)
    if symmetry == True:
        # pyscf.symm.geom.TOLERANCE = 1e-3
        try:
            mol = msym.gen_mol_msym(mol, tol=1e-8)
        # except libmsym.main.Error: # If no symmetry can be found. 
        #     print("WARNING: mysm symmetry not found, therfore leaving symmetry as "+str(mol.topgroup))
        #     sym_mol = mol
        except:
            pass
        # mol.build()
    return mol

def UHF(mol: pyscf.M,charges=None, locs=None)->pyscf.scf.UHF:
    """performs UHF on the molecule

    Args:
        mol (pyscf.M): pySCF Molecule object

    Returns:
        HF (pyscf.scf.UHF): a completed pySCF UHF object
    """
    HF = pyscf.scf.UHF(mol)
    HF.max_cycle = 300
    if charges != None:
        mf = qmmm.mm_charge(HF, locs, charges)
        mf.verbose=0
        mf.kernel()
        grad = mf.nuc_grad_method().kernel()
        return mf, grad
    else:
        HF.kernel()
        try:
            HF.analyze()
        except KeyError:
            print("WARNING Cannot analyze with this symmetry group... Sorry!")
        # print(HF.mo_coeff)
        return HF

def RHF(mol: pyscf.M)->pyscf.scf.RHF:
    """performs RHF on the molecule

    Args:
        mol (pyscf.M): pySCF Molecule object

    Returns:
        HF (pyscf.scf.UHF): a completed pySCF RHF object
    """
    HF = pyscf.scf.RHF(mol)
    HF.kernel()
    try:
        HF.analyze()
    except KeyError:
        print("WARNING Cannot analyze with this symmetry group... Sorry!")    
    return HF

def FCI(HF: pyscf.scf.uhf.UHF, UHF: bool, )-> float:
    """Performs FullCI on a pySCF mean field object

    Args:
        HF (pyscf.scf.UHF or pyscf.scf.RHF): completed pySCF meanfield object

    Returns:
        Energy (float): Energy of the system in Hatree
        rdm1 (list): Reduced 1-body density matrix
        rdm2 (list): Reduces 2-body density matrix
    """
    norb = len(HF.mo_coeff[0])*2
    cisolver = fci.FCI(HF)
    Energy, fcivec = cisolver.kernel()
    occ = HF.mo_occ
    if UHF == True:
        alpha = int(numpy.sum(occ[0]))
        beta = int(numpy.sum(occ[1]))
    else: 
        alpha = int(numpy.sum(occ[0])/2)
        beta = alpha
    print(f"{alpha=}, {beta=}")
    print(f"{norb=}")
    # rdm1, rdm2 = cisolver.make_rdm12(fcivec, int(norb), (alpha, beta)) # Currently not working correctly. #TODO: Fix this
    rdm1 = False
    rdm2 = False
    return Energy, rdm1, rdm2

def genCube(HF, Molecule: pyscf.M, UHF:bool, path:str, basis = "molecular") -> None: # pragma: no cover
    """Generates cube files

    Args:
        HF (HF solver): HF solver containing orbitals
        Molecule (pyscf.M): Molecule information containing basis set
        UHF (bool): Whether Restricted or unrestricted.
        path (str): path to save cube files.
    """
    if basis == "molecular":
        coeffs = HF.mo_coeff
    elif basis == "atomic":
        if UHF == False:
            coeffs = pyscf.lo.nao.nao(Molecule,HF)
        else:
            alpha_coeffs = pyscf.lo.nao.nao(Molecule, HF, mo_coeff=HF.mo_coeff[0])
            beta_coeffs = pyscf.lo.nao.nao(Molecule, HF, mo_coeff=HF.mo_coeff[1])
            coeffs = [alpha_coeffs, beta_coeffs]
    if UHF == False:
        for i in  range(len(coeffs[:,0])):
            # print(mf.mo_coeff)
            try:
                os.mkdir(f"{path}/Orbitals")
            except FileExistsError:
                pass
            cubegen.orbital(Molecule, f"{path}Orbitals/Orbital_{i}-0.cube", coeffs[:,i])
    else:
        for i in range(len(HF.mo_coeff[0][:])):
            try:
                os.mkdir(f"{path}/Orbitals")
            except FileExistsError:
                pass
            cubegen.orbital(Molecule, f"{path}Orbitals/Orbital_{i}-0.cube", HF.mo_coeff[0][i])
            cubegen.orbital(Molecule, f"{path}Orbitals/Orbital_{i}-1.cube", HF.mo_coeff[1][i])
            
def FunctionalChecker(Functionals: list) -> list:
    for Functional in Functionals:
        if Functional not in XC_ALIAS:
            print(f"ERROR: {Functional} Functional not recognised in this code... Known alias functionals are: ")
            print(XC_ALIAS)
            raise Exception("DFT Functional keyword error... Currently only DFT functional alias' are supported... ")
        
    return Functionals

def DFT(Molecule: pyscf.M, XC: str, Dispersion: str, unrestricted: bool, grid:int, GPU:bool, charges=None, locs=None):
    """Performs DFT on a given molecule using pySCF. Allows user to chose either GPU or CPU implementation, however if the GPU implementation is unavailable, it will roll-back to the CPU implementation in pySCF. Default SCF convergence = e-12. Default max SCF cycles = 50

    Args:
        Molecule (pyscf.M): pySCF Molecule object.
        XC (str): Exchange-Correlation functional in pySCF format.
        Dispersion (str): Dispersion correction in pySCF format
        unrestricted (bool): Whether to use unrestricted DFT
        grid (int): DFT integration grid
        GPU (bool): Whether to use gpu4pyscf

    Returns:
        md_DFT: pySCF DFT object that has completed
    """
    if unrestricted == False:
        if Molecule.spin != 0:
            unrestricted = True
    if GPU == True:
        if unrestricted == True:
            mf_DFT = UKSG(Molecule)
        elif unrestricted == False:
            mf_DFT = RKSG(Molecule)
    else:
        if unrestricted == True:
            mf_DFT = pyscf.dft.UKS(Molecule)
        elif unrestricted == False:
            mf_DFT = pyscf.dft.RKS(Molecule)
    mf_DFT.xc = XC
    mf_DFT.grids.level = grid
    if Dispersion != "None":
        mf_DFT.disp = Dispersion
    mf_DFT.conv_tol = 1e-7
    mf_DFT.max_cycle = 300
    mf_DFT.density_fit()
    if charges != None:
        mf = qmmm.mm_charge(mf_DFT, locs, charges, unit="Ang" )
        mf.verbose=0
        # try:
            # mf.init_guess_by_chkfile("tmp.chk")
            # print("Using checkfile")
        # except:
            # print("Not using checkfile")
            # pass
        mf.kernel()
        grad = mf.nuc_grad_method().kernel()
        # mf.dump_chk("tmp.chk")
        return mf, grad
    else:
        mf_DFT.kernel()
        return mf_DFT

def NN_MF(mol: pyscf.gto.Mole,Dispersion, grid, charges=None, locs=None, dm0=None):
    try:
        if charges == None:
            qmmm_calc = False
        else:
            qmmm_calc = True
    except: 
        qmmm_calc = True
    mf = pyscf.scf.UKS(mol)
    mf.conv_tol = 1e-8
    mf.max_cycle = 300
    mf.verbose=0
    mf.density_fit()
    mf.grids.level = grid
    if Dispersion != "None":
        mf.disp = Dispersion
    mf._numint = dm21.NeuralNumInt(dm21.Functional.DM21)
    if qmmm_calc == True:
        qmmm = pyscf.qmmm.mm_charge(mf, locs, charges, unit="Ang")
        qmmm.verbose=0
        qmmm.kernel(dm0=dm0)
        
        return qmmm
    else:
        mf.run()
        # mf.kernel()
        return mf

def CCSD(MF, FrozenCore:bool, Tripples:bool):
    """
    Perform a coupled cluster singles and doubles (CCSD) calculation using PySCF.

    Parameters:
        MF (pyscf.scf.RHF): A Hartree-Fock object initialized with the molecular 
                            information. This is typically obtained from pyscf.
        FrozenCore (bool): If True, include frozen core approximation in the CCSD calculation.
        Tripples (bool): If True, include triples corrections to the CCSD calculation.

    Returns:
        tuple: A tuple containing four elements:
            - myCC (pyscf.cc.CCSD): The coupled cluster object after running the calculation.
            - et (float or None): Energy of the coupled cluster calculation  triple correction, or None if Tripples is False.
            - rdm1 (numpy.ndarray): One-particle reduced density matrix from the CCSD calculation.
            - rdm2 (numpy.ndarray): Two-particle reduced density matrix from the CCSD calculation.
    """
    myCC = cc.CCSD(MF)
    myCC.max_cycle=300
    if FrozenCore == True:
        myCC.set_frozen()
    myCC.run()
    if Tripples == True:
        et = myCC.ccsd_t()
    else: 
        et = 0
    rdm1 = myCC.make_rdm1()
    rdm2 = myCC.make_rdm2()
    return myCC, et, rdm1, rdm2
  
def GeomOpt(mf):
    """
    Perform geometry optimization using PySCF's optimize function.

    Parameters:
        mf (pyscf.scf.RHF or pyscf.scf.UHF): The molecular orbital object from which the geometry will be optimized.

    Returns:
        opt (pyscf.tools.optimize_results): The optimization results containing the optimized geometry and other details.
    """
    opt = optimize(mf, maxsteps=100)
    return opt

def FindTS(mf): # pragma: no cover
    """
    Finds a transition state (TS) using the Quantum State Dynamics (QSD) method.
    
    This function initializes and runs an optimizer for finding the transition state based on the 
    given mean-field object (`mf`). The QSD optimizer is configured to find the TS with specific parameters:
    - hess_update_freq set to 0, which means no Hessian update during optimization.
    - step size set to 0.5 for the optimization steps.
    - hmin set to 1e-3, which specifies the minimum step size for the optimizer.
    
    Parameters:
        mf (object): The mean-field object from pySCF which contains the molecular information and settings.
        
    Returns:
        optimizer (object): An instance of the QSD optimizer configured to find a transition state.
    """
    optimizer = QSD(mf, stationary_point="TS")
    optimizer.kernel(hess_update_freq=0, step=0.5, hmin=1e-3)
    return optimizer

def GetFreq(mf,imaginary_freq = False ):
    """
    Perform harmonic analysis on a molecular system using PySCF and calculate thermodynamic properties.

    Parameters:
        mf (pyscf.dft.rks.RKS or pyscf.scf.hf.RHF): A mean-field object representing the molecular system.
        imaginary_freq (bool, optional): If True, include imaginary frequencies in the analysis. Defaults to False.

    Returns:
        hessian (numpy.ndarray): The Hessian matrix of the molecular system.
        Freq (dict): A dictionary containing frequency information including 'freq_au' (frequencies in atomic units).
        Thermo (pyscf.tools.thermo.ThermoData): An object containing thermodynamic properties.
    """
    hessian = mf.Hessian().kernel()
    Freq = thermo.harmonic_analysis(mf.mol, hessian, imaginary_freq = False)
    Thermo = thermo.thermo(mf, Freq["freq_au"],298.15, 101325)
    return hessian, Freq, Thermo

def get_constant_energy_for_frozen_core_for_uhf(meanfield: pyscf.scf.uhf.UHF, freeze=0)->float:
    """Abhishek's code for calculating the Energy contribution for the cores removed in the frozen core implementation of VQE.

    Args:
        meanfield (pyscf.scf.uhf.UHF | pyscf.scf.rhf.RHF): pySCF HF calculation
        freeze (int, optional): Number of cores to freeze

    Returns:
        energy_core (float):  Nuclear energy + energy of frozen orbitals.
    """
    C = numpy.array(meanfield.mo_coeff)
    frozen_core_mo, active_mo = C[:, :, :freeze], C[:, :, freeze:]
    core_dm_alpha = numpy.dot(frozen_core_mo[0], frozen_core_mo[0].T)
    core_dm_beta = numpy.dot(frozen_core_mo[1], frozen_core_mo[1].T)
    core_dm = numpy.array([core_dm_alpha, core_dm_beta])
    hcore = meanfield.get_hcore()
    energy_core = meanfield.energy_nuc()
    veff = meanfield.get_veff(meanfield.mol, core_dm)

    energy_core += numpy.einsum("ij,ji", core_dm[0], hcore)
    energy_core += numpy.einsum("ij,ji", core_dm[1], hcore)
    energy_core += numpy.einsum("ij,ji", core_dm[0], veff[0]) * 0.5
    energy_core += numpy.einsum("ij,ji", core_dm[1], veff[1]) * 0.5
    return energy_core

def get_frozen(mol: pyscf.gto.Mole)->int:
    atoms = mol._atom
    elements = [atom[0] for atom in atoms]
    NonH = [atom for atom in elements if atom != "H"]
    return int(len(NonH)*2)

def fdiff_forces(atoms:classes.molecule,dispersion, grid, charges, charge_locs, delta, dm0):
    if type(atoms) == pyscf.gto.mole.Mole:
        molecule = classes.molecule()
        molecule.from_gtoMole(atoms)
        atoms = molecule

    forces = [tuple]*atoms.nat
    for i in range(atoms.nat):
        print(f"INFO: atom {atoms.atoms[i].element}")
        atoms.atoms[i].translate_x(delta)
        dmol = atoms.to_gtoMole(False)
        print("INFO: dx+")
        En_p = NN_MF(dmol,dispersion, grid, charges, charge_locs, dm0)
        atoms.atoms[i].translate_x(-2*delta)
        dmol = atoms.to_gtoMole(False)
        print("INFO: dx-")
        En_m = NN_MF(dmol,dispersion, grid, charges, charge_locs, dm0)
        dx = ((En_p.e_tot - En_m.e_tot)/(2*delta))*atoms.bohr2ang
        atoms.atoms[i].translate_x(delta)
        atoms.atoms[i].translate_y(delta)
        dmol = atoms.to_gtoMole(False)
        print("INFO: dy+")
        En_p = NN_MF(dmol,dispersion, grid, charges, charge_locs, dm0)
        atoms.atoms[i].translate_y(-2*delta)
        dmol = atoms.to_gtoMole(False)
        print("INFO: dy-")
        En_m = NN_MF(dmol,dispersion, grid, charges, charge_locs, dm0)
        dy = ((En_p.e_tot - En_m.e_tot)/(2*delta))*atoms.bohr2ang
        atoms.atoms[i].translate_y(delta)
        atoms.atoms[i].translate_z(delta)
        print("INFO: dz+")
        dmol = atoms.to_gtoMole(False)
        En_p = NN_MF(dmol,dispersion, grid, charges, charge_locs, dm0)
        atoms.atoms[i].translate_z(-2*delta)
        dmol = atoms.to_gtoMole(False)
        print("INFO: dz-")
        En_m = NN_MF(dmol,dispersion, grid, charges, charge_locs, dm0)
        dz = ((En_p.e_tot - En_m.e_tot)/(2*delta))*atoms.bohr2ang
        atoms.atoms[i].translate_z(delta)
        forces[i] = [dx, dy, dz]
    return forces