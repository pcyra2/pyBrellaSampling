import os
import time
from pprint import pprint
from pyscf import grad as pygrad
import numpy

import pyscf
from pyscf import lib
from pyscf import qmmm, dft, gto, scf, df
import density_functional_approximation_dm21 as dm21
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
    elif type(atoms) == molecule:
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

def textDump(text, path:str):
    """Prints a text file, 

    Args:
        text: accepts either array of text or string.
        path (str): file path to write to. 
    """
    with open(path, "w")as f:
        if type(text) == str:
            print(text, file=f)
        else:
            for i in text:
                print(i, file=f)
        f.close()

class atom:
    def __init__(self, element, x, y, z):
        self.element = element
        self.x = x
        self.y = y
        self.z = z
    def echo(self):
        return f"{self.element} {self.x} {self.y} {self.z}"
    def translate_x(self, distance):
        self.x += distance
    def translate_y(self, distance):
        self.y += distance
    def translate_z(self, distance):
        self.z += distance

class molecule:
    bohr2ang = 0.529177
    atoms = []
    def __init__(self, ):
        pass
    def from_atoms_list(self, atoms:atom, charge:int, spin:int):
        self.nat = len(atoms)
        self.atoms = atoms
        self.charge = charge
        self.spin = spin
    def from_gtoMole(self, mole:pyscf.gto.Mole):
        atoms = mole._atom
        self.nat = mole.natm
        self.atoms = [atom]*mole.natm
        for i, at in enumerate(atoms):
            self.atoms[i] = atom(at[0], round(at[1][0]*self.bohr2ang,6), round(at[1][1]*self.bohr2ang,6), round(at[1][2]*self.bohr2ang,6))
        self.charge = mole.charge
        self.spin = mole.spin
        self.basis = mole.basis
    def print_coords(self)->str:
        text = ""
        for at in self.atoms:
            text += at.echo()+"\n"
        return text
    def to_gtoMole(self, symmetry:bool):
        """Converts the molecule class into 

        Args:
            basis (str): basis set to describe the molecule
            symmetry (bool): Whether to use symmetry

        Returns:
            mol (pyscf.gto.M): pySCF initialised molecule
        """
        mol = genMol(self, self.charge, self.spin, self.basis, symmetry)
        return mol

def main():
    start = time.perf_counter()

    eh2kcal = 1#627.509469
    grad_fix = 1#-1185.82151
    ang2bohr = 0.529177
    kcal2pN = 1#69.5
    def grad_nuc_mm(qmmm, qm_mol, dm): # Credit to pySCF. This is pulled from v.2.8.0 © Copyright 2025, The PySCF Developers. Pulled purely for compatibility between DM21 and qmmm
            '''Nuclear gradients of the QM-MM nuclear energy
            (in the form of point charge Coulomb interactions)
            with respect to MM atoms.
            '''
            g_mm =  grad_hcore_mm(qmmm, mol, dm)
            if qm_mol is None:
                qm_mol = qm_mol
            mm_mol = qmmm.mm_mol
            coords = mm_mol.atom_coords()
            charges = mm_mol.atom_charges()
            # print(charges)
            # g_mm = numpy.zeros_like(coords)
            for i in range(qm_mol.natm):
                # qm_q = qmcharges[i]
                qm_q = qm_mol.atom_charge(i)
                # print(q1)
                qm_r = qm_mol.atom_coord(i) 
                r = lib.norm(coords - qm_r, axis=1) 
                # print(f"{qm_r=}")
                # print(f"{qm_q=}")
                # print(f"{charges=}")
                # print(f"{coords -qm_r=}")
                # print(f"{r=}")

                UnitV= (coords - qm_r)#/(numpy.abs(coords - qm_r))
                g_mm -= qm_q * numpy.einsum('i,ix,i->ix', charges, UnitV, 1/r**3)
                # q1 = qmcharges[i]*8.8541878176e-12 * 4 * 3.14159265359 * 18897161646.321 
                # q1 = mol.atom_charge(i)
                # # print(q1)
                # r1 = mol.atom_coord(i)
                # r = lib.norm(r1-coords, axis=1)
                # # print(charges)
                # g_mm += q1 * numpy.einsum('i,ix,i->ix', charges, r1-coords, 1/r**3)
            # print(f"{g_mm=}")
            return g_mm

    def get_force(qmmm, qm_mol, qmcharges):
        """
        Calculate the forces on MM charges due to electrostatic interaction with QM system.
        """
        forces = []
        mm_mol = qmmm.mm_mol

        # Loop over all MM charges
        charge_locs = mm_mol.atom_coords()
        charges = mm_mol.atom_charges()
        for i, loc in enumerate(charge_locs):
            mm_loc = loc  # Position of MM charge
            mm_chrg = charges[i]
            force_i = numpy.zeros(3)  # Force on the i-th MM charge
            
            # Calculate the interaction between the QM system and this MM charge
            for j, qm_loc in enumerate(qm_mol.atom_coords()):
                qm_position = qm_loc  # Position of QM particle
                qm_charge = qmcharges[j]    # Charge of QM particle
                
                # Calculate distance vector between QM particle and MM charge
                r_vector = mm_loc - qm_position
                r_magnitude = numpy.linalg.norm(r_vector)
                
                # Compute Coulomb force on MM charge
                if r_magnitude > 1e-6:  # Avoid division by zero
                    force_magnitude = ((qm_charge * mm_chrg)* r_vector) / (r_magnitude ** 3)
                    force_i += force_magnitude 
            
            forces.append(force_i)
        
        return numpy.array(forces)

    def NN(mol, charges, charge_locs, dm0=None):
        mf = dft.UKS(mol)
        mf._numint = dm21.NeuralNumInt(dm21.Functional.DM21mc)
        mf.conv_tol = 1E-6
        mf.conv_tol_grad = 1E-3
        if charges != None:
            qm_mf = qmmm.mm_charge(mf, charge_locs, charges, unit="Ang").run()
        else:
            qm_mf = mf.run()
        return qm_mf

    def fdiff_forces(atoms:molecule, charges, charge_locs, delta, dm0):
        if type(atoms) == pyscf.gto.Mole:
            mol = molecule()
            mol.from_gtoMole(atoms)
            atoms = mol
            
        mol = atoms.to_gtoMole(False)

        start = time.perf_counter()
        forces = [tuple]*atoms.nat
        for i in range(atoms.nat):
            atoms.atoms[i].translate_x(delta)
            dmol = atoms.to_gtoMole(False)
            En_p = NN(dmol, charges, charge_locs, dm0)
            atoms.atoms[i].translate_x(-2*delta)
            dmol = atoms.to_gtoMole(False)
            En_m = NN(dmol, charges, charge_locs, dm0)
            dx = ((En_p.e_tot - En_m.e_tot)/(2*delta))*ang2bohr
            # print(dx)
            atoms.atoms[i].translate_x(delta)
            atoms.atoms[i].translate_y(delta)
            dmol = atoms.to_gtoMole(False)
            En_p = NN(dmol, charges, charge_locs, dm0)
            atoms.atoms[i].translate_y(-2*delta)
            dmol = atoms.to_gtoMole(False)
            En_m = NN(dmol, charges, charge_locs, dm0)
            dy = ((En_p.e_tot - En_m.e_tot)/(2*delta))*ang2bohr
            # print(dy)
            atoms.atoms[i].translate_y(delta)
            atoms.atoms[i].translate_z(delta)
            dmol = atoms.to_gtoMole(False)
            En_p = NN(dmol, charges, charge_locs, dm0)
            atoms.atoms[i].translate_z(-2*delta)
            dmol = atoms.to_gtoMole(False)
            En_m = NN(dmol, charges, charge_locs, dm0)
            dz = ((En_p.e_tot - En_m.e_tot)/(2*delta))*ang2bohr
            # print(dz)
            atoms.atoms[i].translate_z(delta)
            forces[i] = [dx, dy, dz]
            # print(forces[i])
        end = time.perf_counter()

        print(f"INFO: FDiffForces time {end-start}")
        print("Finite difference forces:")
        for i in forces:
            print(i)
        return forces
    
    def grad_hcore_mm(qmmm, qm_mol, dm, mol=None):
            r'''Nuclear gradients of the electronic energy
            with respect to MM atoms:

            ... math::
                g = \sum_{ij} \frac{\partial hcore_{ij}}{\partial R_{I}} P_{ji},

            where I represents MM atoms.

            Args:
                dm : array
                    The QM density matrix.
            '''
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


    # atoms = [str]*nat
    charge = 0
    spin = 0
    method = "PBE"
    basis = "6-31G"
    charges= None#[1]
    charge_locs = None#[[0,1,0] ]
    ncharges = 0#len(charges)
    mol = gto.M(atom="O 0 0 0 ; H 0.758602  0.000000  0.504284; H 0.758602  0.000000  -0.504284 ", basis=basis,unit="Ang", charge=charge, spin=spin)
    nat = mol.natm
    mf = NN(mol, charges, charge_locs)
    qm_time = time.perf_counter()
    dm = mf.make_rdm1()
    grad = fdiff_forces(mol,  charges, charge_locs, 0.01, dm)
    # pc_grad = grad_nuc_mm(mf, mol, dm)
    chg =mf.mulliken_pop()
    result = [str]*(nat+ncharges+1)
    result[0] = f"{mf.e_tot*eh2kcal} {ncharges}"
    for i in range(nat):
        result[i+1] = f"{grad[i][0]*grad_fix} {grad[i][1]*grad_fix} {grad[i][2]*grad_fix} {chg[1][i]}"
    # for i in range(ncharges):
        # result[i+nat+1] = f"{pc_grad[i][0]*kcal2pN} {pc_grad[i][1]*kcal2pN} {pc_grad[i][2]*kcal2pN}"


    print("OUTPUT FILE:")
    for i in result:
        print(i)
    textDump(result, f"qmmm_0.input.result" )

    stop = time.perf_counter()
    print(f"INFO: qm time taken {round(qm_time - start, 2) } s")
    print(f"INFO: Time taken {round(stop - start, 2)} s")


if __name__ == "__main__":
    main()