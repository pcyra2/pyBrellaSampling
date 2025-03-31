from sys import argv as sargv
from sys import exit
import pyBrellaSampling.Tools.utils as utils
import pyBrellaSampling.Tools.pyscf_tools as pyscf_tools

def parse_NAMD(file:str):
    inputFilename = file
    infile = utils.file_read(inputFilename)
    # Gets number of atoms in the Quantum Chemistry region (= QM atoms + Link atoms)
    numQMatms = int(line[0].split()[0])
    # Gets number of point charges
    numPntChr = int(line[0].split()[1].replace("\n",""))

    atoms = [str]*numQMatms
    charges = [float]*numPntChr
    locs = [tuple]*numPntChr

    for i in range(1,len(infile)):
        line = infile[i]

        posx = line.split()[0]
        posy = line.split()[1]
        posz = line.split()[2]

        if i <= numQMatms:
            # ORCA's format requires the fileds to be ordered begining with the
            # atom's element symbol, and followed by the XYZ coordinates.

            element = line.split()[3].replace("\n","")

            atoms[i-1] = f"{element} {posx} {posy} {posz}"

        else:

            # ORCA's format requires the fileds to be ordered begining with the
            # charge, and followed by the XYZ coordinates.

            charges[i-1] = line.split()[3]

            locs[i-1] = (posx,posy,posz)

    return atoms, charges, locs


def main():
    atoms, charges, locs = parse_NAMD(sargv[1])
    mol = pyscf_tools.genMol(atoms,)
    mf = pyscf_tools.UHF(mol,charges,locs)
    
    pass

if __name__ == "__main__":
    main()