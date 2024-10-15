from pyBrellaSampling.Tools.classes import *
import pyBrellaSampling.Tools.utils as utils
import os

def Namd_File(NAMD: NAMDClass,  substep=1, window=0):
    """
    Automates the generation of a NAMD input file. 

    Args:
        NAMD (NAMDClass): Grabs NAMD information.
        substep (int): When breaking long simulations into smaller calcs, use substeps
        window (int): Umbrella window if performing Umbrella sampling

    Raises:
        AttributeError: If shake is not on but timestep large
    
    Returns:
        NAMD_File (str): Formatted NAMD file that can then be ran. 
        
    """
    # NAMD = NAMDClass(Calc, MM, bincoor=StartCoord)
    # NAMD.set_cellvectors(MM.CellVec, MM.CellShape)
    # if QM !=None:
    #     NAMD.set_qm(Calc, QM)
    # if MM.PME != "off":
    #     NAMD.set_pme(MM.PME)
    if NAMD.bincoor == None:
        bincoor = ""
    else:
#         bincoor = f"""bincoordinates      {NAMD.bincoor}
# extendedSystem      {NAMD.bincoor.replace(".coor", ".xsc")}"""
        bincoor = f"""bincoordinates      {NAMD.bincoor}"""

    if NAMD.timestep < 2 and NAMD.rigidBonds != "none":
        return AttributeError, f"Timestep is {NAMD.timestep} but shake is on... This will cause errors."
    if int(NAMD.steps) <= 10000:
        energy = 1
    elif NAMD.runtype == "min":
        energy = 1
    else:
        energy = NAMD.restfreq
    if NAMD.timestep >=2 and NAMD.rigidBonds == "none":
        print(f"WARNING: Timestep is {NAMD.timestep} but shake isn't on... this will slow down simulations. running anyway.")
    NAMD_File = f"""### pybrella {NAMD.outfile} input file 
# File options:
parmfile            {NAMD.parm}
ambercoor           {NAMD.ambercoor}
{bincoor}
DCDfile             {NAMD.outfile}_{substep}.{window}.dcd
DCDfreq             {NAMD.dcdfreq}
restartname         {NAMD.outfile}_{substep}.{window}.restart
restartfreq         {NAMD.restfreq}
outputname          {NAMD.outfile}_{substep}.{window}
outputTiming        {NAMD.timefreq}
outputEnergies      {energy}

# Calculation options:
amber               {NAMD.amber}
switching           {NAMD.switching}
exclude             {NAMD.exclude}
1-4scaling          {NAMD.scaling}
scnb                {NAMD.scnb}
readexclusions      {NAMD.readexclusions}
cutoff              {NAMD.cutoff}
watermodel          {NAMD.watermodel}
pairListdist        {NAMD.pairListDist}
LJcorrection        {NAMD.LJcorrection}
ZeroMomentum        {NAMD.ZeroMomentum}
rigidBonds          {NAMD.rigidBonds}
rigidTolerance      {NAMD.rigidTolerance}
rigidIterations     {NAMD.rigidIterations}
timeStep            {NAMD.timestep}
fullElectFrequency  {NAMD.fullElectFrequency}
nonBondedFreq       {NAMD.nonBondedFreq}
stepspercycle       {NAMD.stepspercycle}


# PME options:
PME                 {NAMD.PME}
PMEGridSizeX        {NAMD.PMEGridSizeX}
PMEGridSizeY        {NAMD.PMEGridSizeY}
PMEGridSizeZ        {NAMD.PMEGridSizeZ}
PMETolerance        {NAMD.PMETolerance}
PMEInterpOrder      {NAMD.PMEInterpOrder}

# Cell options:
cellBasisVector1    {NAMD.cellBasisVector1}
cellBasisVector2    {NAMD.cellBasisVector2}
cellBasisVector3    {NAMD.cellBasisVector3}
cellOrigin          {NAMD.cellOrigin}

# Temperature options:
{NAMD.heating}

# Pressure options:
{NAMD.BrensdenPressure}

# QMMM options:
qmForces            {NAMD.qmForces}
{NAMD.qmLines}
#CUDAFAST

{NAMD.colvarlines}

{NAMD.runtype}                 {NAMD.steps}
"""
    return NAMD_File

def ORCA_Wrapper(QM: QMClass, Calc: CalcClass):
    """
    DEPRECIATED! Attempt at creating a custom wrapper for linking NAMD and ORCA

    This is no longer in use but could be re-introduced. It can be used to help run more efficiently on certain HPC systems.
    Args:
        QM (QMClass): Class containing QM Calc information
        Calc (CalcClass): Class containing Other calc information such as cores

    Returns: 
        file (str): Wrapper script that can be used with the "custom" qm config
    """

    file = f"""#!/usr/bin/python3
# by Marcelo Melo (melomcr@gmail.com)
# Adapted by Ross (ross.amory98@gmail.com)

from sys import argv as sargv
from sys import exit
from os.path import dirname
import subprocess as sp

orcaConfigLines1 = \"\"\"\
!  {QM.Method} {QM.Basis} EnGrad {QM.QMExtra}
%output 
  PrintLevel Mini 
  Print [ P_Mulliken ] 1
  Print [ P_AtCharges_M ] 1
end
%PAL NPROCS {Calc.Threads} END
\"\"\"
orcaConfigLines2 = "%pointcharges \\\"\"

orcaConfigLines3 = \"\"\"\\
%coords
  CTyp xyz
  Charge {QM.Charge}
  Mult {QM.Spin}
  Units Angs
  coords

\"\"\"
inputFilename = sargv[1]
directory = dirname(inputFilename)
orcaInFileName = directory + "/"
orcaInFileName += "qmmm.input"
pcFileName = orcaInFileName
pcFileName += ".pntchrg"

orcaConfigLines2 += pcFileName + "\\"\\n"
orcaOutFileName = orcaInFileName
orcaOutFileName += ".TmpOut"
orcaGradFileName = orcaInFileName
orcaGradFileName += ".engrad"
finalResFileName = inputFilename
finalResFileName += ".result"

infile = open(inputFilename,"r")

line = infile.readline()
numQMatms = int(line.split()[0])
numPntChr = int(line.split()[1].replace("\\n",""))
outLinesQM = []
outLinesPC = []
outLinesPC.append(str(numPntChr) + "\\n")
ident = "  "

lineIndx = 1
for line in infile:

    posx = line.split()[0]
    posy = line.split()[1]
    posz = line.split()[2]

    if lineIndx <= numQMatms:
        element = line.split()[3].replace("\\n","")

        outLinesQM.append(ident + " ".join([element,posx,posy,posz]) + "\\n")
    else:
        charge = line.split()[3]

        outLinesPC.append(" ".join([charge,posx,posy,posz]) + "\\n")

    lineIndx += 1
outLinesQM.append(ident + "end" + "\\n")
outLinesQM.append("end" + "\\n")


infile.close()
with open(orcaInFileName,"w") as outQMFile:

    outQMFile.write(orcaConfigLines1)
    outQMFile.write(orcaConfigLines2)
    outQMFile.write(orcaConfigLines3)

    for line in outLinesQM:
        outQMFile.write(line)

with open(pcFileName,"w") as outPCFile:

    for line in outLinesPC:
        outPCFile.write(line)
cmdline = "cd " + directory + "; "
cmdline += "{QM.QMpath} "
cmdline += orcaInFileName + " '--use-hwthread-cpus --cpu-list CPUS' > " + orcaOutFileName

proc = sp.Popen(args=cmdline, shell=True)
proc.wait()

gradFile = open(orcaGradFileName,"r")
for i in range(3):
    gradFile.readline()

orcaNumQMAtms = int(gradFile.readline().replace("\\n",""))
if orcaNumQMAtms != numQMatms:
    print("ERROR: Expected",numQMatms,"but found",orcaNumQMAtms,"atoms in engrad file!")
    exit(1)
    
for i in range(3):
    gradFile.readline()

finalEnergy = gradFile.readline().replace("\\n","").strip()

print("ORCA energy: ", finalEnergy,"Eh")
finalEnergy = str( float(finalEnergy) * 627.509469 )

print("ORCA energy: ", finalEnergy,"kcal/mol")
for i in range(3):
    gradFile.readline()
grads = []
for i in range(orcaNumQMAtms):

    grads.append( list() )

    for j in range(3):
        gradComp = gradFile.readline().replace("\\n","").strip()
        gradComp = float(gradComp) * -1185.82151
        grads[i].append( str(gradComp) )

gradFile.close()

tmpOutFile = open(orcaOutFileName,"r")

qmCharges = []


chargeSection = False

iterate = True
while iterate:

    line = tmpOutFile.readline()

    if line.find("MULLIKEN ATOMIC CHARGES") != -1:
        chargeSection = True
        line = tmpOutFile.readline()
        continue

    if chargeSection:
        length = len(line.split())
        qmCharges.append(line.split()[length-1].replace("\\n","").strip())
        pass

    if len(qmCharges) == numQMatms:
        break

tmpOutFile.close()

finFile = open(finalResFileName,"w")
finFile.write(finalEnergy + "\\n")
for i in range(numQMatms):

    finFile.write(" ".join(grads[i]) + " " + qmCharges[i] + "\\n")

finFile.close()
exit(0)

"""
    return file

def ORCA_FileGen(Molecule: MolClass, ORCA: ORCAClass):
    """
    Generates an input file for a standalone ORCA calculation

    Args:
        Molecule (MolClass): Molecular information
        ORCA (ORCAClass): ORCA Calculation information
    
    Returns:  
        file (str): ORCA input file
    
    """
    file = f"""! {ORCA.method} {ORCA.basis} {ORCA.dispersion} {ORCA.dificulty} {ORCA.convergence} {ORCA.grid} {ORCA.restart} {ORCA.calculation} {ORCA.extras}
%PAL NPROCS {ORCA.cores} END

*xyz {Molecule.charge} {Molecule.spin}
"""
    for i in range(Molecule.nat):
        file += f"{Molecule.element[i]} {round(float(Molecule.x[i]),4)} {round(float(Molecule.y[i]),4)} {round(float(Molecule.z[i]),4)} \n"
    file += "*"
    return file

