from pyBrellaSampling.classes import NAMDClass
import pyBrellaSampling.utils as utils

def Namd_File(NAMD,  substep=1, window=0):
    # NAMD = NAMDClass(Calc, MM, bincoor=StartCoord)
    # NAMD.set_cellvectors(MM.CellVec, MM.CellShape)
    # if QM !=None:
    #     NAMD.set_qm(Calc, QM)
    # if MM.PME != "off":
    #     NAMD.set_pme(MM.PME)
    if NAMD.bincoor == None:
        bincoor = ""
    else:
        bincoor = f"bincoordinates      {NAMD.bincoor}"
    if NAMD.timestep <= 2 and NAMD.rigidBonds != "none":
        return AttributeError, f"Timestep is {NAMD.timestep} but shake is on... This will cause errors."
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
BerendsenPressure   {NAMD.BrensdenPressure}

# QMMM options:
qmForces            {NAMD.qmForces}
{NAMD.qmLines}

{NAMD.colvarlines}

{NAMD.runtype}                 {NAMD.steps}
"""
    return NAMD_File
