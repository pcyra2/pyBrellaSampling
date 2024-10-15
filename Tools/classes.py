import pandas as pd
import numpy as np
import math as math
from pyBrellaSampling.UserVars.MM_Variables import MM_DefaultVars as MMVars
from pyBrellaSampling.UserVars.SoftwarePaths import NAMD_CPU, NAMD_GPU, ORCA_PATH
# from pyBrellaSampling.Tools.globals import globals.verbosity, WorkDir, globals.DryRun, globals.parmfile
import pyBrellaSampling.Tools.globals as globals

class UmbrellaClass:
    """Class for containing umbrella sampling information.
    
    Attributes:
        args (dict): ArgParse arguments to read in user defined variables
        Min (float): Minimum Umbrella value
        bins (int): Number of Umbrella bins
        Start (float): Current value of the colvar. (taken from the equilibrated start structure)
        Width (float): Width of the Umbrella bins, (Can be negative, if so Min => Max)
        atom1 (int): Atom1 for colective variable (always required)
        atom2 (int): Atom2 for colective variable (always required)
        atom3 (int): Atom3 for colective variable (required if angle or dihedral) else 0
        atom4 (int): Atom4 for colective variable (required if dihedral) else 0
        PullForce (float): Force to use when pulling to each lambda window
        ConstForce (float): Force to use when running constant Umbrella sampling.
        BinVals (list): List of Umbrella bin values
        StartBin (int): Bin closest to the starting value of the collective variable

    """
    def __init__(self, args: dict, Min: float, bins: int, Start:float, Width:float, ):
        """
        Umbrella class Init. 
        Args:
            args (Namespace): ArgParse arguments to read in user defined variables
            Min (float): Minimum Umbrella value
            bins (int): Number of Umbrella bins
            Start (float): Current value of the colvar. (taken from the equilibrated start structure)
            Width (float): Width of the Umbrella bins, (Can be negative, if so Min => Max)
        Attributes:
            atom1 (int): Atom1 for colective variable (always required)
            atom2 (int): Atom2 for colective variable (always required)
            atom3 (int): Atom3 for colective variable (required if angle or dihedral) else 0
            atom4 (int): Atom4 for colective variable (required if dihedral) else 0
            PullForce (float): Force to use when pulling to each lambda window
            ConstForce (float): Force to use when running constant Umbrella sampling.
        """
        self.atom1 = int(args["AtomMask"].split(",")[0])
        self.atom2 = int(args["AtomMask"].split(",")[1])
        self.atom3 = int(args["AtomMask"].split(",")[2])
        self.atom4 = int(args["AtomMask"].split(",")[3])
        self.Min = Min
        self.Bins = bins
        self.Start = Start
        self.Width = Width
        self.PullForce = args["PullForce"]
        self.ConstForce = args["ConstForce"]
    def add_bins(self, BinVals:list):
        """
        Adds values to the Umbrella bins

        Args:
            BinVals (list): List of Umbrella bin values
        """
        self.BinVals = BinVals
    def add_start(self, StartBin:int):
        """
        Defines which bin is the start bin (For use when pulling evenly.)

        Args:
            StartBin (int): Bin closest to the starting value of the collective variable
        """
        self.StartBin = StartBin
    def set_force(self, Force:float): ### To combine Umbrella for use in Standalone calculations
        """
        Can change the force for Steered MD, Not used in Umbrella sampling

        Args:
            Force (float): Force to be applied for SMD
        """
        self.PullForce = Force
        self.ConstForce = Force
class LabelClass:
    """Class for containing bond and dihedral label information (For structure analysis)
    Attributes:
        bond (list[str]): List of bond atom selections in format: ("at1,at2")
        bondName (list[str]): List of bond names
        bondThresh (list[float]): List of standard bond lengths
        dihedral (list[str]): List of dihedral atom selections in format: ("at1,at2,at3,at4")
        dihedralName (list[str]): List of dihedral names
        dihedralTarget1 (list[float]): List of dihedral primary values
        dihedralTarget2 (list[float]): List of dihedral secondary values
        dihedralTarget1Name (list[str]): List of names/ID for primary dihedral values
        dihedralTarget2Name (list[str]): List of names/ID for secondary dihedral values
        parm (str): parameter file to be loaded by vmd
        file (str): coordinate file to be loaded by vmd 
    """
    bond = []
    bondName = []
    bondThresh = []
    dihedral = []
    dihedralName = []
    dihedralTarget1 = []
    dihedralTarget2 = []
    dihedralTarget1Name = []
    dihedralTarget2Name = []
    def __init__(self,  parm: str):
        """
        Initialises a label class, used for structural analysis
        
        Args:
            parm (str): filename of param file for use in vmd.
        """
        self.parm = parm
    def add_bond(self, selection: str, name: str, thresh: float):
        """
        Adds bond labels to the label class.

        Args:
            selection (str): Atomic selection, as a coma delimited string in style "at1,at2"
            name (str): Name/ID of bond
            thresh (float): Standard bond length (bond classed as broken at 1.2x thresh)
        
        Attributes:
            bond (list[str]): List of bond atom selections in format: ("at1,at2")
            bondName (list[str]): List of bond names
            bondThresh (list[float]): List of standard bond lengths
        """
        self.bond.append(selection)
        self.bondName.append(name)
        self.bondThresh.append(thresh)
    def add_dihedral(self, selection: str, name: str, target1: float, t1name:str , target2: float, t2name:str):
        """
        Adds dihedral labels to the label class.
        
        Args:
            selection (str): Atom selection, as a coma delimited string in style "at1, at2, at3, at4"
            name (str): Name/ID of dihedral
            target1 (float): Primary value of dihedral angle
            t1name (str): Name/ID of primary dihedral angle
            target2 (float): Secondary value of dihedral angle
            t2name (str): Name/ID of secondary dihedral angle
        
        Attributes:
            dihedral (list[str]): List of dihedral atom selections in format: ("at1,at2,at3,at4")
            dihedralName (list[str]): List of dihedral names
            dihedralTarget1 (list[float]): List of dihedral primary values
            dihedralTarget2 (list[float]): List of dihedral secondary values
            dihedralTarget1Name (list[str]): List of names/ID for primary dihedral values
            dihedralTarget2Name (list[str]): List of names/ID for secondary dihedral values
        """
        self.dihedral.append(selection)
        self.dihedralName.append(name)
        self.dihedralTarget1.append(target1)
        self.dihedralTarget1Name.append(t1name)
        self.dihedralTarget2.append(target2)
        self.dihedralTarget2Name.append(t2name)
    def file_name(self, name: str):
        """
        Name of coordinate file to be loaded by vmd

        Args:
            name (str): filename of coordinates

        """
        self.file = name
    def clear_Vars(self):
        """
        Clears all label data, whilst maintaining file and param.
        """
        self.bond = []
        self.bondName = []
        self.bondThresh = []
        self.dihedral = []
        self.dihedralName = []
        self.dihedralTarget1 = []
        self.dihedralTarget2 = []
        self.dihedralTarget1Name = []
        self.dihedralTarget2Name = []
class DataClass:
    """Class for containing generic analysis data.
    Attributes:
        dat: List of data DataFrames
    """
    dat = []
    def __init__(self, Name: str,):
        """
        Data class init. for use in structure/trajectory analysis

        Args:
            Name (str): Name of data
        """
        self.name = Name
    def add_data(self, name:str, window:int, data:float):
        """
        Adds data to the data class
        Args:
            name (str): name/ID of analysis variable 
            window (int): window
            data (float): data point

        """
        self.dat.append(pd.DataFrame(
            data={"Name": name, "Window": window, "Data": data}))
    def __repr__(self):
        return f"{self.name}: \n{self.dat}"
class BondClass:
    """Class for containing bond data
    Attributes:
        atom1 (int): Atom1 in bond
        atom2 (int): Atom2 in bond
        name (str): Name/ID of bond
        threshold (float): Standard value of bond, (Broken if 1.2x bigger)
    """
    def __init__(self, atom1: int, atom2: int, name: str, threshold: float):
        """
        Bond class init.

        Args:
            atom1 (int): Atom1 in bond
            atom2 (int): Atom2 in bond
            name (str): Name/ID of bond
            threshold (float): Standard value of bond, (Broken if 1.2x bigger)

        """
        self.at1 = int(atom1)
        self.at2 = int(atom2)
        self.name = name
        self.thresh = threshold
    def __repr__(self):
        return f"Atoms: {self.at1} {self.at2}, name: {self.name}, threshold: {self.thresh}"
class DihedralClass:
    """Class for containing dihedral data
    Attributes:
        atom1 (int): Atom1 in dihedral
        atom2 (int): Atom2 in dihedral
        atom3 (int): Atom3 in dihedral
        atom4 (int): Atom4 in dihedral
        name (str): Name of dihedral angle
        target1 (float): Value of primary dihedral state
        t1name (str): Name or ID of primary dihedral state
        target2 (float): Value of secondary dihedral state
        t2name (str): Name or ID of secondary dihedral state
    """
    def __init__(self, atom1: int, atom2: int, atom3: int, atom4: int, name: str, target1: float, t1name: str, target2: float, t2name: str):
        """
        Dihedral class Init.

        Args:
            atom1 (int): Atom1 in dihedral
            atom2 (int): Atom2 in dihedral
            atom3 (int): Atom3 in dihedral
            atom4 (int): Atom4 in dihedral
            name (str): Name of dihedral angle
            target1 (float): Value of primary dihedral state
            t1name (str): Name or ID of primary dihedral state
            target2 (float): Value of secondary dihedral state
            t2name (str): Name or ID of secondary dihedral state

        """
        self.at1 = int(atom1)
        self.at2 = int(atom2)
        self.at3 = int(atom3)
        self.at4 = int(atom4)
        self.name = name
        self.target1 = target1
        self.target1Name = t1name
        self.target2 = target2
        self.target2Name = t2name
    def __repr__(self):
        return f"Atoms: {self.at1} {self.at2} {self.at3} {self.at4}, Name: {self.name}, {self.target1Name} = {self.target1}, {self.target2Name} = {self.target2}"

class WhamClass:
    """Class for controlling WHAM calculation information
    Attributes:
        Name: Name of wham calculation
        Force: Force for biasing windows
        Type: discrete or continuous
    """
    def __init__(self, Name: str, Force: float, type="discrete"):
        """
        Wham Class initialisation. Used for containing information about the wham calc. 
        Args:
            Name (str): Name of wham calculation
            Force (float): Force for biasing windows
            type (str): discrete or continuous
        """
        self.Name = Name
        self.Force = Force
        assert type == "discrete"  or type == "periodic", f"ERROR: {type} Wham type not supported. Should be either discrete or continuous"
        self.Type = type

class QMClass:
    """ Class containing QM information 
    Attributes:
        SelFile (str): PDB file containing QM-Zone infrmation, (Should be generated by vmd)
        QMpath (str): Path to the QM executable
        QMSel (str): string containing vmd selection algebra of qm zone
        Charge (int): Net charge of QM zone
        Spin (int): Total spin of QM zone
        Method (str): QM method supported by ORCA
        Basis (str): QM basis set supported by ORCA
        QMExtra (str): Any extra commands to provide to ORCA
        SelFile (str): Selection file
    """
    SelFile="../syst-qm.pdb"
    def __init__(self, args: dict):
        """
        QM Class init.

        Args:
            args (Namespace): ArgParse formatted user inputs.

        Attributes:
            QMpath (str): Path to the QM executable
            QMSel (str): string containing vmd selection algebra of qm zone
            Charge (int): Net charge of QM zone
            Spin (int): Total spin of QM zone
            Method (str): QM method supported by ORCA
            Basis (str): QM basis set supported by ORCA
            QMExtra (str): Any extra commands to provide to ORCA
        """
        self.QMpath = args["QmPath"]
        self.QMSel = args["QmSelection"]
        self.Charge = args["QmCharge"]
        self.Spin = args["QmSpin"]
        self.Method = args["QmMethod"]
        self.Basis = args["QmBasis"]
        self.QMExtra = args["QmArgs"]
    def set_selfile(self, File: str):
        """
        Sets the selection file

        Args:
            File (str): Name of file

        Attribute:
            SelFile (str): Selection file
        """
        self.SelFile = File

class MMClass:
    """Class for containing molecular mechanics calculation information
    
    Attributes:
        CPUNamd (str): Path to CPU version of NAMD (For QMMM)
        GPUNamd (str): Path to GPU version of NAMD (For MD)
        Shake (str): NAMD shake settings, usually all or none
        TimeOut (int): Frequency to print out TIMING and PERFORMANCE data
        RestOut (int): Frequency to update the restart files
        TrajOut (int): Frequency to add frames to the trajectory file
        Temp (float): Temperature of simulation.
        Ensemble (str): min, heat, NVT or NVP
        Steps (int): Number of steps in the simulation
        TimeStep (float): Time step in fs
        CellVec (float): Cell vector in Angstrom
        Force (float): Force constant used by the colvar module
        parmfile (str): file of the amber parameters (.parm7)
        ambercoor (str): file of the amber coordinates (.rst7)
        PME (str): Defines the PME interactions (off or on)
        CutOff (float): Sets the long range cutoff distance
    """
    CellVec = 150.2205127 # Default although should change. 
    CellShape = "oct"
    Temp = MMVars["temperature"]
    CutOff = MMVars["cutoff"]
    PME = "off"
    parmfile = f"../{globals.parmfile}"
    ambercoor = "../start.rst7"
    Shake = "none"
    def __init__(self, ):
        """
        Contains MD parameters

        Attributes:
            CPUNamd (str): Path to CPU version of NAMD (For QMMM)
            GPUNamd (str): Path to GPU version of NAMD (For MD)
        """
        self.CPUNamd = NAMD_CPU
        self.GPUNamd = NAMD_GPU
    def Set_Shake(self, Shake: str):
        """
        Defines the shake/rattle state. (For restraining bonds during large time steps.)

        Args:
            Shake (str): NAMD shake settings, usually all, water or none
        """
        self.Shake = Shake
    def Set_Outputs(self, TimeOut: int, RestOut: int, TrajOut: int):
        """
        Defines how often to print simulation data.

        Args:
            TimeOut (int): Frequency to print out TIMING and PERFORMANCE data
            RestOut (int): Frequency to update the restart files
            TrajOut (int): Frequency to add frames to the trajectory file
        """
        self.TimeOut = TimeOut
        self.RestOut = RestOut
        self.TrajOut = TrajOut
    def Set_Temp(self, Temp: float):
        """
        Sets/Changes the temperature of the simulation

        Args:
            Temp (float): Temperature of simulation.
        """
        self.Temp = Temp
    def Set_Ensemble(self, Ensemble: str):
        """
        Sets the ensemble of the simulation.
        Args:
            Ensemble (str): min, heat, NVT or NVP

        """
        self.Ensemble = Ensemble
    def Set_Length(self, Steps: int, TimeStep: float=0.05):
        """
        Defines simulation length 

        Args:
            Steps (int): Number of steps in the simulation
            TimeStep (float=0.05): Time step in fs
        """
        self.Steps = int(Steps)
        self.TimeStep = float(TimeStep)
    def Change_Cell(self, CellVec: float):
        """
        Changes the cell vector, (used when reading in from the parm file)

        Args:
            CellVec (float): Cell vector in Angstrom

        """
        self.CellVec = CellVec
    def Set_Force(self, Force: float):
        """
        Sets the force when running steered MD or umbrella sampling

        Args:
            Force (float): Force constant used by the colvar module
        """
        self.Force = Force
    def Set_Files(self, parm: str, ambercoor: str):
        """
        Defines amber files that are required when using amber parameters

        Args:
            parm (str): file of the amber parameters (.parm7)
            ambercoor (str): file of the amber coordinates (.rst7)
        """
        self.parmfile = parm
        self.ambercoor = ambercoor

class CalcClass:
    """Class for containing generic calculation information
    
    Attributes:
        QM (str): "off" or "on" for toggling QMMM calculations.
        Threads (int): Number of CPU threads per individual calculation (i.e. per umbrella window)
        Memory (int): Ammount of RAM per individual calculation (i.e. per umbrella window)
        MaxSteps (int): Maximum number of steps per individual calculation (For breaking long simulations to fit on short HPC queues)
        Name (str): Name of Job
        OutFile (str): file name
        Id (int): SLURM jobID
        QM (str): off or on
    """
    QM = "off"
    def __init__(self, args: dict):
        """
        Initiates Calculation class, to set some standard calculation variables. 
        Args:
            args (dict): ArgParse user input

        Attributes:
            Threads (int): Number of CPU threads per individual calculation (i.e. per umbrella window)
            Memory (int): Ammount of RAM per individual calculation (i.e. per umbrella window)
            MaxSteps (int): Maximum number of steps per individual calculation (For breaking long simulations to fit on short HPC queues)
        """
        self.Threads = args["CoresPerJob"]
        self.Memory = args["MemoryPerJob"]
        try:
            self.MaxSteps = args["MaxStepsPerCalc"]
        except KeyError:
            self.MaxSteps = 0
    def Job_Name(self, Name: str,):
        """
        Sets the job name 

        Args:
            Name (str): Name of Job

        """  
        self.Name = Name
    def Set_OutFile(self, OutFile: str):
        """
        sets the output file name

        Args:
            OutFile (str): file name
        """
        self.OutFile = OutFile
    def Set_Id(self, Id: int):
        """
        Gets the slurm ID of the calculation         

        Args:
            Id (int): SLURM jobID
        """
        self.Id = Id
    def Set_QM(self, QM: str):
        """
        Sets the QM state of the calculation
        Args:
            QM (str): off or on
        """
        self.QM = QM

class NAMDClass:
    """Class for containing NAMD Variables
    
    Attributes:
        amber (str): Whether to use amber coordinates and parameters
        switching (str): NAMD Variable (refer to manual)
        exclude (str): NAMD Variable (refer to manual)
        scaling (float): NAMD Variable (refer to manual)
        scnb (float): NAMD Variable (refer to manual)
        readexclusions (str): NAMD Variable (refer to manual)
        watermodel (str): Water model to be used (TIP3P)
        pairListDist (float): NAMD Variable (refer to manual)
        LJcorrection (str): NAMD Variable (refer to manual)
        ZeroMomentum (str): NAMD Variable (refer to manual)
        rigidTolerance (str): NAMD Variable (refer to manual)
        rigidIterations (int): NAMD Variable (refer to manual)
        fullElectFrequency (float): NAMD Variable (refer to manual)
        nonBondedFreq (int): NAMD Variable (refer to manual)
        stepspercycle (int): NAMD Variable (refer to manual)
        PME (str): Toggle for PME
        PMEGridSizeX (int):  NAMD Variable (refer to manual)
        PMEGridSizeY (int): NAMD Variable (refer to manual)
        PMEGridSizeZ (int): NAMD Variable (refer to manual)
        PMETolerance (str): NAMD Variable (refer to manual)
        PMEInterpOrder (int): NAMD Variable (refer to manual)
        qmForces (str): Whether to use QMMM
        qmLines (str): QM config lines
        colvarlines (str): Colvar config lines for SMD and Umbrella sampling
        parm (str): parameter file
        ambercoor (str): Amber coordinate file
        outfile (str): output filename
        dcdfreq (int): frequency to update trajectories
        restfreq (int): frequency to update restart file
        timefreq (int): frequency to print timing and performance information
        cutoff (float): non-bonded cutoff
        timestep (float): timestep for MD
        rigidBonds (str): whether to use shake
        steps (int): Number of steps 
        BrensdenPressure (str): Whether to control pressure of system
        heating (str): Heating information
        runtype (str): whether dynamics or minimisation
        cellBasisVector1 (str): NAMD Variable (refer to manual)
        cellBasisVector2 (str): NAMD Variable (refer to manual)
        cellBasisVector3 (str): NAMD Variable (refer to manual)
        cellOrigin (str): NAMD Variable (refer to manual)
        bincoor (str): coordinates to actully be used as start coodinates
        ambercoor (str): coordinates related to amber parameters (Must be in amber format)
        parm (str): AMBER parameter file 
    """
    amber = MMVars["amber"]
    switching = MMVars["switching"]
    exclude = MMVars["exclude"]
    scaling = MMVars["1-4scaling"]
    scnb = MMVars["scnb"]
    readexclusions = MMVars["readexclusions"]
    watermodel = MMVars["watermodel"]
    pairListDist = MMVars["pairListdist"]
    LJcorrection = MMVars["LJcorrection"]
    ZeroMomentum = MMVars["ZeroMomentum"]
    rigidTolerance = MMVars["rigidTolerance"]
    rigidIterations = MMVars["rigidIterations"]
    fullElectFrequency = MMVars["fullElectFrequency"]
    nonBondedFreq = MMVars["nonBondedFreq"]
    stepspercycle = MMVars["stepspercycle"]
    PME = MMVars["PME"]
    PMEGridSizeX = MMVars["PMEGridSizeX"]
    PMEGridSizeY = MMVars["PMEGridSizeY"]
    PMEGridSizeZ = MMVars["PMEGridSizeZ"]
    PMETolerance = MMVars["PMETolerance"]
    PMEInterpOrder = MMVars["PMEInterpOrder"]
    qmForces = MMVars["qmForces"]
    qmLines = f""
    colvarlines = ""
    def __init__(self, Calc: CalcClass, MM: MMClass,):
        """
        NAMD Class initialise. Requres some MD Parameters 
        Args:
            Calc (CalcClass): for job name
            MM (MMClass): for MD parameters.

        Raises:
            AttributeError: if unknown ensemble is provided (either min, heat, NVT or NVP)
        
        Attributes:
        parm (str): parameter file
        ambercoor (str): Amber coordinate file
        outfile (str): output filename
        dcdfreq (int): frequency to update trajectories
        restfreq (int): frequency to update restart file
        timefreq (int): frequency to print timing and performance information
        cutoff (float): non-bonded cutoff
        timestep (float): timestep for MD
        rigidBonds (str): whether to use shake
        steps (int): Number of steps 
        BrensdenPressure (str): Whether to control pressure of system
        heating (str): Heating information
        runtype (str): whether dynamics or minimisation
        """
        self.parm = MM.parmfile
        self.ambercoor = MM.ambercoor
        # self.bincoor = bincoor
        self.outfile = Calc.Name
        self.dcdfreq = MM.TrajOut
        self.restfreq = MM.RestOut
        self.timefreq = MM.TimeOut
        self.cutoff = MM.CutOff
        self.timestep = MM.TimeStep
        self.rigidBonds = MM.Shake
        self.steps = MM.Steps
        if MM.Ensemble == "NVT": # Formats NVT calculation
            self.BrensdenPressure = ""
            self.heating = f"""langevin            on
langevinDamping     5
langevinTemp        {MM.Temp}
langevinHydrogen    off
temperature         {MM.Temp}
"""
            self.runtype = "run"
        elif MM.Ensemble == "NPT": # Formats NPT calculation
            self.BrensdenPressure = f"""useGroupPressure     yes 
useFlexibleCell     no
useConstantArea     no

langevinPiston on
langevinPistonTarget 1.01325
langevinPistonPeriod 100
langevinPistonDecay 50
langevinPistonTemp {MM.Temp}

"""
            self.heating = f"""langevin            on
langevinDamping     5
langevinTemp        {MM.Temp}
langevinHydrogen    off
temperature         {MM.Temp}
"""
            self.runtype = "run"
        elif "heat" in MM.Ensemble.casefold(): # Heat is classed as an ensemble 
            self.BrensdenPressure = "langevin       off"
            self.heating = f"""temperature         0
reassignFreq        {math.floor(MM.Steps/(MM.Temp/0.2))}
reassignIncr        0.2
reassignHold        {MM.Temp}
"""
            self.runtype = "run"
        elif "min" in MM.Ensemble.casefold(): # Min is classed as an ensemble
            self.BrensdenPressure = "langevin         off"
            self.heating = """temperature             0"""
            self.runtype = "minimize"
        else:
            return AttributeError, f"ERROR: {MM.Ensemble} not configured."
    def set_pme(self, val="off"):
        """
        Whether to use PME
        Args:
            val (undefined): off or on
        Attributes:
            PME (str): off or on
        """
        self.PME = val
    def set_cellvectors(self, CellVector: float, CellShape="oct"): ### Currently only supports truncated Octahedron cells
        """
        Sorts out NAMD's cell shapes as it does not natively support octahedral.
        Args:
            CellVector (float): first number from the .parm7 file
            CellShape (str, optional): Shape of cell. Only supports octahedral. Defaults to "oct".
        
        TODO:
            CellShape: Support other cell vectors

        Attributes:
            cellBasisVector1 (str): NAMD Variable (refer to manual)
            cellBasisVector2 (str): NAMD Variable (refer to manual)
            cellBasisVector3 (str): NAMD Variable (refer to manual)
            cellOrigin (str): NAMD Variable (refer to manual)

        Raises:
            AttributeError: If any unit cell provided that is not curently supported (Only Oct at the moment.)
        """
        if CellShape.casefold() == "oct":
            self.cellBasisVector1 = f"{CellVector} 0.0 0.0"
            self.cellBasisVector2 = f"{(-1/3)*CellVector} {(2/3)*np.sqrt(2)*CellVector} 0.0"
            self.cellBasisVector3 = f"{(-1/3)*CellVector} {(-1/3)*np.sqrt(2)*CellVector} {(-1/3)*np.sqrt(6)*CellVector}"
            self.cellOrigin = "0 0 0"
        else:
            return AttributeError, "Cell shape not supported currently, please add this functionality."
    def set_qm(self, Calc: CalcClass, QM: QMClass, index=0):
        """
        Initialises the QM section of the NAMD input file. 

        Args:
            Calc (CalcClass): Required for number of cores
            QM (QMClass): Contains QM information
            index (int): For umbrella sampling. 

        Attributes:
            qmForces (str): Toggles QMMM
            qmLines (str): Variables required for QMMM
        """
        if Calc.QM == "False":
            self.qmForces = "off"
            self.qmLines = f""
        else:
            self.qmForces = "on"
#             self.qmLines = f"""
# qmParamPDB              "../syst-qm.pdb"
# qmColumn                "beta"
# qmBondColumn            "occ"
# QMsimsPerNode           1
# QMElecEmbed             on
# QMSwitching             on
# QMSwitchingType         shift
# QMPointChargeScheme     round
# QMBondScheme            "cs"
# qmBaseDir               "/dev/shm/NAMD_{index}"
# qmConfigLine            "! {QM.Method} {QM.Basis} EnGrad {QM.QMExtra}"
# qmConfigLine            "%%output PrintLevel Mini Print\[ P_Mulliken \] 1 Print\[P_AtCharges_M\] 1 end"
# qmConfigLine            "%PAL NPROCS {Calc.Threads} END"
# qmMult                  "1 {QM.Spin}"
# qmCharge                "1 {QM.Charge}"
# qmSoftware              "orca"
# qmExecPath              "{QM.QMpath}"
# QMOutStride             1
# qmEnergyStride          1
# QMPositionOutStride     1
#
# """
            self.qmLines = f"""
            qmParamPDB              "{QM.SelFile}"
            qmColumn                "beta"
            qmBondColumn            "occ"
            QMsimsPerNode           1
            QMElecEmbed             on
            QMSwitching             on
            QMSwitchingType         shift
            QMPointChargeScheme     round
            QMBondScheme            "cs"
            qmBaseDir               "/dev/shm/RUNDIR_{index}"
            qmConfigLine            "! {QM.Method} {QM.Basis} EnGrad {QM.QMExtra}"
            qmConfigLine            "%%output PrintLevel Mini Print\[ P_Mulliken \] 1 Print\[P_AtCharges_M\] 1 end"
            qmConfigLine            "%PAL NPROCS {Calc.Threads} END"
            qmMult                  "1 {QM.Spin}"
            qmCharge                "1 {QM.Charge}"
            qmSoftware              "orca"
            qmExecPath              "{QM.QMpath}"
            QMOutStride             1
            qmEnergyStride          1
            QMPositionOutStride     1

            """
    def set_colvars(self, file: str, toggle="on"):
        """
        Sorts out collective variables for SMD and Umbrella Sampling

        Args:
            file (str): filename and location of the colvars config file.
            toggle (str): on or off depending on whether used. 

        Attributes:
            colvarlines (str): Config lines for colvars.
        """
        self.colvarlines = f"""# Colvar options:
colvars         {toggle}
colvarsConfig   {file}
"""
    def set_startcoords(self, bincoor:str, ambercoor:str ="../start.rst7", parm:str =f"../{globals.parmfile}"):
        """
        Controls the coordinate files
        Args:
            bincoor (str): coordinates to actully be used as start coodinates
            ambercoor (str): coordinates related to amber parameters (Must be in amber format)
            parm (str): AMBER parameter file 

        """
        self.bincoor = bincoor
        self.ambercoor = ambercoor
        self.parm = parm

class SLURMClass:
    """Class for containing SLURM variables
    Attributes:
        HostName (str): Hostname returned by socket.gethostname()
        WallTime (int): Maximum calculation length hours
        Cores (int): Number of cores
        Memory (int): Ammount of memory
        Partition (str): SLURM partition
        NodeCores (int): Max cores per node
        account (str): SLURM billing account
        qos (str): SLURM QoS
        Software (list): Lines containing software loading commands such as modules
        dependency (str): SLURM Dependency line
    """
    account = ""
    qos = ""
    dependency = ""
    def __init__(self, args: dict):
        """inits slurm class

        Args:
            args (dict): User defined variables

        Attributes:
            HostName (str): Hostname returned by socket.gethostname()
            WallTime (int): Maximum calculation length hours
            Cores (int): Number of cores
            Memory (int): Ammount of memory
            Partition (str): SLURM partition
            NodeCores (int): Max cores per node
            account (str): SLURM billing account
            qos (str): SLURM QoS
            Software (list): Lines containing software loading commands such as modules

        """
        self.HostName = args["HostName"]
        self.WallTime = args["MaxWallTime"]
        self.Cores = args["CoresPerJob"]
        self.Memory = args["MemoryPerJob"]
        # self.k = args["HostName"]
        self.Partition = args["Partition"]
        self.NodeCores = args["MaxCores"]
        if args["Account"].casefold() != "none":
            acc = args["Account"]
            self.account = f"#SBATCH --account={acc}"
        if args["QualityofService"].casefold() != "none":
            qos = args["QualityofService"]
            self.qos = f"#SBATCH --qos={qos}"
        self.Software = args["SoftwareLines"]
    def set_dependancy(self, ID:int):
        """
        Can be used to chain slurm Calculations together.

        Args:
            ID (int): SLURM ID

        Attributes:
            dependency (str): SLURM Dependency line
        """
        self.dependency = f"#SBATCH --dependency=afterok:{ID}"
    def set_arrayJob(self, ArrayJob:str, Length:int):
        """
        Hacky work around for HPC systems with Node exclusivity so you can run multiple jobs per node in a single calculation.

        Args:
            ArrayJob (str): yes or no
            Length (int): Size of array
        
        Attributes:
            ArrayLength (int): Reduced array size
            JobsPerNode (int): Number of mini jobs to be run as a single job.

        """
        self.ArrayJob = ArrayJob
        if "archer" in self.HostName.casefold():
            sublength = math.ceil(Length/(self.NodeCores/self.Cores))
            self.ArrayLength = sublength
            self.JobsPerNode = math.floor(self.NodeCores/self.Cores)
        else:
            self.ArrayLength = Length
            self.JobsPerNode = 1    # If Nodes arent exclusive, this can be set to 1 and it will run a normal array job.
    # def set_software(self, Software:str):
    #     self.Software = Software
    def set_accountInfo(self, Qos:str, Account:str):
        """
        sorts out the SLURM user account info. 

        Args:
            Qos (str): SLURM QoS
            Account (str): SLURM Account

        """
        self.qos = Qos
        self.account = Account
    def set_IDNumber(self, ):
        """_summary_
        Attribute:
            ID_Number (int): 
        """
        if "archer" in self.HostName.casefold():
            Number = 1
        elif "sulis" in self.HostName.casefold():
            Number = 4
        elif "ada" in self.HostName.casefold():
            Number = 4
        else:
            Number = 1
            print("WARNING, The HPC that you are using is not in my database... ")
            print("The runnner.sh script may not work as intended and so you may have to tune the sed command to grab the slurm ID")
        self.ID_Number = Number

class MolClass:
    """Class for containing molecular information including coordinates and charges. This is used for benchmarking only
    
    Attributes:
        name (str): name of molecule
        charge (int): net charge of molecule
        spin (int): total spin of molecule
        element (list): list of elements ( in order of coordinates)
        x (list): list of x coordinates
        y (list): list of y coordinates
        z (list): list of z coordinates
        nat (int): Total number of atoms
    """
    def __init__(self, name:str):
        """
        Gets the name of the molecule

        Args:
            name (str): Name of the molecule

        """
        self.name = name
    def set_charge(self, charge: int, spin: int):
        """
        Sets the charge and spin of the molecule

        Args:
            charge (int): Net charge of the molecule (Only whole charge supported)
            spin (int): Spin state of the molecule as defined by ORCA

        """
        assert type(charge) is int, "Charge must be an integer, partial charges are not supported"
        self.charge = charge
        assert type(spin) is int, "Spin must be an integer, partial spins are not supported"
        self.spin = spin
    def add_coordinates(self, at: list, x: list, y: list, z: list, nat: int):
        """
        Adds all structural information to the molecule
        Args:
            at (list): list of elements ( in order of coordinates)
            x (list): list of x coordinates
            y (list): list of y coordinates
            z (list): list of z coordinates
            nat (int): Total number of atoms

        """
        assert len(x) == len(y) == len(z) == len(at) == nat, "Coordinates are broken... x,y&z coordinates are different lengths"
        self.element = at
        self.x = x
        self.y = y
        self.z = z
        self.nat = nat

class ORCAClass:
    """Class for ORCA input information

    Attributes:
        dispersion (str): Any dispersion correction supported by ORCA
        restart (str): Whether to start from the .gbw restart file with ORCA
        path (str): path to the executable
        method (str): (DFT functional, or CCSD or HF)
        basis (str): Any basis set supported by ORCA
        cores (int): Number of CPU cores to use.
        convergence (str): Convergence threshold Threshold
        grid (str): ORCA grid size
        dificulty (str): NormalConv should work for most jobs. 
        calculation (str): QM claculation type. 
        extras (str): Any other commands to parse to the ORCA input.
    """
    dispersion = ""
    restart = "NOAUTOSTART"
    def __init__(self, path:str):
        """
        Sets the path to ORCA
        Args:
            path (str): path to the executable

        """
        self.path = path
    def set_method(self, method:str):
        """
        sets the ORCA QM method (DFT functional, or CCSD or HF)

        Args:
            method (str): (DFT functional, or CCSD or HF)

        """
        self.method = method
    def set_basis(self, basis:str):
        """
        Sets the basis set to be used

        Args:
            basis (str): Any basis set supported by ORCA

        """
        self.basis = basis
    def set_dispersion(self, disp: str):
        """
        Sets the dispersion correction

        Args:
            disp (str): Any dispersion correction supported by ORCA

        """
        self.dispersion = disp
    def set_cores(self, cores: int):
        """
        Sorts out parallel calculations in ORCA
        
        Args:
            cores (int): Number of CPU cores to use.

        """
        self.cores = cores
    def set_convergence(self, conv: str):
        """Set the Convergence threshold

        Args:
            conv (str): Threshold
        """
        self.convergence = conv
    def change_autostart(self, restart: str):
        """Whether to use .gbw restart files

        Args:
            restart (str): AUTOSTART or NOAUTOSTART

        """
        self.restart = restart
    def set_grid(self, grid: str):
        """SCF grid to use

        Args:
            grid (str): ORCA grid size

        """
        self.grid = grid
    def set_dificulty(self, diff:str):
        """Sets how aggressive the SCF algorithm should be. 

        Args:
            diff (str): NormalConv should work for most jobs. 

        """
        self.dificulty = diff
    def set_type(self, caltype:str):
        """QM calculation type (SP, opt, force ect. )

        Args:
            caltype (str): QM claculation type. 

        """
        self.calculation = caltype
    def set_extras(self, extras:str):
        """Any other things to parse to ORCA

        Args:
            extras (str): Any other commands to parse to the ORCA input.

        """
        self.extras = extras

class QMCalcClass:
    """QM ORCA class for running standalone QM calculations in benchmarking, and containing its energies

    Attributes:
        SCFEnergy (float): ORCA SCF energy
        vdw (float): ORCA VdW energy
        TotalEnergy (float): ORCA Total energy
        molecule (str): Name of molecule
        method (str): QM Method
        basis (str): QM basis set
        dispersion (str): QM dispersion correction
        path (str): Path to the calculation
        time (float): Sets total calculation time
        reason (str): Reason for failed calculation
        runline (str): Line containing command to run the calculation.
    """
    SCFEnergy = 0
    vdw = 0
    TotalEnergy = 0
    def __init__(self, molecule: str, method: str, basis: str, dispersion: str):
        """Initialises a calculation

        Args:
            molecule (str): Name of molecule
            method (str): QM Method
            basis (str): QM basis set
            dispersion (str): QM dispersion correction

        """
        self.molecule = molecule
        self.functional = method
        self.basis = basis
        self.dispersion = dispersion
    def set_path(self, path: str): 
        """Path to the actual calculation

        Args:
            path (str): Path to the calculation
        """
        self.path = path
    def set_time(self, time:float):
        """Sets Calculation time

        Args:
            time (float): Sets total calculation time
        """
        self.time = time
    def set_scfenergy(self, scfenergy: float):
        """Gets the ORCA SCF energy

        Args:
            scfenergy (float): ORCA SCF energy
        """
        self.SCFEnergy = scfenergy
    def set_vdw(self, vdwenergy: float):
        """ORCA VdW energy

        Args:
            vdwenergy (str): ORCA VdW Energy
        """
        self.vdw = vdwenergy
    def set_TotalEnergy(self, totenergy:float):
        """Sets Total energy of calculation

        Args:
            totenergy (float): ORCA total energy calculation
        """
        self.TotalEnergy = totenergy
    def set_reason(self, reason:str):
        """Gets the reason for calculation failure

        Args:
            reason (str): Reason for failed calculation

        """
        self.reason = reason
    def set_runline(self, line: str):
        """Line that runs the calculation

        Args:
            line (str): Line containing command to run the calculation.

        """
        self.runline = line

class ReactionClass:
    """ Class containing reaction data for two step processes. 
    Attributes:
        Reactant (str): Name of molecule 1
        Reactant_Energy (float): Total energy corresponding to molecule 1
        Product (str): Name of molecule 2
        Product_Energy (float): Total energy corresponding to molecule 2
        functional (str): QM Method
        basis (str): QM Basis Set
        dispersion (str): QM Dispersion correction
        deltaE (float): energy2-energy1
        Timing (float): Total time to calculate both molecules
    """
    def __init__(self, name1:str, energy1:float, name2:str, energy2:float, method:str, basis:str, dispersion:str):
        """Initialises with energies

        Args:
            name1 (str): Name of molecule 1
            energy1 (float): Total energy corresponding to molecule 1
            name2 (str): Name of molecule 2
            energy2 (float): Total energy corresponding to molecule 2
            method (str): QM Method
            basis (str): QM Basis Set
            dispersion (str): QM Dispersion correction
        
        Attributes:
            deltaE (float): energy2-energy1
        """
        self.Reactant = name1
        self.Reactant_Energy = energy1
        self.Product = name2
        self.Product_Energy = energy2
        self.functional = method
        self.basis = basis
        self.dispersion = dispersion
        self.deltaE = energy2 - energy1
    def add_timings(self, Time1:float, Time2:float):
        """ Calculates the total time for the reaction calculation

        Args:
            Time1 (float): Time to calculate molecule 1
            Time2 (float): Time to calculate molecule 2

        Attributes:
            Timing (float): Total time to calculate both molecules
        """
        self.Timing = Time1 + Time2
    def add_Error(self, err:float):
        """
        Adds the absolute error of the calculation to the data class

        Args:
            err (float): Absolute error.

        """
        self.error = err
    
# class JobClass:
#     """Class for containing job specific information
#     Attributes:
#         WorkDir (str): Path to location of calculation
#         JobType (str): Type of calculation to perform
#         globals.verbosity (int): Level of globals.verbosity  (0: Errors, 1: Warnings, 2: Info, 3: Debug)
#         Stage (str): Stage of calculation.

#     TODO:
#         Need to depreciate this class..     
#     """
#     def __init__(self, args: Namespace):
#         """
#         JobClass init.

#         Args:
#             args (Namespace): ArgParse user inputs to contain the Job information

#         Attributes:
#             WorkDir (str): Path to location of calculation
#             JobType (str): Type of calculation to perform
#             globals.verbosity (int): Level of globals.verbosity  (0: Errors, 1: Warnings, 2: Info, 3: Debug)
#             Stage (str): Stage of calculation.
#         """
#         self.WorkDir = args.WorkDir
#         self.JobType = args.JobType.casefold()
#         self.globals.verbosity = args.globals.verbosity
#         self.Stage = args.Stage.casefold()