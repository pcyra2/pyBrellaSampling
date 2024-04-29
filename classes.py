import pandas as pd
import numpy as np
import math as math
from pyBrellaSampling.utils import MM_DefaultVars as MMVars
import os


class UmbrellaClass:
    """Class for containing umbrella sampling information."""
    def __init__(self, args, Min, bins, Start, Width, ):
        self.atom1 = int(args.AtomMask.split(",")[0])
        self.atom2 = int(args.AtomMask.split(",")[1])
        self.atom3 = int(args.AtomMask.split(",")[2])
        self.atom4 = int(args.AtomMask.split(",")[3])
        self.Min = Min
        self.Bins = bins
        self.Start = Start
        self.Width = Width
        self.PullForce = args.PullForce
        self.ConstForce = args.ConstForce
    def add_bins(self, BinVals):
        self.BinVals = BinVals
    def add_start(self, StartBin):
        self.StartBin = StartBin
    def set_force(self, Force): ### To combine Umbrella for use in Standalone calculations
        self.PullForce = Force
        self.ConstForce = Force

class JobClass:
    """Class for containing job specific information"""
    def __init__(self, args):
        self.WorkDir = args.WorkDir
        self.JobType = args.JobType.casefold()
        self.Verbosity = args.Verbosity
        self.Stage = args.Stage.casefold()

class LabelClass:
    """Class for containing bond and dihedral label information"""
    bond = []
    bondName = []
    bondThresh = []
    dihedral = []
    dihedralName = []
    dihedralTarget1 = []
    dihedralTarget2 = []
    dihedralTarget1Name = []
    dihedralTarget2Name = []
    def __init__(self,  parm):
        self.parm = parm
    def add_bond(self, selection, name, thresh):
        self.bond.append(selection)
        self.bondName.append(name)
        self.bondThresh.append(thresh)
    def add_dihedral(self, selection, name, target1, t1name, target2, t2name):
        self.dihedral.append(selection)
        self.dihedralName.append(name)
        self.dihedralTarget1.append(target1)
        self.dihedralTarget1Name.append(t1name)
        self.dihedralTarget2.append(target2)
        self.dihedralTarget2Name.append(t2name)
    def file_name(self, name):
        self.file = name
    def clear_Vars(self):
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
    """Class for containing generic analysis data."""
    dat = []
    def __init__(self, Name, ):
        self.name = Name
    def add_data(self, name, window, data):
        self.dat.append(pd.DataFrame(
            data={"Name": name, "Window": window, "Data": data}))
    def __repr__(self):
        return f"{self.name}: \n{self.dat}"

class BondClass:
    """Class for containing bond data"""
    def __init__(self, atom1, atom2, name, threshold):
        self.at1 = int(atom1)
        self.at2 = int(atom2)
        self.name = name
        self.thresh = threshold
    def __repr__(self):
        return f"Atoms: {self.at1} {self.at2}, name: {self.name}, threshold: {self.thresh}"

class DihedralClass:
    """Class for containing dihedral data"""
    def __init__(self, atom1, atom2, atom3, atom4, name, target1, t1name, target2, t2name):
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
    """Class for contining wham data"""
    def __init__(self, Name, Force, type="discrete"):
        self.Name = Name
        self.Force = Force
        self.Type = type

class QMClass:
    SelFile="../syst-qm.pdb"
    """Class for containing quantum calculation information"""
    def __init__(self, args):
        self.QMpath = args.QmPath
        self.QMSel = args.QmSelection
        self.Charge = args.QmCharge
        self.Spin = args.QmSpin
        self.Method = args.QmMethod
        self.Basis = args.QmBasis
        self.QMExtra = args.QmArgs
    def set_selfile(self, File):
        self.SelFile = File

class MMClass:
    """Class for containing molecular mechanics calculation information"""
    CellVec = 150.2205127
    CellShape = "oct"
    Temp = MMVars["temperature"]
    CutOff = MMVars["cutoff"]
    PME = "off"
    parmfile = "../complex.parm7"
    ambercoor = "../start.rst7"
    Shake = "none"
    def __init__(self, args):
        self.NamdPath = args.MDCPUPath
        self.GPUNamd = args.MDGPUPath
    def Set_Shake(self, Shake):
        self.Shake = Shake
    def Set_Outputs(self, TimeOut, RestOut, TrajOut):
        self.TimeOut = TimeOut
        self.RestOut = RestOut
        self.TrajOut = TrajOut
    def Set_Temp(self, Temp):
        self.Temp = Temp
    def Set_Ensemble(self, Ensemble):
        self.Ensemble = Ensemble
    def Set_Length(self, Steps, TimeStep=0.05):
        self.Steps = Steps
        self.TimeStep = TimeStep
    def Change_Cell(self, CellVec):
        self.CellVec = CellVec
    def Set_Force(self, Force):
        self.Force = Force
    def Set_Files(self, parm, ambercoor):
        self.parmfile = parm
        self.ambercoor = ambercoor

class CalcClass:
    """Class for containing generic calculation information"""
    QM = "off"
    def __init__(self, args):
        self.Threads = args.CoresPerJob
        self.Memory = args.MemoryPerJob
        self.MaxSteps = args.MaxStepsPerCalc
    def Job_Name(self, Name, ):
        self.Name = Name
    def Set_OutFile(self, OutFile):
        self.OutFile = OutFile
    def Set_Id(self, Id):
        self.Id = Id
    def Set_QM(self, val):
        self.QM = val


class NAMDClass:
    """Class for containing NAMD Variables"""
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
    def __init__(self, Calc, MM,):
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
        if MM.Ensemble == "NVT":
            self.BrensdenPressure = "off"
            self.heating = f"""langevin            on
langevinDamping     5
langevinTemp        {MM.Temp}
langevinHydrogen    off
temperature         {MM.Temp}
"""
            self.runtype = "run"
        elif MM.Ensemble == "NPT":
            self.BrensdenPressure = "on"
            self.heating = f"""langevin            on
langevinDamping     5
langevinTemp        {MM.Temp}
langevinHydrogen    off
temperature         {MM.Temp}
"""
            self.runtype = "run"
        elif "heat" in MM.Ensemble.casefold():
            self.BrensdenPressure = "off"
            self.heating = f"""temperature         0
reassignFreq        {math.floor(MM.Steps/(MM.Temp/0.2))}
reassignIncr        0.2
reassignHold        {MM.Temp}
"""
            self.runtype = "run"
        elif "min" in MM.Ensemble.casefold():
            self.BrensdenPressure = "off"
            self.heating = """temperature             0
langevin            off"""
            self.runtype = "minimize"
        else:
            return AttributeError, f"{MM.Ensemble} not configured."
    def set_pme(self, val="off"):
        self.PME = val
    def set_cellvectors(self, CellVector, CellShape="oct"): ### Currently only supports truncated Octahedron cells
        if CellShape.casefold() == "oct":
            self.cellBasisVector1 = f"{CellVector} 0.0 0.0"
            self.cellBasisVector2 = f"{(-1/3)*CellVector} {(2/3)*np.sqrt(2)*CellVector} 0.0"
            self.cellBasisVector3 = f"{(-1/3)*CellVector} {(-1/3)*np.sqrt(2)*CellVector} {(-1/3)*np.sqrt(6)*CellVector}"
            self.cellOrigin = "0 0 0"
        else:
            return AttributeError, "Cell shape not supported currently, please add this functionality."
    def set_qm(self, Calc, QM, index=0):
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
            qmBaseDir               "/dev/shm/NAMD_{index}"
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

    def set_colvars(self, file, toggle="on"):
        self.colvarlines = f"""# Colvar options:
colvars         {toggle}
colvarsConfig   {file}
"""
    def set_startcoords(self, bincoor, ambercoor="../start.rst7", parm="../complex.parm7"):
        self.bincoor = bincoor
        self.ambercoor = ambercoor
        self.parm = parm

class SLURMClass:
    """Class for containing SLURM variables"""
    account = ""
    qos = ""
    dependency = ""
    def __init__(self, args):
        self.HostName = args.HostName
        self.WallTime = args.MaxWallTime
        self.Cores = args.CoresPerJob
        self.Memory = args.MemoryPerJob
        self.k = args.HostName
        self.Partition = args.Partition
        self.NodeCores = args.MaxCores
        if args.Account.casefold() != "none":
            self.account = f"#SBATCH --account={args.Account}"
        if args.QualityofService.casefold() != "none":
            self.qos = f"#SBATCH --qos={args.QualityofService}"
        self.Software = args.SoftwareLines
    def set_dependancy(self, ID):
        self.dependency = f"#SBATCH --dependency=afterok:{ID}"
    def set_arrayJob(self, ArrayJob, Length):
        self.ArrayJob = ArrayJob
        if "archer" in self.HostName.casefold():
            sublength = math.ceil(Length/(self.NodeCores/self.Cores))
            self.ArrayLength = sublength
            self.JobsPerNode = math.floor(self.NodeCores/self.Cores)
        else:
            self.ArrayLength = Length
            self.JobsPerNode = 1    # If Nodes arent exclusive, this can be set to 1 and it will run a normal array job.
    def set_software(self, Software):
        self.Software = Software
    def set_accountInfo(self, Qos, Account):
        self.qos = Qos
        self.account = Account
    def set_IDNumber(self, ):
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
    """Class for containing molecular information including coordinates and charges"""
    def __init__(self, name):
        self.name = name
    def set_charge(self, charge, spin):
        assert type(charge) is int, "Charge must be an integer, partial charges are not supported"
        self.charge = charge
        assert type(spin) is int, "Spin must be an integer, partial spins are not supported"
        self.spin = spin
    def add_coordinates(self, at, x, y, z, nat):
        assert len(x) == len(y) == len(z) == len(at) == nat, "Coordinates are broken... x,y&z coordinates are different lengths"
        self.element = at
        self.x = x
        self.y = y
        self.z = z
        self.nat = nat


class ORCAClass:
    """Class for ORCA input information"""
    dispersion = ""
    restart = "NOAUTOSTART"
    def __init__(self, path):
        self.path = path
    def set_method(self, method):
        self.method = method
    def set_basis(self, basis):
        self.basis = basis
    def set_dispersion(self, disp):
        self.dispersion = disp
    def set_cores(self, cores):
        self.cores = cores
    def set_convergence(self, conv):
        self.convergence = conv
    def change_autostart(self, restart):
        self.restart = restart
    def set_grid(self, grid):
        self.grid = grid
    def set_dificulty(self, diff):
        self.dificulty = diff
    def set_type(self, caltype):
        self.calculation = caltype
    def set_extras(self, extras):
        self.extras = extras

class QMCalcClass:
    SCFEnergy = 0
    vdw = 0
    TotalEnergy = 0
    def __init__(self, molecule, method, basis, dispersion):
        self.molecule = molecule
        self.functional = method
        self.basis = basis
        self.dispersion = dispersion
    def set_path(self, path): 
        self.path = path
    def set_time(self, time):
        self.time = time
    def set_scfenergy(self, scfenergy):
        self.SCFEnergy = scfenergy
    def set_vdw(self, vdwenergy):
        self.vdw = vdwenergy
    def set_TotalEnergy(self, totenergy):
        self.TotalEnergy = totenergy
    def set_reason(self, reason):
        self.reason = reason
    def set_runline(self, line):
        self.runline = line

class ReactionClass:
    def __init__(self, name1, energy1, name2, energy2, method, basis, dispersion):
        self.Reactant = name1
        self.Reactant_Energy = energy1
        self.Product = name2
        self.Product_Energy = energy2
        self.functional = method
        self.basis = basis
        self.dispersion = dispersion
        self.deltaE = energy2 - energy1
    def add_timings(self, Time1, Time2):
        self.Timing = Time1 + Time2
    def add_Error(self, err):
        self.error = err
    
