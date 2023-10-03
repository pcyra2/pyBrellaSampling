import pandas as pd
import numpy as np
import math as math
from pyBrellaSampling.utils import MM_DefaultVars as MMVars
import os

class UmbrellaClass:
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
    def __init__(self, args):
        self.WorkDir = args.WorkDir
        self.JobType = args.JobType.casefold()
        self.Verbosity = args.Verbosity
        self.Stage = args.Stage.casefold()

class LabelClass:
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

class DataClass:
    dat = []
    def __init__(self, Name, ):
        self.name = Name
    def add_data(self, name, window, data):
        self.dat.append(pd.DataFrame(
            data={"Name": name, "Window": window, "Data": data}))
    def __repr__(self):
        return f"{self.name}: \n{self.dat}"

class BondClass:
    def __init__(self, atom1, atom2, name, threshold):
        self.at1 = int(atom1)
        self.at2 = int(atom2)
        self.name = name
        self.thresh = threshold
    def __repr__(self):
        return f"Atoms: {self.at1} {self.at2}, name: {self.name}, threshold: {self.thresh}"

class DihedralClass:
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
    def __init__(self, Name, Force, type="discrete"):
        self.Name = Name
        self.Force = Force
        self.Type = type

class QMClass:
    def __init__(self, args):
        self.QMpath = args.QmPath
        self.QMSel = args.QmSelection
        self.Charge = args.QmCharge
        self.Spin = args.QmSpin
        self.Method = args.QmMethod
        self.Basis = args.QmBasis
        self.QMExtra = args.QmArgs

class MMClass:
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
            qmParamPDB              "../syst-qm.pdb"
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
            qmSoftware              "custom"
            qmExecPath              "WorkDir/{index}/runORCA.py"
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
        if args.Account != "None":
            self.account = f"#SBATCH --account={args.Account}"
        if args.QualityofService != "None":
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

