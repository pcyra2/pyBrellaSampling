import pandas as pd


class UmbrellaClass:
    def __init__(self, args, Min, bins, Start, Width, ):
        self.atom1 = args.AtomMask.split(",")[0]
        self.atom2 = args.AtomMask.split(",")[1]
        self.atom3 = args.AtomMask.split(",")[2]
        self.atom4 = args.AtomMask.split(",")[3]
        self.Min = Min
        self.Bins = bins
        self.Start = Start
        self.Width = Width
    def add_bins(self, BinVals):
        self.BinVals = BinVals
    def add_start(self, StartBin):
        self.StartBin = StartBin


class JobClass:
    def __init__(self, args):
        self.WorkDir = args.WorkDir
        self.JobType = args.JobType
        self.Verbosity = args.Verbosity

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
        self.at1 = atom1
        self.at2 = atom2
        self.name = name
        self.thresh = threshold
    def __repr__(self):
        return f"Atoms: {self.at1} {self.at2}, name: {self.name}, threshold: {self.thresh}"

class DihedralClass:
    def __init__(self, atom1, atom2, atom3, atom4, name, target1, t1name, target2, t2name):
        self.at1 = atom1
        self.at2 = atom2
        self.at3 = atom3
        self.at4 = atom4
        self.name = name
        self.target1 = target1
        self.target1Name = t1name
        self.target2 = target2
        self.target2Name = t2name
    def __repr__(self):
        return f"Atoms: {self.at1} {self.at2} {self.at3} {self.at4}, Name: {self.name}, {self.target1Name} = {self.target1}, {self.target2Name} = {self.target2}"

class WhamClass:
    def __init__(self, Name, Force):
        self.Name = Name
        self.Force = Force

class CalcClass:
    CellVec = 150.2205127
    Temp = 300
    QM = "off"
    CutOff = 8.0
    def __init__(self,args):
        self.NamdPath = args.MDCPUPath
        self.GPUNamd = args.MDGPUPath
        self.QMpath = args.QmPath
        self.QMSel = args.QmSelection
        self.Charge = args.QmCharge
        self.Spin = args.QmSpin
        self.Method = args.QmMethod
        self.Basis = args.QmBasis
        self.Threads = args.CoresPerJob
        self.Memory = args.MemoryPerJob
    def Job_Name(self, Name,):
        self.Name = Name
    def Set_OutFile(self, OutFile):
        self.OutFile = OutFile
    def Set_Force(self, Force):
        self.Force = Force
    def Set_Id(self, Id):
        self.Id = Id
    def Change_Cell(self, CellVec):
        self.CellVec = CellVec
    def Set_QM(self, val):
        self.QM = val
    def Set_Ensemble(self, Ensemble):
        self.Ensemble = Ensemble
    def Set_Temp(self, Temp):
        self.Temp = Temp
    def Set_Length(self, Steps, TimeStep):
        self.Steps = Steps
        self.TimeStep = TimeStep
    def Set_Outputs(self, TimeOut, RestOut, TrajOut):
        self.TimeOut = TimeOut
        self.RestOut = RestOut
        self.TrajOut = TrajOut
    def Set_Shake(self, Shake):
        self.Shake = Shake