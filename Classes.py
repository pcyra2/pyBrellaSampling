class CalcClass:
    NamdPath = args.NamdPath
    GPUNamd = "/home/pcyra2/Downloads/NAMD_Git-2021-09-30_Linux-x86_64-multicore-CUDA/namd2"
    QMpath = args.QmPath
    QMSel = args.QmSelection
    Charge = args.QmCharge
    Spin = args.QmSpin
    Method = args.QmMethod
    Basis = args.QmBasis
    Threads = args.NumCores
    CellVec = 150.2205127
    Temp = 300
    QM = "off"
    CutOff = 8.0

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
class UmbrellaClass:
    atom1 = args.AtomMask.split(",")[0]
    atom2 = args.AtomMask.split(",")[1]
    atom3 = args.AtomMask.split(",")[2]
    atom4 = args.AtomMask.split(",")[3]
    def __init__(self, Min, bins, Start, Width,):
        self.Min = Min
        self.Bins = bins
        self.Start = Start
        self.Width = Width
    def add_bins(self, BinVals):
        self.BinVals = BinVals
    def add_start(self, StartBin):
        self.StartBin = StartBin
class JobClass:
    WorkDir = args.WorkDir
    JobType = args.JobType
    Verbosity = args.verbose
class WhamClass:
    def __init__(self, Name, Force):
        self.Name = Name
        self.Force = Force
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