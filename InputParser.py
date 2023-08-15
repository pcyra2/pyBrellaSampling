import pyBrellaSampling.utils as utils
import argparse as ap

### Aim to use case insensitive variables for the end user, especially in the input files.
### Use ".casefold()" when comparing strings and ensure new variables are not differing by case...
### To add a new variable, add the variable name and default value to the relevant array first,
### Then add the variable to the arg_parse function, specifying the default to the dictionary.

### Argument priority: Commandline > Input file > default
### Final inputs are all combined into the "args" variable which is a argparse namespace.

def VariableParser(sysargs):
    JobDict = JobInput("./Job.conf")
    WorkDir = JobDict["WorkDir"]
    ComputeDict = ComputeInput(f"{WorkDir}Compute.conf")
    MMDict = MMInput(f"{WorkDir}MM.conf")
    QMDict = QMInput(f"{WorkDir}QM.conf")
    UmbrellaDict = UmbrellaInput(f"{WorkDir}Umbrella.conf")
    FileDict = {**JobDict, **ComputeDict, **MMDict, **QMDict, **UmbrellaDict}
    args = arg_parse(FileDict,sysargs)
    if args.JobType.casefold() == "inpfile":
        InputFileGen(args)
    return args

def InputFileGen(args):
    argsDict = vars(args)
    JobDict = JobInput(f"{args.WorkDir}Job.conf")
    with open(f"{args.WorkDir}Job.conf", "w") as f:
        for i in JobDict.keys():
            print(f"{i}={argsDict[i]}", file=f)
    ComputeDict = ComputeInput(f"{args.WorkDir}Compute.conf")
    with open(f"{args.WorkDir}Compute.conf", "w") as f:
        for i in ComputeDict.keys():
            print(f"{i}={argsDict[i]}", file=f)
    MMDict = MMInput(f"{args.WorkDir}MM.conf")
    with open(f"{args.WorkDir}MM.conf", "w") as f:
        for i in MMDict.keys():
            print(f"{i}={argsDict[i]}", file=f)
    QMDict = QMInput(f"{args.WorkDir}QM.conf")
    with open(f"{args.WorkDir}QM.conf", "w") as f:
        for i in QMDict.keys():
            print(f"{i}={argsDict[i]}", file=f)
    UmbrellaDict = UmbrellaInput(f"{args.WorkDir}Umbrella.conf")
    with open(f"{args.WorkDir}Umbrella.conf", "w") as f:
        for i in UmbrellaDict.keys():
            print(f"{i}={argsDict[i]}", file=f)

def arg_parse(dict, sysargs):
    parser = ap.ArgumentParser(description="Commandline arguments. This method of calculation input is being deprecated. Please do not use.")
    ### Core Job arguments
    Core = parser.add_argument_group("Core Job Arguments")
    Core.add_argument('-wd', '--WorkDir', type=str,
                        help="Home location for the calculations", default=dict["WorkDir"])
    Core.add_argument('-jt', '--JobType', type=str,
                        help="Type of calculation to run", default=dict["JobType"])
    Core.add_argument('-v', '--Verbosity', type=int,
                        help="Verbosity: 0 = none, 1 = info", default=dict["Verbosity"])
    Core.add_argument('-dr', '--DryRun', type=str,
                        help="Indicates whether programs are executed or not", default=dict["DryRun"])

    Compute = parser.add_argument_group("Compute Arguments")
    ### Compute Arguments
    Compute.add_argument('-cores', '--CoresPerJob', type=int,
                        help="Number of cores per individual calculation", default=dict["CoresPerJob"])
    Compute.add_argument('-mem','--MemoryPerJob', type=int,
                        help="Gb of memory per individual calculation", default=dict["MemoryPerJob"])
    Compute.add_argument('-MaxCalc', '--MaxStepsPerCalc', type=int,
                         help="The maximum number of steps per calculation. splits jobs into sub-steps. useful for short wall times. 0 == No cap.",
                         default=dict["MaxStepsPerCalc"])

    ### MM Arguments
    MM = parser.add_argument_group("Molecular Dynamics Arguments")
    MM.add_argument('-MDcpu', '--MDCPUPath', type=str,
                        help="Path to NAMD CPU executable", default=dict["MDCPUPath"])
    MM.add_argument('-MDgpu', '--MDGPUPath', type=str,
                        help="Path to NAMD GPU executable", default=dict["MDGPUPath"])

    ### QM Arguments
    QM = parser.add_argument_group("QM Arguments")
    QM.add_argument('-qp', '--QmPath', type=str,
                        help="Path to QM software", default=dict["QmPath"])
    QM.add_argument('-qsel', '--QmSelection', type=str,
                        help="Selection algebra for QM atoms", default=dict["QmSelection"])
    QM.add_argument('-qc', '--QmCharge', type=int,
                        help="Charge of QM region", default=dict["QmCharge"])
    QM.add_argument('-qspin', '--QmSpin', type=int,
                        help="Spin of QM region", default=dict["QmSpin"])
    QM.add_argument('-qm', '--QmMethod', type=str,
                        help="Qm method", default=dict["QmMethod"])
    QM.add_argument('-qb', '--QmBasis', type=str,
                        help="QM basis set", default=dict["QmBasis"])
    QM.add_argument('-qargs', '--QmArgs', type=str, help="Extra arguments for ORCA calculation", default=dict["QmArgs"])

    ### Umbrella Arguments
    Umbrella = parser.add_argument_group("Umbrella Sampling arguments")
    Umbrella.add_argument('-min', '--UmbrellaMin', type=float,
                        help="Minimum Umbrella distance", default=dict["UmbrellaMin"])
    Umbrella.add_argument('-width', '--UmbrellaWidth', type=float,
                        help="Umbrella bin width in Angstroms or degrees", default=dict["UmbrellaWidth"])
    Umbrella.add_argument('-bins', '--UmbrellaBins', type=int,
                        help="Number of umbrella bins", default=dict["UmbrellaBins"])
    Umbrella.add_argument('-pf', '--PullForce', type=float,
                        help="Force for pulls in KCal A-2", default=dict["PullForce"])
    Umbrella.add_argument('-f', '--ConstForce', type=float,
                        help="Force for standard Umbrella runs", default=dict["ConstForce"]) ### NAMD uses 1/2 k rather than just k
    Umbrella.add_argument('-sd', '--StartDistance', type=float,
                        help="Distance of initial simulation", default=dict["StartDistance"])
    Umbrella.add_argument('-mask', '--AtomMask', type=str,
                        help="Mask for the restrained atoms.", default=dict["AtomMask"])
    Umbrella.add_argument('-stg', '--Stage', type=str,
                        help="Stage of ummbrella simulation", default=dict["Stage"])
    Umbrella.add_argument('-wf', '--WhamFile', type=str,
                        help="Name prefix of wham data.(XXX.i.colvars.traj", default=dict["WhamFile"])
    ### Parse commandline arguments
    args = parser.parse_args(sysargs)
    return args



def JobInput(path):
    InpVars = ["WorkDir", "JobType", "Verbosity", "DryRun"]
    InpValues = ["./", None, 0, "True"]
    assert len(InpVars) == len(InpValues)
    try:
        lines = utils.file_read(path)
    except FileNotFoundError:
        print("WARNING, No config found for Job input, This is a bad idea... Using defaults.")
        Dict = {}
        for i in range(len(InpVars)):
            Dict[InpVars[i]] = InpValues[i]
        return Dict
    for line in lines:
        words = line.split("=")
        for i in range(len(InpVars)):
            if words[0].casefold() == InpVars[i].casefold():
                InpValues[i] = words[1].replace("\n","")
    Dict = {}
    for i in range(len(InpVars)):
        Dict[InpVars[i]] = InpValues[i]
    return Dict

def ComputeInput(path):
    InpVars = ["CoresPerJob", "MemoryPerJob", "MaxStepsPerCalc"]
    InpValues = [10, 10, 0]
    assert len(InpVars) == len(InpValues)
    try:
        lines = utils.file_read(path)
    except FileNotFoundError:
        print("WARNING, No config found for Compute input, Using defaults.")
        Dict = {}
        for i in range(len(InpVars)):
            Dict[InpVars[i]] = InpValues[i]
        return Dict
    for line in lines:
        words = line.split("=")
        for i in range(len(InpVars)):
            if words[0].casefold() == InpVars[i].casefold():
                InpValues[i] = words[1].replace("\n","")
    Dict = {}
    for i in range(len(InpVars)):
        Dict[InpVars[i]] = InpValues[i]
    return Dict

def MMInput(path):
    InpVars = ["MDCPUPath", "MDGPUPath"]
    InpValues = ["/gpfs01/software/NAMD_2.13_Linux-x86_64-multicore/namd2", "/home/pcyra2/Downloads/NAMD_Git-2021-09-30_Linux-x86_64-multicore-CUDA/namd2"]
    assert len(InpVars) == len(InpValues)
    try:
        lines = utils.file_read(path)
    except FileNotFoundError:
        print("WARNING, No config found for MM input, Using defaults.")
        Dict = {}
        for i in range(len(InpVars)):
            Dict[InpVars[i]] = InpValues[i]
        return Dict
    for line in lines:
        words = line.split("=")
        for i in range(len(InpVars)):
            if words[0].casefold() == InpVars[i].casefold():
                InpValues[i] = words[1].replace("\n","")
    Dict = {}
    for i in range(len(InpVars)):
        Dict[InpVars[i]] = InpValues[i]
    return Dict

def QMInput(path):
    InpVars = ["QmPath", "QmSelection", "QmCharge", "QmSpin", "QmMethod", "QmBasis", "QmArgs"]
    InpValues = ["/gpfs01/home/pcyra2/Software/ORCA/orca", "resname CTN POP MG", 1, 1, "PBE", "6-31G*", "D3BJ TightSCF CFLOAT NormalSCF"]
    assert len(InpVars) == len(InpValues)
    try:
        lines = utils.file_read(path)
    except FileNotFoundError:
        print("WARNING, No config found for QM input, Using defaults.")
        Dict = {}
        for i in range(len(InpVars)):
            Dict[InpVars[i]] = InpValues[i]
        return Dict
    for line in lines:
        words = line.split("=")
        for i in range(len(InpVars)):
            if words[0].casefold() == InpVars[i].casefold():
                InpValues[i] = words[1].replace("\n","")
    Dict = {}
    for i in range(len(InpVars)):
        Dict[InpVars[i]] = InpValues[i]
    return Dict

def UmbrellaInput(path):
    InpVars = ["UmbrellaMin", "UmbrellaWidth", "UmbrellaBins", "PullForce", "ConstForce", "StartDistance", "AtomMask", "Stage", "WhamFile"]
    InpValues = [1.3, 0.05, 54, 5000, 150, 1.4, "13716,13731,0,0", "Setup", "prod"]
    assert len(InpVars) == len(InpValues)
    try:
        lines = utils.file_read(path)
    except FileNotFoundError:
        print("WARNING, No config found for Umbrella input, Using defaults.")
        Dict = {}
        for i in range(len(InpVars)):
            Dict[InpVars[i]] = InpValues[i]
        return Dict
    for line in lines:
        words = line.split("=")
        for i in range(len(InpVars)):
            if words[0].casefold() == InpVars[i].casefold():
                InpValues[i] = words[1].replace("\n","")
    Dict = {}
    for i in range(len(InpVars)):
        Dict[InpVars[i]] = InpValues[i]
    return Dict
