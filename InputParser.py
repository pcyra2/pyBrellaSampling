import pyBrellaSampling.utils as utils
import argparse as ap

### Aim to use case insensitive variables for the end user, especially in the input files.
### Use ".casefold()" when comparing strings and ensure new variables are not differing by case...
### To add a new variable, add the variable name and default value to the relevant array first,
### Then add the variable to the arg_parse function, specifying the default to the dictionary.

### Argument priority: Commandline > Input file > default
### Final inputs are all combined into the "args" variable which is a argparse namespace.

def VariableParser(sysargs, JT="Umbrella"):
    if JT == "Umbrella":
        JobDict = JobInput("./Job.conf")
        WorkDir = JobDict["WorkDir"]
        ComputeDict = ComputeInput(f"{WorkDir}Compute.conf")
        MMDict = MMInput(f"{WorkDir}MM.conf")
        QMDict = QMInput(f"{WorkDir}QM.conf")
        UmbrellaDict = UmbrellaInput(f"{WorkDir}Umbrella.conf")
        HPCDict = HPCInput(f"{WorkDir}HPC.conf")
        StandaloneDict = StandaloneJobInput(f"{WorkDir}Standalone.conf")
        if JobDict["JobType"].casefold() == "umbrella":
            FileDict = {**JobDict, **ComputeDict, **MMDict, **QMDict,
                        **HPCDict, **StandaloneDict,**UmbrellaDict}
        elif JobDict["JobType"].casefold() == "Inpfile":
            FileDict = {**JobDict, **ComputeDict, **MMDict, **QMDict, **UmbrellaDict, **HPCDict, **StandaloneDict}
        elif JobDict["JobType"].casefold() == "mm" or JobDict["JobType"].casefold() == "qmmm":
            FileDict = {**JobDict, **ComputeDict, **MMDict, **QMDict, **UmbrellaDict, **HPCDict,**StandaloneDict}
        else:
            FileDict = {**JobDict, **ComputeDict, **MMDict, **QMDict, **UmbrellaDict, **HPCDict,**StandaloneDict}
        args = arg_parse_Umbrella(FileDict,sysargs)
        if args.JobType.casefold() == "inpfile":
            InputFileGen(args, JobType="Umbrella")
        if args.Verbosity >= 2:
            print(vars(args))
        return args
    elif JT == "Standalone":
        JobDict = JobInput("./Job.conf")
        WorkDir = JobDict["WorkDir"]
        ComputeDict = ComputeInput(f"{WorkDir}Compute.conf")
        MMDict = MMInput(f"{WorkDir}MM.conf")
        QMDict = QMInput(f"{WorkDir}QM.conf")
        HPCDict = HPCInput(f"{WorkDir}HPC.conf")
        StandaloneDict = StandaloneJobInput(f"{WorkDir}Standalone.conf")
        FileDict = {**JobDict, **ComputeDict, **MMDict, **QMDict, **HPCDict, **StandaloneDict}
        args = arg_parse_Standalone(FileDict,sysargs)
        if args.JobType.casefold() == "inpfile":
            InputFileGen(args, JobType="Standalone")
        if args.Verbosity >= 2:
            print(vars(args))
        return args

def InputFileGen(args, JobType="Umbrella"):
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
    if JobType == "Umbrella":
        UmbrellaDict = UmbrellaInput(f"{args.WorkDir}Umbrella.conf")
        with open(f"{args.WorkDir}Umbrella.conf", "w") as f:
            for i in UmbrellaDict.keys():
                print(f"{i}={argsDict[i]}", file=f)
    HPCDict = HPCInput(f"{args.WorkDir}HPC.conf")
    with open(f"{args.WorkDir}HPC.conf", "w") as f:
        for i in HPCDict.keys():
            print(f"{i}={argsDict[i]}", file=f)
    if JobType == "Standalone":
        StandaloneDict = StandaloneJobInput(f"{args.WorkDir}Standalone.conf")
        with open(f"{args.WorkDir}Standalone.conf", "w") as f:
            for i in StandaloneDict.keys():
                print(f"{i}={argsDict[i]}", file=f)

def arg_parse_Umbrella(dict, sysargs):
    parser = ap.ArgumentParser(description=f"""Commandline arguments. This method of calculation input is being deprecated. Please do not use.
It is recommended to use -jt inpfile to generate input file templates with default values that you can then edit.""")
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
    Umbrella.add_argument("--StartFile", default=dict["StartFile"], type=str, help="Initial coordinate file if not starting from \"start.rst7\"")

    ### HPC Arguments
    HPC = parser.add_argument_group("HPC/SLURM arguments")
    HPC.add_argument("-MaxTime", "--MaxWallTime", type=int,
                     help="Maximum wall time (Hours) for your jobs (either leave as node max, or set as job length)",
                     default=dict["MaxWallTime"])
    HPC.add_argument("-Host", "--HostName", type=str,
                     help="HostName of the HPC", default=dict["HostName"])
    HPC.add_argument("--Partition", type=str, help="Calculation partition name",
                     default=dict["Partition"])
    HPC.add_argument("--MaxCores", type=int,
                     help="Maximum number of cores available to a node (For array splitting)", default=dict["MaxCores"])
    HPC.add_argument("-QoS", "--QualityofService", type=str,
                     help="Slurm QoS, set to None if not relevant.", default=dict["QualityofService"])
    HPC.add_argument("--Account", type=str,
                     help="Slurm account, (Not username), Set to None if not relevant", default=dict["Account"])
    HPC.add_argument("-Software", "--SoftwareLines", type=str,
                    help="List of commands like \"module load XXX\" to load software. Keep each line surrounded by quotes.",
                    default=dict["SoftwareLines"], nargs="*")

    # Standalone = parser.add_argument_group("Standalone Job arguments")
    # Standalone.add_argument("--Name", type=str,
    #                         default=dict["Name"], help="Name for the calculation")
    # Standalone.add_argument("--Ensemble", type=str,
    #                         choices=["min", "heat", "NVT", "NPT"],
    #                         help="Ensemble for Calculation", default=dict["Ensemble"])
    # Standalone.add_argument("--QM", type=str, choices=["True", "False"],
    #                         default=dict["QM"], help="Whether this is a QMMM calculation or not.")
    # Standalone.add_argument("-st", "--Steps", type=int,
    #                         default=dict["Steps"], help="Number of simulation steps.")
    # Standalone.add_argument("-dt", "--TimeStep", type=float,
    #                         default=dict["TimeStep"], help="Time step for the simulation. We recommend 2 for MM, 0.5 for QMMM")
    # Standalone.add_argument("--ParmFile", type=str,
    #                         default=dict["ParmFile"], help="Parameter file name")
    # Standalone.add_argument("--AmberCoordinates", type=str,
    #                         default=dict["AmberCoordinates"], help="Amber coordinate file name that relates to the parameter file")
    # # Standalone.add_argument("--StartFile", type=str, default=dict["StartFile"], help="Either Amber coordinates or NAMD coordinates. These are the coordinates that it starts from.")
    # Standalone.add_argument("--RestartOut", type=int, default=dict["RestartOut"], help="Frequency to generate a restart file")
    # Standalone.add_argument("--TrajOut", type=int, default=dict["TrajOut"], help="Frequency to add to the trajectory file")
    # Standalone.add_argument("--SMD", type=str, choices=["off", "on"], default=dict["SMD"], help="Wheter to use steered molecular dynamics")
    # Standalone.add_argument("--Force", type=float, default=dict["Force"], help="Force for Steered MD")
    # Standalone.add_argument("--StartValue", type=float, default=dict["StartValue"], help="Start value for SMD")
    # Standalone.add_argument("--EndValue", type=float, default=dict["EndValue"], help="End value for SMD. MAKE == Start if wanting constant.")

    # Standalone.add_argument("", type=, default=dict[""], help="")

    ### Parse commandline arguments
    args = parser.parse_args(sysargs)
    # print(vars(args))
    return args


def arg_parse_Standalone(dict, sysargs):
    parser = ap.ArgumentParser(description=f"""Commandline arguments. This method of calculation input is being deprecated. Please do not use.
It is recommended to use -jt inpfile to generate input file templates with default values that you can then edit.""")
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

    ### HPC Arguments
    HPC = parser.add_argument_group("HPC/SLURM arguments")
    HPC.add_argument("-MaxTime", "--MaxWallTime", type=int,
                     help="Maximum wall time (Hours) for your jobs (either leave as node max, or set as job length)",
                     default=dict["MaxWallTime"])
    HPC.add_argument("-Host", "--HostName", type=str,
                     help="HostName of the HPC", default=dict["HostName"])
    HPC.add_argument("--Partition", type=str, help="Calculation partition name",
                     default=dict["Partition"])
    HPC.add_argument("--MaxCores", type=int,
                     help="Maximum number of cores available to a node (For array splitting)", default=dict["MaxCores"])
    HPC.add_argument("-QoS", "--QualityofService", type=str,
                     help="Slurm QoS, set to None if not relevant.", default=dict["QualityofService"])
    HPC.add_argument("--Account", type=str,
                     help="Slurm account, (Not username), Set to None if not relevant", default=dict["Account"])
    HPC.add_argument("-Software", "--SoftwareLines", type=str,
                    help="List of commands like \"module load XXX\" to load software. Keep each line surrounded by quotes.",
                    default=dict["SoftwareLines"], nargs="*")

    Standalone = parser.add_argument_group("Standalone Job arguments")
    Standalone.add_argument("--Name", type=str,
                            default=dict["Name"], help="Name for the calculation")
    Standalone.add_argument("--Ensemble", type=str,
                            choices=["min", "heat", "NVT", "NPT"],
                            help="Ensemble for Calculation", default=dict["Ensemble"])
    Standalone.add_argument("--QM", type=str, choices=["True", "False"],
                            default=dict["QM"], help="Whether this is a QMMM calculation or not.")
    Standalone.add_argument("-st", "--Steps", type=int,
                            default=dict["Steps"], help="Number of simulation steps.")
    Standalone.add_argument("-dt", "--TimeStep", type=float,
                            default=dict["TimeStep"], help="Time step for the simulation. We recommend 2 for MM, 0.5 for QMMM")
    Standalone.add_argument("--ParmFile", type=str,
                            default=dict["ParmFile"], help="Parameter file name")
    Standalone.add_argument("--AmberCoordinates", type=str,
                            default=dict["AmberCoordinates"], help="Amber coordinate file name that relates to the parameter file")
    Standalone.add_argument("--StartFile", type=str, default=dict["StartFile"], help="Either Amber coordinates or NAMD coordinates. These are the coordinates that it starts from.")
    Standalone.add_argument("--RestartOut", type=int, default=dict["RestartOut"], help="Frequency to generate a restart file")
    Standalone.add_argument("--TrajOut", type=int, default=dict["TrajOut"], help="Frequency to add to the trajectory file")
    Standalone.add_argument("--SMD", type=str, choices=["off", "on"], default=dict["SMD"], help="Wheter to use steered molecular dynamics")
    Standalone.add_argument("--Force", type=float, default=dict["Force"], help="Force for Steered MD")
    Standalone.add_argument("--StartValue", type=float, default=dict["StartValue"], help="Start value for SMD")
    Standalone.add_argument("--EndValue", type=float, default=dict["EndValue"], help="End value for SMD. MAKE == Start if wanting constant.")
    Standalone.add_argument('-mask', '--AtomMask', type=str,
                        help="Mask for the restrained atoms.", default=dict["AtomMask"])
    ### Parse commandline arguments
    args = parser.parse_args(sysargs)
    # print(vars(args))
    return args

def JobInput(path):
    InpVars = ["WorkDir", "JobType", "Verbosity", "DryRun"]
    InpValues = ["./", "inpfile", 0, "True"]
    assert len(InpVars) == len(InpValues)
    try:
        lines = utils.file_read(path)
    except FileNotFoundError:
        print("WARNING, No config found for Job input, This is a bad idea... Using defaults.")
        Dict = {}
        for i in range(len(InpVars)):
            Dict[InpVars[i]] = str(InpValues[i])
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
    InpValues = [10, 10, 1000]
    assert len(InpVars) == len(InpValues)
    try:
        lines = utils.file_read(path)
    except FileNotFoundError:
        print("WARNING, No config found for Compute input, Using defaults.")
        Dict = {}
        for i in range(len(InpVars)):
            Dict[InpVars[i]] = str(InpValues[i])
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
            Dict[InpVars[i]] = str(InpValues[i])
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
            Dict[InpVars[i]] = str(InpValues[i])
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
    InpVars = ["UmbrellaMin", "UmbrellaWidth", "UmbrellaBins", "PullForce", "ConstForce", "StartDistance", "AtomMask", "Stage", "WhamFile", "StartFile"]
    InpValues = [1.3, 0.05, 54, 5000, 150, 1.4, "0,0,0,0", "Setup", "prod", "start.rst7"]
    assert len(InpVars) == len(InpValues)
    try:
        lines = utils.file_read(path)
    except FileNotFoundError:
        print("WARNING, No config found for Umbrella input, Using defaults.")
        Dict = {}
        for i in range(len(InpVars)):
            Dict[InpVars[i]] = str(InpValues[i])
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

def HPCInput(path):
    InpVars = ["MaxWallTime", "HostName","Partition", "MaxCores", "QualityofService", "Account", "SoftwareLines"]
    InpValues = [24, "login.archer2.ac.uk", "standard", 128, None, None, "module load ORCA", ]
    InpVars2 = []
    InpValues2 = []
    assert len(InpVars) == len(InpValues)
    try:
        lines = utils.file_read(path)
    except FileNotFoundError:
        print("WARNING, No config found for HPC input, Using defaults.")
        Dict = {}
        for i in range(len(InpVars)):
            Dict[InpVars[i]] = str(InpValues[i])
        return Dict
    for i in range(len(lines)):
        words = lines[i].split("=",1)
        for j in range(len(InpVars)):
            if words[0].casefold() == InpVars[j].casefold():
                InpValues2.append(words[1].replace("\n",""))
                InpVars2.append(words[0])
    Dict = {}
    for i in range(len(InpVars2)):
        if Dict.get(InpVars2[i]) == None:
            Dict[InpVars2[i]] = InpValues2[i]
            # print(Dict.get(InpVars2[i]))
        else:       ### This allows for multiple iof the same word keyword (i.e. Multiple software lines)
            Vals = Dict.get(InpVars2[i])
            if type(Vals) == str:
                Dict[InpVars2[i]] = [Vals, InpValues2[i]]
            else:
                Vals.append(InpValues2[i])
                Dict[InpVars2[i]] = Vals
            # print(Vals)
    return Dict

def StandaloneJobInput(path):
    InpVars = ["Name", "ParmFile", "AmberCoordinates", "StartFile", "Ensemble",
               "QM", "Steps", "TimeStep", "RestartOut", "TrajOut", "SMD", "Force", "StartValue", "EndValue", "AtomMask"]
    InpValues = ["QMMM_Job","complex.parm7", "start.rst7", "Start.rst7" ,"min",
                 "true", 1000, 0.05, 10, 50, "off", 1, 1, 2, "0,0,0,0"]
    assert len(InpVars) == len(InpValues)
    try:
        lines = utils.file_read(path)
    except FileNotFoundError:
        print("WARNING, No config found for Standalone Job input, Using defaults.")
        Dict = {}
        for i in range(len(InpVars)):
            Dict[InpVars[i]] = str(InpValues[i])
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

def BondsInput(path, Labels):
    try:
        data = utils.file_read(path)
    except FileNotFoundError:
        print("No bond information file found...")
        return Labels
    for lines in data:
        variables = lines.split()
        if "name" in variables[0].casefold():
            pass
        else:
            Labels.add_bond(selection=f"{variables[1]},{variables[2]}", name=variables[0], thresh=float(variables[3]))
    return Labels

def DihedralInput(path, Labels):
    try:
        data = utils.file_read(path)
    except FileNotFoundError:
        print("No dihedral information file found...")
        return Labels
    for lines in data:
        variables = lines.split()
        if "name" in variables[0].casefold():
            pass
        else:
            Labels.add_dihedral(selection=f"{variables[1]},{variables[2]},{variables[3]},{variables[4]}",
                        name=variables[0], target1=float(variables[6]), t1name=variables[5],
                        target2=float(variables[8]), t2name=variables[7])
    return Labels