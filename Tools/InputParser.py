import pyBrellaSampling.Tools.utils as utils
import pyBrellaSampling.UserVars.Defaultinputs as UserVars
from pyBrellaSampling.UserVars.HPC_Config import HPC_Config
from pyBrellaSampling.UserVars.SoftwarePaths import ORCA_PATH, NAMD_CPU, NAMD_GPU
from pyBrellaSampling.Tools.classes import LabelClass
import argparse as argparse
import socket

### Aim to use case insensitive variables for the end user, especially in the input files.
### Use ".casefold()" when comparing strings and ensure new variables are not differing by case...
### To add a new variable, add the variable name and default value to the relevant array first,
### Then add the variable to the arg_parse function, specifying the default to the dictionary.

### Argument priority: Input file > Commandline > default
### Final inputs are all combined into the "args" variable which is a argparse namespace.


def UmbrellaInput(sysargs):
    """
    Description of UmbrellaInput

    Args:
        sysargs (list): CLI user variables

    Returns:
        arg_dict (dict): User variables parsed as a dictionary.

    """
    HPC = "None"
    HostName = socket.gethostname()
    for alias, data in HPC_Config.items():
        if HostName == data["HostName"]:
            HPC = alias
            continue
    defaults = UserVars.Umbrella_Inp
    HPC_Conf = HPC_Config[HPC]
    parser = argparse.ArgumentParser(description=f"""Commandline arguments. This method of calculation input is being deprecated. Please do not use.
It is recommended to use --Stage inpfile to generate input file templates with default values that you can then edit.""")
    ### Core Job arguments
    Core = parser.add_argument_group("Core Job Arguments")
    Core.add_argument('-wd', '--WorkDir', type=str,
                        help="Home location for the calculations", default=defaults["WorkDir"])
    Core.add_argument('-i', '--Input', type=str, default="None",
                      help="Global input file, This will overwrite all other variables and input files.")
    Core.add_argument('-v', '--Verbosity', type=int,
                        help="Verbosity: 0 = none, 1 = info", default=defaults["Verbosity"])
    Core.add_argument('-dr', '--DryRun', type=str,
                        help="Indicates whether programs are executed or not", default=defaults["DryRun"])

    Compute = parser.add_argument_group("Compute Arguments")
    ### Compute Arguments
    Compute.add_argument('-cores', '--CoresPerJob', type=int,
                        help="Number of cores per individual calculation", default=defaults["CoresPerJob"])
    Compute.add_argument('-mem','--MemoryPerJob', type=int,
                        help="Gb of memory per individual calculation", default=defaults["MemoryPerJob"])
    Compute.add_argument('-MaxCalc', '--MaxStepsPerCalc', type=int,
                         help="The maximum number of steps per calculation. splits jobs into sub-steps. useful for short wall times. 0 == No cap.",
                         default=defaults["MaxStepsPerCalc"])

    ### MM Arguments
    MM = parser.add_argument_group("Molecular Dynamics Arguments")
    MM.add_argument('-MDcpu', '--MDCPUPath', type=str,
                        help="Path to NAMD CPU executable", default=NAMD_CPU)
    MM.add_argument('-MDgpu', '--MDGPUPath', type=str,
                        help="Path to NAMD GPU executable", default=NAMD_GPU)

    ### QM Arguments
    QM = parser.add_argument_group("QM Arguments")
    QM.add_argument("--QM", type=str, choices=["True", "False"], default=defaults["QM"],
                    help="Whether to use QMMM Umbrella Sampling")
    QM.add_argument("-qf", "--QmFile", type=str, default=defaults["QmFile"],
                        help="Name of file containing QM information.")
    QM.add_argument('-qp', '--QmPath', type=str,
                        help="Path to QM software", default=ORCA_PATH)
    QM.add_argument('-qsel', '--QmSelection', type=str,
                        help="Selection algebra for QM atoms", default=defaults["QmSelection"])
    QM.add_argument('-qc', '--QmCharge', type=int,
                        help="Charge of QM region", default=defaults["QmCharge"])
    QM.add_argument('-qspin', '--QmSpin', type=int,
                        help="Spin of QM region", default=defaults["QmSpin"])
    QM.add_argument('-qm', '--QmMethod', type=str,
                        help="Qm method", default=defaults["QmMethod"])
    QM.add_argument('-qb', '--QmBasis', type=str,
                        help="QM basis set", default=defaults["QmBasis"])
    QM.add_argument('-qargs', '--QmArgs', type=str, help="Extra arguments for ORCA calculation", default=defaults["QmArgs"])

    ### Umbrella Arguments
    Umbrella = parser.add_argument_group("Umbrella Sampling arguments")
    Umbrella.add_argument("--UmbrellaFile", type=str, default=defaults["UmbrellaFile"],
                        help="Name of file containing Umbrella information.")
    Umbrella.add_argument('-min', '--UmbrellaMin', type=float,
                        help="Minimum Umbrella distance", default=defaults["UmbrellaMin"])
    Umbrella.add_argument('-width', '--UmbrellaWidth', type=float,
                        help="Umbrella bin width in Angstroms or degrees", default=defaults["UmbrellaWidth"])
    Umbrella.add_argument('-bins', '--UmbrellaBins', type=int,
                        help="Number of umbrella bins", default=defaults["UmbrellaBins"])
    Umbrella.add_argument('-pf', '--PullForce', type=float,
                        help="Force for pulls in KCal A-2", default=defaults["PullForce"])
    Umbrella.add_argument('-f', '--ConstForce', type=float,
                        help="Force for standard Umbrella runs", default=defaults["ConstForce"]) ### NAMD uses 1/2 k rather than just k
    Umbrella.add_argument('-sd', '--StartDistance', type=float,
                        help="Distance of initial simulation", default=defaults["StartDistance"])
    Umbrella.add_argument('-mask', '--AtomMask', type=str,
                        help="Mask for the restrained atoms.", default=defaults["AtomMask"])
    Umbrella.add_argument('-stg', '--Stage', type=str,
                        help="Stage of umbrella simulation", default=defaults["Stage"])
    Umbrella.add_argument('-af', '--AnalysisFile', type=str,
                        help="Name prefix to perform custom analysis.(XXX.i.colvars.traj", default=defaults["AnalysisFile"])
    Umbrella.add_argument("--StartFile", default=defaults["StartFile"], type=str, help="Initial coordinate file if not starting from \"start.rst7\"")
    Umbrella.add_argument("--ParmFile", default=defaults["ParmFile"], type=str, 
                              help="Amber parameter file")
    Umbrella.add_argument("--EquilLength", default=defaults["EquilLength"], type=int, 
                          help="Length of equilibration in ps (per window)")
    Umbrella.add_argument("--ProdLength", type=int, default=defaults["ProdLength"],
                          help="Length of production umbrellasampling per window in ps.")
    if HPC != "None":
    ### HPC Arguments
        HPC = parser.add_argument_group("HPC/SLURM arguments")
        HPC.add_argument("-MaxTime", "--MaxWallTime", type=int,
                        help="Maximum wall time (Hours) for your jobs (either leave as node max, or set as job length)",
                        default=HPC_Conf["MaxWallTime"])
        HPC.add_argument("-Host", "--HostName", type=str,
                        help="HostName of the HPC", default=HPC_Conf["HostName"])
        HPC.add_argument("--Partition", type=str, help="Calculation partition name",
                        default=HPC_Conf["Partition"])
        HPC.add_argument("--MaxCores", type=int,
                        help="Maximum number of cores available to a node (For array splitting)", default=HPC_Conf["MaxCores"])
        HPC.add_argument("-QoS", "--QualityofService", type=str,
                        help="Slurm QoS, set to None if not relevant.", default=HPC_Conf["QualityofService"])
        HPC.add_argument("--Account", type=str,
                        help="Slurm account, (Not username), Set to None if not relevant", default=HPC_Conf["Account"])
        HPC.add_argument("-Software", "--SoftwareLines", type=str,
                        help="List of commands like \"module load XXX\" to load software. Keep each line surrounded by quotes.",
                        default=HPC_Conf["SoftwareLines"], nargs="*")
    args = parser.parse_args(sysargs)
    arg_dict = vars(args)
    if HPC == "None":
        arg_dict = arg_dict | HPC_Conf
    workdir = arg_dict["WorkDir"]
    qmfile = arg_dict["QmFile"]
    umbfile = arg_dict["UmbrellaFile"]
    inpfile = arg_dict["Input"]
    if arg_dict["Stage"].casefold() == "inpfile": # Generates an Umbrella.inp file in the work directory using default values
        utils.file_2dwrite(f"{workdir}Umbrella.inp", x=list(defaults.keys()), y=list(defaults.values()), delim="=" )
    if arg_dict["QmFile"] != "None":
        qm_input = utils.file_read(f"{workdir}{qmfile}")
        for i in qm_input:
            if "#" in i: # ignore comment lines.
                continue
            var, val = i.split("=")
            var = var.replace(" ", "") # removes random spaces.
            val = val.replace("\n","") # removes newline errors.
            if var in arg_dict:
                arg_dict[var] = val
            else:
                raise ValueError(f"ERROR: Unknown variable: {var} provided in {qmfile}")
    if arg_dict["UmbrellaFile"] != "None":
        umb_input = utils.file_read(f"{workdir}{umbfile}")
        for i in umb_input:
            if "#" in i: # ignore comment lines.
                continue
            var, val = i.split("=")
            var = var.replace(" ", "")
            val = val.replace("\n","")
            if var in arg_dict:
                arg_dict[var] = val
            else:
                raise ValueError(f"ERROR: Unknown variable: {var} provided in {umbfile}")
    if arg_dict["Input"] != "None":
        inp = utils.file_read(f"{workdir}{inpfile}")
        for i in inp:
            if "#" in i: # ignore comment lines.
                continue
            var, val = i.split("=")
            var = var.replace(" ", "")
            val = val.replace("\n","")
            if var in arg_dict:
                arg_dict[var] = val
            else:
                raise ValueError(f"ERROR: Unknown variable: {var} provided in {inpfile}")
    return arg_dict

def BondsInput(path: str, Labels: LabelClass):
    """Reads in Bond information and updates the Label class

    Args:
        path (str): Path of file containing bond information
        Labels (LabelClass): Label class containing labels

    Returns:
        Labels (LabelClass): Updated label class
    """
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

def DihedralInput(path: str, Labels: LabelClass):
    """
    Reads in Dihedral bond information
    Args:
        path (str): Path to file containing dihedral information
        Labels (LabelClass): Label class containing Label information

    Returns:
        Labels (LabelClass): Updated label class with new dihedral information

    """
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

def StandaloneInput(sysargs: list):
    """
    Controls the inputs of a standalone calculation. 

    Args:
        sysargs (list): CLI inputs (see attributes)
    
    Attributes:
        --WorkDir (str): Working Directory
        --Verbosity (int): Verbosity: 0 = none, 1 = info
        --DryRun (str):  Indicates whether programs are executed or not
        --CoresPerJob (int): Number of cores per individual calculation
        --MemoryPerJob (int): Gb of memory per individual calculation
        --MDCPUPath (str): Path to NAMD CPU executable
        --MDGPUPath (str): Path to NAMD GPU executable
        --QmFile (str): Name of file containing QM information.
        --QmPath (str): Path to QM software
        --QmSelection (str): Selection algebra for QM atoms
        --QmCharge (int): Charge of QM region
        --QmSpin (int): Spin of QM region
        --QmMethod (str): Qm method
        --QmBasis (str): QM basis set
        --QmArgs (str): Extra arguments for ORCA calculation
        --HPC (bool): Whether to run on a HPC
        --MaxWallTime (int): Maximum wall time (Hours) for your jobs (either leave as node max, or set as job length)
        --HostName (str): HostName of the HPC
        --Partition (str): Calculation partition name
        --MaxCores (int): Maximum number of cores available to a node (For array splitting)
        --QualityofService (str): Slurm QoS, set to None if not relevant.
        --Account (str): Slurm account, (Not username), Set to None if not relevant
        --SoftwareLines (list): List of commands like "module load XXX" to load software. Keep each line surrounded by quotes.
        --MMFile (str): Name of input file containing MD information
        --Name (str): Name for the calculation
        --Ensemble (str): Ensemble for Calculation
        --QM (bool): Whether this is a QMMM calculation or not.
        --Steps (int): Number of simulation steps.
        --TimeStep (float): Time step for the simulation. We recommend 2 for MM, 0.5 for QMMM
        --ParmFile (str): Parameter file name
        --AmberCoordinates (str): Amber coordinate file name that relates to the parameter file
        --StartFile (str): Either Amber coordinates or NAMD coordinates. These are the coordinates that it starts from.
        --RestartOut (int): Frequency to generate a restart file
        --TrajOut (int): Frequency to add to the trajectory file
        --SMD (bool): Whether to use steered molecular dynamics
        --Force (int): Force for Steered MD
        --StartValue (int): Start value for SMD
        --EndValue (int): End value for SMD. MAKE == Start if wanting constant.
        --AtomMask (str): Mask for the restrained atoms.

    Raises:
        ValueError: If an unknown variable is parsed
    
    Returns:
        arg_dict (dict): Dictionary of user variables. 
    """
    HPC = "None"
    HostName = socket.gethostname()
    for alias, data in HPC_Config.items():
        if HostName == data["HostName"]:
            HPC = alias
            continue
    defaults = UserVars.Standalone_Inp
    HPC_Conf = HPC_Config[HPC]
    parser = argparse.ArgumentParser(description=f"""Commandline arguments. This method of calculation input is being deprecated. Please do not use.
It is recommended to use -jt inpfile to generate input file templates with default values that you can then edit.""")
    ### Core Job arguments
    Core = parser.add_argument_group("Core Job Arguments")
    Core.add_argument('-wd', '--WorkDir', type=str,
                        help="Home location for the calculations", default=defaults["WorkDir"])
    Core.add_argument('-v', '--Verbosity', type=int,
                        help="Verbosity: 0 = none, 1 = info", default=defaults["Verbosity"])
    Core.add_argument('-dr', '--DryRun', type=str,
                        help="Indicates whether programs are executed or not", default=defaults["DryRun"])

    Compute = parser.add_argument_group("Compute Arguments")
    ### Compute Arguments
    Compute.add_argument('-c', '--CoresPerJob', type=int,
                        help="Number of cores per individual calculation", default=defaults["CoresPerJob"])
    Compute.add_argument('-m','--MemoryPerJob', type=int,
                        help="Gb of memory per individual calculation", default=defaults["MemoryPerJob"])

    ### MM Arguments
    MM = parser.add_argument_group("Molecular Dynamics Arguments")
    MM.add_argument('-MDcpu', '--MDCPUPath', type=str,
                        help="Path to NAMD CPU executable", default=NAMD_CPU)
    MM.add_argument('-MDgpu', '--MDGPUPath', type=str,
                        help="Path to NAMD GPU executable", default=NAMD_GPU)

    ### QM Arguments
    QM = parser.add_argument_group("QM Arguments")
    QM.add_argument("-qf", "--QmFile", type=str, default=defaults["QmFile"],
                    help="Name of file containing QM information.")
    QM.add_argument('-qp', '--QmPath', type=str,
                        help="Path to QM software", default=ORCA_PATH)
    QM.add_argument('-qsel', '--QmSelection', type=str,
                        help="Selection algebra for QM atoms", default=defaults["QmSelection"])
    QM.add_argument('-qc', '--QmCharge', type=int,
                        help="Charge of QM region", default=defaults["QmCharge"])
    QM.add_argument('-qspin', '--QmSpin', type=int,
                        help="Spin of QM region", default=defaults["QmSpin"])
    QM.add_argument('-qm', '--QmMethod', type=str,
                        help="Qm method", default=defaults["QmMethod"])
    QM.add_argument('-qb', '--QmBasis', type=str,
                        help="QM basis set", default=defaults["QmBasis"])
    QM.add_argument('-qargs', '--QmArgs', type=str, 
                    help="Extra arguments for ORCA calculation", default=defaults["QmArgs"])

    ### HPC Arguments
    HPC = parser.add_argument_group("HPC/SLURM arguments")
    HPC.add_argument("--HPC", type=bool, default=defaults["HPC"], help="Whether to run on a HPC")
    HPC.add_argument("-MaxTime", "--MaxWallTime", type=int,
                     help="Maximum wall time (Hours) for your jobs (either leave as node max, or set as job length)",
                     default=HPC_Conf["MaxWallTime"])
    HPC.add_argument("-Host", "--HostName", type=str,
                     help="HostName of the HPC", default=HPC_Conf["HostName"])
    HPC.add_argument("--Partition", type=str, help="Calculation partition name",
                     default=HPC_Conf["Partition"])
    HPC.add_argument("--MaxCores", type=int,
                     help="Maximum number of cores available to a node (For array splitting)", default=HPC_Conf["MaxCores"])
    HPC.add_argument("-QoS", "--QualityofService", type=str,
                     help="Slurm QoS, set to None if not relevant.", default=HPC_Conf["QualityofService"])
    HPC.add_argument("--Account", type=str,
                     help="Slurm account, (Not username), Set to None if not relevant", default=HPC_Conf["Account"])
    HPC.add_argument("-Software", "--SoftwareLines", type=str,
                    help="List of commands like \"module load XXX\" to load software. Keep each line surrounded by quotes.",
                    default=HPC_Conf["SoftwareLines"], nargs="*")

    Standalone = parser.add_argument_group("Standalone Job arguments")
    Standalone.add_argument("-i", "--MMFile", type=str, default=defaults["MMFile"],
                            help="Name of input file containing MD information")
    Standalone.add_argument("--Name", type=str,
                            default=defaults["Name"], help="Name for the calculation")
    Standalone.add_argument("--Ensemble", type=str,
                            choices=["min", "heat", "NVT", "NPT"],
                            help="Ensemble for Calculation", default=defaults["Ensemble"])
    Standalone.add_argument("--QM", type=str, choices=["True", "False"],
                            default=defaults["QM"], help="Whether this is a QMMM calculation or not.")
    Standalone.add_argument("-st", "--Steps", type=int,
                            default=defaults["Steps"], help="Number of simulation steps.")
    Standalone.add_argument("-dt", "--TimeStep", type=float,
                            default=float(defaults["TimeStep"]), help="Time step for the simulation. We recommend 2 for MM, 0.5 for QMMM")
    Standalone.add_argument("--ParmFile", type=str,
                            default=defaults["ParmFile"], help="Parameter file name")
    Standalone.add_argument("--AmberCoordinates", type=str,
                            default=defaults["AmberCoordinates"], help="Amber coordinate file name that relates to the parameter file")
    Standalone.add_argument("--StartFile", type=str, default=defaults["StartFile"], 
                            help="Either Amber coordinates or NAMD coordinates. These are the coordinates that it starts from.")
    Standalone.add_argument("--RestartOut", type=int, default=defaults["RestartOut"], 
                            help="Frequency to generate a restart file")
    Standalone.add_argument("--TrajOut", type=int, default=defaults["TrajOut"], 
                            help="Frequency to add to the trajectory file")
    Standalone.add_argument("--SMD", type=str, choices=["True", "False"], default=defaults["SMD"], 
                            help="Wheter to use steered molecular dynamics")
    Standalone.add_argument("--Force", type=float, default=defaults["Force"], 
                            help="Force for Steered MD")
    Standalone.add_argument("--StartValue", type=float, default=defaults["StartValue"], 
                            help="Start value for SMD")
    Standalone.add_argument("--EndValue", type=float, default=defaults["EndValue"], 
                            help="End value for SMD. MAKE == Start if wanting constant.")
    Standalone.add_argument('-mask', '--AtomMask', type=str,
                        help="Mask for the restrained atoms.", default=defaults["AtomMask"])
    args = parser.parse_args(sysargs)
    arg_dict = vars(args)
    workdir = arg_dict["WorkDir"]
    qmfile = arg_dict["QmFile"]
    mmfile = arg_dict["MMFile"]
    # print(arg_dict["DryRun"])
    if arg_dict["QmFile"] != "None":
        qm_input = utils.file_read(f"{workdir}{qmfile}")
        for i in qm_input:
            if "#" in i: # ignore comment lines.
                continue
            var, val = i.split("=")
            var = var.replace(" ", "") # removes random spaces.
            val = val.replace("\n","") # removes newline errors/
            if var in arg_dict:
                arg_dict[var] = val
            else:
                raise ValueError(f"ERROR: Unknown variable: {var} provided in {qmfile}")
    if arg_dict["MMFile"] != "None":
        mm_input = utils.file_read(f"{workdir}{mmfile}")
        for i in mm_input:
            if "#" in i: # ignore comment lines.
                continue
            var, val = i.split("=")
            var = var.replace(" ", "")
            val = val.replace("\n","")
            if var in arg_dict:
                arg_dict[var] = val
            else:
                raise ValueError(f"ERROR: Unknown variable: {var} provided in {mmfile}")
    # utils.print_attributes(arg_dict) # Prints out all vars for use in docstring generation
    return arg_dict
    
def BenchmarkInput(sysargs: list):
    """
    Parses command line arguments for a benchmark calculation

    Args:
        sysargs (list): Commandline variables to be parsed to ArgParse
    
    Returns:
        args (dict): Dictionary containing all user-defined variables.
    """
    Benchmark_Inp = UserVars.Benchmark_Inp
    parser = argparse.ArgumentParser(description=f"""CLI interface for running a QM benchmark in ORCA""")
    gen = parser.add_argument_group("Generic Inputs")
    gen.add_argument("-v", "--verbosity", type=int, help="Control the verbosity of the job", choices=[0,1,2,3 ], default=Benchmark_Inp["Verbosity"])
    gen.add_argument("-wd","--WorkDir", type=str, help="Working directory path", default=Benchmark_Inp["WorkDir"])
    gen.add_argument("-c", "--Cores", type=int, help="The number of CPU cores to give each ORCA Job.", default=Benchmark_Inp["Cores"])
    gen.add_argument("-dr", "--DryRun", type=str, help="Should the calculation be performed?", default=Benchmark_Inp["DryRun"])
    gen.add_argument("-stg", "--Stage", type=str, help="The stage of the calculation",choices=["init", "calc", "analysis",], default=Benchmark_Inp["Stage"])
    gen.add_argument("-i", "--input", type=str, help="Optinonal input file rather than using the CLI interface.", default=None)
    bench = parser.add_argument_group("Benchmark Inputs")
    bench.add_argument("-type", "--BenchmarkType", type=str, help="Define the type of benchmark to perform", choices=["Energy", "Gradient", "Structure"], default=Benchmark_Inp["BenchmarkType"])
    bench.add_argument("-cd", "--CoordinateLoc", type=str, help="Location for the structures to perform the benchmark on", default=Benchmark_Inp["CoordinateLoc"])
    bench.add_argument("-reactions", "--ReactionList", type=str, 
                       help="Name of the file containing the list of reaction steps. This must be line delimited single step reactions with the name of the \
                       structure being identical of that to the file in the coordinates folder.", default=Benchmark_Inp["ReactionList"])
    qm = parser.add_argument_group("ORCA specific inputs")
    qm.add_argument("--Path", type=str, help="Path to ORCA executable", default=ORCA_PATH)
    qm.add_argument("--SCF",type=str, help="The orca scf convergence setting", choices=["NORMALSCF", "TIGHTSCF", "VERYTIGHTSCF", "EXTREMESCF"], default=Benchmark_Inp["SCF"])
    qm.add_argument("--Grid", type=str, help="Specify the orca integration grid", choices=["DEFGRID1", "DEFGRID2", "DEFGRID3"], default=Benchmark_Inp["Grid"])
    qm.add_argument("--Restart", type=str, help="Togle whether to use a restart file or not", choices=["AUTOSTART", "NOAUTOSTART"], default=Benchmark_Inp["Restart"])
    qm.add_argument("--Convergence", type=str, help="Chose th orca convergence stratergy", choices=["EasyConv", "NormalConv", "SlowConv", "VerySlowConv", "ForceConv"], default=Benchmark_Inp["Convergence"])
    qm.add_argument("--Extras", type=str, help="Extra ORCA commands to use for all calculations (default is \"MINIPRINT\")", default=Benchmark_Inp["Extras"])
    args = parser.parse_args(sysargs)
    arg_dict = vars(args)
    if arg_dict["input"] != None:
        print("ERROR: Not implemented yet.. ")
    return arg_dict


# def argStandalone(sysargs: list):
#     """Parses the cli variables for a standalone calculation
    
#     Args: 
#         sysargs (list):
        
#     Returns:
#         parser (argparse.ArgumentParser):
#     """
#     HPC = "None"
#     HostName = socket.gethostname()
#     for alias, data in HPC_Config.items():
#         if HostName == data["HostName"]:
#             HPC = alias
#             continue
#     defaults = UserVars.Standalone_Inp
#     HPC_Conf = HPC_Config[HPC]
#     parser = argparse.ArgumentParser(description=f"""Commandline arguments. This method of calculation input is being deprecated. Please do not use.
# It is recommended to use -jt inpfile to generate input file templates with default values that you can then edit.""")
#     ### Core Job arguments
#     Core = parser.add_argument_group("Core Job Arguments")
#     Core.add_argument('-wd', '--WorkDir', type=str,
#                         help="Home location for the calculations", default=defaults["WorkDir"])
#     Core.add_argument('-v', '--Verbosity', type=int,
#                         help="Verbosity: 0 = none, 1 = info", default=defaults["Verbosity"])
#     Core.add_argument('-dr', '--DryRun', type=str,
#                         help="Indicates whether programs are executed or not", default=defaults["DryRun"])

#     Compute = parser.add_argument_group("Compute Arguments")
#     ### Compute Arguments
#     Compute.add_argument('-c', '--CoresPerJob', type=int,
#                         help="Number of cores per individual calculation", default=defaults["CoresPerJob"])
#     Compute.add_argument('-m','--MemoryPerJob', type=int,
#                         help="Gb of memory per individual calculation", default=defaults["MemoryPerJob"])

#     ### MM Arguments
#     MM = parser.add_argument_group("Molecular Dynamics Arguments")
#     MM.add_argument('-MDcpu', '--MDCPUPath', type=str,
#                         help="Path to NAMD CPU executable", default=NAMD_CPU)
#     MM.add_argument('-MDgpu', '--MDGPUPath', type=str,
#                         help="Path to NAMD GPU executable", default=NAMD_GPU)

#     ### QM Arguments
#     QM = parser.add_argument_group("QM Arguments")
#     QM.add_argument("-qf", "--QmFile", type=str, default=defaults["QmFile"],
#                     help="Name of file containing QM information.")
#     QM.add_argument('-qp', '--QmPath', type=str,
#                         help="Path to QM software", default=ORCA_PATH)
#     QM.add_argument('-qsel', '--QmSelection', type=str,
#                         help="Selection algebra for QM atoms", default=defaults["QmSelection"])
#     QM.add_argument('-qc', '--QmCharge', type=int,
#                         help="Charge of QM region", default=defaults["QmCharge"])
#     QM.add_argument('-qspin', '--QmSpin', type=int,
#                         help="Spin of QM region", default=defaults["QmSpin"])
#     QM.add_argument('-qm', '--QmMethod', type=str,
#                         help="Qm method", default=defaults["QmMethod"])
#     QM.add_argument('-qb', '--QmBasis', type=str,
#                         help="QM basis set", default=defaults["QmBasis"])
#     QM.add_argument('-qargs', '--QmArgs', type=str, 
#                     help="Extra arguments for ORCA calculation", default=defaults["QmArgs"])

#     ### HPC Arguments
#     HPC = parser.add_argument_group("HPC/SLURM arguments")
#     HPC.add_argument("--HPC", type=bool, default=defaults["HPC"], help="Whether to run on a HPC")
#     HPC.add_argument("-MaxTime", "--MaxWallTime", type=int,
#                      help="Maximum wall time (Hours) for your jobs (either leave as node max, or set as job length)",
#                      default=HPC_Conf["MaxWallTime"])
#     HPC.add_argument("-Host", "--HostName", type=str,
#                      help="HostName of the HPC", default=HPC_Conf["HostName"])
#     HPC.add_argument("--Partition", type=str, help="Calculation partition name",
#                      default=HPC_Conf["Partition"])
#     HPC.add_argument("--MaxCores", type=int,
#                      help="Maximum number of cores available to a node (For array splitting)", default=HPC_Conf["MaxCores"])
#     HPC.add_argument("-QoS", "--QualityofService", type=str,
#                      help="Slurm QoS, set to None if not relevant.", default=HPC_Conf["QualityofService"])
#     HPC.add_argument("--Account", type=str,
#                      help="Slurm account, (Not username), Set to None if not relevant", default=HPC_Conf["Account"])
#     HPC.add_argument("-Software", "--SoftwareLines", type=str,
#                     help="List of commands like \"module load XXX\" to load software. Keep each line surrounded by quotes.",
#                     default=HPC_Conf["SoftwareLines"], nargs="*")

#     Standalone = parser.add_argument_group("Standalone Job arguments")
#     Standalone.add_argument("-i", "--MMFile", type=str, default=defaults["MMFile"],
#                             help="Name of input file containing MD information")
#     Standalone.add_argument("--Name", type=str,
#                             default=defaults["Name"], help="Name for the calculation")
#     Standalone.add_argument("--Ensemble", type=str,
#                             choices=["min", "heat", "NVT", "NPT"],
#                             help="Ensemble for Calculation", default=defaults["Ensemble"])
#     Standalone.add_argument("--QM", type=str, choices=["True", "False"],
#                             default=defaults["QM"], help="Whether this is a QMMM calculation or not.")
#     Standalone.add_argument("-st", "--Steps", type=int,
#                             default=defaults["Steps"], help="Number of simulation steps.")
#     Standalone.add_argument("-dt", "--TimeStep", type=float,
#                             default=float(defaults["TimeStep"]), help="Time step for the simulation. We recommend 2 for MM, 0.5 for QMMM")
#     Standalone.add_argument("--ParmFile", type=str,
#                             default=defaults["ParmFile"], help="Parameter file name")
#     Standalone.add_argument("--AmberCoordinates", type=str,
#                             default=defaults["AmberCoordinates"], help="Amber coordinate file name that relates to the parameter file")
#     Standalone.add_argument("--StartFile", type=str, default=defaults["StartFile"], 
#                             help="Either Amber coordinates or NAMD coordinates. These are the coordinates that it starts from.")
#     Standalone.add_argument("--RestartOut", type=int, default=defaults["RestartOut"], 
#                             help="Frequency to generate a restart file")
#     Standalone.add_argument("--TrajOut", type=int, default=defaults["TrajOut"], 
#                             help="Frequency to add to the trajectory file")
#     Standalone.add_argument("--SMD", type=bool, choices=[True, False], default=defaults["SMD"], 
#                             help="Wheter to use steered molecular dynamics")
#     Standalone.add_argument("--Force", type=float, default=defaults["Force"], 
#                             help="Force for Steered MD")
#     Standalone.add_argument("--StartValue", type=float, default=defaults["StartValue"], 
#                             help="Start value for SMD")
#     Standalone.add_argument("--EndValue", type=float, default=defaults["EndValue"], 
#                             help="End value for SMD. MAKE == Start if wanting constant.")
#     Standalone.add_argument('-mask', '--AtomMask', type=str,
#                         help="Mask for the restrained atoms.", default=defaults["AtomMask"])
#     return parser
    
# def VariableParser(sysargs, JT="Umbrella"):
#     if JT == "Umbrella":
#         JobDict = JobInput("./Job.conf")
#         WorkDir = JobDict["WorkDir"]
#         ComputeDict = ComputeInput(f"{WorkDir}Compute.conf")
#         MMDict = MMInput(f"{WorkDir}MM.conf")
#         QMDict = QMInput(f"{WorkDir}QM.conf")
#         UmbrellaDict = UmbrellaInput(f"{WorkDir}Umbrella.conf")
#         HPCDict = HPCInput(f"{WorkDir}HPC.conf")
#         StandaloneDict = StandaloneJobInput(f"{WorkDir}Standalone.conf")
#         if JobDict["JobType"].casefold() == "umbrella":
#             FileDict = {**JobDict, **ComputeDict, **MMDict, **QMDict,
#                         **HPCDict, **StandaloneDict,**UmbrellaDict}
#         elif JobDict["JobType"].casefold() == "Inpfile":
#             FileDict = {**JobDict, **ComputeDict, **MMDict, **QMDict, **UmbrellaDict, **HPCDict, **StandaloneDict}
#         elif JobDict["JobType"].casefold() == "mm" or JobDict["JobType"].casefold() == "qmmm":
#             FileDict = {**JobDict, **ComputeDict, **MMDict, **QMDict, **UmbrellaDict, **HPCDict,**StandaloneDict}
#         else:
#             FileDict = {**JobDict, **ComputeDict, **MMDict, **QMDict, **UmbrellaDict, **HPCDict,**StandaloneDict}
#         args = arg_parse_Umbrella(FileDict,sysargs)
#         if args.JobType.casefold() == "inpfile":
#             InputFileGen(args, JobType="Umbrella")
#         if args.Verbosity >= 2:
#             print(vars(args))
#         return args
#     elif JT == "Standalone":
#         JobDict = JobInput("./Job.conf")
#         WorkDir = JobDict["WorkDir"]
#         ComputeDict = ComputeInput(f"{WorkDir}Compute.conf")
#         MMDict = MMInput(f"{WorkDir}MM.conf")
#         QMDict = QMInput(f"{WorkDir}QM.conf")
#         HPCDict = HPCInput(f"{WorkDir}HPC.conf")
#         StandaloneDict = StandaloneJobInput(f"{WorkDir}Standalone.conf")
#         FileDict = {**JobDict, **ComputeDict, **MMDict, **QMDict, **HPCDict, **StandaloneDict}
#         args = arg_parse_Standalone(FileDict,sysargs)
#         if args.JobType.casefold() == "inpfile":
#             InputFileGen(args, JobType="Standalone")
#         if args.Verbosity >= 2:
#             print(vars(args))
#         return args
        
# def InputFileGen(args, JobType="Umbrella"):
#     argsDict = vars(args)
#     JobDict = JobInput(f"{args.WorkDir}Job.conf")
#     with open(f"{args.WorkDir}Job.conf", "w") as f:
#         for i in JobDict.keys():
#             print(f"{i}={argsDict[i]}", file=f)
#     ComputeDict = ComputeInput(f"{args.WorkDir}Compute.conf")
#     with open(f"{args.WorkDir}Compute.conf", "w") as f:
#         for i in ComputeDict.keys():
#             print(f"{i}={argsDict[i]}", file=f)
#     MMDict = MMInput(f"{args.WorkDir}MM.conf")
#     with open(f"{args.WorkDir}MM.conf", "w") as f:
#         for i in MMDict.keys():
#             print(f"{i}={argsDict[i]}", file=f)
#     QMDict = QMInput(f"{args.WorkDir}QM.conf")
#     with open(f"{args.WorkDir}QM.conf", "w") as f:
#         for i in QMDict.keys():
#             print(f"{i}={argsDict[i]}", file=f)
#     if JobType == "Umbrella":
#         UmbrellaDict = UmbrellaInput(f"{args.WorkDir}Umbrella.conf")
#         with open(f"{args.WorkDir}Umbrella.conf", "w") as f:
#             for i in UmbrellaDict.keys():
#                 print(f"{i}={argsDict[i]}", file=f)
#     HPCDict = HPCInput(f"{args.WorkDir}HPC.conf")
#     with open(f"{args.WorkDir}HPC.conf", "w") as f:
#         for i in HPCDict.keys():
#             print(f"{i}={argsDict[i]}", file=f)
#     if JobType == "Standalone":
#         StandaloneDict = StandaloneJobInput(f"{args.WorkDir}Standalone.conf")
#         with open(f"{args.WorkDir}Standalone.conf", "w") as f:
#             for i in StandaloneDict.keys():
#                 print(f"{i}={argsDict[i]}", file=f)

# def arg_parse_Umbrella(dict, sysargs):
#     parser = ap.ArgumentParser(description=f"""Commandline arguments. This method of calculation input is being deprecated. Please do not use.
# It is recommended to use -jt inpfile to generate input file templates with default values that you can then edit.""")
#     ### Core Job arguments
#     Core = parser.add_argument_group("Core Job Arguments")
#     Core.add_argument('-wd', '--WorkDir', type=str,
#                         help="Home location for the calculations", default=dict["WorkDir"])
#     Core.add_argument('-jt', '--JobType', type=str,
#                         help="Type of calculation to run", default=dict["JobType"])
#     Core.add_argument('-v', '--Verbosity', type=int,
#                         help="Verbosity: 0 = none, 1 = info", default=dict["Verbosity"])
#     Core.add_argument('-dr', '--DryRun', type=str,
#                         help="Indicates whether programs are executed or not", default=dict["DryRun"])

#     Compute = parser.add_argument_group("Compute Arguments")
#     ### Compute Arguments
#     Compute.add_argument('-cores', '--CoresPerJob', type=int,
#                         help="Number of cores per individual calculation", default=dict["CoresPerJob"])
#     Compute.add_argument('-mem','--MemoryPerJob', type=int,
#                         help="Gb of memory per individual calculation", default=dict["MemoryPerJob"])
#     Compute.add_argument('-MaxCalc', '--MaxStepsPerCalc', type=int,
#                          help="The maximum number of steps per calculation. splits jobs into sub-steps. useful for short wall times. 0 == No cap.",
#                          default=dict["MaxStepsPerCalc"])

#     ### MM Arguments
#     MM = parser.add_argument_group("Molecular Dynamics Arguments")
#     MM.add_argument('-MDcpu', '--MDCPUPath', type=str,
#                         help="Path to NAMD CPU executable", default=dict["MDCPUPath"])
#     MM.add_argument('-MDgpu', '--MDGPUPath', type=str,
#                         help="Path to NAMD GPU executable", default=dict["MDGPUPath"])

#     ### QM Arguments
#     QM = parser.add_argument_group("QM Arguments")
#     QM.add_argument('-qp', '--QmPath', type=str,
#                         help="Path to QM software", default=dict["QmPath"])
#     QM.add_argument('-qsel', '--QmSelection', type=str,
#                         help="Selection algebra for QM atoms", default=dict["QmSelection"])
#     QM.add_argument('-qc', '--QmCharge', type=int,
#                         help="Charge of QM region", default=dict["QmCharge"])
#     QM.add_argument('-qspin', '--QmSpin', type=int,
#                         help="Spin of QM region", default=dict["QmSpin"])
#     QM.add_argument('-qm', '--QmMethod', type=str,
#                         help="Qm method", default=dict["QmMethod"])
#     QM.add_argument('-qb', '--QmBasis', type=str,
#                         help="QM basis set", default=dict["QmBasis"])
#     QM.add_argument('-qargs', '--QmArgs', type=str, help="Extra arguments for ORCA calculation", default=dict["QmArgs"])

#     ### Umbrella Arguments
#     Umbrella = parser.add_argument_group("Umbrella Sampling arguments")
#     Umbrella.add_argument('-min', '--UmbrellaMin', type=float,
#                         help="Minimum Umbrella distance", default=dict["UmbrellaMin"])
#     Umbrella.add_argument('-width', '--UmbrellaWidth', type=float,
#                         help="Umbrella bin width in Angstroms or degrees", default=dict["UmbrellaWidth"])
#     Umbrella.add_argument('-bins', '--UmbrellaBins', type=int,
#                         help="Number of umbrella bins", default=dict["UmbrellaBins"])
#     Umbrella.add_argument('-pf', '--PullForce', type=float,
#                         help="Force for pulls in KCal A-2", default=dict["PullForce"])
#     Umbrella.add_argument('-f', '--ConstForce', type=float,
#                         help="Force for standard Umbrella runs", default=dict["ConstForce"]) ### NAMD uses 1/2 k rather than just k
#     Umbrella.add_argument('-sd', '--StartDistance', type=float,
#                         help="Distance of initial simulation", default=dict["StartDistance"])
#     Umbrella.add_argument('-mask', '--AtomMask', type=str,
#                         help="Mask for the restrained atoms.", default=dict["AtomMask"])
#     Umbrella.add_argument('-stg', '--Stage', type=str,
#                         help="Stage of ummbrella simulation", default=dict["Stage"])
#     Umbrella.add_argument('-wf', '--WhamFile', type=str,
#                         help="Name prefix of wham data.(XXX.i.colvars.traj", default=dict["WhamFile"])
#     Umbrella.add_argument("--StartFile", default=dict["StartFile"], type=str, help="Initial coordinate file if not starting from \"start.rst7\"")
#     Umbrella.add_argument("-exclude", "--WhamExclude", type=str, help="Comma delimited list of umbrella windows to exclude from Wham calculations", 
#                           default=dict["WhamExclude"])

#     ### HPC Arguments
#     HPC = parser.add_argument_group("HPC/SLURM arguments")
#     HPC.add_argument("-MaxTime", "--MaxWallTime", type=int,
#                      help="Maximum wall time (Hours) for your jobs (either leave as node max, or set as job length)",
#                      default=dict["MaxWallTime"])
#     HPC.add_argument("-Host", "--HostName", type=str,
#                      help="HostName of the HPC", default=dict["HostName"])
#     HPC.add_argument("--Partition", type=str, help="Calculation partition name",
#                      default=dict["Partition"])
#     HPC.add_argument("--MaxCores", type=int,
#                      help="Maximum number of cores available to a node (For array splitting)", default=dict["MaxCores"])
#     HPC.add_argument("-QoS", "--QualityofService", type=str,
#                      help="Slurm QoS, set to None if not relevant.", default=dict["QualityofService"])
#     HPC.add_argument("--Account", type=str,
#                      help="Slurm account, (Not username), Set to None if not relevant", default=dict["Account"])
#     HPC.add_argument("-Software", "--SoftwareLines", type=str,
#                     help="List of commands like \"module load XXX\" to load software. Keep each line surrounded by quotes.",
#                     default=dict["SoftwareLines"], nargs="*")

#     # Standalone = parser.add_argument_group("Standalone Job arguments")
#     # Standalone.add_argument("--Name", type=str,
#     #                         default=dict["Name"], help="Name for the calculation")
#     # Standalone.add_argument("--Ensemble", type=str,
#     #                         choices=["min", "heat", "NVT", "NPT"],
#     #                         help="Ensemble for Calculation", default=dict["Ensemble"])
#     # Standalone.add_argument("--QM", type=str, choices=["True", "False"],
#     #                         default=dict["QM"], help="Whether this is a QMMM calculation or not.")
#     # Standalone.add_argument("-st", "--Steps", type=int,
#     #                         default=dict["Steps"], help="Number of simulation steps.")
#     # Standalone.add_argument("-dt", "--TimeStep", type=float,
#     #                         default=dict["TimeStep"], help="Time step for the simulation. We recommend 2 for MM, 0.5 for QMMM")
#     # Standalone.add_argument("--ParmFile", type=str,
#     #                         default=dict["ParmFile"], help="Parameter file name")
#     # Standalone.add_argument("--AmberCoordinates", type=str,
#     #                         default=dict["AmberCoordinates"], help="Amber coordinate file name that relates to the parameter file")
#     # # Standalone.add_argument("--StartFile", type=str, default=dict["StartFile"], help="Either Amber coordinates or NAMD coordinates. These are the coordinates that it starts from.")
#     # Standalone.add_argument("--RestartOut", type=int, default=dict["RestartOut"], help="Frequency to generate a restart file")
#     # Standalone.add_argument("--TrajOut", type=int, default=dict["TrajOut"], help="Frequency to add to the trajectory file")
#     # Standalone.add_argument("--SMD", type=str, choices=["off", "on"], default=dict["SMD"], help="Wheter to use steered molecular dynamics")
#     # Standalone.add_argument("--Force", type=float, default=dict["Force"], help="Force for Steered MD")
#     # Standalone.add_argument("--StartValue", type=float, default=dict["StartValue"], help="Start value for SMD")
#     # Standalone.add_argument("--EndValue", type=float, default=dict["EndValue"], help="End value for SMD. MAKE == Start if wanting constant.")

#     # Standalone.add_argument("", type=, default=dict[""], help="")

#     ### Parse commandline arguments
#     args = parser.parse_args(sysargs)
#     # print(vars(args))
#     return args

# def JobInput(path):
#     InpVars = ["WorkDir", "JobType", "Verbosity", "DryRun"]
#     InpValues = ["./", "inpfile", 0, "True"]
#     assert len(InpVars) == len(InpValues)
#     try:
#         lines = utils.file_read(path)
#     except FileNotFoundError:
#         print("WARNING, No config found for Job input, This is a bad idea... Using defaults.")
#         Dict = {}
#         for i in range(len(InpVars)):
#             Dict[InpVars[i]] = str(InpValues[i])
#         return Dict
#     for line in lines:
#         words = line.split("=")
#         for i in range(len(InpVars)):
#             if words[0].casefold() == InpVars[i].casefold():
#                 InpValues[i] = words[1].replace("\n","")
#     Dict = {}
#     for i in range(len(InpVars)):
#         Dict[InpVars[i]] = InpValues[i]
#     return Dict

# def ComputeInput(path):
#     InpVars = ["CoresPerJob", "MemoryPerJob", "MaxStepsPerCalc"]
#     InpValues = [10, 10, 1000]
#     assert len(InpVars) == len(InpValues)
#     try:
#         lines = utils.file_read(path)
#     except FileNotFoundError:
#         print("WARNING, No config found for Compute input, Using defaults.")
#         Dict = {}
#         for i in range(len(InpVars)):
#             Dict[InpVars[i]] = str(InpValues[i])
#         return Dict
#     for line in lines:
#         words = line.split("=")
#         for i in range(len(InpVars)):
#             if words[0].casefold() == InpVars[i].casefold():
#                 InpValues[i] = words[1].replace("\n","")
#     Dict = {}
#     for i in range(len(InpVars)):
#         Dict[InpVars[i]] = InpValues[i]
#     return Dict

# def MMInput(path):
#     InpVars = ["MDCPUPath", "MDGPUPath"]
#     InpValues = ["/work/e280/e280-Hirst/pcyra2/Software/NAMD_3.0b3_Linux-x86_64-multicore/namd3", "/home/pcyra2/Software/NAMD/NAMD_3.0b4_Linux-x86_64-multicore-CUDA/namd3"]
#     assert len(InpVars) == len(InpValues)
#     try:
#         lines = utils.file_read(path)
#     except FileNotFoundError:
#         print("WARNING, No config found for MM input, Using defaults.")
#         Dict = {}
#         for i in range(len(InpVars)):
#             Dict[InpVars[i]] = str(InpValues[i])
#         return Dict
#     for line in lines:
#         words = line.split("=")
#         for i in range(len(InpVars)):
#             if words[0].casefold() == InpVars[i].casefold():
#                 InpValues[i] = words[1].replace("\n","")
#     Dict = {}
#     for i in range(len(InpVars)):
#         Dict[InpVars[i]] = InpValues[i]
#     return Dict

# def QMInput(path):
#     InpVars = ["QmPath", "QmSelection", "QmCharge", "QmSpin", "QmMethod", "QmBasis", "QmArgs"]
#     InpValues = ["/work/y07/shared/apps/core/orca/5.0.3/orca", "resname CTN POP MG", 3, 1, "PBE", "6-31G*", "D3BJ TightSCF CFLOAT "]
#     assert len(InpVars) == len(InpValues)
#     try:
#         lines = utils.file_read(path)
#     except FileNotFoundError:
#         print("WARNING, No config found for QM input, Using defaults.")
#         Dict = {}
#         for i in range(len(InpVars)):
#             Dict[InpVars[i]] = str(InpValues[i])
#         return Dict
#     for line in lines:
#         words = line.split("=")
#         for i in range(len(InpVars)):
#             if words[0].casefold() == InpVars[i].casefold():
#                 InpValues[i] = words[1].replace("\n","")
#     Dict = {}
#     for i in range(len(InpVars)):
#         Dict[InpVars[i]] = InpValues[i]
#     return Dict

# def UmbrellaInput(path):
#     InpVars = ["UmbrellaMin", "UmbrellaWidth", "UmbrellaBins", "PullForce", "ConstForce", "StartDistance", "AtomMask", "Stage", "WhamFile", "StartFile", "WhamExclude"]
#     InpValues = [1.3, 0.05, 54, 5000, 300, 1.4, "0,0,0,0", "Setup", "prod", "start.rst7",""]
#     assert len(InpVars) == len(InpValues)
#     try:
#         lines = utils.file_read(path)
#     except FileNotFoundError:
#         print("WARNING, No config found for Umbrella input, Using defaults.")
#         Dict = {}
#         for i in range(len(InpVars)):
#             Dict[InpVars[i]] = str(InpValues[i])
#         return Dict
#     for line in lines:
#         words = line.split("=")
#         for i in range(len(InpVars)):
#             if words[0].casefold() == InpVars[i].casefold():
#                 InpValues[i] = words[1].replace("\n","")
#     Dict = {}
#     for i in range(len(InpVars)):
#         Dict[InpVars[i]] = InpValues[i]
#     return Dict

# def HPCInput(path):
#     InpVars = ["MaxWallTime", "HostName","Partition", "MaxCores", "QualityofService", "Account", "SoftwareLines"]
#     InpValues = [24, "login.archer2.ac.uk", "standard", 128, None, None, "module load ORCA", ]
#     InpVars2 = []
#     InpValues2 = []
#     assert len(InpVars) == len(InpValues)
#     try:
#         lines = utils.file_read(path)
#     except FileNotFoundError:
#         print("WARNING, No config found for HPC input, Using defaults.")
#         Dict = {}
#         for i in range(len(InpVars)):
#             Dict[InpVars[i]] = str(InpValues[i])
#         return Dict
#     for i in range(len(lines)):
#         words = lines[i].split("=",1)
#         for j in range(len(InpVars)):
#             if words[0].casefold() == InpVars[j].casefold():
#                 InpValues2.append(words[1].replace("\n",""))
#                 InpVars2.append(words[0])
#     Dict = {}
#     for i in range(len(InpVars2)):
#         if Dict.get(InpVars2[i]) == None:
#             Dict[InpVars2[i]] = InpValues2[i]
#             # print(Dict.get(InpVars2[i]))
#         else:       ### This allows for multiple iof the same word keyword (i.e. Multiple software lines)
#             Vals = Dict.get(InpVars2[i])
#             if type(Vals) == str:
#                 Dict[InpVars2[i]] = [Vals, InpValues2[i]]
#             else:
#                 Vals.append(InpValues2[i])
#                 Dict[InpVars2[i]] = Vals
#             # print(Vals)
#     return Dict

# def StandaloneJobInput(path):
#     InpVars = ["Name", "ParmFile", "AmberCoordinates", "StartFile", "Ensemble",
#                "QM", "Steps", "TimeStep", "RestartOut", "TrajOut", "SMD", "Force", "StartValue", "EndValue", "AtomMask"]
#     InpValues = ["QMMM_Job","complex.parm7", "start.rst7", "Start.rst7" ,"min",
#                  "true", 1000, 0.05, 10, 50, "off", 1, 1, 2, "0,0,0,0"]
#     assert len(InpVars) == len(InpValues)
#     try:
#         lines = utils.file_read(path)
#     except FileNotFoundError:
#         print("WARNING, No config found for Standalone Job input, Using defaults.")
#         Dict = {}
#         for i in range(len(InpVars)):
#             Dict[InpVars[i]] = str(InpValues[i])
#         return Dict
#     for line in lines:
#         words = line.split("=")
#         for i in range(len(InpVars)):
#             if words[0].casefold() == InpVars[i].casefold():
#                 InpValues[i] = words[1].replace("\n","")
#     Dict = {}
#     for i in range(len(InpVars)):
#         Dict[InpVars[i]] = InpValues[i]
#     return Dict

