import pyBrellaSampling.Tools.InputParser as input
import pyBrellaSampling.Tools.utils as utils
import pyBrellaSampling.Tools.FileGen as FileGen
from pyBrellaSampling.Tools.classes import *
import pyBrellaSampling.Tools.globals as globals
import subprocess
import sys
import time
import os.path as path

def main_cli():
    starttime = time.time()
    args = input.StandaloneInput(sys.argv[1:],)
    verbosity = args["Verbosity"]
    WorkDir = args["WorkDir"]
    dr = args["DryRun"]
    if dr == "True" or dr == True:
        DryRun = True
    elif dr == "False" or dr == False:
        DryRun = False
    parmfile = args["ParmFile"]
    assert path.isfile(f"{WorkDir}{parmfile}"), f"ERROR: {parmfile} does not exist in the working directory!"
    globals.init(v=verbosity,wd=WorkDir,dr=DryRun, parm=parmfile)
    Calc, MM, QM = calc_setup(args)
    if globals.DryRun == False:
        calc_run(Calc=Calc, MM=MM, QM=QM, )
    endtime = time.time()
    print(f"Total time is {endtime - starttime}")

        
def Class_init(args: dict):
    """Initialises some of the main classes used by the calculation. 

    Args:
        args (dict): User defined variables.

    Returns:
        Calc (CalcClass): Generic calculation info
        MM (MMClass): MM Variables
        QM (QMClas): QM variables, if == False : No QMMM
        Umbrella (UmbrellaClass): Umbrella variables for restraint information.
    """
    Calc = CalcClass(args=args)
    MM = MMClass()
    Calc.Job_Name(Name=args["Name"])
    Calc.Set_OutFile(OutFile=args["Name"])
    Calc.Set_Id(0)
    MM.Set_Length(Steps=args["Steps"], TimeStep=args["TimeStep"])
    MM.Set_Ensemble(Ensemble=args["Ensemble"])
    MM.Set_Files(parm=args["ParmFile"], ambercoor=args["AmberCoordinates"])
    MM.Set_Outputs(TimeOut=args["TrajOut"], RestOut=args["RestartOut"], 
                   TrajOut=args["TrajOut"]) ### Output timings when trajectory is printed.
    args["PullForce"] = args["Force"]
    args["ConstForce"] = args["Force"]
    Umbrella = UmbrellaClass(args, Min=float(args["StartValue"]), bins=0,
                            Start=float(args["StartValue"]),
                            Width=(float(args["EndValue"]) - float(args["StartValue"])))
    Umbrella.add_start(0)
    Umbrella.set_force(args["Force"])
    Umbrella.add_bins([str(args["EndValue"])])
    bincoor = args["StartFile"]
    if bincoor == "None" or bincoor == "" or bincoor == MM.ambercoor:
                bincoor = None
    if ".ncrst" in str(bincoor) or ".rst7" in str(bincoor):
        print("WARNING: Your Binary Coordinates look like amber coordinates... They should be NAMD coodinates (.coor). This may cause issues...\n" if globals.verbosity >= 1 else "", end="")
    if args["QM"] == True or args["QM"] == "True":
        print("INFO: Running a QMMM calculation \n" if globals.verbosity >= 2 else "", end="")
        QM = QMClass(args=args)
        QM.set_selfile("syst-qm.pdb")
    else:
        QM = False
    if MM.TimeStep > 1:
        print("WARNING: TimeStep is greater than 1 fs. Setting Rattle to True\n" if globals.verbosity >=1 else "", end="")
        MM.Set_Shake("all")
    NAMD = NAMDClass(Calc=Calc, MM=MM)
    return Calc, MM, QM, Umbrella, NAMD

def calc_setup(args: dict):
    """
    Generates the input file, including any syst-*.pdb files.
    
    Args:
        args (dict): User defined variables

    Returns:
        Calc (CalcClass): Calculation variables
        MM (MMClass): Required for knowing how to run the calculation (With or without GPU ect. )
        QM (QMClass): Again required for knowing how to run the calculation.
    """
    Calc, MM, QM, Umbrella, NAMD = Class_init(args) # inits classes
    MM.CellVec = utils.get_cellVec(MM) # Obtains the cell vectors from the parm file. 
    NAMD.set_cellvectors(MM.CellVec)
    NAMD.set_startcoords(args["StartFile"], ambercoor=MM.ambercoor, parm=MM.parmfile)
    if QM != False:
        print("INFO: Setting up a QMMM calculation\n" if globals.verbosity >= 2 else "", end="")
        NAMD.set_qm(Calc=Calc, QM=QM, index="QMMM")
        utils.QM_Gen(QM.QMSel, globals.WorkDir)
        if globals.DryRun == False:
            print("INFO: Setting up the QM pdb file.\n" if globals.verbosity >=2 else "", end="")
            logfile = subprocess.run(["vmd", "-dispdev", "text", "-e", "qm_prep.tcl"],
                                     text = True, capture_output = True) # Generates sys-qm.pdb using vmd
            with open(f"{globals.WorkDir}tcl-qm.log","w") as f:
                print(logfile, file=f)
    else:
        NAMD.set_pme("on")
        # pass
    if args["SMD"].casefold() == "true":
        print(f"WARNING: SMD method is currently untested... " if globals.verbosity >=1 else "", end="")
        NAMD = init_SMD(NAMD=NAMD,Umbrella=Umbrella )
        utils.ColVarPDB_Gen(Umbrella)
        if globals.DryRun == False:
            print("INFO: Setting up the Colvar pdb file.\n" if globals.verbosity >=2 else "", end="")
            logfile = subprocess.run(["vmd", "-dispdev", "text", "-e", "Colvar_prep.tcl"],
                                     text = True, capture_output = True) # Generates sys-col.pdb using vmd. 
            with open(f"{globals.WorkDir}tcl-colvar.log","w") as f:
                print(logfile, file=f)
    print("INFO: Setting up the conf file\n" if globals.verbosity >=2 else "", end="")
    file = FileGen.Namd_File(NAMD)
    utils.file_write(f"{globals.WorkDir}{Calc.Name}.conf", [file]) # Outputs the calculation file for namd. 
    return Calc, MM, QM

def calc_run(Calc: CalcClass, MM: MMClass, QM: QMClass):
    """
    Handles the running of the calculations. Choosing the right version of NAMD
    Args:
        Calc (CalcClass): Gets the name of the calculation
        MM (MMClass): Gets the path to executables
        QM (QMClass): Gets wheter to use the GPU.

    """
    print(f"INFO: Running the {Calc.Name} Calculation.\n" if globals.verbosity >=2 else "", end="")
    if  QM == False:
        print(f"INFO: {MM.GPUNamd} +oneWthPerCore +setcpuaffinity +devices 0 +cs {globals.WorkDir}{Calc.Name}.conf > {globals.WorkDir}{Calc.Name}_1.0.out\n" 
              if globals.verbosity >=2 else "", end="")
        subprocess.run([MM.GPUNamd +
                        f" +oneWthPerCore +setcpuaffinity +devices 0 +cs {globals.WorkDir}{Calc.Name}.conf > {globals.WorkDir}{Calc.Name}_1.0.out"],
                       shell=True, capture_output=True)
    else:
        subprocess.run(["mkdir /dev/shm/NAMD_QMMM"], shell=True)
        subprocess.run([MM.CPUNamd +
                        f" +p1 {globals.WorkDir}{Calc.Name}.conf > {globals.WorkDir}{Calc.Name}_1.0.out"],
                       shell=True, capture_output=True)
        subprocess.run(["rm -r /dev/shm/NAMD_QMMM"], shell=True)

def init_SMD(NAMD: NAMDClass, Umbrella: UmbrellaClass):
    """
    Initialises Steered MD, Warning this isnt tested properly yet... 
    Args:
        NAMD (NAMDClass): Carries colvar variables to input file
        Umbrella (UmbrellaClass): Contains original colvars

    Returns: 
        NAMD (NAMDClass): Updated NAMD class with colvars. 

    """
    if Umbrella.Width != 0:
        Type = "pull"
    else:
        Type = "constant"
    NAMD.set_colvars(file="colvars.conf", toggle="on")
    colvarfile = utils.colvar_gen(Umbrella, i=0,type=Type, force=Umbrella.ConstForce, relLoc="./")
    utils.file_write(path=f"{globals.WorkDir}colvars.conf", lines=[colvarfile])
    return NAMD

# def InputFileRun(Path,args,Run=True):
#     print(utils.MM_DefaultVars.keys())
#     CalcStarts = []
#     CalcEnds = []
#     data = utils.file_read(Path)
#     if CalcStarts[0] > 0:
#         for i in range(CalcStarts[0]):
#             words = data[i].split()
#             if "threads" in words[0].casefold():
#                 args.threads = int(words[1])
#     for i in range(len(data)):
#         if "$JOB" in data[i]:
#             CalcStarts.append(i)
#         if "$END" in data[i]:
#             CalcEnds.append(i)
#     assert len(CalcStarts) == len(CalcEnds), "Job input file is not terminated... "
#     for i in range(len(CalcStarts)):
#         MMVars = utils.MM_DefaultVars
#         for j in range(CalcStarts[i], CalcEnds[i]):
#             words = data[j].split()
#             if "name" in words[0].casefold():
#                 JobName = words[1]
#             for keys in MMVars.keys():
#                 if words[0].casefold() == keys.casefold():
#                     MMVars[keys] = words[1]
#             Calc = CalcClass(args)
#             Calc.Set_Shake(MMVars["rigidBonds"])

#     return FileNotFoundError
