import sys
import time
import os
import os.path as path
import subprocess
import math

import pyBrellaSampling.Tools.InputParser as input
from pyBrellaSampling.Tools.classes import *
import pyBrellaSampling.Tools.globals as globals # Global Variables
import pyBrellaSampling.Tools.utils as utils
import pyBrellaSampling.Tools.FileGen as FileGen
import pyBrellaSampling.Umbrella.wham as Wham
import pyBrellaSampling.Tools.analysis as Anal

def main():
    """Main Umbrella sampling executable. This deals with everything and can be run from the command line. 
    """
    starttime = time.time()
    args = input.UmbrellaInput(sys.argv[1:],) 
    verbosity = args["Verbosity"]
    assert type(verbosity) is int, f"ERROR: verbosity statement must be an integer, not {verbosity}"
    WorkDir = args["WorkDir"]
    assert path.isdir(f"{WorkDir}"), f"ERROR: User defined working directory: {WorkDir} Does not exist!"
    dr = args["DryRun"]   
    if dr == "True" or dr == True:
        DryRun = True
    elif dr == "False" or dr == False: # pragma: no cover
        DryRun = False
    assert type(DryRun) is bool, f"ERROR: DryRun must be a boolean, not {DryRun}!"
    parmfile = args["ParmFile"]
    assert path.isfile(f"{WorkDir}{parmfile}"), f"ERROR: {parmfile} does not exist in the working directory!"
    globals.init(v=verbosity, wd=WorkDir,dr=DryRun, parm=parmfile)
    Calc, MM, QM, Umbrella, SLURM = class_init(args)
    equil_length = int(args["EquilLength"])*2000 # Number of steps to equilibrate each window. converts from ps to steps in 0.5 fs timsteps
    prod_length = int(args["ProdLength"])*2000 # Number of steps to run production for each window.
    assert type(equil_length) is int, f"ERROR: equil length must be an integer. Cant do half steps... {equil_length}"
    assert type(prod_length) is int, f"ERROR: prod length must be an integer, cant do half steps {prod_length}"
    AnalIgnore = [] # Bins to ignore for analysis... Use with caution.
    if args["Stage"].casefold() == "setup":
        setup(QM, Umbrella)
    if args["Stage"].casefold() == "min":
        min(MM, Calc, args["StartFile"])
    if args["Stage"].casefold() == "heat":
        heat(MM, Calc)
    if args["Stage"].casefold() == "pull":
        pull(Umbrella, MM, QM, Calc)
    if args["Stage"].casefold() == "equil":
        MM.Set_Ensemble("NVT")
        Calc.Job_Name("equil")
        MM.Set_Length(equil_length, 0.5)  # 2000 steps at 0.5 fs = 1 ps ~ 1 day
        MM.Set_Outputs(100, 100, 80)
        const(MM, QM, Calc, Umbrella, "pull_1")
        SLURM.set_arrayJob("equil_1.txt", Umbrella.Bins)
        utils.slurm_gen("NAME", SLURM, "sh array_job.sh", globals.WorkDir)
        utils.batch_sub(5,20, )
    if args["Stage"].casefold() == "prod":
        Calc.Job_Name("prod")
        MM.Set_Ensemble("NVT")
        MM.Set_Length(prod_length, 0.5)  # 8000 steps at 0.5 fs = 4 ps, ~ 3.5 days
        MM.Set_Outputs(100, 100, 80)  # Timings, Restart, Trajectory
        if Calc.MaxSteps == 0:
            num = 1
        else:
            num = math.ceil(equil_length/Calc.MaxSteps)
        const(MM, QM, Calc, Umbrella, f"equil_{num}")
    if args["Stage"].casefold() == "wham":
        if Calc.MaxSteps == 0:
            NumJobs = 1
        else:
            if "prod" in args["AnalysisFile"].casefold():
                steps = prod_length
            elif "equil" in args["AnalysisFile"].casefold():
                steps = equil_length
            else:
                steps = 0
            NumJobs = math.ceil(steps / Calc.MaxSteps)
        if "_" not in args["AnalysisFile"].casefold():
            if globals.verbosity >= 1:
                print(f"Number of steps to glue together is {NumJobs}")
            Anal.glue_stick(Umbrella, NumJobs=NumJobs, file=args["AnalysisFile"])
        if Umbrella.atom3 != 0:
            periodicity = "periodic"
        else:
            periodicity = "discrete"
        wham = WhamClass(args["AnalysisFile"], Umbrella.ConstForce, periodicity)
        Wham.Init_Wham(Umbrella, wham, WhamIgnore=AnalIgnore)
        Wham.Run_Wham(Umbrella, WhamIgnore=AnalIgnore)
    if args["Stage"].casefold() == "analysis": # pragma: no cover
        Anal.analysis(Umbrella)
    if args["Stage"].casefold() == "convergence":
        Wham.convergence(Calc, args["AnalysisFile"], equil_length, prod_length, Umbrella)
    if args["Stage"].casefold() == "vis":
        VisInit( Umbrella, args["AnalysisFile"])
        if globals.DryRun == "False":
            VisLoad(args["AnalysisFile"])
    endtime = time.time()
    print(f"Total time is {endtime - starttime}")

def class_init(args: dict):
    """ Initialises some key dictionary

    Args:
        args (dict): User inputs and variables

    Returns:
        Calc (CalcClass): Class containing calculation variables
        MM (MMClass): Class containing MD variables
        QM (QMClass): Class containing QM variables
        Umbrella (UmbrellaClass): Class containing Umbrella variables
        SLURM (SLURMClass): Class containing Slurm variables for HPC usage
    """
    Umbrella = UmbrellaClass(args,float(args["UmbrellaMin"]), int(args["UmbrellaBins"]),
                             float(args["StartDistance"]), float(args["UmbrellaWidth"]),)
    Calc = CalcClass(args)
    MM = MMClass()
    MM.CellVec = utils.get_cellVec(MM)
    QM = QMClass(args)
    bins = utils.init_bins(Umbrella.Bins, Umbrella.Width, Umbrella.Min)
    Umbrella.add_bins(bins)
    SLURM = SLURMClass(args)
    SLURM.set_IDNumber()
    return Calc, MM, QM, Umbrella, SLURM

def setup(QM: QMClass, Umbrella: UmbrellaClass):
    """
    Generates a setup directory and creates the syst-*.pdb files
    
    Args:
        QM (QMClass): QM class containing QM section
        Umbrella (UmbrellaClass): Umbrella Class containing colvar information
    """
    if path.isdir(f"{globals.WorkDir}setup") == False: ### Make setup directory. Not essential but makes dirs cleaner
        print("INFO: Making the setup directory" if globals.verbosity >= 2 else "", end="")
        os.mkdir(f"{globals.WorkDir}setup")
    utils.QM_Gen(QM.QMSel)
    if globals.DryRun == False: 
        print("INFO: Setting up the syst-qm.pdb file." if globals.verbosity >= 2 else "", end="")
        logfile = subprocess.run(["vmd", "-dispdev", "text", "-e", "qm_prep.tcl"],
                                     text = True, capture_output = True)
        with open(f"{globals.WorkDir}setup/tcl-qm.log","w") as f:
            print(logfile, file=f)
    utils.ColVarPDB_Gen(Umbrella)
    if globals.DryRun == False:
        print("INFO: Setting up the Colvar pdb file." if globals.verbosity >= 2 else "", end="")
        logfile = subprocess.run(["vmd", "-dispdev", "text", "-e", "Colvar_prep.tcl"],
                                    text = True, capture_output = True)
        print(logfile)
        with open(f"{globals.WorkDir}setup/tcl-colvar.log","w") as f:
            print(logfile, file=f)

def min(MM: MMClass, Calc: CalcClass, StartFile: str):    
    """sets up and runs the initial minimisation run (Pure MD)

    Args:
        MM (MMClass): Class containing dynamics information
        Calc (CalcClass): Class containing calculation information
        StartFile (str): Start coordinates pre minimisation.
    """
    print(globals.WorkDir)
    MM.Set_Ensemble("min")
    MM.Set_Outputs(1000, 100, 0)
    MM.Set_Length(10000)
    Calc.Job_Name("min")
    NAMD = NAMDClass(Calc, MM)
    NAMD.set_pme("on")
    if StartFile == "start.rst7":
        NAMD.set_startcoords(None, ambercoor="start.rst7", parm=globals.parmfile)
    else:
        NAMD.set_startcoords(bincoor=StartFile, ambercoor="start.rst7", parm=globals.parmfile)
    NAMD.set_cellvectors(MM.CellVec)
    file = FileGen.Namd_File(NAMD)
    utils.file_write(f"{globals.WorkDir}min.conf", [file])
    print(f"{globals.WorkDir}min.conf")
    if globals.DryRun == False: # pragma: no cover
        print("INFO: Running the minimisation script" if globals.verbosity >=2 else"", end="" )
        subprocess.run([MM.GPUNamd + f" +oneWthPerCore +setcpuaffinity +devices 0 min.conf > min_1.0.out"],
                       shell=True, capture_output=True)
        print("INFO: Minimisation complete, cleaning up directory" if globals.verbosity >=2 else"", end="" )
        output = subprocess.run(["mv min* ./setup", "cp ./setup/min_1.0.restart.coor .",
                    "cp ./setup/min_1.0.out .", "cp ./setup/min.conf ."],
                   shell=True, capture_output=True)
        output = subprocess.run(["cp ./setup/min_1.0.restart.coor .", ],
                   shell=True, capture_output=True)
        output = subprocess.run(["cp ./setup/min_1.0.out ."],
                   shell=True, capture_output=True)
        output = subprocess.run(["cp ./setup/min.conf ."],
                   shell=True, capture_output=True)
        
def heat(MM: MMClass, Calc: CalcClass):
    """
    Sets up and runs the heating calculation
    Args:
        MM (MMClass): Class containing dynamics information
        Calc (CalcClass): Class containing calculation information
    """
    MM.Set_Ensemble("heat")
    MM.Set_Outputs(200, 10, 100)
    MM.Set_Length(10000)
    Calc.Job_Name("heat")
    NAMD = NAMDClass(Calc, MM)
    NAMD.set_pme("on")
    NAMD.set_startcoords("min_1.0.restart.coor", ambercoor="start.rst7", parm=globals.parmfile)
    NAMD.set_cellvectors(MM.CellVec)
    file = FileGen.Namd_File(NAMD)
    utils.file_write(f"{globals.WorkDir}heat.conf", [file])
    if globals.DryRun == False: # pragma: no cover
        print("INFO: Running the heat calculation" if globals.verbosity >=2 else "", end="")
        subprocess.run([MM.GPUNamd + f" +oneWthPerCore +setcpuaffinity +devices 0 heat.conf > heat_1.0.out"],
                       shell=True, capture_output=True)
        subprocess.run(["mv heat* ./setup", ],
                   shell=True, capture_output=True)
        subprocess.run(["cp ./setup/heat_1.0.restart.coor .", ],
                    shell=True, capture_output=True)
        subprocess.run(["cp ./setup/heat_1.0.out .", ],
                    shell=True, capture_output=True)
        subprocess.run(["cp ./setup/heat.conf ."],
                    shell=True, capture_output=True)

def pull(Umbrella: UmbrellaClass, MM: MMClass, QM: QMClass, Calc: CalcClass):
    make_umbrellaDirs(Umbrella)
    print("INFO: Setting up pulls" if globals.verbosity >=2 else "", end="")
    MM.Set_Ensemble("NVT")
    MM.Set_Length(50,0.5)
    MM.Set_Outputs(5,1,10)
    MM.Set_Shake("none")
    MM.Set_Force(Umbrella.PullForce)
    Calc.Job_Name("pull")
    Calc.Set_OutFile("pull")
    NAMD = NAMDClass(Calc, MM)
    NAMD.set_cellvectors(MM.CellVec)
    Joblist = [None]*Umbrella.Bins
    for i in range(Umbrella.Bins):
       if abs(Umbrella.BinVals[i] - Umbrella.Start) < abs(Umbrella.Width * 0.5):
           Umbrella.add_start(i)
    for i in range(Umbrella.Bins):
        if i < Umbrella.StartBin:
            prevPull = f"../{i + 1}/pull_1.{i + 1}.restart.coor"
        elif i > Umbrella.StartBin:
            prevPull = f"../{i-1}/pull_1.{i-1}.restart.coor"
        elif i == Umbrella.StartBin:
            prevPull = f"../heat_1.0.restart.coor"#
            print(f"Pull starts from directory {i}")
        NAMD.set_qm(Calc, QM, i)
        NAMD.set_startcoords(prevPull)
        NAMD.set_colvars("colvars.pull.conf")
        file = FileGen.Namd_File(NAMD, window=i)
        utils.file_write(f"{globals.WorkDir}{i}/pull.conf", [file])
        colvarfile = utils.colvar_gen(Umbrella, i, "pull", Umbrella.PullForce )
        utils.file_write(f"{globals.WorkDir}{i}/colvars.pull.conf", [colvarfile])
        Joblist[i] = f"mkdir /dev/shm/RUNDIR_{i} ; cd ./{i} ; {MM.CPUNamd} pull.conf > pull_1.{i}.out ; cd ../ ; rm -r /dev/shm/RUNDIR_{i}"
    utils.file_write(f"{globals.WorkDir}pull.txt",Joblist)
    make_runfile(Umbrella, Joblist)
    if globals.DryRun == False: # pragma: no cover
        print("INFO: Running serial pulls." if globals.verbosity >=2 else "", end="")
        run_out = subprocess.run([f"sh {globals.WorkDir}pull.sh"], shell=True, capture_output=True)

def const(MM: MMClass, QM: QMClass, Calc: CalcClass, Umbrella: UmbrellaClass, PreviousJob:str):
    MM.Set_Force(Umbrella.ConstForce)
    if Calc.MaxSteps == 0:
        NumJobs = 1
        MM.Set_Length(MM.Steps, MM.TimeStep)
    else:
        NumJobs = math.ceil(MM.Steps/Calc.MaxSteps)
        MM.Set_Length(Calc.MaxSteps, MM.TimeStep)
    NAMD = NAMDClass(Calc, MM)
    NAMD.set_pme("off")
    NAMD.set_cellvectors(MM.CellVec)
    NAMD.set_colvars("colvars.const.conf")
    JobList = [None] * Umbrella.Bins * NumJobs
    for i in range(Umbrella.Bins):
        colvar = utils.colvar_gen(Umbrella, i, "constant", Umbrella.ConstForce)
        utils.file_write(f"{globals.WorkDir}{i}/colvars.const.conf",[colvar])
        for j in range(NumJobs):
            place = j * Umbrella.Bins + i
            JobList[place] = f"sleep 2 ; ( mkdir /dev/shm/RUNDIR_{i} ; cd ./{i} ; {MM.CPUNamd} {Calc.Name}_{j+1}.conf > {Calc.Name}_{j+1}.{i}.out ; cd ../ ; rm -r /dev/shm/RUNDIR_{i} ) &"
            NAMD.set_qm(Calc, QM, i)
            if j == 0:
                NAMD.set_startcoords(f"{PreviousJob}.{i}.restart.coor")
                file = FileGen.Namd_File(NAMD, j+1, i)
                utils.file_write(f"{globals.WorkDir}{i}/{Calc.Name}_{j + 1}.conf",
                                 [file])
            else:
                NAMD.set_startcoords(f"{Calc.Name}_{j}.{i}.restart.coor")
                file = FileGen.Namd_File(NAMD, j + 1, i)
                utils.file_write(f"{globals.WorkDir}{i}/{Calc.Name}_{j + 1}.conf",
                                 [file])
    for i in range(NumJobs):
        lines = i * Umbrella.Bins
        utils.file_write(f"{globals.WorkDir}{Calc.Name}_{i + 1}.txt",
                         JobList[lines : lines + Umbrella.Bins])
    if globals.DryRun == False:
        print("ERROR: Local umbrella run is not currently supported. We recommend you submit to a HPC")

def make_umbrellaDirs(Umbrella: UmbrellaClass):
    """Makes the directory structure for umbrella calculation

    Args:
        Umbrella (UmbrellaClass): Contains information about bins.
    """
    print("INFO: Making umbrella directories" if globals.verbosity >= 2 else "", end="")
    for i in range(Umbrella.Bins):
        dir_path = str(globals.WorkDir) + str(i)
        if path.exists(dir_path):
            print(f"INFO: {str(i)} exists. Deleting!" if globals.verbosity >= 2 else "", end="")
            try:
                os.rmdir(dir_path)
            except OSError:
                print(f"INFO: {i} directory not empty, deletion failed..." if globals.verbosity >= 2 else "", end="")
                pass
        try:
            os.mkdir(dir_path)
        except FileExistsError:
            print(f"INFO: {i} directory exists, Skipping making new directory" if globals.verbosity >= 2 else "", end="")
            pass

def make_runfile(Umbrella: UmbrellaClass, Joblist: list):
    """Generates the pull.sh file to be run either locally or on the HPC

    Args:
        Umbrella (UmbrellaClass): Contains bin information
        Joblist (list): List of jobs to run.
    """
    print("INFO: Generating pull.sh script." if globals.verbosity >=2 else "", end="")
    with open(f"{globals.WorkDir}pull.sh", 'w') as f:
        print("#!/bin/bash", file=f)
        for i in range(Umbrella.StartBin, Umbrella.Bins):
            print(Joblist[i], file=f)
        for i in range(0, Umbrella.StartBin):
            print(Joblist[Umbrella.StartBin - i - 1], file=f)

def VisInit(Umbrella: UmbrellaClass, File: str, extension="restart.coor"):
    """Initialises a tcl script to allow the opening of a stage of the calculation as a single trajectory..

    Args:
        Umbrella (UmbrellaClass): UmbrellaClass containing directory structure information
        File (str): name of stage eg equil_1
        extension (str, optional): File extenstion of the structures. Loads NAMD restart files but could be adapted. Defaults to "restart.coor".
    """
    Lines = [None] * (Umbrella.Bins + 1)
    Lines[0] = f"mol new {globals.parmfile}"
    for i in range(Umbrella.Bins):
        Lines[i+1] = f"mol addfile ./{i}/{File}.{i}.{extension}"
    utils.file_write(f"{globals.WorkDir}/{File}_load.tcl", lines=Lines)

def VisLoad(File:str): # pragma: no cover
    """Loads up a tcl script into vmd

    Args:
        File (str): file name of step trying to visualise eg. equil_1
    """
    subprocess.run([f"vmd -e {File}_load.tcl"], shell=True, capture_output=True)