from pyBrellaSampling.Tools import FileGen, utils
from pyBrellaSampling.UserVars.QM_Methods import * # End User can change this. 
# from pyBrellaSampling.UserVars.Defaultinputs import Benchmark_Inp
# from pyBrellaSampling.UserVars.SoftwarePaths import ORCA_PATH
from pyBrellaSampling.Tools.classes import *
from pyBrellaSampling.Tools.InputParser import BenchmarkInput as Input_Parser
import sys
import time
import numpy
import subprocess
import pandas
import os
from pyBrellaSampling.Tools.utils import kcal 


def main_cli():
    starttime = time.time()
    args = Input_Parser(sys.argv[1:])
    global verbosity
    verbosity = args["verbosity"]
    global WorkDir
    WorkDir = args["WorkDir"]
    global DryRun
    DryRun = args["DryRun"]
    Stage=args["Stage"]
    print("Welcome to the Benchmark suite \n" if verbosity >= 2 else "",end="")
    ORCA = ORCA_init(args)
    Structures = Get_Structures(args)
    calcs = init_benchmark(Structures=Structures, ORCA=ORCA, )
    if Stage == "run" and DryRun == False:
        Finished, Failed = get_energy(calcs)
        if len(Finished) == 0: ### No jobs have started therefore a clean run.
            run_benchmark("jobfile.dat")
        else: ## Some jobs are complete, therefore a re-run to fix calculations or add to the dataset.
            print(f"WARNING: Some jobs have already completed. If you wish to overwrite them, please delete the outputs. \n" 
                  if verbosity >= 1 else "", end="")
            print(f"WARNING: Running jobs that have not completed \n" if verbosity >= 1 else "", end = "")
            time.sleep(10)
            for i in Failed:
                fix_script = []
                fix_script.append(i.runline)
            utils.file_write(f"job_fix.dat", fix_script)
            run_benchmark("job_fix.dat")
    elif Stage == "analysis":
        Finished, Failed = get_energy(calcs)
        timing = []
        for i in Finished:
            # print(i.scf_energy)
            timing.append(i.time)
        print(f"REPORT: Total number of completed jobs is {len(Finished)}, {len(Failed)} have failed... \n")
        fix_script = []
        for i in Failed:
            print(f"REPORT: Failed Job is {i.path}\n" if verbosity >= 2 else "", end="")
            print(f"REPORT: Reason is {i.reason}\n" if verbosity >= 2 else "", end="")
            fix_script.append(i.runline)
        Reactions = Benchmark_Energy(Finished, args["ReactionList"])
        df = Error_Gen(Reactions=Reactions)
        Error_Anal(df)
    endtime = time.time()
    print(f"INFO: Calculation time is {endtime - starttime} \n" if verbosity >=2 else "", end="")
    
def ORCA_init(args: dict):
    """
    Initis ORCA Class

    Args:
        args (dict): User Variables

    Raises:
        ValueError: When Benchmark type is not supported.

    Returns:
        ORCA (ORCAClass): Initialised ORCA class
    """
    print("DEBUG: ORCA init \n" if verbosity >= 3 else "", end="")
    BenchType = args["BenchmarkType"]
    ORCA = ORCAClass(args["Path"])
    if BenchType.casefold() == "energy":
        ORCA.set_type("SP")
    elif BenchType.casefold() == "gradient":
        ORCA.set_type("Engrad")
    elif BenchType.casefold() == "structure":
        ORCA.set_type("OPT")
    else:
        print("ERROR: {BenchType} Calculation not supported... Ending setup. ")
        return ValueError
    ORCA.set_convergence(args["SCF"])
    ORCA.set_dificulty(args["Convergence"])
    ORCA.change_autostart(args["Restart"])
    ORCA.set_grid(args["Grid"])
    ORCA.set_cores(args["Cores"])
    ORCA.set_extras(args["Extras"])
    return ORCA

def Get_Structures(args: dict):
    """
    Reads the reaction file to obtain a list of structures. Then reads in their coordinates ect.
    Args:
        args (dict): User Variables, specifically "ReactionList" and "CoordLoc"

    Returns:
        Structures (list): List of Stuctures and their coordinates. 

    """
    print("DEBUG: Get_Structures \n" if verbosity >= 3 else "", end="")
    mols = get_molecules(args["ReactionList"])
    Structures = []
    coordloc=args["CoordinateLoc"]
    for i in range(len(mols)):
        molecule = MolClass(str(mols[i]))
        data = utils.file_read(f"{WorkDir}/{coordloc}/{mols[i]}.xyz")
        nat = int(data[0])
        info = data[1].split() ## Charge, Spin
        molecule.set_charge(int(info[0]), int(info[1]))
        a = numpy.empty([nat], dtype=str)
        x = numpy.zeros(nat)
        y = numpy.zeros(nat)
        z = numpy.zeros(nat)
        print(f"DEBUG: Number of atoms for {mols[i]} is {nat} \n" if verbosity >= 3 else "",end="")
        for j in range(2,len(data)):
            cols = data[j].split()
            # print(f"DEBUG: Atom line is {cols}\n" if verbosity >= 3 else "",end="")
            if len(cols) != 4:
                print("WARNING: There is more than 4 columns in this coordinate file, you should check that they are correct." if verbosity >= 1 else "", end="")
            a[j-2] = str(cols[0])
            x[j-2] = float(cols[1])
            y[j-2] = float(cols[2])
            z[j-2] = float(cols[3])
        molecule.add_coordinates(a, x, y, z, nat)
        Structures.append(molecule)
    return Structures
       

def get_molecules(file: str):
    """
    Obtains the names of all structures to be used in the benchmark

    Args:
        file (str): file containing a list of 2 step reactions, each column is a structure name 

    Returns:
        coordinates (list[str]): List of unique coordinates that then need to be found for benchmarking. 
    """
    coordinates = []
    reactionfile = utils.file_read(f"{file}")
    for i in reactionfile:
        molecules = i.split() 
        for j in range(2):
            if molecules[j] not in coordinates: # Ensures you dont duplicate structures.
                coordinates.append(str(molecules[j]))
    print(f"DEBUG: coordinates are: {coordinates}")
    return coordinates

def init_benchmark(Structures: list, ORCA: ORCAClass):
    """
    Sets up the folder structure and generates the input files and job lists.

    Args:
        Structures (list): List of Structures and their information
        ORCA (ORCAClass): Initialised ORCA class

    Raises:
        FileExistsError: When directory exists, returns a warning.

    Returns:
        calculations (list[QMCalcClass]): Returns a list of calculations
    """
    print("INFO: Setting up benchmark\n" if verbosity >= 1 else "", end="")    
    jobfile=[]
    calculations = []
    for i in Functionals:
        try:
            os.mkdir(f"{WorkDir}/{i}")
        except FileExistsError:
            print(f"WARNING: {i} directory exists, skipping \n" if verbosity >= 1 else "", end="")
        ORCA.set_method(i)
        if i in DispersionCorrFunc:
            print(f"WARNING: {i} contains a dispersion correction, therefore skipping other corrections. \n" if verbosity > 1 else "", end="")
            dispersion = [""]
        else:
            dispersion = Dispersion_Corrections
        for j in Basis_Sets:
            ORCA.set_basis(j)
            try:
                os.mkdir(f"""{WorkDir}/{i}/{j.replace("*","s")}""")
            except FileExistsError:
                print(f"WARNING: {i} {j} directory exists, skipping \n" if verbosity >= 1 else "", end="")
            for k in dispersion:
                if i == "HF" and k == "D3ZERO":
                    continue
                ORCA.set_dispersion(k)
                if k == "":
                    directory = "NONE"
                else:
                    directory = k
                try:
                    os.mkdir(f"""{WorkDir}/{i}/{j.replace("*","s")}/{directory}""")
                except FileExistsError:
                    print(f"WARNING: {i} {j} {directory} directory exists, skipping \n" if verbosity >= 1 else "", end="")
                for l in range(len(Structures)):
                    path = f"""{WorkDir}/{i}/{j.replace("*","s")}/{directory}/{Structures[l].name}"""
                    try:
                        os.mkdir(path)
                    except FileExistsError:
                        print(f"WARNING: {i} {j} {directory} {Structures[l].name} directory exists, skipping \n" if verbosity > 1 else "", end="")
                    file = FileGen.ORCA_FileGen(Structures[l], ORCA) 
                    if DryRun == False:
                        utils.file_write(f"{path}/ORCA.inp", [file])
                    jobfile.append(f"cd {path} ; WD=$PWD ; cp * /dev/shm/QM ; cd /dev/shm/QM ;  {ORCA.path} ORCA.inp '--use-hwthread-cpus --bind-to hwthread' > ORCA.out ; mv * $WD ; cd $WD ; cd ../../../../")
                    calculation = QMCalcClass(Structures[l].name, i, j, directory)
                    calculation.set_path(path)
                    calculation.set_runline(f"cd {path} ; WD=$PWD ; cp * /dev/shm/QM ; cd /dev/shm/QM ;  {ORCA.path} ORCA.inp '--use-hwthread-cpus --bind-to hwthread' > ORCA.out ; mv * $WD ; cd $WD ; cd ../../../../")
                    calculations.append(calculation)
    utils.file_write("jobfile.dat", jobfile)
    return calculations

def run_benchmark(file:str = "jobfile.dat"):
    """
    Runs each line in the jobfile 
    Args:
        file (str): File containing list of commands to run. 

    """

    jobs = utils.file_read(f"{WorkDir}{file}")
    print(f"INFO: Running Jobs \n" if verbosity >= 2 else "", end="")
    for lines in jobs:
        print(f"INFO: Running {lines} \n" if verbosity >= 2 else "", end="")
        runout = subprocess.run([f"{lines}"], shell=True, capture_output=True)
        print(f"{runout.stdout}\n " if verbosity >= 2 else "", end="")

def get_energy(calculations: list):
    """
    Checks a list of calculations, if complete -> extracts data, if incomplete -> logs.    
    Args:
        calculations (list): List of calculations to check (Usually all calculations. )
    
    Returns:
        Finished (list): List of finished calculations and relevant information.
        Failed (list): List of failed calculations and an estimate as to why they failed. 

    """
    Finished = []
    Failed = []
    for calc in calculations:
        started = 1
        try:
            lines = utils.file_read(f"{calc.path}/ORCA.out")
        except FileNotFoundError:
            calc.set_reason("Job never started")
            started = 0
            lines = []
        completed = 0
        if len(lines) > 0:
            line=lines[len(lines)-1]
        else:
            line=""
        if "TOTAL RUN TIME" in line:
            completed = 1
            words = line.split()
            days = int(words[3])
            hours = int(words[5]) + 24*days
            mins = int(words[7]) + 60*hours
            secs = int(words[9]) + 60*mins
            totaltime = secs
            calc.set_time(totaltime)
        if started == 1:
            try:
                lines = utils.file_read(f"{calc.path}/ORCA_property.txt")
            except FileNotFoundError:
                lines=[]
                completed=0
            for line in lines:
                if "SCF Energy:" in line:
                    words=line.split()
                    calc.set_scfenergy(float(words[2]) * kcal)
                if "Van der Waals Correction:" in line:
                    words=line.split()
                    calc.set_vdw(float(words[4]) * kcal)
                if "Total Energy" in line:
                    words=line.split()
                    calc.set_TotalEnergy(float(words[2]) * kcal)
        if completed == 1:
            Finished.append(calc)
        else:
            if started == 1:
                ending = subprocess.run([f"tail -n 2 {calc.path}/ORCA.out"], shell=True, capture_output=True)
                calc.set_reason(ending.stdout)
            Failed.append(calc)
    return Finished, Failed

def Benchmark_Energy(Finished: list, ReactionList: str):
    """
    Generates a dataframe/table of all calculations including energies, timings, ect.
    Args:
        Finished (list): List of finished calculations (Excludes failed calculations)
        ReactionList (str): File name of reaction list as this contains reference energies and reactions.

    Returns:
        Reactions (list): List of reactions with energy values. 
    """
    DataFrame=pandas.DataFrame([F.__dict__ for F in Finished ], columns=["molecule", "functional", "basis", "dispersion", "time", "TotalEnergy", "SCFEnergy", "vdw"])
    DataFrame.to_csv("./DataFrame", sep="\t")
    data = utils.file_read(ReactionList)
    Reactions = []
    for i in data:
        cols = i.split()
        reaction = ReactionClass(str(cols[0]), float(cols[2]), str(cols[1]), float(cols[3]), "df-CCSD(T)", "cc-pVTZ", "NONE" )
        # Benchmark = reaction
        Reactions.append(reaction)
        BenchDelta = reaction.deltaE
        for j in range(len(Finished)):
            for k in range(j,len(Finished)):
                a = Finished[j]
                b = Finished[k]
                if a.functional != b.functional or a.basis != b.basis or a.dispersion != b.dispersion:
                    continue
                if a.molecule == str(cols[0]) and b.molecule == str(cols[1]):
                    reaction = ReactionClass(a.molecule, a.TotalEnergy, b.molecule, b.TotalEnergy, a.functional, a.basis, a.dispersion)
                    err = numpy.absolute(reaction.deltaE - BenchDelta)
                    reaction.add_Error(err)
                    reaction.add_timings(a.time, b.time)
                if reaction in Reactions:
                    continue
                else:
                    Reactions.append(reaction)
    ReactionsDF = pandas.DataFrame([R.__dict__ for R in Reactions])
    ReactionsDF.to_csv("./ReactionDF.dat", sep="\t")
    return Reactions
    
def Error_Gen(Reactions: list):
    """
    Calculates the MAE from the standard energies from reactionlist file. 
    Args:
        Reactions (list): List of reactions, provided by the reactionlist input file. 
    
    Returns:
        df (pandas.DataFrame): Dataframe containing an overview of the benchmark.
    """
    MAE_list = []
    MaxE_list = []
    mList = []
    bList = []
    dList = []
    TimeList = []
    for i in Functionals:
        for j in Basis_Sets:
            for k in Dispersion_Corrections:
                if k == "":
                    k = "NONE"
                Errors = []
                timer = []
                for l in Reactions:
                    if l.functional == i and l.basis == j and l.dispersion == k:
                        Errors.append(numpy.absolute(l.error))
                        timer.append(l.Timing)
                if len(Errors) > 0:
                    MAE = numpy.average(Errors)
                    MaxE = numpy.max(Errors)
                    MAE_list.append(MAE)
                    MaxE_list.append(MaxE)
                    mList.append(i)
                    bList.append(j)
                    dList.append(k)
                    TimeList.append(numpy.average(timer))
    dictionary = {"Method" : mList, "Basis" : bList, "Dispersion" : dList, "MAE" : MAE_list, "MaxError" : MaxE_list, "Average Time" : TimeList}
    df = pandas.DataFrame(dictionary)
    df.to_csv("Summary.dat", sep="\t")
    return df

def Error_Anal(df: pandas.DataFrame):
    """
    Performs basic error analysis on the summary dataframe. Score = AverageError x Average Time

    Args:
        df (pandas.DataFrame): Summary dataframe

    """
    for i in Functionals:
        Errors = numpy.average(df[df["Method"] == i]["MAE"])
        Time = numpy.average(df[df["Method"] == i]["Average Time"])
        print(f"{i} average MAE is {Errors}")
        print(f"{i} average Time is {Time}")
        Score = Errors * Time
        print(f"{i} average Score is {Score}")
        print("\n")
    for i in Basis_Sets:
        Errors = numpy.average(df[df["Basis"] == i]["MAE"])
        Time = numpy.average(df[df["Basis"] == i]["Average Time"])
        print(f"{i} average MAE is {Errors}")
        print(f"{i} average Time is {Time}")
        Score = Errors * Time
        print(f"{i} average Score is {Score}")
        print("\n")
    for i in Dispersion_Corrections:
        print(str(i))
        Errors = numpy.average(df[df["Dispersion"] == str(i)]["MAE"])
        Time = numpy.average(df[df["Dispersion"] == str(i)]["Average Time"])
        print(f"{i} average MAE is {Errors}")
        print(f"{i} average Time is {Time}")
        Score = Errors * Time
        print(f"{i} average Score is {Score}")
        print("\n")
    i = "NONE"
    print("No dispersion correction")
    Errors = numpy.average(df[df["Dispersion"] == str(i)]["MAE"])
    Time = numpy.average(df[df["Dispersion"] == str(i)]["Average Time"])
    print(f"No dispersion correction average MAE is {Errors}")
    print(f"No dispersion correction average Time is {Time}")
    Score = Errors * Time
    print(f"No dispersion correction average Score is {Score}")
    print("\n")
    df.sort_values("MAE").to_csv("Errors_Ranked.dat", sep="\t")