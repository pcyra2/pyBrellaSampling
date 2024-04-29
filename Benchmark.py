import pyBrellaSampling.InputParser as input
import pyBrellaSampling.FileGen as FileGen
import pyBrellaSampling.utils as utils
from pyBrellaSampling.UserVars.QM_Methods import *
from pyBrellaSampling.classes import *
import sys
import time
import numpy
import subprocess
import pandas
from pyBrellaSampling.utils import kcal 

def main_cli():
    starttime = time.time()
    args = input.Benchmark_input(sys.argv[1:])
    Benchmark(args)
    endtime = time.time()

def Benchmark(args):
    verbosity = args.verbosity
    print("Welcome to the Benchmark suite \n" if verbosity > 1 else "",end="")
    WorkDir = args.WorkDir
    BenchType = args.BenchmarkType
    Stage=args.Stage
    ORCA = ORCAClass(args.Path)
    if BenchType.casefold() == "energy":
        ORCA.set_type("SP")
    elif BenchType.casefold() == "gradient":
        ORCA.set_type("Engrad")
    elif BenchType.casefold() == "structure":
        ORCA.set_type("OPT")
    else:
        print("ERROR: {BenchType} Calculation not supported... Ending setup. ")
        return ValueError
    ORCA.set_convergence(args.SCF)
    ORCA.set_dificulty(args.Convergence)
    ORCA.change_autostart(args.Restart)
    ORCA.set_grid(args.Grid)
    ORCA.set_cores(args.Cores)
    ORCA.set_extras(args.Extras)
    mols = get_molecules(args.ReactionList)
    Structures = []
    for i in range(len(mols)):
        molecule = MolClass(str(mols[i]))
        data = utils.file_read(f"{WorkDir}/{args.CoordinateLoc}/{mols[i]}.xyz")
        nat = int(data[0])
        info = data[1].split() ## Charge, Spin
        molecule.set_charge(int(info[0]), int(info[1]))
        a = numpy.empty([nat], dtype=str)
        x = numpy.zeros(nat)
        y = numpy.zeros(nat)
        z = numpy.zeros(nat)
        print(f"DEBUG: Number of atoms for {mols[i]} is {nat} \n" if verbosity > 1 else "",end="")
        for j in range(2,len(data)):
            cols = data[j].split()
            print(f"DEBUG: Atom line is {cols}\n" if verbosity > 1 else "",end="")
            if len(cols) != 4:
                print("WARNING: There is more than 4 columns in this coordinate file, you should check that they are correct." if verbosity > 0 else "", )
            a[j-2] = str(cols[0])
            x[j-2] = float(cols[1])
            y[j-2] = float(cols[2])
            z[j-2] = float(cols[3])
        molecule.add_coordinates(a, x, y, z, nat)
        Structures.append(molecule)
    if args.DryRun == "False" and Stage =="init":
        calcs = init_benchmark(Structures=Structures, Methods=Functionals, Basis_sets=Basis_Sets,
                               Dispersion_list=Dispersion_Corrections, ORCA=ORCA, verbosity=verbosity, WorkDir=WorkDir, run=True)
    elif args.DryRun == "False" and Stage =="calc":
        calcs = init_benchmark(Structures=Structures, Methods=Functionals, Basis_sets=Basis_Sets,
                               Dispersion_list=Dispersion_Corrections, ORCA=ORCA, verbosity=verbosity, WorkDir=WorkDir, run=True)
        run_benchmark(WorkDir = WorkDir, verbosity=verbosity)
    else:
        calcs = init_benchmark(Structures=Structures, Methods=Functionals, Basis_sets=Basis_Sets,
                               Dispersion_list=Dispersion_Corrections, ORCA=ORCA, verbosity=verbosity, WorkDir=WorkDir, run=False)
#    if args.DryRun == "False":
#        run_benchmark(WorkDir = WorkDir, verbosity=verbosity)
    if Stage == "analysis":
        Finished, Failed = get_energy(calcs)
        timing = []
        for i in Finished:
            # print(i.scf_energy)
            timing.append(i.time)
        print(f"REPORT: Total number of completed jobs is {len(Finished)}, {len(Failed)} have failed...")
        fix_script = []
        for i in Failed:
            print(f"Failed Job is {i.path}\n " if verbosity > 0 else "", end="")
            print(f"Reason is {i.reason}\n" if verbosity > 0 else "", end="")
            fix_script.append(i.runline)
        utils.file_write(f"job_fix.dat", fix_script)
        if args.DryRun == "False":
            print("Running fix")
            run_benchmark(WorkDir = WorkDir, verbosity=verbosity, file="job_fix.dat")
        Reactions = Benchmark_Energy(Finished, args.ReactionList)
        df = Error_Gen(Method_list=Functionals, Basis_list=Basis_Sets, Dispersion_list=Dispersion_Corrections, Reactions=Reactions)
        Error_Anal(df,Method_list=Functionals, Basis_list=Basis_Sets, Dispersion_list=Dispersion_Corrections,)
        



def get_molecules(file):
    """Obtains the names of all structures to be used in the benchmark"""
    coordinates = []
    reactionfile = utils.file_read(f"{file}")
    for i in reactionfile:
        molecules = i.split()
        for j in range(2):
            if molecules[j] not in coordinates:
                coordinates.append(str(molecules[j]))
    print(coordinates)
    return coordinates


def init_benchmark(Structures, Methods, Basis_sets, Dispersion_list, ORCA, verbosity=1, WorkDir="./", run=False):
    """Sets up the folder structure and generates the input files and job lists."""
    print("INFO: Setting up benchmark\n" if verbosity > 1 else "", end="")    
    jobfile=[]
    calculations = []
    for i in Methods:
        try:
            os.mkdir(f"{WorkDir}/{i}")
        except FileExistsError:
            print(f"WARNING: {i} directory exists, skipping \n" if verbosity > 1 else "", end="")
        ORCA.set_method(i)
        if i in utils.DispersionCorrFunc:
            print(f"WARNING: {i} contains a dispersion correction, therefore skipping other corrections. \n" if verbosity > 1 else "", end="")
            dispersion = [""]
        else:
            dispersion = Dispersion_list
        for j in Basis_sets:
            ORCA.set_basis(j)
            try:
                os.mkdir(f"""{WorkDir}/{i}/{j.replace("*","s")}""")
            except FileExistsError:
                print(f"WARNING: {i} {j} directory exists, skipping \n" if verbosity > 1 else "", end="")
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
                    print(f"WARNING: {i} {j} {directory} directory exists, skipping \n" if verbosity > 1 else "", end="")
                for l in range(len(Structures)):
                    path = f"""{WorkDir}/{i}/{j.replace("*","s")}/{directory}/{Structures[l].name}"""
                    try:
                        os.mkdir(path)
                    except FileExistsError:
                        print(f"WARNING: {i} {j} {directory} {Structures[l].name} directory exists, skipping \n" if verbosity > 1 else "", end="")
                    file = FileGen.ORCA_FileGen(Structures[l], ORCA) 
                    if run == True:
                        utils.file_write(f"{path}/ORCA.inp", [file])
                    jobfile.append(f"cd {path} ; WD=$PWD ; cp * /dev/shm/QM ; cd /dev/shm/QM ;  {ORCA.path} ORCA.inp '--use-hwthread-cpus --bind-to hwthread' > ORCA.out ; mv * $WD ; cd $WD ; cd ../../../../")
                    calculation = QMCalcClass(Structures[l].name, i, j, directory)
                    calculation.set_path(path)
                    calculation.set_runline(f"cd {path} ; WD=$PWD ; cp * /dev/shm/QM ; cd /dev/shm/QM ;  {ORCA.path} ORCA.inp '--use-hwthread-cpus --bind-to hwthread' > ORCA.out ; mv * $WD ; cd $WD ; cd ../../../../")
                    calculations.append(calculation)
    if run == True:
        utils.file_write("jobfile.dat", jobfile)
    return calculations

def run_benchmark(WorkDir, verbosity, file="jobfile.dat"):
    print("Running")
    jobs = utils.file_read(f"{WorkDir}{file}")
    print(f"Running Jobs \n" if verbosity > 0 else "", end="")
    for lines in jobs:
        print(f"Running {lines} \n" if verbosity > 1 else "", end="")
        runout = subprocess.run([f"{lines}"], shell=True, capture_output=True)
        print(f"{runout.stdout}\n " if verbosity > 1 else "", end="")

def get_energy(calculations):
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
       # for line in lines:
       #     if "FINAL SINGLE POINT ENERGY" in line:
       #         completed = 1
       #         words = line.split()
       #         calc.set_scfenergy(float(words[4]))
       #     if "TOTAL RUN TIME" in line:
       #         words = line.split()
       #         days = int(words[3])
       #         hours = int(words[5]) + 24*days
       #         mins = int(words[7]) + 60*hours
       #         secs = int(words[9]) + 60*mins
       #         totaltime = secs
       #         calc.set_time(totaltime)
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

def Benchmark_Energy(Finished, ReactionList):
    DataFrame=pandas.DataFrame([F.__dict__ for F in Finished ], columns=["molecule", "functional", "basis", "dispersion", "time", "TotalEnergy", "SCFEnergy", "vdw"])
    DataFrame.to_csv("./DataFrame", sep="\t")
    data = utils.file_read(ReactionList)
    Reactions = []
    for i in data:
        cols = i.split()
        reaction = ReactionClass(str(cols[0]), float(cols[2]), str(cols[1]), float(cols[3]), "df-CCSD(T)", "cc-pVTZ", "NONE" )
        Benchmark = reaction
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
    

def Error_Gen(Method_list, Basis_list, Dispersion_list, Reactions):
    MAE_list = []
    MaxE_list = []
    mList = []
    bList = []
    dList = []
    TimeList = []
    for i in Method_list:
        for j in Basis_list:
            for k in Dispersion_list:
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
    df.to_csv("Errors.dat", sep="\t")
    return df

def Error_Anal(df,Method_list, Basis_list, Dispersion_list,):
    for i in Method_list:
        Errors = numpy.average(df[df["Method"] == i]["MAE"])
        Time = numpy.average(df[df["Method"] == i]["Average Time"])
        print(f"{i} average MAE is {Errors}")
        print(f"{i} average Time is {Time}")
        Score = Errors * Time
        print(f"{i} average Score is {Score}")
        print("\n")
    for i in Basis_list:
        Errors = numpy.average(df[df["Basis"] == i]["MAE"])
        Time = numpy.average(df[df["Basis"] == i]["Average Time"])
        print(f"{i} average MAE is {Errors}")
        print(f"{i} average Time is {Time}")
        Score = Errors * Time
        print(f"{i} average Score is {Score}")
        print("\n")
    for i in Dispersion_list:
        print(str(i))
        Errors = numpy.average(df[df["Dispersion"] == str(i)]["MAE"])
        Time = numpy.average(df[df["Dispersion"] == str(i)]["Average Time"])
        print(f"{i} average MAE is {Errors}")
        print(f"{i} average Time is {Time}")
        Score = Errors * Time
        print(f"{i} average Score is {Score}")
        print("\n")
    i = "NONE"
    print(str(i))
    Errors = numpy.average(df[df["Dispersion"] == str(i)]["MAE"])
    Time = numpy.average(df[df["Dispersion"] == str(i)]["Average Time"])
    print(f"{i} average MAE is {Errors}")
    print(f"{i} average Time is {Time}")
    Score = Errors * Time
    print(f"{i} average Score is {Score}")
    print("\n")
    df.sort_values("MAE").to_csv("Errors_Ranked.dat", sep="\t")