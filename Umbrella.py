import pyBrellaSampling.generate as generate
import pyBrellaSampling.utils
import pyBrellaSampling.utils as utils
import pyBrellaSampling.FileGen as FileGen
import pyBrellaSampling.wham as Wham
import pyBrellaSampling.analysis as Anal
from pyBrellaSampling.classes import *
import pyBrellaSampling.InputParser as input

import subprocess
import matplotlib.pyplot as plt
import os
import os.path as path
import math
def main(args):
    Umbrella = UmbrellaClass(args,args.UmbrellaMin, args.UmbrellaBins,
                             args.StartDistance, args.UmbrellaWidth,)
    Job = JobClass(args)
    Calc = CalcClass(args)
    MM = MMClass(args)
    QM = QMClass(args)
    bins = utils.init_bins(Umbrella.Bins, Umbrella.Width, Umbrella.Min)
    Umbrella.add_bins(bins)
    # print(Umbrella.atom1, Umbrella.atom2, Umbrella.atom3, Umbrella.atom4)
    equil_length = 2000
    prod_length = 8000
    SLURM = SLURMClass(args)
    # UmbrellaCalculation()
    if path.isdir(f"{Job.WorkDir}setup") == False: ### Make setup directory. Not essential but makes dirs cleaner
        os.mkdir(f"{Job.WorkDir}setup")
    if Job.Stage == "init" or Job.Stage == "full":
        utils.QM_Gen(QM.QMSel, Job.WorkDir)
        if args.DryRun == "False":
            if Job.Verbosity >= 1:
                print("Setting up the QM pdb file.")
            logfile = subprocess.run(["vmd", "-dispdev", "text", "-e", "qm_prep.tcl"],
                                     text = True, capture_output = True)
            with open(f"{Job.WorkDir}tcl-qm.log","w") as f:
                print(logfile, file=f)
        utils.ColVarPDB_Gen(Umbrella, Job)
        if args.DryRun == "False":
            if Job.Verbosity >= 1:
                print("Setting up the Colvar pdb file.")
            logfile = subprocess.run(["vmd", "-dispdev", "text", "-e", "Colvar_prep.tcl"],
                                     text = True, capture_output = True)
            with open(f"{Job.WorkDir}tcl-colvar.log","w") as f:
                print(logfile, file=f)
    if Job.Stage == "full" or Job.Stage == "min":
        # generate.Min_Setup(Calc, Job, MM, QM)
        min_setup(MM, Calc, Job, args.StartFile)
        if args.DryRun == "False":
            min_run(MM,Job)
    if Job.Stage == "full" or Job.Stage == "heat":
        heat_setup(MM, Calc, Job)
        if args.DryRun == "False":
            heat_run(MM,Job)
    if Job.Stage == "setup" or Job.Stage == "full" or Job.Stage == "pull":
        make_umbrellaDirs(Umbrella, Job)
        pull_setup(Umbrella, MM, QM, Calc, Job)
        if args.DryRun == "False":
            output = run_pullScript()
            if Job.Verbosity >= 1:
                print(output)
    if Job.Stage == "equil" or Job.Stage == "full":
        MM.Set_Ensemble("NVT")
        Calc.Job_Name("equil")
        MM.Set_Length(equil_length, 0.5)  # 2000 steps at 0.5 fs = 1 ps ~ 1 day
        MM.Set_Outputs(100, 100, 80)  # Timings, Restart, Trajectory
        equil_setup(MM, QM, Job, Calc, Umbrella, "pull_1")
        SLURM.set_arrayJob("equil_1.txt", 54)
        print("gen slurm")
        utils.slurm_gen("equil", SLURM, "sh array_job.sh", Job.WorkDir)
        if args.DryRun == "False":
            equil_run(Job)
            # subprocess.run(["sh equil.sh"], shell=True, capture_output=True)
    if Job.Stage == "prod" or Job.Stage == "full":
        Calc.Job_Name("prod")
        MM.Set_Ensemble("NVT")
        MM.Set_Length(prod_length, 0.5)  # 8000 steps at 0.5 fs = 4 ps, ~ 3.5 days
        MM.Set_Outputs(100, 100, 80)  # Timings, Restart, Trajectory
        equil_setup(MM, QM, Job, Calc, Umbrella, f"equil_{math.ceil(equil_length/Calc.MaxSteps)}")
        if args.DryRun == "False":
            equil_run(Job)
    if Job.Stage == "wham" or Job.Stage == "full":
        if Calc.MaxSteps == 0:
            NumJobs = 1
        else:
            if "prod" in args.WhamFile.casefold():
                steps = prod_length
            elif "equil" in args.WhamFile.casefold():
                steps = equil_length
            NumJobs = math.ceil(steps / Calc.MaxSteps)
        if Job.Verbosity >= 1:
            print(f"Number of steps to glue together is {NumJobs}")
        Anal.glue_stick(Umbrella, Job, NumJobs=NumJobs, file=args.WhamFile)
        if Umbrella.atom3 != 0:
            periodicity = "Periodic"
        else:
            periodicity = "discrete"
        wham = WhamClass(args.WhamFile, Umbrella.ConstForce, periodicity)
        Wham.Init_Wham(Job, Umbrella, wham)
        Wham.Run_Wham(Job, Umbrella)
    if Job.Stage == "analysis" or Job.Stage == "full":
        Labels = LabelClass("../complex.parm7")
        Labels = input.BondsInput(f"{Job.WorkDir}Bonds.dat", Labels)
        # Labels.add_bond("13722,13724", "C4-Hup", 1.1) # Atom index
        # Labels.add_bond("13722,13723", "C4-Hdown", 1.1) # Atom index
        # Labels.add_bond("13730,13715", "C2-C7", 1.6)
        # Labels.add_bond("13717,13730", "C3-C7", 2.5)
        # Labels.add_bond("13722,13708", "C4-O", 2.8)
        Labels = input.DihedralInput(f"{Job.WorkDir}Dihedral.dat", Labels)
        # Labels.add_dihedral("13741,13738,13735,13730", "SobVert", 50, "Sob", -60, "Vert")
        # Labels.add_dihedral("13761,13756,13730,13731", "C7Methyl", 30, "up", 180, "down")
        # Labels.add_dihedral("13761,13756,13717,13718", "C3Methyl", 50, "up", 180, "down")
        # Labels.add_dihedral("13715,13730,13728,13725", "c6Ring", -55, "Chair", -20, "Boat")
        core_load = []
        core_load.append("mol new complex.parm7")
        dataframe = DataClass("Production")
        for i in range(Umbrella.Bins):
            # Labels.file_name([f"pull_1.{i}.dcd"])
            # Labels.file_name([f"equil_1.{i}.dcd", f"equil_2.{i}.dcd"])
            Labels.file_name([f"prod_1.{i}.dcd", f"prod_2.{i}.dcd", f"prod_3.{i}.dcd", f"prod_4.{i}.dcd", f"prod_5.{i}.dcd", f"prod_6.{i}.dcd", f"prod_7.{i}.dcd", f"prod_8.{i}.dcd"])
            # Labels.file_name([f"prod-NoRI_1.{i}.dcd", f"prod-NoRI_2.{i}.dcd", f"prod-NoRI_3.{i}.dcd", f"prod-NoRI_4.{i}.dcd", f"prod-NoRI_5.{i}.dcd", f"prod-NoRI_6.{i}.dcd", f"prod-NoRI_7.{i}.dcd",f"prod-NoRI_8.{i}.dcd",])
            # Labels.file_name(f"prod_PBE_1_{i}.nc")
            Bonds, Dihedrals = Anal.Label_Maker(Labels, f"{Job.WorkDir}{i}/label_maker.tcl")
            if args.DryRun == "False":
                subprocess.run([f"cd {Job.WorkDir}{i} ; vmd -dispdev text -e label_maker.tcl ; cd ../"], shell=True, capture_output=True)
            dataframe = Anal.Labal_Analysis(Bonds, Dihedrals, f"{Job.WorkDir}{i}/", i, dataframe)
            core_load.append(f"mol addfile ./{i}/equil.{i}.restart.coor ")
        with open(f"{Job.WorkDir}prod_load.tcl",'w') as f:
            for i in range(len(core_load)):
                print(core_load[i], file=f)
        # print(dataframe.dat)
        df = pd.concat(dataframe.dat)
        print(df.loc[df["Name"] == "c6Ring", "Data"].shape)
        # reactioncoordinate = "C2-C7"
        print(type(Umbrella.atom1), Umbrella.atom2, Umbrella.atom3)
        for bond in Bonds:
            if ((bond.at1 == Umbrella.atom1) or (bond.at2 == Umbrella.atom1)) and ((bond.at1 == Umbrella.atom2) or (bond.at2 == Umbrella.atom2)) and (Umbrella.atom3 == 0):
                reactioncoordinate = bond.name
                print(f"Reaction coordinate is {reactioncoordinate}")
            # elif bond.at1 == Umbrella.atom1 and bond.at2 == Umbrella.atom2:
            #     reactioncoordinate = bond.name
            #     print(f"Reaction coordinate is {reactioncoordinate}")
            # elif bond.at1 == Umbrella.atom2 and bond.at2 == Umbrella.atom1:#and Umbrella.atom3 == 0:
            #     reactioncoordinate = bond.name
            #     print(f"Reaction coordinate is {reactioncoordinate}")
            # print(type(bond.at1), bond.at2)
        for dihed in Dihedrals:
            if ((dihed.at1 == Umbrella.atom1) or (dihed.at4 == Umbrella.atom1)) and ((dihed.at2 == Umbrella.atom2) or (dihed.at3 == Umbrella.atom2)) and ((dihed.at3 == Umbrella.atom3) or (dihed.at2 == Umbrella.atom3)) and ((dihed.at4 == Umbrella.atom4) or (dihed.at1 == Umbrella.atom4)):
                reactioncoordinate = dihed.name
                print(f"Reaction coordinate is {reactioncoordinate}")
            # print(dihed.at1, dihed.at2, dihed.at3, dihed.at4)
            # print(Umbrella.atom1, Umbrella.atom2, Umbrella.atom3, Umbrella.atom4)
        try:
            os.mkdir(f"{Job.WorkDir}Figures/")
        except:
            if Job.Verbosity >= 1:
                print("Figures directory already exists")
            else:
                pass
        for bond in Bonds:
            plt.hist(df.loc[df["Name"] == bond.name, "Data",], 100,color="black")
            plt.title(f"{bond.name} bond")
            plt.xlabel("Distance")
            plt.ylabel("Count")
            plt.savefig(f"{Job.WorkDir}Figures/{bond.name}.eps",transparent=True)
            if args.Verbosity >= 1:
                plt.show()
            else:
                plt.clf()
            plt.hist2d(df.loc[df["Name"] == reactioncoordinate, "Data"], df.loc[df["Name"] == bond.name, "Data"],
                       100, cmap="binary")
            # for i in Umbrella.BinVals:
            #     plt.axvline(i, color="red", linestyle="dashed", linewidth=0.2)
            plt.title(f"Reaction coordinate vs. {bond.name} bond")
            plt.ylabel(f"{bond.name} bond distance")
            plt.xlabel("Reaction coordinate")
            plt.savefig(f"{Job.WorkDir}Figures/{bond.name}_2d.eps",transparent=True)
            if args.Verbosity >= 1:
                plt.show()
            else:
                plt.clf()
        for dihed in Dihedrals:
            plt.hist(df.loc[df["Name"] == dihed.name, "Data"], 100, color="black")
            plt.title(f"{dihed.name} dihedral")
            plt.xlabel("Angle")
            plt.ylabel("Count")
            plt.savefig(f"{Job.WorkDir}Figures/{dihed.name}.eps",transparent=True)
            if args.Verbosity >= 1:
                plt.show()
            else:
                plt.clf()
            plt.hist2d(df.loc[df["Name"] == reactioncoordinate, "Data"], df.loc[df["Name"] == dihed.name, "Data"],100, cmap="binary")
            # for i in Umbrella.BinVals:
            #     plt.axvline(i, color="red", linestyle="dashed", linewidth=0.2)
            plt.title(f"Reaction coordinate vs. {dihed.name} dihedral")
            plt.ylabel(f"{dihed.name} dihedral angle")
            plt.xlabel("Reaction coordinate")
            plt.savefig(f"{Job.WorkDir}Figures/{dihed.name}_2d.eps",transparent=True)
            if args.Verbosity >= 1:
                plt.show()
            else:
                plt.clf()

        df.to_csv(f"{Job.WorkDir}Figures/Data.csv")
    if Job.Stage == "convergence":
        if Calc.MaxSteps == 0:
            NumJobs = 1
        else:
            if "prod" in args.WhamFile.casefold():
                steps = prod_length
            elif "equil" in args.WhamFile.casefold():
                steps = equil_length
            NumJobs = math.ceil(steps / Calc.MaxSteps)
        if Job.Verbosity >= 1:
            print(f"Number of steps to glue together is {NumJobs}")
        try:
            os.mkdir(f"{Job.WorkDir}WHAM/Conv")
        except:
            if Job.Verbosity >= 1:
                print("Figures directory already exists")
            else:
                pass
        for i in range(1,NumJobs+1):
            Anal.glue_stick(Umbrella, Job, NumJobs=i, file=args.WhamFile)
            if Umbrella.atom3 != 0:
                periodicity = "Periodic"
            else:
                periodicity = "discrete"
            wham = WhamClass(args.WhamFile, Umbrella.ConstForce, periodicity)
            Wham.Init_Wham(Job, Umbrella, wham)
            Wham.Run_Wham(Job, Umbrella)
            subprocess.run(f"head -n {Umbrella.Bins} {Job.WorkDir}WHAM/out.pmf > {Job.WorkDir}WHAM/Conv/{i}.pmf", shell=True)
            subprocess.run(f"mv {Job.WorkDir}WHAM/PMF.eps {Job.WorkDir}WHAM/Conv/{i}.eps", shell=True)
            subprocess.run(f"sed -i \"0,/+\/-/s/+\/-/Err1/\" {Job.WorkDir}WHAM/Conv/{i}.pmf " , shell=True)
            subprocess.run(
                f"sed -i \"0,/+\/-/s/+\/-/Err2/\" {Job.WorkDir}WHAM/Conv/{i}.pmf ",
                shell=True)
            subprocess.run(
                f"sed -i \"s/#Coor/Coor/g\" {Job.WorkDir}WHAM/Conv/{i}.pmf",
                shell=True)
    if Job.Stage == "vis":
        VisInit(Job, Umbrella, args.WhamFile)
        if args.DryRun == "False":
            VisLoad(args.WhamFile)
    if Job.Stage == "test":
        file = FileGen.ORCA_Wrapper(QM, Calc)
        utils.file_write(f"{Job.WorkDir}runORCA.py", [file])
        SLURM.set_arrayJob("equil_1.txt", 54)
        utils.mpi_gen(SLURM,Umbrella)

def min_setup(MM, Calc, Job, startfile):
    MM.Set_Ensemble("min")
    MM.Set_Outputs(1000, 100, 0)
    MM.Set_Length(10000)
    Calc.Job_Name("min")
    NAMD = NAMDClass(Calc, MM)
    NAMD.set_pme("on")
    if startfile == "start.rst7":
        NAMD.set_startcoords(None, ambercoor="start.rst7", parm="complex.parm7")
    else:
        NAMD.set_startcoords(bincoor=startfile, ambercoor="start.rst7", parm="complex.parm7")
    NAMD.set_cellvectors(MM.CellVec)
    file = FileGen.Namd_File(NAMD)
    utils.file_write(f"{Job.WorkDir}min.conf", [file])

def min_run(MM, Job, GPU=True):
    if Job.Verbosity >= 1:
        print("Running the minimisation script without checking!")
    print("Running the minimisation calculation")
    if GPU == False:
        subprocess.run([f"{MM.NamdPath} ++autoProvision +setcpuaffinity min.conf > min_1.0.out"],shell = True,
                       capture_output = True)
    else:
        if Job.Verbosity >= 1:
            print("Running locally on GPU!")
        subprocess.run([
                           MM.GPUNamd + f" +oneWthPerCore +setcpuaffinity +devices 0 min.conf > min_1.0.out"],
                       shell=True, capture_output=True)
    if Job.Verbosity >= 1:
        print("Cleaning up directory!")
    output = subprocess.run(["mv min* ./setup", "cp ./setup/min_1.0.restart.coor .",
                    "cp ./setup/min_1.0.out .", "cp ./setup/min.conf ."],
                   shell=True, capture_output=True)
    if Job.Verbosity >= 1:
        print(output)
    output = subprocess.run(["cp ./setup/min_1.0.restart.coor .", ],
                   shell=True, capture_output=True)
    if Job.Verbosity >= 1:
        print(output)
    output = subprocess.run(["cp ./setup/min_1.0.out ."],
                   shell=True, capture_output=True)
    if Job.Verbosity >= 1:
        print(output)
    output = subprocess.run(["cp ./setup/min.conf ."],
                   shell=True, capture_output=True)
    if Job.Verbosity >= 1:
        print(output)
def heat_setup(MM, Calc, Job):
    MM.Set_Ensemble("heat")
    MM.Set_Outputs(200, 10, 100)
    MM.Set_Length(10000)
    Calc.Job_Name("heat")
    NAMD = NAMDClass(Calc, MM)
    NAMD.set_pme("on")
    NAMD.set_startcoords("min_1.0.restart.coor", ambercoor="start.rst7", parm="complex.parm7")
    NAMD.set_cellvectors(MM.CellVec)
    file = FileGen.Namd_File(NAMD)
    utils.file_write(f"{Job.WorkDir}heat.conf", [file])

def heat_run(MM, Job,GPU=True):
    if Job.Verbosity >= 1:
        print("Running the Heat script without checking!")
    print("Running the Heating Calculation")
    if GPU == False:
        subprocess.run([MM.NamdPath + f" +autoProvision +setcpuaffinity heat.conf > heat_1.0.out"],shell = True,
                       capture_output = True)
    else:
        if Job.Verbosity >= 1:
            print("Running locally on GPU!")
        subprocess.run([
                           MM.GPUNamd + f" +oneWthPerCore +setcpuaffinity +devices 0 heat.conf > heat_1.0.out"],
                       shell=True, capture_output=True)
    subprocess.run(["mv heat* ./setup", ],
                   shell=True, capture_output=True)
    subprocess.run(["cp ./setup/heat_1.0.restart.coor .", ],
                   shell=True, capture_output=True)
    subprocess.run(["cp ./setup/heat_1.0.out .", ],
                   shell=True, capture_output=True)
    subprocess.run(["cp ./setup/heat.conf ."],
                   shell=True, capture_output=True)

def pull_setup(Umbrella, MM, QM, Calc, Job, DryRun="False"):
    if Job.Verbosity >= 1:
        print("Setting up pulls")
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
       # print(Umbrella.BinVals[i] - Umbrella.Start)
       if abs(Umbrella.BinVals[i] - Umbrella.Start) < abs(Umbrella.Width * 0.5):
           Umbrella.add_start(i)
       # if Umbrella.Width > 0:
       #      if Umbrella.BinVals[i] == Umbrella.Start:
       #          Umbrella.add_start(i)
       # else:
       #      if Umbrella.BinVals[i] == Umbrella.Start:
       #          Umbrella.add_start(i)
    for i in range(Umbrella.Bins):
        # if Umbrella.Width > 0:
        #     if Umbrella.BinVals[i] > Umbrella.Start:
        #         prevPull = f"../{i-1}/pull_1.{i-1}.restart.coor"
        #     elif Umbrella.BinVals[i] < Umbrella.Start:
        #         prevPull = f"../{i + 1}/pull_1.{i + 1}.restart.coor"
        #     elif Umbrella.BinVals[i] == Umbrella.Start:
        #         # Umbrella.add_start(i)
        #         prevPull = f"../heat_1.0.restart.coor"
        #         if Job.Verbosity >= 1:
        #             print(f"Pull starts from directory {i}")
        # else:
        #     if Umbrella.BinVals[i] < Umbrella.Start:
        #         prevPull = f"../{i-1}/pull_1.{i-1}.restart.coor"
        #     elif Umbrella.BinVals[i] > Umbrella.Start:
        #         prevPull = f"../{i + 1}/pull_1.{i + 1}.restart.coor"
        #     elif Umbrella.BinVals[i] == Umbrella.Start:
        #         # Umbrella.add_start(i)
        #         prevPull = f"../heat_1.0.restart.coor"
        #         if Job.Verbosity >= 1:
        #             print(f"Pull starts from directory {i}")
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
        utils.file_write(f"{Job.WorkDir}{i}/pull.conf", [file])
        colvarfile = utils.colvar_gen(Umbrella, i, "pull", Umbrella.PullForce )
        utils.file_write(f"{Job.WorkDir}{i}/colvars.pull.conf", [colvarfile])

        Joblist[i] = f"mkdir /dev/shm/NAMD_{i} ; cd ./{i} ; {MM.NamdPath} pull.conf > pull_1.{i}.out ; cd ../ ; rm -r /dev/shm/NAMD_{i}"
    utils.file_write(f"{Job.WorkDir}pull.txt",Joblist)
    make_runfile(Job, Umbrella, Joblist)
# def run_setup(args, Umbrella, Calc, Job, MM, QM, DryRun=False):
#     if Job.Verbosity >= 1:
#         print("Setting up pulls")
#     make_umbrellaDirs(args)
#     MM.Set_Ensemble("NVT")
#
#     # Umbrella, pullJobs = generate.Pull_Setup(Umbrella, Calc, Job, MM, QM)
#     # make_runfile(Job, Umbrella, pullJobs)
#     if DryRun == False:
#         run_pullScript()
#     return "Umbrella pull has completed setup"

def run_pullScript(loc="./"):
    print("Running pull command")
    run_out = subprocess.run([f"sh {loc}pull.sh"], shell=True, capture_output=True)
    return run_out

def equil_setup(MM, QM, Job, Calc, Umbrella, PreviousJob):
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
        utils.file_write(f"{Job.WorkDir}{i}/colvars.const.conf",[colvar])
        for j in range(NumJobs):
            place = j * Umbrella.Bins + i
            JobList[
                place] = f"( mkdir /dev/shm/NAMD_{i} ; cd ./{i} ; {MM.NamdPath} {Calc.Name}_{j+1}.conf > {Calc.Name}_{j+1}.{i}.out ; cd ../ ; rm -r /dev/shm/NAMD_{i} ) &"
            NAMD.set_qm(Calc, QM, i)
            if j == 0:
                NAMD.set_startcoords(f"{PreviousJob}.{i}.restart.coor")
                file = FileGen.Namd_File(NAMD, j+1, i)
                utils.file_write(f"{Job.WorkDir}{i}/{Calc.Name}_{j + 1}.conf",
                                 [file])
            else:
                NAMD.set_startcoords(f"{Calc.Name}_{j}.{i}.restart.coor")
                file = FileGen.Namd_File(NAMD, j + 1, i)
                utils.file_write(f"{Job.WorkDir}{i}/{Calc.Name}_{j + 1}.conf",
                                 [file])
    for i in range(NumJobs):
        lines = i * Umbrella.Bins
        utils.file_write(f"{Job.WorkDir}{Calc.Name}_{i + 1}.txt",
                         JobList[lines : lines + Umbrella.Bins])


def equil_run(Job):
    if Job.Verbosity >= 1:
        print("Running Equil command locally, This is not recommended...")
    print("WANING Running the Equilibration calculation... THIS DOESNT WORK CURRENTLY")
    raise Exception("Local equil run is not currently supported.")

    # equilJobs = generate.Prod_Setup(Umbrella, Calc, Job, MM, QM, "pull")
    # with open(f"{Job.WorkDir}equil.sh", 'w') as f:
    #     for i in range(Umbrella.Bins):
    #         print(equilJobs[i], file=f)

def make_umbrellaDirs(Umbrella, Job):
    if Job.Verbosity >= 1:
        print("Making umbrella directories")
    # generate.make_dirs(args.UmbrellaBins, args.WorkDir
    for i in range(Umbrella.Bins):
        dir_path = str(Job.WorkDir) + str(i)
        if path.exists(dir_path):
            if Job.Verbosity >= 1:
                print(f"{str(i)} exists. Deleting!")
            try:
                os.rmdir(dir_path)
            except OSError:
                if Job.Verbosity >= 1:
                    print(f"{i} directory not empty, deletion failed...")
                pass
        try:
            os.mkdir(dir_path)
        except FileExistsError:
            if Job.Verbosity >= 1:
                print(f"{i} directory exists, Skipping making new directory")
            pass

# def setup_pulls(Umbrella, Calc, Job, MM, QM):
#     if Job.Verbosity >= 1:
#         print("Setting up pulls")
#     Umbrella, pullJobs = generate.Pull_Setup(Umbrella, Calc, Job, MM, QM)
#     return Umbrella, pullJobs

def make_runfile(Job, Umbrella, pullJobs):
    if Job.Verbosity >= 1:
        print("Generating pull.sh script.")
    with open(f"{Job.WorkDir}pull.sh", 'w') as f:
        print("#!/bin/bash", file=f)
        for i in range(Umbrella.StartBin, Umbrella.Bins):
            print(pullJobs[i], file=f)
        for i in range(0, Umbrella.StartBin):
            print(pullJobs[Umbrella.StartBin - i - 1], file=f)

def VisInit(Job, Umbrella, File, extension="restart.coor"):
    Lines = [None] * (Umbrella.Bins + 1)
    Lines[0] = f"mol new complex.parm7"
    for i in range(Umbrella.Bins):
        Lines[i+1] = f"mol addfile ./{i}/{File}.{i}.{extension}"
    utils.file_write(f"{Job.WorkDir}/{File}_load.tcl", lines=Lines)

def VisLoad(File):
    subprocess.run([f"vmd -e {File}_load.tcl"], shell=True, capture_output=True)

