import pyBrellaSampling.generate as generate
import pyBrellaSampling.wham as Wham
import pyBrellaSampling.analysis as Anal
from pyBrellaSampling.classes import *

import subprocess
import matplotlib.pyplot as plt

def main(args):
    Umbrella = UmbrellaClass(args,args.UmbrellaMin, args.UmbrellaBins,
                             args.StartDistance, args.UmbrellaWidth,)
    Job = JobClass(args)
    Calc = CalcClass(args)
    bins = generate.init_bins(Umbrella.Bins, Umbrella.Width, Umbrella.Min)
    Umbrella.add_bins(bins)
    # UmbrellaCalculation()
    if args.Stage.casefold() == "init" or args.Stage.casefold() == "full":
        generate.QM_Gen(Calc.QMSel, Job.WorkDir)
        if args.DryRun == False:
            if Job.Verbosity >= 1:
                print("Setting up the QM pdb file.")
            logfile = subprocess.run(["vmd", "-dispdev", "text", "-e", "qm_prep.tcl"],
                                     text = True, capture_output = True)
            with open(f"{Job.WorkDir}tcl.log","w") as f:
                print(logfile, file=f)
    if args.Stage.casefold() == "init" or args.Stage.casefold() == "full" or args.Stage.casefold() == "min":
        generate.Min_Setup(Calc, Job)
        if args.DryRun == False:
            if Job.Verbosity >= 1:
                print("Running the minimisation script without checking!")
            print("Running the minimisation calculation")
            # subprocess.run([Calc.NamdPath + f" +p{2*Calc.Threads} min.conf > min.out"],shell = True,
            #                capture_output = True)
            if Job.Verbosity >= 1:
                print("Running locally on GPU!")
            subprocess.run([Calc.GPUNamd + f" +p{2*Calc.Threads} +setcpuaffinity +devices 0 min.conf > min.out"],
                           shell = True, capture_output = True)
            if Job.Verbosity >= 1:
                print("Cleaning up directory!")
            subprocess.run(["mv min.* ./setup", "cp ./setup/min.restart.coor .",
                            "cp ./setup/min.out .", "cp ./setup/min.conf ."],
                           shell = True, capture_output = True)
            subprocess.run(["cp ./setup/min.restart.coor .",],
                           shell=True, capture_output=True)
            subprocess.run(["cp ./setup/min.out ."],
                           shell=True, capture_output=True)
            subprocess.run(["cp ./setup/min.conf ."],
                           shell=True, capture_output=True)
    if args.Stage.casefold() == "init" or args.Stage.casefold() == "full" or args.Stage.casefold() == "heat":
        generate.Heat_Setup(Calc, Job)
        if args.DryRun == False:
            if Job.Verbosity >= 1:
                print("Running the Heat script without checking!")
            print("Running the Heating Calculation")
            # subprocess.run([Calc.NamdPath + f" +p{2*Calc.Threads} heat.conf > heat.out"],shell = True,
            #                capture_output = True)
            if Job.Verbosity >= 1:
                print("Running locally on GPU!")
            subprocess.run([Calc.GPUNamd + f" +p{2*Calc.Threads} +setcpuaffinity +devices 0 heat.conf > heat.out"],
                           shell=True, capture_output=True)
            subprocess.run(["mv heat.* ./setup",],
                           shell=True, capture_output=True)
            subprocess.run(["cp ./setup/heat.restart.coor .",],
                           shell=True, capture_output=True)
            subprocess.run(["cp ./setup/heat.out .",],
                           shell=True, capture_output=True)
            subprocess.run(["cp ./setup/heat.conf ."],
                           shell=True, capture_output=True)
    if args.Stage.casefold() == "setup" or args.Stage.casefold() == "full" or args.Stage.casefold() == "pull":
        run_setup(args, Umbrella, Calc, Job, args.DryRun)
    if args.Stage.casefold() == "equil" or "full":
        Calc.Set_Force(args.ConstForce)
        Calc.Set_Ensemble("NVT")
        Calc.Job_Name("equil")
        Calc.Set_Length(2000, 0.5)  # 2000 steps at 0.5 fs = 1 ps ~ 1 day
        Calc.Set_Outputs(100, 100, 80)  # Timings, Restart, Trajectory
        equilJobs = generate.Prod_Setup(Umbrella, Calc, Job, "pull")
        with open(f"{Job.WorkDir}equil.sh",'w') as f:
            for i in range(Umbrella.Bins):
                print(equilJobs[i], file=f)
        if args.DryRun == False:
            if Job.Verbosity >= 1:
                print("Running Equil command locally, This is not recommended...")
            print("Running the Equilibration calculation...")
            subprocess.run(["sh equil.sh"], shell=True, capture_output=True)
    if args.Stage.casefold() == "prod" or args.Stage.casefold() == "full":
        Calc.Job_Name("prod-NoRI")
        Calc.Set_Length(8000, 0.5)  # 8000 steps at 0.5 fs = 4 ps, ~ 3.5 days
        Calc.Set_Outputs(100, 100, 80)  # Timings, Restart, Trajectory
        equilJobs = generate.Prod_Setup(Umbrella, Calc, Job, "equil")
        with open(f"{Job.WorkDir}prod.sh", 'w') as f:
            for i in range(Umbrella.Bins):
                print(equilJobs[i], file=f)
        if args.DryRun == False:
            if Job.Verbosity >= 1:
                print(
                "Running prod command locally, This is not recommended...")
            print("Running the Production calculation...")
            subprocess.run(["sh prod.sh"], shell=True, capture_output=True)
    if args.Stage.casefold() == "wham" or args.Stage.casefold() == "full":
        wham = WhamClass(args.WhamFile, 300)
        Wham.Init_Wham(Job, Umbrella, wham)
    if args.Stage.casefold() == "analysis" or "full":
        Labels = LabelClass("../complex.parm7")
        Labels.add_bond("13722,13724", "C4-Hup", 1.1) # Atom index
        Labels.add_bond("13722,13723", "C4-Hdown", 1.1) # Atom index
        Labels.add_bond("13730,13715", "C2-C7", 1.6)
        Labels.add_bond("13717,13730", "C3-C7", 2.5)
        Labels.add_bond("13722,13708", "C4-O", 2.8)
        Labels.add_dihedral("13741,13738,13735,13730", "SobVert", 50, "Sob", -60, "Vert")
        Labels.add_dihedral("13761,13756,13730,13731", "C7Methyl", 30, "up", 180, "down")
        Labels.add_dihedral("13761,13756,13717,13718", "C3Methyl", 50, "up", 180, "down")
        Labels.add_dihedral("13715,13730,13728,13725", "c6Ring", -55, "Chair", -20, "Boat")
        core_load = []
        core_load.append("mol new complex.parm7")
        dataframe = DataClass("Production")
        for i in range(Umbrella.Bins):
            Labels.file_name(f"equil.{i}.dcd")
            # Labels.file_name(f"prod_PBE_1_{i}.nc")
            Bonds, Dihedrals = Anal.Label_Maker(Labels, f"{Job.WorkDir}{i}/label_maker.tcl")
            if args.DryRun == False:
                subprocess.run([f"cd {Job.WorkDir}{i} ; vmd -dispdev text -e label_maker.tcl ; cd ../"], shell=True, capture_output=True)
            dataframe = Anal.Labal_Analysis(Bonds, Dihedrals, f"{Job.WorkDir}{i}/", i, dataframe)
            core_load.append(f"mol addfile ./{i}/equil.{i}.restart.coor ")
        with open(f"{Job.WorkDir}prod_load.tcl",'w') as f:
            for i in range(len(core_load)):
                print(core_load[i], file=f)
        # print(dataframe.dat)
        df = pd.concat(dataframe.dat)
        print(df.loc[df["Name"] == "c6Ring", "Data"].shape)
        for bond in Bonds:
            plt.hist(df.loc[df["Name"] == bond.name, "Data"], 100)
            plt.title(f"{bond.name} bond")
            plt.xlabel("Distance")
            plt.ylabel("Count")
            plt.savefig(f"{Job.WorkDir}{bond.name}.eps")
            plt.show()
            plt.hist2d(df.loc[df["Name"] == bond.name, "Data"],
                       df.loc[df["Name"] == "C2-C7", "Data"],100, cmap="binary")
            plt.title(f"Reaction coordinate vs. {bond.name} bond")
            plt.xlabel(f"{bond.name} bond distance")
            plt.ylabel("Reaction coordinate")
            plt.savefig(f"{Job.WorkDir}{bond.name}_2d.eps")
            plt.show()
        for dihed in Dihedrals:
            plt.hist(df.loc[df["Name"] == dihed.name, "Data"], 100)
            plt.title(f"{dihed.name} dihedral")
            plt.xlabel("Angle")
            plt.ylabel("Count")
            plt.savefig(f"{Job.WorkDir}{dihed.name}.eps")
            plt.show()
            plt.hist2d(df.loc[df["Name"] == dihed.name, "Data"],
                       df.loc[df["Name"] == "C2-C7", "Data"], 100, cmap="binary")
            plt.title(f"Reaction coordinate vs. {dihed.name} dihedral")
            plt.xlabel(f"{dihed.name} dihedral angle")
            plt.ylabel("Reaction coordinate")
            plt.savefig(f"{Job.WorkDir}{dihed.name}_2d.eps")
            plt.show()

        df.to_csv(f"{Job.WorkDir}/Data.csv")

def run_setup(args, Umbrella, Calc, Job, DryRun=False):
    if Job.Verbosity >= 1:
        print("Setting up pulls")
    make_umbrellaDirs(args)
    Calc.Set_Ensemble("NVT")
    Umbrella, pullJobs = generate.Pull_Setup(Umbrella, Calc, Job)
    make_runfile(Job, Umbrella, pullJobs)
    if DryRun == False:
        run_pullScript()
    return "Umbrella pull has completed setup"

def make_umbrellaDirs(args):
    if args.Verbosity >= 1:
        print("Making umbrella directories")
    generate.make_dirs(args.UmbrellaBins, args.WorkDir)

def setup_pulls(Umbrella, Calc, Job):
    if Job.Verbosity >= 1:
        print("Setting up pulls")
    Umbrella, pullJobs = generate.Pull_Setup(Umbrella, Calc, Job)
    return Umbrella, pullJobs

def make_runfile(Job, Umbrella, pullJobs):
    if Job.Verbosity >= 1:
        print("Generating pull.sh script.")
    with open(f"{Job.WorkDir}pull.sh", 'w') as f:
        print("#!/bin/bash", file=f)
        for i in range(Umbrella.StartBin, Umbrella.Bins):
            print(pullJobs[i], file=f)
        for i in range(0, Umbrella.StartBin):
            print(pullJobs[Umbrella.StartBin - i - 1], file=f)

def run_pullScript(loc="./"):
    print("Running pull command")
    run_out = subprocess.run([f"sh {loc}pull.sh"], shell=True, capture_output=True)
    return run_out