import sys
import time


import pyBrellaSampling.InputParser as input
import pyBrellaSampling.Umbrella as UmbrellaRun
import pyBrellaSampling.QMMM as QMMMRun


def main():
    starttime = time.time()

    ### Parse user inputs.
    # args = argParser()
    args = input.VariableParser(sys.argv[1:])
    ### Calculation vars.
    # CalcVars = Calculation_Parse(args.JobType)
    # if CalcVars[0] == "Umbrella":
    #     init = CalcVars[1]
    #     setup = CalcVars[2]
    #     run = CalcVars[3]
    #     analysis = CalcVars[4]
    #     wham = CalcVars[5]
    # bins = Generate.init_bins(Umbrella.Bins, Umbrella.Width, Umbrella.Min)
    # Umbrella.add_bins(bins)
    if args.JobType.casefold() == "umbrella":   #JobType is case insensitive
        UmbrellaRun.main(args)
    if args.JobType.casefold() == "QMMM":
        QMMMRun.main(args)
        # UmbrellaCalculation(args, init, setup, run, wham, analysis)
    # if setup == True:
    #     Setup.run_setup(args, Umbrella, Calc, Job, args.DryRun)
        # Generate.make_dirs(args.UmbrellaBins, args.WorkDir)
        # Calc.Set_Ensemble("NVT")
        # Umbrella, pullJobs = Generate.Pull_Setup(Umbrella, Calc, Job)
        # log.info(f"Pull Start Bin = {Umbrella.StartBin}")
        # with open(f"{Job.WorkDir}pull.sh",'w') as f:
        #     print("#!/bin/bash", file=f)
        #     for i in range(Umbrella.StartBin, Umbrella.Bins):
        #         print(pullJobs[i], file=f)
        #     for i in range(0, Umbrella.StartBin):
        #         print(pullJobs[Umbrella.StartBin -i -1],file=f)
        # if args.DryRun == False:
        #     log.warning("Running Pull command")
        #     subprocess.run(["sh pull.sh"], shell=True, capture_output=True)

    # if run == True:
    #     Calc.Set_Force(args.Force)
    #     Calc.Set_Ensemble("NVT")
    #     Calc.Job_Name("equil")
    #     Calc.Set_Length(2000, 0.5)  # 2000 steps at 0.5 fs = 1 ps ~ 1 day
    #     Calc.Set_Outputs(100, 100, 80)  # Timings, Restart, Trajectory
    #     equilJobs = Generate.Prod_Setup(Umbrella, Calc, Job, "pull")
    #     with open(f"{Job.WorkDir}equil.sh",'w') as f:
    #         for i in range(Umbrella.Bins):
    #             print(equilJobs[i], file=f)
    #     if args.DryRun == False:
    #         log.warning("Running Equil command locally, This is not recommended...")
    #         print("Running the Equilibration calculation...")
    #         subprocess.run(["sh equil.sh"], shell=True, capture_output=True)
    #     Calc.Job_Name("prod-NoRI")
    #     Calc.Set_Length(8000, 0.5)  # 8000 steps at 0.5 fs = 4 ps, ~ 3.5 days
    #     Calc.Set_Outputs(100, 100, 80)  # Timings, Restart, Trajectory
    #     equilJobs = Generate.Prod_Setup(Umbrella, Calc, Job, "equil")
    #     with open(f"{Job.WorkDir}prod.sh", 'w') as f:
    #         for i in range(Umbrella.Bins):
    #             print(equilJobs[i], file=f)
    #     if args.DryRun == False:
    #         log.warning(
    #             "Running prod command locally, This is not recommended...")
    #         print("Running the Production calculation...")
    #         subprocess.run(["sh prod.sh"], shell=True, capture_output=True)

    # if wham == True:
    #     wham = WhamClass(args.WhamFile, 300)
    #     Wham.Init_Wham(Job, Umbrella, wham)

    # if analysis == True:
    #     Labels = LabelClass("../complex.parm7")
    #     Labels.add_bond("13722,13724", "C4-Hup", 1.1) # Atom index
    #     Labels.add_bond("13722,13723", "C4-Hdown", 1.1) # Atom index
    #     Labels.add_bond("13730,13715", "C2-C7", 1.6)
    #     Labels.add_bond("13717,13730", "C3-C7", 2.5)
    #     Labels.add_bond("13722,13708", "C4-O", 2.8)
    #     Labels.add_dihedral("13741,13738,13735,13730", "SobVert", 50, "Sob", -60, "Vert")
    #     Labels.add_dihedral("13761,13756,13730,13731", "C7Methyl", 30, "up", 180, "down")
    #     Labels.add_dihedral("13761,13756,13717,13718", "C3Methyl", 50, "up", 180, "down")
    #     Labels.add_dihedral("13715,13730,13728,13725", "c6Ring", -55, "Chair", -20, "Boat")
    #     core_load = []
    #     core_load.append("mol new complex.parm7")
    #     dataframe = DataClass("Production")
    #     for i in range(Umbrella.Bins):
    #         Labels.file_name(f"prod.{i}.dcd")
    #         # Labels.file_name(f"prod_PBE_1_{i}.nc")
    #         Bonds, Dihedrals = Anal.Label_Maker(Labels, f"{Job.WorkDir}{i}/label_maker.tcl")
    #         if args.DryRun == False:
    #             subprocess.run([f"cd {Job.WorkDir}{i} ; vmd -dispdev text -e label_maker.tcl ; cd ../"], shell=True, capture_output=True)
    #         dataframe = Anal.Labal_Analysis(Bonds, Dihedrals, f"{Job.WorkDir}{i}/", i, dataframe)
    #         core_load.append(f"mol addfile ./{i}/prod.{i}.restart.coor ")
    #     with open(f"{Job.WorkDir}prod_load.tcl",'w') as f:
    #         for i in range(len(core_load)):
    #             print(core_load[i], file=f)
    #     # print(dataframe.dat)
    #     df = pd.concat(dataframe.dat)
    #     print(df.loc[df["Name"] == "c6Ring", "Data"].shape)
    #     for bond in Bonds:
    #         plt.hist(df.loc[df["Name"] == bond.name, "Data"], 100)
    #         plt.title(f"{bond.name} bond")
    #         plt.xlabel("Distance")
    #         plt.ylabel("Count")
    #         plt.savefig(f"{Job.WorkDir}{bond.name}.eps")
    #         plt.show()
    #         plt.hist2d(df.loc[df["Name"] == bond.name, "Data"],
    #                    df.loc[df["Name"] == "C2-C7", "Data"],100, cmap="binary")
    #         plt.title(f"Reaction coordinate vs. {bond.name} bond")
    #         plt.xlabel(f"{bond.name} bond distance")
    #         plt.ylabel("Reaction coordinate")
    #         plt.savefig(f"{Job.WorkDir}{bond.name}_2d.eps")
    #         plt.show()
    #     for dihed in Dihedrals:
    #         plt.hist(df.loc[df["Name"] == dihed.name, "Data"], 100)
    #         plt.title(f"{dihed.name} dihedral")
    #         plt.xlabel("Angle")
    #         plt.ylabel("Count")
    #         plt.savefig(f"{Job.WorkDir}{dihed.name}.eps")
    #         plt.show()
    #         plt.hist2d(df.loc[df["Name"] == dihed.name, "Data"],
    #                    df.loc[df["Name"] == "C2-C7", "Data"], 100, cmap="binary")
    #         plt.title(f"Reaction coordinate vs. {dihed.name} dihedral")
    #         plt.xlabel(f"{dihed.name} dihedral angle")
    #         plt.ylabel("Reaction coordinate")
    #         plt.savefig(f"{Job.WorkDir}{dihed.name}_2d.eps")
    #         plt.show()
    #
    #     df.to_csv(f"{Job.WorkDir}/Data.csv")
    endtime = time.time()
    print(f"Total time is {endtime - starttime}")

# def argParser():
#     parser = ap.ArgumentParser(description="Commandline arguments. This method of calculation input is being deprecated. Please do not use.")
#     # parser.add_argument('', '', type=, help=, default=, )
#     parser.add_argument('-i', '--input', type=str, help="Input file containing default parameters for this script. Must be in JSON format...", default="input.in")
#     ### Core Calculation bits
#     parser.add_argument('-wd', '--WorkDir', type=str,
#                         help="Home location for the calculations", default="./", )
#     parser.add_argument('-NPATH', '--NamdPath', type=str,
#                         help="Path to NAMD executable",
#                         default="/gpfs01/software/NAMD_2.13_Linux-x86_64-multicore/namd2",)
#     parser.add_argument('-cores', '--NumCores', type=int,
#                         help="Number of cores per individual calculation", default=10, )
#     parser.add_argument('-jt', '--JobType', type=str,
#                         help="Type of calculation to run", default="NONE", )
#     parser.add_argument('-v', '--verbose', type=int,
#                         help="Verbosity: 0 = none, 1 = info", default=0, )
#     parser.add_argument('-dr', '--DryRun', type=bool,
#                         help="Indicates whether programs are executed or not", default=False, )
#     ### Umbrella Bits
#     parser.add_argument('-min', '--UmbrellaMin', type=float,
#                         help="Minimum Umbrella distance", default=1.30, )
#     parser.add_argument('-width', '--UmbrellaWidth', type=float,
#                         help="Umbrella bin width in Angstroms or degrees", default=0.05, )
#     parser.add_argument('-bins', '--UmbrellaBins', type=int,
#                         help="Number of umbrella bins", default=54, )
#     parser.add_argument('-pf', '--PullForce', type=int,
#                         help="Force for pulls in KCal A-2", default=5000, )
#     parser.add_argument('-f', '--Force', type=int,
#                         help="Force for standard Umbrella runs", default=600, ) ### NAMD uses 1/2 k rather than just k
#     parser.add_argument('-sd', '--StartDistance', type=float, help="Distance of initial simulation", default=1.6, )
#     parser.add_argument('-mask', '--AtomMask', type=str, help="Mask for the restrained atoms.", default="13716,13731,0,0", )
#     ### QMMM Bits
#     parser.add_argument('-qsel', '--QmSelection', type=str,
#                         help="Selection algebra for QM atoms",
#                         default="resname CTN POP MG", )
#     parser.add_argument('-qc', '--QmCharge', type=int,
#                         help="Charge of QM region", default=3, )
#     parser.add_argument('-qspin', '--QmSpin', type=int,
#                         help="Spin of QM region", default=1, )
#     parser.add_argument('-qm', '--QmMethod', type=str,
#                         help="Qm method", default="PBE", )
#     parser.add_argument('-qb', '--QmBasis', type=str,
#                         help="QM basis set", default="6-31G*", )
#     parser.add_argument('-qp', '--QmPath', type=str, help="Path to QM software",
#                         default="/gpfs01/home/pcyra2/Software/ORCA/orca", )
#     ### Wham Bits
#     parser.add_argument('-wf', '--WhamFile', type=str,
#                         help="Name prefix of wham data.(XXX.i.colvars.traj", default="prod")
#     ### Parse commandline arguments
#     args = parser.parse_args()
#     args = Utils.input_parser(args)
#     return args

# def Calculation_Parse(JobType):
#     init = False
#     setup = False
#     run = False
#     analysis = False
#     wham = False
#     if JobType == "Init" or JobType == "Full":
#         init = True
#     if JobType == "Setup" or JobType == "Full":
#         setup = True
#     if JobType == "Run" or JobType == "Full":
#         run = True
#     if JobType == "Analysis" or JobType == "Full":
#         analysis = True
#     if JobType == "Wham" or JobType == "Full":
#         wham = True
#     return ["Umbrella", init, setup, run, analysis, wham]

# def UmbrellaCalculation(args, init=False, setup=False, run=False, wham=False, analysis=False):
    # Umbrella = UmbrellaClass(args,args.UmbrellaMin, args.UmbrellaBins,
    #                          args.StartDistance, args.UmbrellaWidth,)
    # Job = JobClass(args)
    # Calc = CalcClass(args)
    # bins = Generate.init_bins(Umbrella.Bins, Umbrella.Width, Umbrella.Min)
    # Umbrella.add_bins(bins)
    # # UmbrellaCalculation()
    # if init == True:
    #     Generate.QM_Gen(Calc.QMSel, Job.WorkDir)
    #     if args.DryRun == False:
    #         log.warning("Setting up the QM pdb file.")
    #         logfile = subprocess.run(["vmd", "-dispdev", "text", "-e", "qm_prep.tcl"],
    #                                  text = True, capture_output = True)
    #         with open(f"{Job.WorkDir}tcl.log","w") as f:
    #             print(logfile, file=f)
    #     Generate.Min_Setup(Calc, Job)
    #     if args.DryRun == False:
    #         log.warning("Running the minimisation script without checking!")
    #         print("Running the minimisation calculation")
    #         # subprocess.run([Calc.NamdPath + f" +p{2*Calc.Threads} min.conf > min.out"],shell = True,
    #         #                capture_output = True)
    #         log.warning("Running locally on GPU!")
    #         subprocess.run([Calc.GPUNamd + f" +p{2*Calc.Threads} +setcpuaffinity +devices 0 min.conf > min.out"],
    #                        shell = True, capture_output = True)
    #         log.info("Cleaning up directory!")
    #         subprocess.run(["mv min.* ./setup", "cp ./setup/min.restart.coor .",
    #                         "cp ./setup/min.out .", "cp ./setup/min.conf ."],
    #                        shell = True, capture_output = True)
    #         subprocess.run(["cp ./setup/min.restart.coor .",],
    #                        shell=True, capture_output=True)
    #         subprocess.run(["cp ./setup/min.out ."],
    #                        shell=True, capture_output=True)
    #         subprocess.run(["cp ./setup/min.conf ."],
    #                        shell=True, capture_output=True)
    #     Generate.Heat_Setup(Calc, Job)
    #     if args.DryRun == False:
    #         log.warning("Running the Heat script without checking!")
    #         print("Running the Heating Calculation")
    #         # subprocess.run([Calc.NamdPath + f" +p{2*Calc.Threads} heat.conf > heat.out"],shell = True,
    #         #                capture_output = True)
    #         log.warning("Running locally on GPU!")
    #         subprocess.run([Calc.GPUNamd + f" +p{2*Calc.Threads} +setcpuaffinity +devices 0 heat.conf > heat.out"],
    #                        shell=True, capture_output=True)
    #         subprocess.run(["mv heat.* ./setup",],
    #                        shell=True, capture_output=True)
    #         subprocess.run(["cp ./setup/heat.restart.coor .",],
    #                        shell=True, capture_output=True)
    #         subprocess.run(["cp ./setup/heat.out .",],
    #                        shell=True, capture_output=True)
    #         subprocess.run(["cp ./setup/heat.conf ."],
    #                        shell=True, capture_output=True)
    # if setup == True:
    #     UmbrellaRun.run_setup(args, Umbrella, Calc, Job, args.DryRun)
    # if run == True:
    #     Calc.Set_Force(args.Force)
    #     Calc.Set_Ensemble("NVT")
    #     Calc.Job_Name("equil")
    #     Calc.Set_Length(2000, 0.5)  # 2000 steps at 0.5 fs = 1 ps ~ 1 day
    #     Calc.Set_Outputs(100, 100, 80)  # Timings, Restart, Trajectory
    #     equilJobs = Generate.Prod_Setup(Umbrella, Calc, Job, "pull")
    #     with open(f"{Job.WorkDir}equil.sh",'w') as f:
    #         for i in range(Umbrella.Bins):
    #             print(equilJobs[i], file=f)
    #     if args.DryRun == False:
    #         log.warning("Running Equil command locally, This is not recommended...")
    #         print("Running the Equilibration calculation...")
    #         subprocess.run(["sh equil.sh"], shell=True, capture_output=True)
    #     Calc.Job_Name("prod-NoRI")
    #     Calc.Set_Length(8000, 0.5)  # 8000 steps at 0.5 fs = 4 ps, ~ 3.5 days
    #     Calc.Set_Outputs(100, 100, 80)  # Timings, Restart, Trajectory
    #     equilJobs = Generate.Prod_Setup(Umbrella, Calc, Job, "equil")
    #     with open(f"{Job.WorkDir}prod.sh", 'w') as f:
    #         for i in range(Umbrella.Bins):
    #             print(equilJobs[i], file=f)
    #     if args.DryRun == False:
    #         log.warning(
    #             "Running prod command locally, This is not recommended...")
    #         print("Running the Production calculation...")
    #         subprocess.run(["sh prod.sh"], shell=True, capture_output=True)
    # if wham == True:
    #     wham = WhamClass(args.WhamFile, 300)
    #     Wham.Init_Wham(Job, Umbrella, wham)
    # if analysis == True:
    #     Labels = LabelClass("../complex.parm7")
    #     Labels.add_bond("13722,13724", "C4-Hup", 1.1) # Atom index
    #     Labels.add_bond("13722,13723", "C4-Hdown", 1.1) # Atom index
    #     Labels.add_bond("13730,13715", "C2-C7", 1.6)
    #     Labels.add_bond("13717,13730", "C3-C7", 2.5)
    #     Labels.add_bond("13722,13708", "C4-O", 2.8)
    #     Labels.add_dihedral("13741,13738,13735,13730", "SobVert", 50, "Sob", -60, "Vert")
    #     Labels.add_dihedral("13761,13756,13730,13731", "C7Methyl", 30, "up", 180, "down")
    #     Labels.add_dihedral("13761,13756,13717,13718", "C3Methyl", 50, "up", 180, "down")
    #     Labels.add_dihedral("13715,13730,13728,13725", "c6Ring", -55, "Chair", -20, "Boat")
    #     core_load = []
    #     core_load.append("mol new complex.parm7")
    #     dataframe = DataClass("Production")
    #     for i in range(Umbrella.Bins):
    #         Labels.file_name(f"equil.{i}.dcd")
    #         # Labels.file_name(f"prod_PBE_1_{i}.nc")
    #         Bonds, Dihedrals = Anal.Label_Maker(Labels, f"{Job.WorkDir}{i}/label_maker.tcl")
    #         if args.DryRun == False:
    #             subprocess.run([f"cd {Job.WorkDir}{i} ; vmd -dispdev text -e label_maker.tcl ; cd ../"], shell=True, capture_output=True)
    #         dataframe = Anal.Labal_Analysis(Bonds, Dihedrals, f"{Job.WorkDir}{i}/", i, dataframe)
    #         core_load.append(f"mol addfile ./{i}/equil.{i}.restart.coor ")
    #     with open(f"{Job.WorkDir}prod_load.tcl",'w') as f:
    #         for i in range(len(core_load)):
    #             print(core_load[i], file=f)
    #     # print(dataframe.dat)
    #     df = pd.concat(dataframe.dat)
    #     print(df.loc[df["Name"] == "c6Ring", "Data"].shape)
    #     for bond in Bonds:
    #         plt.hist(df.loc[df["Name"] == bond.name, "Data"], 100)
    #         plt.title(f"{bond.name} bond")
    #         plt.xlabel("Distance")
    #         plt.ylabel("Count")
    #         plt.savefig(f"{Job.WorkDir}{bond.name}.eps")
    #         plt.show()
    #         plt.hist2d(df.loc[df["Name"] == bond.name, "Data"],
    #                    df.loc[df["Name"] == "C2-C7", "Data"],100, cmap="binary")
    #         plt.title(f"Reaction coordinate vs. {bond.name} bond")
    #         plt.xlabel(f"{bond.name} bond distance")
    #         plt.ylabel("Reaction coordinate")
    #         plt.savefig(f"{Job.WorkDir}{bond.name}_2d.eps")
    #         plt.show()
    #     for dihed in Dihedrals:
    #         plt.hist(df.loc[df["Name"] == dihed.name, "Data"], 100)
    #         plt.title(f"{dihed.name} dihedral")
    #         plt.xlabel("Angle")
    #         plt.ylabel("Count")
    #         plt.savefig(f"{Job.WorkDir}{dihed.name}.eps")
    #         plt.show()
    #         plt.hist2d(df.loc[df["Name"] == dihed.name, "Data"],
    #                    df.loc[df["Name"] == "C2-C7", "Data"], 100, cmap="binary")
    #         plt.title(f"Reaction coordinate vs. {dihed.name} dihedral")
    #         plt.xlabel(f"{dihed.name} dihedral angle")
    #         plt.ylabel("Reaction coordinate")
    #         plt.savefig(f"{Job.WorkDir}{dihed.name}_2d.eps")
    #         plt.show()
    #
    #     df.to_csv(f"{Job.WorkDir}/Data.csv")