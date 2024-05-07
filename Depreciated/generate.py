# import logging as log
# import numpy as np
#
# from pyBrellaSampling.utils import colvar_gen
#
#
# def Namd_File(Calc, parm, coord1, coord2, qmvars, i, heating="off", colfile = "const"):
#     if Calc.Ensemble == "NVP":
#         pressure = "on"
#     else:
#         pressure = "off"
#     if heating == "on":
#         heat_lines = f"""# Heating Control
# temperature     0
# reassignFreq    50
# reassignIncr    0.2
# reassignHold    300
#
# """
#     else:
#         heat_lines = f"""# Temp Control
# langevin            on
# langevinDamping     5
# langevinTemp        {Calc.Temp}
# langevinHydrogen    off
# temperature         {Calc.Temp}
# """
#     try:
#         if Calc.Force > 1000:
#             Rest_Lines = f"""colvars on
# colvarsConfig colvars.{colfile}.conf"""
#         elif Calc.Force < 1000:
#             Rest_Lines = \
# f"""colvars on
# colvarsConfig colvars.{colfile}.conf"""
#     except AttributeError:
#         Rest_Lines = "# No Restraints applied"
#     file = f"""## Umbrella {Calc.OutFile} input file
# # File Options
# parmfile        {parm}
# ambercoor       {coord1}
# bincoordinates  {coord2}
# DCDfile         {Calc.OutFile}.dcd
# DCDfreq         {Calc.TrajOut}
# restartname     {Calc.OutFile}.restart
# restartfreq     {Calc.RestOut}
# outputname      {Calc.OutFile}
# outputTiming    {Calc.TimeOut}
#
#
# # Calculation Options
# amber           on
# switching       off
# exclude         scaled1-4
# 1-4scaling      0.833333333
# scnb            2.0
# readexclusions  yes
# cutoff          {Calc.CutOff}
# watermodel      tip3
# pairListdist        11
# LJcorrection    on
# ZeroMomentum    on
# rigidBonds      {Calc.Shake}
# rigidTolerance  1.0e-8
# rigidIterations 100
# timeStep        {Calc.TimeStep}
#
# fullElectFrequency  1
# nonBondedFreq       1
# stepspercycle       10
#
# # PME Vars
# PME             off
# PMEGridSizeX    300
# PMEGridSizeY    300
# PMEGridSizeZ    300
# PMETolerance    1.0e-6
# PMEInterpOrder  4
#
# cellBasisVector1    {Calc.CellVec} 0.0 0.0
# cellBasisVector2    {(-1/3)*Calc.CellVec} {(2/3)*np.sqrt(2)*Calc.CellVec} 0.0
# cellBasisVector3    {(-1/3)*Calc.CellVec} {(-1/3)*np.sqrt(2)*Calc.CellVec} {(-1/3)*np.sqrt(6)*Calc.CellVec}
# cellOrigin          0 0 0
#
# {heat_lines}
#
# # Pressure Control
# BerendsenPressure   {pressure}
#
# # QMMM settings
# qmForces                {Calc.QM}
# qmParamPDB              "../syst-qm.pdb"
# qmColumn                "beta"
# qmBondColumn            "occ"
# QMSimsPerNode           1
# QMElecEmbed             on
# QMSwitching             on
# QMSwitchingType         shift
# QMPointChargeScheme     round
# QMBondScheme            "cs"
# qmBaseDir               "/dev/shm/NAMD_{i}"
# qmConfigLine            "{qmvars}"
# qmConfigLine            "%%output PrintLevel Mini Print\[ P_Mulliken \] 1 Print\[P_AtCharges_M\] 1 end"
# qmConfigLine            "%PAL NPROCS {Calc.Threads} END"
# qmMult                  "1 {Calc.Spin}"
# qmCharge                "1 {Calc.Charge}"
# qmSoftware              "orca"
# qmExecPath              "{Calc.QMpath}"
# QMOutStride             1
# qmEnergyStride          1
# QMPositionOutStride     1
#
# {Rest_Lines}
#
# run {Calc.Steps}
# """
#     return file
#
#
# # def make_dirs(num, wdir):
# #     for i in range(num):
# #         dir_path = str(wdir) + str(i)
# #         if path.exists(dir_path):
# #             log.info(str(i) + " exists. Deleting!")
# #             try:
# #                 os.rmdir(dir_path)
# #             except OSError:
# #                 log.warning(f"{i} directory not empty, deletion failed...")
# #                 pass
# #         try:
# #             os.mkdir(dir_path)
# #         except FileExistsError:
# #             log.warning(f"{i} directory exists, Skipping making new directory")
# #             pass
#         # pull_file = Namd_File("")
#
# def Min_Setup(Calc, Job):
#     qm_config = f"! {Calc.Method} {Calc.Basis} D3BJ EnGrad TightSCF CFLOAT EasyConv"
#     log.info(f"QM config line is {qm_config}")
#     file = f"""## Minimisation file:
# # File Options
# parmfile        complex.parm7
# ambercoor       start.rst7
# restartname     min.restart
# restartfreq     100
# outputname      min
# outputTiming    1000
#
# # Calculation Options
# amber           on
# switching       off
# exclude         scaled1-4
# 1-4scaling      0.833333333
# scnb            2.0
# readexclusions  yes
# cutoff          8.0
# watermodel      tip3
# pairListdist    11
# LJcorrection    on
# ZeroMomentum    on
#
# fullElectFrequency  1
# nonBondedFreq       1
# stepspercycle       10
#
# # PME Vars
# PME             on
# PMEGridSizeX    300
# PMEGridSizeY    300
# PMEGridSizeZ    300
# PMETolerance    1.0e-6
# PMEInterpOrder  4
#
# cellBasisVector1    {Calc.CellVec} 0.0 0.0
# cellBasisVector2    {(-1/3)*Calc.CellVec} {(2/3)*np.sqrt(2)*Calc.CellVec} 0.0
# cellBasisVector3    {(-1/3)*Calc.CellVec} {(-1/3)*np.sqrt(2)*Calc.CellVec} {(-1/3)*np.sqrt(6)*Calc.CellVec}
# cellOrigin          0 0 0
#
# # Temp Control
# langevin            off
# temperature         0
#
# # Pressure Control
# BerendsenPressure   off
#
#
# minimize                10000
# """
#     with open(f"{Job.WorkDir}min.conf",'w') as f:
#         print(file, file=f)
#
# def Heat_Setup(Calc, Job):
#     qm_config = f"! {Calc.Method} {Calc.Basis} D3BJ EnGrad TightSCF CFLOAT EasyConv"
#     log.info(f"QM config line is {qm_config}")
#     Calc.Set_Length(100000, 2) # 500 steps at 2 fs = 1 ps
#     Calc.Set_Outputs(200,10,100) # Timings, Restart, Trajectory
#     Calc.Set_Shake("all") # restrain all bonds involving Hydrogen
#     Calc.Job_Name("heat")
#     Calc.Set_QM("off")
#     Calc.Set_OutFile(f"{Calc.Name}")
#     Calc.Set_Ensemble("Heating")
#     file = Namd_File(Calc, "complex.parm7", "start.rst7", "min.restart.coor",qm_config, 0, heating="on")
#     with open(f"{Job.WorkDir}heat.conf", 'w') as f:
#         print(file, file=f)
#
# def Pull_Setup(Umbrella, Calc, Job):
#     qm_config = f"! {Calc.Method} {Calc.Basis} D3BJ EnGrad TightSCF CFLOAT NormalSCF"  # Removed EasyConv for convergence assistance
#     log.info(f"QM config line is {qm_config}")
#     Calc.Set_Length(50, 0.5) # 100 steps at 2 fs = 200 fs
#     Calc.Set_Outputs(5,1,10) # Timings, Restart, Trajectory
#     Calc.Set_Shake("all") # restrain all bonds involving Hydrogen
#     Calc.Set_Force(10)
#     Calc.Set_QM("on")
#     if Calc.Ensemble == "NVT":
#         pressure = "off"
#     Joblist = [None]*Umbrella.Bins
#     for i in range(Umbrella.Bins):
#        if Umbrella.Width > 0:
#             if Umbrella.BinVals[i] == Umbrella.Start:
#                 Umbrella.add_start(i)
#        else:
#             if Umbrella.BinVals[i] == Umbrella.Start:
#                 Umbrella.add_start(i)
#     for i in range(Umbrella.Bins):
#         if Umbrella.Width > 0:
#             if Umbrella.BinVals[i] > Umbrella.Start:
#                 prevPull = f"../{i-1}/pull.{i-1}.restart.coor"
#             elif Umbrella.BinVals[i] < Umbrella.Start:
#                 prevPull = f"../{i + 1}/pull.{i + 1}.restart.coor"
#             elif Umbrella.BinVals[i] == Umbrella.Start:
#                 Umbrella.add_start(i)
#                 prevPull = f"../heat.restart.coor"
#         else:
#             if Umbrella.BinVals[i] < Umbrella.Start:
#                 prevPull = f"../{i-1}/pull.{i-1}.restart.coor"
#             elif Umbrella.BinVals[i] > Umbrella.Start:
#                 prevPull = f"../{i + 1}/pull.{i + 1}.restart.coor"
#             elif Umbrella.BinVals[i] == Umbrella.Start:
#                 Umbrella.add_start(i)
#                 prevPull = f"../heat.restart.coor"
#         Calc.Job_Name("pull")
#         Calc.Set_OutFile(f"{Calc.Name}.{i}")
#         file = Namd_File(Calc, "../complex.parm7", "../start.rst7",prevPull, qm_config, i, colfile="pull")
#         with open(f"{Job.WorkDir}{i}/pull.conf", 'w') as f:
#             print(file, file=f)
#         file = colvar_gen(Umbrella, i, "pull", Calc.Force)
#         with open(f"{Job.WorkDir}{i}/colvars.pull.conf", 'w') as f:
#             print(file, file=f)
#         Joblist[i] = f"mkdir /dev/shm/NAMD_{i} ; cd ./{i} ; {Calc.NamdPath} pull.conf > pull.{i}.out ; cd ../ ; rm -r /dev/shm/NAMD_{i}"
#     return Umbrella, Joblist
#
# def Prod_Setup(Umbrella, Calc, Job, startCoord):
#     qm_config = f"! {Calc.Method} {Calc.Basis} D3BJ EnGrad TightSCF CFLOAT SlowConv"
#     log.info(f"QM config line is {qm_config}")
#     Calc.Set_Shake("none")  # restrain all bonds involving Hydrogen
#     Calc.Set_QM("on")
#     if Calc.Ensemble == "NVT":
#         pressure = "off"
#     Joblist = [None] * Umbrella.Bins
#     for i in range(Umbrella.Bins):
#         Calc.Set_OutFile(f"{Calc.Name}.{i}")
#         file = Namd_File(Calc, "../complex.parm7", "../start.rst7", f"{startCoord}.{i}.restart.coor",
#                          qm_config, i)
#         with open(f"{Job.WorkDir}{i}/{Calc.Name}.conf", 'w') as f:
#             print(file, file=f)
#         file = colvar_gen(Umbrella, i, "const", Calc.Force)
#         with open(f"{Job.WorkDir}{i}/colvars.const.conf", 'w') as f:
#             print(file, file=f)
#         Joblist[
#             i] = f"mkdir /dev/shm/NAMD_{i} ; cd ./{i} ; {Calc.NamdPath} {Calc.Name}.conf > {Calc.Name}.{i}.out ; cd ../ ; rm -r /dev/shm/NAMD_{i}"
#     return Joblist