import pyBrellaSampling.utils as utils
import pyBrellaSampling.FileGen as FileGen
from pyBrellaSampling.classes import *
import subprocess

def main(args):
    Job = JobClass(args=args)
    if args.DryRun == "False":
        Run = True
    else:
        Run = False
    try:
        InputFileRun(Path=f"{Job.WorkDir}Standalone.inp", args=args, Run=Run)
    except FileNotFoundError:
        Calc = CalcClass(args=args)
        MM = MMClass(args=args)
        Calc.Job_Name(Name=args.Name)
        Calc.Set_OutFile(OutFile=args.Name)
        Calc.Set_Id(0)
        print(vars(Calc))
        MM.Set_Length(Steps=args.Steps, TimeStep=args.TimeStep)
        MM.Set_Ensemble(Ensemble=args.Ensemble)
        MM.Set_Files(parm=args.ParmFile, ambercoor=args.AmberCoordinates)
        MM.Set_Outputs(TimeOut=args.TrajOut, RestOut=args.RestartOut, TrajOut=args.TrajOut) ### Output timings when trajectory is printed.
        print(vars(MM))
        Umbrella = UmbrellaClass(args, Min=args.StartValue, bins=0,
                                 Start=args.StartValue,
                                 Width=(args.EndValue - args.StartValue))
        Umbrella.add_start(0)
        Umbrella.set_force(args.Force)
        bincoor = args.StartFile
        if bincoor == "None":
            bincoor = None
        if args.QM.casefold() == "true":
            QM = QMClass(args=args)
        else:
            QM = False
        NAMD = calc_setup(Job=Job, Calc=Calc, MM=MM, Umbrella=Umbrella, QM=QM, SMD=args.SMD, Run=args.DryRun, bc=bincoor)
        if args.DryRun == "False":
            if QM==False:
                if MM.Ensemble != "min" or MM.Ensemble != "heat":
                    calc_run(Job=Job, MM=MM, Calc=Calc, CUDA="FAST")
                else:
                    calc_run(Job=Job, MM=MM, Calc=Calc)
            else:
                calc_run(Job=Job, MM=MM, Calc=Calc, GPU=False)

def calc_setup(Job, Calc, MM, Umbrella,  QM=False, SMD="off", Run=True, bc="start.rst7"):
    if MM.TimeStep > 1:
        if Job.Verbosity >= 1:
            print("TimeStep is greater than 1 fs. Setting Rattle to True")
        MM.Set_Shake("all")
    NAMD = NAMDClass(Calc=Calc, MM=MM)
    print("#####")
    print(NAMD.amber)
    MM.CellVec = utils.get_cellVec(Job, MM)
    NAMD.set_cellvectors(MM.CellVec)
    NAMD.set_startcoords(bc, ambercoor=MM.ambercoor, parm=MM.parmfile)
    if QM != False:
        if Job.Verbosity >= 1:
            print("Setting up a QMMM calculation")
        NAMD.set_qm(Calc=Calc, QM=QM)
        utils.QM_Gen(QM.QMSel, Job.WorkDir)
        if Run == "True":
            if Job.Verbosity >= 1:
                print("Setting up the QM pdb file.")
            logfile = subprocess.run(["vmd", "-dispdev", "text", "-e", "qm_prep.tcl"],
                                     text = True, capture_output = True)
            with open(f"{Job.WorkDir}tcl-qm.log","w") as f:
                print(logfile, file=f)
    else:
        NAMD.set_pme("on")
        # pass
    if SMD != "off":
        NAMD = init_SMD(Job=Job, NAMD=NAMD,Umbrella=Umbrella )
        utils.ColVarPDB_Gen(Umbrella, Job)
        if Run == "True":
            if Job.Verbosity >= 1:
                print("Setting up the Colvar pdb file.")
            logfile = subprocess.run(["vmd", "-dispdev", "text", "-e", "Colvar_prep.tcl"],
                                     text = True, capture_output = True)
            with open(f"{Job.WorkDir}tcl-colvar.log","w") as f:
                print(logfile, file=f)
    if Job.Verbosity >= 1:
        print("Setting up the conf file")
    file = FileGen.Namd_File(NAMD)
    utils.file_write(f"{Job.WorkDir}{Calc.Name}.conf", [file])
    return NAMD

def calc_run(Job, MM, Calc, GPU=True, CUDA="SLOW"):
    if Job.Verbosity >=1:
        print(f"Running the {Job.Name} Calculation.")
    if CUDA == "FAST":
        subprocess.run([f"sed -i \"s/#CUDAFAST/CUDASOAintegrate on/g\" {Job.WorkDir}{Calc.Name}.conf "],
                       shell=True, capture_output=True)
    if GPU == True:
        print(f"{MM.GPUNamd} +oneWthPerCore +setcpuaffinity +devices 0 +cs {Job.WorkDir}{Calc.Name}.conf > {Job.WorkDir}{Calc.Name}_1.0.out")
        subprocess.run([MM.GPUNamd +
                        f" +oneWthPerCore +setcpuaffinity +devices 0 +cs {Job.WorkDir}{Calc.Name}.conf > {Job.WorkDir}{Calc.Name}_1.0.out"],
                       shell=True, capture_output=True)
    else:
        subprocess.run([MM.CPUNamd +
                        f" +p1 {Job.WorkDir}{Calc.Name}.conf > {Job.WorkDir}{Calc.Name}_1.0.out"],
                       shell=True, capture_output=True)

def init_SMD(Job, NAMD, Umbrella):
    if Umbrella.width != 0:
        Type = "pull"
    else:
        Type = "constant"
    NAMD.set_colvars(file="colvars.conf", toggle="on")
    colvarfile = utils.colvar_gen(Umbrella, i=0,type=Type, force=Umbrella.ConstForce)
    utils.file_write(path=f"{Job.WorkDir}colvars.conf", lines=[colvarfile])
    return NAMD

def InputFileRun(Path,args,Run=True):
    print(utils.MM_DefaultVars.keys())
    CalcStarts = []
    CalcEnds = []
    data = utils.file_read(Path)
    if CalcStarts[0] > 0:
        for i in range(CalcStarts[0]):
            words = data[i].split()
            if "threads" in words[0].casefold():
                args.threads = int(words[1])
    for i in range(len(data)):
        if "$JOB" in data[i]:
            CalcStarts.append(i)
        if "$END" in data[i]:
            CalcEnds.append(i)
    assert len(CalcStarts) == len(CalcEnds), "Job input file is not terminated... "
    for i in range(len(CalcStarts)):
        MMVars = utils.MM_DefaultVars
        for j in range(CalcStarts[i], CalcEnds[i]):
            words = data[j].split()
            if "name" in words[0].casefold():
                JobName = words[1]
            for keys in MMVars.keys():
                if words[0].casefold() == keys.casefold():
                    MMVars[keys] = words[1]
            Calc = CalcClass(args)
            Calc.Set_Shake(MMVars["rigidBonds"])

    return FileNotFoundError



