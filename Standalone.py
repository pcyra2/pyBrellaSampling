import pyBrellaSampling.utils as utils
import pyBrellaSampling.FileGen as FileGen
from pyBrellaSampling.classes import *


def main(args):
    Job = JobClass(args=args)
    Calc = CalcClass(args=args)
    MM = MMClass(args=args)
    QM = QMClass(args=args)
    MM.Set_Length(Steps=args.Steps, TimeStep=args.TimeStep)
    MM.Set_Ensemble(Ensemble=args.Ensemble)
    MM.Set_Files(parm=args.ParmFile, ambercoor=args.AmberCoordinates)
    MM.Set_Outputs(TimeOut=args.TrajOut, RestOut=args.RestartOut, TrajOut=args.TrajOut) ### Output timings when trajectory is printed.
    if MM.Steps > 1:
        if Job.Verbosity >= 1:
            print("TimeStep is greater than 1 fs. Setting Rattle to True")
        MM.Set_Shake("all")
    NAMD = NAMDClass(Calc=Calc, MM=MM, )
    if args.QM.casefold() == "true":
        if Job.Verbosity >= 1:
            print("Setting up a QMMM calculation")
        NAMD.set_qm(Calc=Calc, QM=QM)

    if args.SMD != "off":
        Umbrella = UmbrellaClass(args, Min=args.StartValue, bins=0, Start=args.StartValue, Width=(args.EndValue-args.StartValue))
        Umbrella.add_start(0)
        NAMD = init_SMD(NAMD)
# def calc_setup(MM):
#     MM.

def init_SMD(Job, NAMD, Type, force):
    NAMD.set_colvars(file="colvars.conf", toggle="on")
    colvarfile = utils.colvar_gen(type=Type, force=force)
    utils.file_write(path=f"{Job.WorkDir}colvars.conf", lines=[colvarfile])
    return NAMD