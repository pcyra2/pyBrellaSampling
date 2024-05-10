from pyBrellaSampling.Tools.classes import UmbrellaClass, CalcClass
import pyBrellaSampling.Tools.utils as utils
import pyBrellaSampling.pyBrella as umbrella
import pyBrellaSampling.Tools.InputParser as input
import pickle
import os
import shutil
from pyBrellaSampling.Tools.globals import globals.verbosity, WorkDir, globals.DryRun, globals.parmfile
from testfixtures import compare

Test_Dir = "./Testing/TestFiles/Umbrella/"
ClassExamples_Dir = "./Testing/TestFiles/Classes/"
WorkDir = Test_Dir
globals.verbosity = 0
globals.DryRun = True
globals.parmfile = "complex.parm7"
Arguments = [f"-wd={Test_Dir}", "-v=0","-dr=True","-cores=1",
                 "-mem=1", "-MaxCalc=0", "-MDcpu=NAMDPATHCPU", "-MDgpu=NAMDPATHGPU",
                 "--QmPath=ORCAPATH", "-qsel=ATOMSEL", "-qc=1", "-qspin=1",
                 "-qm=FUNCTIONAL", "-qb=BASIS", "-qargs=D3BJ", "-min=1.0", "-width=1", "-bins=10",
                 "-pf=1", "-f=150", "-sd=1","-mask=Comma,Separated,Atom,Index", "-stg=full", "-af=WHAM"]
args = input.UmbrellaInput(Arguments)

def class_load(name):
    with open(f"{ClassExamples_Dir}{name}.pickle", 'rb') as f:
        obj = pickle.load(f)
    return obj

Umbrella = class_load("UmbrellaClass")
Calc = class_load("CalcClass")
QM = class_load("QMClass")
MM = class_load("MMClass")


# def test_makeumbrelladirs():
#     # for i in range(args.UmbrellaBins):
#     #     os.rmdir(f"{Test_Dir}{i}")
#     umbrella.make_umbrellaDirs(Umbrella,)
#     for i in range(args["UmbrellaBins"]):
#         assert os.path.isdir(f"{Test_Dir}{i}"), "Problem with generating paths"
#         os.rmdir(f"{Test_Dir}{i}")


# def test_pullsetup():
#     umbrella.make_umbrellaDirs(Umbrella, )
#     umbrella.pull_setup(Umbrella, MM, QM, Calc, )
#     GenFile = utils.file_read(f"{Test_Dir}pull.txt")
#     TestFile = utils.file_read(f"{Test_Dir}pullArray.example")
#     assert GenFile == TestFile, "Array file doesnt work"
#     os.remove(f"{Test_Dir}pull.txt")
#     GenFile = utils.file_read(f"{Test_Dir}pull.sh")
#     TestFile = utils.file_read(f"{Test_Dir}pullScript.example")
#     assert GenFile == TestFile, "Script file doesnt work"
#     os.remove(f"{Test_Dir}pull.sh")
#     GenFile = utils.file_read(f"{Test_Dir}0/pull.conf")
#     TestFile = utils.file_read(f"{Test_Dir}pull.example")
#     assert GenFile == TestFile, "pull.conf file doesnt work"
#     GenFile = utils.file_read(f"{Test_Dir}0/colvars.pull.conf")
#     TestFile = utils.file_read(f"{Test_Dir}colvars.pull.example")
#     assert GenFile == TestFile, "colvars.pull.conf file doesnt work"
#     for i in range(args.UmbrellaBins):
#         shutil.rmtree(f"{Test_Dir}{i}")

# def test_minsetup():
#     umbrella.min_setup(MM, Calc, Job, startfile="start.rst7")
#     GenFile = utils.file_read(f"{Test_Dir}min.conf")
#     TestFile = utils.file_read(f"{Test_Dir}min.example")
#     assert GenFile == TestFile, "Min file doesnt work"
#     os.remove(f"{Test_Dir}min.conf")

# def test_heatsetup():
#     umbrella.heat_setup(MM, Calc, Job)
#     GenFile = utils.file_read(f"{Test_Dir}heat.conf")
#     TestFile = utils.file_read(f"{Test_Dir}heat.example")
#     assert GenFile == TestFile, "Heat file doesnt work"
#     os.remove(f"{Test_Dir}heat.conf")

# def test_equilsetup():
#     umbrella.make_umbrellaDirs(Umbrella, Job)
#     Calc.Name = "equil"
#     umbrella.equil_setup(MM, QM, Job, Calc, Umbrella, "prevJob")
#     GenFile = utils.file_read(f"{Test_Dir}equil_1.txt")
#     TestFile = utils.file_read(f"{Test_Dir}equil_1.example")
#     assert GenFile == TestFile, "Equil array file not working"
#     os.remove(f"{Test_Dir}equil_1.txt")
#     GenFile = utils.file_read(f"{Test_Dir}0/equil_1.conf")
#     TestFile = utils.file_read(f"{Test_Dir}equil.example")
#     assert GenFile == TestFile, "Equil conf file not working"
#     GenFile = utils.file_read(f"{Test_Dir}0/colvars.const.conf")
#     TestFile = utils.file_read(f"{Test_Dir}colvars.const.example")
#     assert GenFile == TestFile, "colvars.equil.conf file doesnt work"
#     for i in range(args.UmbrellaBins):
#         shutil.rmtree(f"{Test_Dir}{i}")