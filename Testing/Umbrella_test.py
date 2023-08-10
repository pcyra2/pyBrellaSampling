from pyBrellaSampling.classes import UmbrellaClass, CalcClass, JobClass
import pyBrellaSampling.utils as utils
import pyBrellaSampling.Umbrella as setup
import pyBrellaSampling.InputParser as input
import pickle
import os
import shutil
from testfixtures import compare

Test_Dir = "./Testing/TestFiles/Umbrella/"
ClassExamples_Dir = "./Testing/TestFiles/Classes/"

Arguments = [f"-wd={Test_Dir}", "-jt=umbrella", "-v=0","-dr=True","-cores=1",
                 "-mem=1", "-MDcpu=NAMDPATHCPU", "-MDgpu=NAMDPATHGPU",
                 "--QmPath=ORCAPATH", "-qsel=ATOMSEL", "-qc=1", "-qspin=1",
                 "-qm=FUNCTIONAL", "-qb=BASIS", "-min=1.0", "-width=1", "-bins=10",
                 "-pf=1", "-sd=1","-mask=Comma,Separated,Atom,Index", "-stg=full", "-wf=WHAM"]
args = input.VariableParser(Arguments)


def class_load(name):
    with open(f"{ClassExamples_Dir}{name}.pickle", 'rb') as f:
        obj = pickle.load(f)
    return obj

Umbrella = class_load("UmbrellaClass")
Calc = class_load("CalcClass")
Job = class_load("JobClass")

def test_makeumbrelladirs():
    setup.make_umbrellaDirs(args)
    for i in range(args.UmbrellaBins):
        assert os.path.isdir(f"{Test_Dir}{i}"), "Problem with generating paths"
        os.rmdir(f"{Test_Dir}{i}")

def test_SetupPulls():
    Job.WorkDir = Test_Dir
    setup.make_umbrellaDirs(args)
    GenUmbrella, GenPullJobs = setup.setup_pulls(Umbrella, Calc, Job)
    for i in range(args.UmbrellaBins):
        shutil.rmtree(f"{Test_Dir}{i}")
    with open(f"{Test_Dir}Umbrella.pickle", 'rb') as f:
        # pickle.dump(GenUmbrella,f)
        TestUmbrella=pickle.load(f)
    # utils.file_write(f"{Test_Dir}pullJobs.example", GenPullJobs)
    TestJobs = utils.file_read(f"{Test_Dir}pullJobs.example")
    TestJobs = [line.replace("\n","") for line in TestJobs] #When File is re-read in, linebreaks are added.
    assert GenPullJobs == TestJobs, "JobFile generation failed"
    assert compare(GenUmbrella, TestUmbrella) == None, "Updated umbrellaclass failed"
    Job.WorkDir = f"{Test_Dir}example_"
    setup.make_runfile(Job, TestUmbrella, TestJobs)
    test_runscript = utils.file_read(f"{Test_Dir}pullsh.example")
    gen_runscript = utils.file_read(f"{Test_Dir}example_pull.sh")
    os.remove(f"{Test_Dir}example_pull.sh")
    assert test_runscript == gen_runscript, "Run Script generation failed"


def test_PullRun():
    runout = setup.run_pullScript(f"{Test_Dir}")
    assert runout != None, "Cannot find Namd executable. Is the module loaded"

def test_FullRunSetup():
    Job.WorkDir = Test_Dir
    setup.run_setup(args, Umbrella, Calc, Job, DryRun=True)
    for i in range(args.UmbrellaBins):
        shutil.rmtree(f"{Test_Dir}{i}")