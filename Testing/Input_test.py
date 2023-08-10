import pickle
import pyBrellaSampling.InputParser as input

Test_Dir = "./Testing/TestFiles/Input/"

with open(f"{Test_Dir}UmbrellaDefaults.pickle", "rb") as f:
    # pickle.dump(DefInput,f)
    UmbrellaDefaults = pickle.load(f)
with open(f"{Test_Dir}QMDefaults.pickle", "rb") as f:
    # pickle.dump(DefInput,f)
    QMDefaults = pickle.load(f)
with open(f"{Test_Dir}MMDefaults.pickle", "rb") as f:
    # pickle.dump(DefInput,f)
    MMDefaults = pickle.load(f)
with open(f"{Test_Dir}ComputeDefaults.pickle", "rb") as f:
    # pickle.dump(DefInput,f)
    ComputeDefaults = pickle.load(f)
with open(f"{Test_Dir}JobDefaults.pickle", "rb") as f:
    # pickle.dump(DefInput,f)
    JobDefaults = pickle.load(f)

def test_JobInput():
    GenInput = input.JobInput(f"{Test_Dir}Job.example")
    with open(f"{Test_Dir}Job.pickle", "rb") as f:
        # pickle.dump(GenInput, f)
        TestInput = pickle.load(f)
    assert GenInput == TestInput, "JobInput doesnt work"
    DefInput = input.JobInput(f"{Test_Dir}Job.conf")
    assert JobDefaults == DefInput, "Job defaults arent working, Have the defaults changed?"

def test_ComputeInput():
    GenInput = input.ComputeInput(f"{Test_Dir}Compute.example")
    with open(f"{Test_Dir}Compute.pickle", "rb") as f:
        # pickle.dump(GenInput, f)
        TestInput = pickle.load(f)
    assert GenInput == TestInput, "ComputeInput doesnt work"
    DefInput = input.ComputeInput(f"{Test_Dir}Compute.conf")
    assert ComputeDefaults == DefInput, "Compute defaults arent working, Have the defaults changed?"

def test_MMInput():
    GenInput = input.MMInput(f"{Test_Dir}MM.example")
    with open(f"{Test_Dir}MM.pickle", "rb") as f:
        # pickle.dump(GenInput, f)
        TestInput = pickle.load(f)
    assert GenInput == TestInput, "MMInput doesnt work"
    DefInput = input.MMInput(f"{Test_Dir}MM.conf")
    assert MMDefaults == DefInput, "MM defaults arent working, Have the defaults changed?"

def test_QMInput():
    GenInput = input.QMInput(f"{Test_Dir}QM.example")
    with open(f"{Test_Dir}QM.pickle", "rb") as f:
        # pickle.dump(GenInput, f)
        TestInput = pickle.load(f)
    assert GenInput == TestInput, "ComputeInput doesnt work"
    DefInput = input.QMInput(f"{Test_Dir}QM.conf")
    assert QMDefaults == DefInput, "QM defaults arent working, Have the defaults changed?"

def test_UmbrellaInput():
    GenInput = input.UmbrellaInput(f"{Test_Dir}Umbrella.example")
    with open(f"{Test_Dir}Umbrella.pickle", "rb") as f:
        # pickle.dump(GenInput, f)
        TestInput = pickle.load(f)
    assert GenInput == TestInput, "ComputeInput doesnt work"
    DefInput = input.UmbrellaInput(f"{Test_Dir}Umbrella.conf")
    assert UmbrellaDefaults == DefInput, "Umbrella defaults arent working, Have the defaults changed?"

def test_ArgParser():
    Arguments = ["-wd=WORKDIR", "-jt=JOBTYPE", "-v=0","-dr=True","-cores=1",
                 "-mem=1", "-MDcpu=NAMDPATHCPU", "-MDgpu=NAMDPATHGPU",
                 "--QmPath=ORCAPATH", "-qsel=ATOMSEL", "-qc=1", "-qspin=1",
                 "-qm=FUNCTIONAL", "-qb=BASIS", "-min=1.0", "-width=1", "-bins=1",
                 "-pf=1", "-sd=1","-mask=1,2,3,4", "-stg=TESTING", "-wf=WHAM"]
    parsedArgs = input.VariableParser(Arguments)
    with open(f"{Test_Dir}Arguments.pickle","rb") as f:
        # pickle.dump(parsedArgs,f)
        TestArgs = pickle.load(f)
    assert parsedArgs == TestArgs, "Argument parser isnt working"
