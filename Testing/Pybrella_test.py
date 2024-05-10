import pyBrellaSampling.pyBrella as pybrella
from pyBrellaSampling.Tools.classes import *



import pyBrellaSampling.Tools.utils as utils

from unittest.mock import patch
import os
import sys
import shutil
import numpy
from testfixtures import compare


Test_Dir = "./Testing/TestFiles/Umbrella/"

# globals.init(wd=Test_Dir)

def test_setup():
    Arguments = [f"pyBrella",f"--WorkDir={Test_Dir}", "-v=2","-dr=False","-cores=1",
                 "-mem=1", "-MaxCalc=0", "-MDcpu=NAMDPATHCPU", "-MDgpu=NAMDPATHGPU",
                 "--QmPath=ORCAPATH", "-qsel=ATOMSEL", "-qc=1", "-qspin=1",
                 "-qm=FUNCTIONAL", "-qb=BASIS", "-qargs=D3BJ", "-min=1.0", "-width=1", "-bins=10",
                 "-pf=1", "-f=150", "-sd=1","-mask=0,1,2,3", "-stg=setup", "-af=WHAM"]
    with patch.object(sys, "argv", Arguments):
        pybrella.main()
    GenFile = utils.file_read(f"{Test_Dir}Colvar_prep.tcl")
    TestFile = utils.file_read(f"{Test_Dir}Colvar_prep.example")
    assert GenFile == TestFile, "VMD Colvar script not working"
    os.remove(f"{Test_Dir}Colvar_prep.tcl")
    GenFile = utils.file_read(f"{Test_Dir}qm_prep.tcl")
    TestFile = utils.file_read(f"{Test_Dir}qm_prep.example")
    assert GenFile == TestFile, "VMD QM script not working"
    os.remove(f"{Test_Dir}qm_prep.tcl")
    assert os.path.isfile(f"{Test_Dir}setup/tcl-qm.log"), "vmd QM script not working"
    assert os.path.isfile(f"{Test_Dir}setup/tcl-colvar.log"), "vmd colvar script not working"
    shutil.rmtree(f"{Test_Dir}setup")

def test_min():
    # pybrella.min(MM, Calc, StartFile="start.rst7")
    Arguments = [f"pyBrella",f"--WorkDir={Test_Dir}", "-v=0","-dr=True","-cores=1",
                 "-mem=1", "-MaxCalc=0", "-MDcpu=NAMDPATHCPU", "-MDgpu=NAMDPATHGPU",
                 "--QmPath=ORCAPATH", "-qsel=ATOMSEL", "-qc=1", "-qspin=1",
                 "-qm=FUNCTIONAL", "-qb=BASIS", "-qargs=D3BJ", "-min=1.0", "-width=1", "-bins=10",
                 "-pf=1", "-f=150", "-sd=1","-mask=0,1,2,3", "-stg=min", "-af=WHAM"]
    with patch.object(sys, "argv", Arguments):
        pybrella.main()
    GenFile = utils.file_read(f"{Test_Dir}min.conf")
    TestFile = utils.file_read(f"{Test_Dir}min.example")
    assert GenFile == TestFile, "Min file doesnt work"
    os.remove(f"{Test_Dir}min.conf")
    Arguments = [f"pyBrella",f"--WorkDir={Test_Dir}", "-v=0","-dr=True","-cores=1",
                 "-mem=1", "-MaxCalc=0", "-MDcpu=NAMDPATHCPU", "-MDgpu=NAMDPATHGPU",
                 "--QmPath=ORCAPATH", "-qsel=ATOMSEL", "-qc=1", "-qspin=1",
                 "-qm=FUNCTIONAL", "-qb=BASIS", "-qargs=D3BJ", "-min=1.0", "-width=1", "-bins=10",
                 "-pf=1", "-f=150", "-sd=1","-mask=0,1,2,3", "-stg=min", "-af=WHAM", "--StartFile=alt"]
    with patch.object(sys, "argv", Arguments):
        pybrella.main()
    GenFile = utils.file_read(f"{Test_Dir}min.conf")
    TestFile = utils.file_read(f"{Test_Dir}min_alt.example")
    assert GenFile == TestFile, "Failed to change start file..."
    os.remove(f"{Test_Dir}min.conf")
    
def test_heat():
    Arguments = [f"pyBrella",f"--WorkDir={Test_Dir}", "-v=0","-dr=True","-cores=1",
                 "-mem=1", "-MaxCalc=0", "-MDcpu=NAMDPATHCPU", "-MDgpu=NAMDPATHGPU",
                 "--QmPath=ORCAPATH", "-qsel=ATOMSEL", "-qc=1", "-qspin=1",
                 "-qm=FUNCTIONAL", "-qb=BASIS", "-qargs=D3BJ", "-min=1.0", "-width=1", "-bins=10",
                 "-pf=1", "-f=150", "-sd=1","-mask=0,1,2,3", "-stg=heat", "-af=WHAM"]
    with patch.object(sys, "argv", Arguments):
        pybrella.main()
    GenFile = utils.file_read(f"{Test_Dir}heat.conf")
    TestFile = utils.file_read(f"{Test_Dir}heat.example")
    assert GenFile == TestFile, "Heat file doesnt work"
    os.remove(f"{Test_Dir}heat.conf")
    
def test_pull():
    Bins = 10
    Arguments = [f"pyBrella",f"--WorkDir={Test_Dir}", "-v=0","-dr=True","-cores=1",
                 "-mem=1", "-MaxCalc=0", "-MDcpu=NAMDPATHCPU", "-MDgpu=NAMDPATHGPU",
                 "--QmPath=ORCAPATH", "-qsel=ATOMSEL", "-qc=1", "-qspin=1",
                 "-qm=FUNCTIONAL", "-qb=BASIS", "-qargs=D3BJ", "-min=1.0", "-width=1", f"-bins={Bins}",
                 "-pf=1", "-f=150", "-sd=2","-mask=0,1,2,3", "-stg=pull", "-af=WHAM"]
    with patch.object(sys, "argv", Arguments):
        pybrella.main()
    GenFile = utils.file_read(f"{Test_Dir}pull.txt")
    TestFile = utils.file_read(f"{Test_Dir}pullArray.example")
    assert GenFile == TestFile, "JobList not working"
    os.remove(f"{Test_Dir}pull.txt")
    GenFile = utils.file_read(f"{Test_Dir}pull.sh")
    TestFile = utils.file_read(f"{Test_Dir}pullScript.example")
    assert GenFile == TestFile, "Pull.sh not working"
    os.remove(f"{Test_Dir}pull.sh")
    for i in range(Bins):
        assert os.path.isfile(f"{Test_Dir}{i}/colvars.pull.conf"), "Problem with generation of pull files.."
        assert os.path.isfile(f"{Test_Dir}{i}/pull.conf"), "Problem with generation of pull files.."
        if i == 0:
            GenFile = utils.file_read(f"{Test_Dir}{i}/pull.conf")
            TestFile = utils.file_read(f"{Test_Dir}pull{i}.example")
            assert GenFile == TestFile, "Pulling from big to little not working"
        if i == 1:
            GenFile = utils.file_read(f"{Test_Dir}{i}/pull.conf")
            TestFile = utils.file_read(f"{Test_Dir}pull{i}.example")
            assert GenFile == TestFile, "Pulling from heat_1.0 not working"
        if i == 2:
            GenFile = utils.file_read(f"{Test_Dir}{i}/pull.conf")
            TestFile = utils.file_read(f"{Test_Dir}pull{i}.example")
            assert GenFile == TestFile, "Pulling from little to big not working"
        shutil.rmtree(f"{Test_Dir}{i}")

def test_equil():
    Bins = 10
    Arguments = [f"pyBrella",f"--WorkDir={Test_Dir}", "-v=0","-dr=True","-cores=1",
                 "-mem=1", "-MaxCalc=0", "-MDcpu=NAMDPATHCPU", "-MDgpu=NAMDPATHGPU",
                 "--QmPath=ORCAPATH", "-qsel=ATOMSEL", "-qc=1", "-qspin=1",
                 "-qm=FUNCTIONAL", "-qb=BASIS", "-qargs=D3BJ", "-min=1.0", "-width=1", f"-bins={Bins}",
                 "-pf=1", "-f=150", "-sd=2","-mask=0,1,2,3", "-stg=pull", "-af=WHAM"]
    with patch.object(sys, "argv", Arguments):
        pybrella.main()
    Arguments = [f"pyBrella",f"--WorkDir={Test_Dir}", "-v=0","-dr=True","-cores=1",
                 "-mem=1", "-MaxCalc=0", "-MDcpu=NAMDPATHCPU", "-MDgpu=NAMDPATHGPU",
                 "--QmPath=ORCAPATH", "-qsel=ATOMSEL", "-qc=1", "-qspin=1",
                 "-qm=FUNCTIONAL", "-qb=BASIS", "-qargs=D3BJ", "-min=1.0", "-width=1", f"-bins={Bins}",
                 "-pf=1", "-f=150", "-sd=2","-mask=0,1,2,3", "-stg=equil", "-af=WHAM"]
    with patch.object(sys, "argv", Arguments):
        pybrella.main()
    os.remove(f"{Test_Dir}pull.sh")
    os.remove(f"{Test_Dir}pull.txt")
    GenFile = utils.file_read(f"{Test_Dir}equil_1.txt")
    TestFile = utils.file_read(f"{Test_Dir}equilArray.example")
    assert GenFile == TestFile, "Equil Array job not working."
    os.remove(f"{Test_Dir}equil_1.txt")
    GenFile = utils.file_read(f"{Test_Dir}0/equil_1.conf")
    TestFile = utils.file_read(f"{Test_Dir}equil.example")
    assert GenFile == TestFile, "Equil input file not working."
    GenFile = utils.file_read(f"{Test_Dir}0/colvars.pull.conf")
    TestFile = utils.file_read(f"{Test_Dir}colvars.pull.example")
    assert GenFile == TestFile, "pull colvars not working."
    GenFile = utils.file_read(f"{Test_Dir}0/colvars.const.conf")
    TestFile = utils.file_read(f"{Test_Dir}colvars.const.example")
    assert GenFile == TestFile, "const colvars not working."
    for i in range(Bins):
        shutil.rmtree(f"{Test_Dir}{i}")
    
def test_prod():
    Bins = 10
    Arguments = [f"pyBrella",f"--WorkDir={Test_Dir}", "-v=0","-dr=True","-cores=1",
                 "-mem=1", "-MaxCalc=0", "-MDcpu=NAMDPATHCPU", "-MDgpu=NAMDPATHGPU",
                 "--QmPath=ORCAPATH", "-qsel=ATOMSEL", "-qc=1", "-qspin=1",
                 "-qm=FUNCTIONAL", "-qb=BASIS", "-qargs=D3BJ", "-min=1.0", "-width=1", f"-bins={Bins}",
                 "-pf=1", "-f=150", "-sd=2","-mask=0,1,2,3", "-stg=pull", "-af=WHAM"]
    with patch.object(sys, "argv", Arguments):
        pybrella.main()
    Arguments = [f"pyBrella",f"--WorkDir={Test_Dir}", "-v=0","-dr=True","-cores=1",
                 "-mem=1", "-MaxCalc=0", "-MDcpu=NAMDPATHCPU", "-MDgpu=NAMDPATHGPU",
                 "--QmPath=ORCAPATH", "-qsel=ATOMSEL", "-qc=1", "-qspin=1",
                 "-qm=FUNCTIONAL", "-qb=BASIS", "-qargs=D3BJ", "-min=1.0", "-width=1", f"-bins={Bins}",
                 "-pf=1", "-f=150", "-sd=2","-mask=0,1,2,3", "-stg=prod", "-af=WHAM"]
    with patch.object(sys, "argv", Arguments):
        pybrella.main()
    os.remove(f"{Test_Dir}pull.sh")
    os.remove(f"{Test_Dir}pull.txt")
    GenFile = utils.file_read(f"{Test_Dir}prod_1.txt")
    TestFile = utils.file_read(f"{Test_Dir}prodArray.example")
    assert GenFile == TestFile, "Equil Array job not working."
    GenFile = utils.file_read(f"{Test_Dir}0/prod_1.conf")
    TestFile = utils.file_read(f"{Test_Dir}prod.example")
    assert GenFile == TestFile, "Equil input file not working."
    for i in range(Bins):
        shutil.rmtree(f"{Test_Dir}{i}")
        
def test_wham():
    Test_Dir = "./Testing/TestFiles/Wham/"
    Bins=54
    Arguments = [f"pyBrella",f"--WorkDir={Test_Dir}", "-v=0","-dr=False","-cores=1",
                 "-mem=1", "-MaxCalc=0", "-MDcpu=NAMDPATHCPU", "-MDgpu=NAMDPATHGPU",
                 "--QmPath=ORCAPATH", "-qsel=ATOMSEL", "-qc=1", "-qspin=1",
                 "-qm=FUNCTIONAL", "-qb=BASIS", "-qargs=D3BJ", "-min=1.3", "-width=0.05", f"-bins={Bins}",
                 "-pf=1", "-f=300", "-sd=2","-mask=0,1,2,3", "-stg=wham", "-af=prod_1"]
    with patch.object(sys, "argv", Arguments):
        pybrella.main()
    TestFile = utils.file_read(f"{Test_Dir}prod_1.example")
    GenFile = utils.file_read(f"{Test_Dir}WHAM/plot_free_energy.dat")
    TestEn = numpy.zeros(Bins)
    GenEn = numpy.zeros(Bins)
    for i in range(len(GenFile)):
        if "#" in GenFile[i]: # Break after the # window line
            continue
        words = GenFile[i].split()
        GenEn[i] = float(words[1])
    for i in range(len(TestFile)):
        if "#" in TestFile[i]: # Break after the # window line
            continue
        words = TestFile[i].split()
        TestEn[i] = float(words[1])
    assert np.array_equiv(GenEn, TestEn), "PMF energies are different to the example.."
    file_list = ["PMF.dat", "PMF.eps", "prod_1metadata.dat", "prod_1UmbrellaHist.dat", "wham.sh", "out.pmf"]
    for i in file_list:
        assert os.path.isfile(f"{Test_Dir}/WHAM/{i}"), f"{i} file not generated by the wham script.."
    shutil.rmtree(f"{Test_Dir}WHAM")
    
def test_convergence():
    Test_Dir = "./Testing/TestFiles/Wham/"
    Bins=54
    Arguments = [f"pyBrella",f"--WorkDir={Test_Dir}", "-v=0","-dr=False","-cores=1",
                 "-mem=1", "-MaxCalc=4000", "-MDcpu=NAMDPATHCPU", "-MDgpu=NAMDPATHGPU",
                 "--QmPath=ORCAPATH", "-qsel=ATOMSEL", "-qc=1", "-qspin=1",
                 "-qm=FUNCTIONAL", "-qb=BASIS", "-qargs=D3BJ", "-min=1.3", "-width=0.05", f"-bins={Bins}",
                 "-pf=1", "-f=300", "-sd=2","-mask=0,1,2,3", "-stg=convergence", "-af=prod"]
    with patch.object(sys, "argv", Arguments):
        pybrella.main()
    for i in range(1,2):
        assert os.path.isfile(f"{Test_Dir}WHAM/Conv/{i}.eps"), f"Convergence file {i}.eps not generated"
        assert os.path.isfile(f"{Test_Dir}WHAM/Conv/{i}.pmf"), f"Convergence file {i}.pmf not generated"
    shutil.rmtree(f"{Test_Dir}WHAM")

def test_vis():
    Arguments = [f"pyBrella",f"--WorkDir={Test_Dir}", "-v=0","-dr=True","-cores=1",
                 "-mem=1", "-MaxCalc=4000", "-MDcpu=NAMDPATHCPU", "-MDgpu=NAMDPATHGPU",
                 "--QmPath=ORCAPATH", "-qsel=ATOMSEL", "-qc=1", "-qspin=1",
                 "-qm=FUNCTIONAL", "-qb=BASIS", "-qargs=D3BJ", "-min=1.3", "-width=0.05", f"-bins=10",
                 "-pf=1", "-f=300", "-sd=2","-mask=0,1,2,3", "-stg=vis", "-af=pull_1"]
    with patch.object(sys, "argv", Arguments):
        pybrella.main()
    GenFile=utils.file_read(f"{Test_Dir}pull_1_load.tcl")
    TestFile=utils.file_read(f"{Test_Dir}pull_1_load.example")
    assert GenFile == TestFile, "Visualisation file not working."