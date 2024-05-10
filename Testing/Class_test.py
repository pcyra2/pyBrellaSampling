from pyBrellaSampling.Tools.classes import *
import pyBrellaSampling.Tools.InputParser as input
from pyBrellaSampling.Tools.globals import *

import pickle as pickle
from testfixtures import compare
# import argparse as ap

Test_Dir = "./Testing/TestFiles/Classes/"
init(wd=Test_Dir) # inits global vars.

# with open("./Testing/TestFiles/Input/Arguments.pickle", 'rb') as f:
#     args = pickle.load(f)
arguments = [f"-wd=./Testing/TestFiles/Input/", "-v=0","-dr=True","-cores=1",
                 "-mem=1", "-MaxCalc=0", "-MDcpu=NAMDPATHCPU", "-MDgpu=NAMDPATHGPU",
                 "--QmPath=ORCAPATH", "-qsel=ATOMSEL", "-qc=1", "-qspin=1",
                 "-qm=FUNCTIONAL", "-qb=BASIS", "-qargs=D3BJ", "-min=1.0", "-width=1", "-bins=10",
                 "-pf=1","-f=150", "-sd=1","-mask=1,2,3,4", "-stg=full", "-af=WHAM"]
args = input.UmbrellaInput(arguments)

test_array = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]



def test_BondClass():
    GenBond = BondClass(1, 2, "Bond", 1.2)
    # with open(f"{Test_Dir}BondClass.pickle", 'wb') as f:
    #     pickle.dump(GenBond,f)
    with open(f"{Test_Dir}BondClass.pickle", "rb") as f:
        TestBond = pickle.load(f)
    assert compare(GenBond, TestBond) == None, "Bond class is incorrect."

def test_DihedralClass():
    GenDihedral = DihedralClass(1, 2, 3, 4,
                                         "Dihedral", 0,
                                         "Min", 180,"Max")
    # with open(f"{Test_Dir}DihedralClass.pickle", 'wb') as f:
    #     pickle.dump(GenDihedral,f)
    with open(f"{Test_Dir}DihedralClass.pickle", "rb") as f:
        TestDihedral = pickle.load(f)
    assert compare(GenDihedral, TestDihedral) == None, "Dihedral class is incorrect."

def test_LabelClass():
    GenLabelClass = LabelClass("param.file")
    GenLabelClass.file_name("coord.file")
    GenLabelClass.add_bond("1,2", "testbond", 1.2)
    GenLabelClass.add_dihedral("1,2,3,4", "testdihedral", 0,
                            "Min", 180, "Max")
    # with open(f"{Test_Dir}LabelClass.pickle", 'wb') as f:
    #     pickle.dump(GenLabelClass,f)
    with open(f"{Test_Dir}LabelClass.pickle", "rb") as f:
        # pickle.dump(GenLabelClass,f)
        TestLabelClass = pickle.load(f)
    assert compare(GenLabelClass, TestLabelClass) == None, "Label class is incorrect."

def test_UmbrellaClass():
    GenUmbrellaClass = UmbrellaClass(args, 1,10,1,1)
    GenUmbrellaClass.add_start(4)
    GenUmbrellaClass.add_bins(test_array)
    # with open(f"{Test_Dir}UmbrellaClass.pickle", 'wb') as f:
    #     pickle.dump(GenUmbrellaClass,f)
    with open(f"{Test_Dir}UmbrellaClass.pickle", "rb") as f:
        # pickle.dump(GenUmbrellaClass,f)
        TestUmbrellaClass = pickle.load(f)
    assert compare(GenUmbrellaClass,TestUmbrellaClass) == None, "Umbrella class is incorrect"

# def test_JobClass():
#     GenJobClass = JobClass(args)
#     with open(f"{Test_Dir}JobClass.pickle", "rb") as f:
#         # pickle.dump(GenJobClass, f)
#         TestJobClass = pickle.load(f)
#     assert compare(TestJobClass,GenJobClass) == None, "Job class is incorrect"

def test_DataClass():
    GenDataClass = DataClass("DataClass")
    GenDataClass.add_data("TestData", 1, test_array)
    # with open(f"{Test_Dir}DataClass.pickle", 'wb') as f:
    #     pickle.dump(GenDataClass,f)
    with open(f"{Test_Dir}DataClass.pickle","rb") as f:
        TestDataClass = pickle.load(f)
    assert compare(TestDataClass,GenDataClass) == None, "Dataframe class is incorrect"

def test_CalcClass():
    GenCalcClass = CalcClass(args)
    GenCalcClass.Job_Name("TestCalc")
    GenCalcClass.Set_OutFile("TestCalcOutFile")
    GenCalcClass.Set_Id(1)
    GenCalcClass.Set_QM(True)
    # with open(f"{Test_Dir}CalcClass.pickle", 'wb') as f:
    #     pickle.dump(GenCalcClass,f)
    with open(f"{Test_Dir}CalcClass.pickle", 'rb') as f:
        TestCalcClass = pickle.load(f)
    assert compare(TestCalcClass, GenCalcClass) == None, "CalcClass is incorrect"

def test_WhamClass():
    GenWhamClass = WhamClass("Wham", 100)
    # with open(f"{Test_Dir}WhamClass.pickle", 'wb') as f:
    #     pickle.dump(GenWhamClass,f)
    with open(f"{Test_Dir}WhamClass.pickle", "rb") as f:
        TestWhamClass = pickle.load(f)
    assert compare(GenWhamClass, TestWhamClass) == None, "WhamClass is incorrect"


def test_QMClass():
    GenQMClass = QMClass(args)
    # with open(f"{Test_Dir}QMClass.pickle", 'wb') as f:
    #     pickle.dump(GenQMClass,f)
    with open(f"{Test_Dir}QMClass.pickle", "rb") as f:
        TestQMClass = pickle.load(f)
    assert compare(GenQMClass, TestQMClass) == None, "QMClass is incorrect"

def test_MMClass():
    GenMMClass = MMClass()
    GenMMClass.Set_Shake("bonds")
    GenMMClass.Set_Outputs(1,1,1,)
    GenMMClass.Set_Temp(300)
    GenMMClass.Set_Ensemble("NVT")
    GenMMClass.Set_Length(1000, 0.05)
    GenMMClass.Change_Cell(10)
    GenMMClass.Set_Force(args["ConstForce"])
    # with open(f"{Test_Dir}MMClass.pickle", 'wb') as f:
    #     pickle.dump(GenMMClass,f)
    with open(f"{Test_Dir}MMClass.pickle", "rb") as f:
        TestMMClass = pickle.load(f)
    assert compare(GenMMClass, TestMMClass ) == None, "MMClass is incorrect"

def test_NAMDClass():
    with open(f"{Test_Dir}CalcClass.pickle", "rb") as f:
        Calc = pickle.load(f)
    with open(f"{Test_Dir}MMClass.pickle", "rb") as f:
        MM = pickle.load(f)
    with open(f"{Test_Dir}QMClass.pickle", "rb") as f:
        QM = pickle.load(f)
    print(Calc.QM)
    MM.Ensemble = "heat"
    GenNAMDClass = NAMDClass(Calc,MM)
    assert GenNAMDClass.BrensdenPressure == "off", "heat not working"
    MM.Ensemble = "NPT"
    GenNAMDClass = NAMDClass(Calc, MM)
    assert GenNAMDClass.BrensdenPressure == "on", "NPT not working"
    MM.Ensemble = "min"
    Calc.QM = "False"
    GenNAMDClass = NAMDClass(Calc, MM)
    GenNAMDClass.set_qm(Calc, QM,)
    assert GenNAMDClass.qmForces == "off", "turn qm off doesnt work"
    assert GenNAMDClass.runtype == "minimize", "minimize not working"
    Calc.QM = True
    MM.Ensemble = "NVT"
    GenNAMDClass = NAMDClass(Calc, MM)
    GenNAMDClass.set_pme("on")
    GenNAMDClass.set_cellvectors(10)
    GenNAMDClass.set_qm(Calc,QM, 0)
    assert GenNAMDClass.colvarlines == "", "Colvars not turned off by default"
    GenNAMDClass.set_colvars("ColvarFile")
    GenNAMDClass.set_startcoords("bincoor", "ambercoor", "parm")
    # with open(f"{Test_Dir}NAMDClass.pickle", 'wb') as f:
    #     pickle.dump(GenNAMDClass,f)
    with open(f"{Test_Dir}NAMDClass.pickle", "rb") as f:
        # pickle.dump(GenNAMDClass, f)
        TestNAMDClass = pickle.load(f)
    assert compare(GenNAMDClass, TestNAMDClass) == None, "NAMDClass is incorrect"

