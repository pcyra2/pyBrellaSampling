from pyBrellaSampling.classes import *
import pyBrellaSampling.utils as utils
import pickle as pickle
from testfixtures import compare
import argparse as ap


Test_Dir = "./Testing/TestFiles/Classes/"

test_input = utils.dict_read(f"./Testing/TestFiles/Utils/input.example")
args = ap.Namespace(**test_input)

test_array = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]


def test_BondClass():
    GenBond = BondClass(1, 2, "Bond", 1.2)
    with open(f"{Test_Dir}BondClass.pickle", "rb") as f:
        TestBond = pickle.load(f)
    assert compare(GenBond, TestBond) == None, "Bond class is incorrect."

def test_DihedralClass():
    GenDihedral = DihedralClass(1, 2, 3, 4,
                                         "Dihedral", 0,
                                         "Min", 180,"Max")
    with open(f"{Test_Dir}DihedralClass.pickle", "rb") as f:
        TestDihedral = pickle.load(f)
    assert compare(GenDihedral, TestDihedral) == None, "Bond class is incorrect."

def test_LabelClass():
    GenLabelClass = LabelClass("param.file")
    GenLabelClass.file_name("coord.file")
    GenLabelClass.add_bond("1,2", "testbond", 1.2)
    GenLabelClass.add_dihedral("1,2,3,4", "testdihedral", 0,
                            "Min", 180, "Max")
    with open(f"{Test_Dir}LabelClass.pickle", "rb") as f:
        TestLabelClass = pickle.load(f)
    assert compare(GenLabelClass, TestLabelClass) == None, "Label class is incorrect."

def test_UmbrellaClass():
    GenUmbrellaClass = UmbrellaClass(args, 1,10,1,1)
    GenUmbrellaClass.add_start(4)
    GenUmbrellaClass.add_bins(test_array)
    with open(f"{Test_Dir}UmbrellaClass.pickle", "rb") as f:
        TestUmbrellaClass = pickle.load(f)
    assert compare(GenUmbrellaClass,TestUmbrellaClass) == None, "Umbrella class is incorrect"

def test_JobClass():
    GenJobClass = JobClass(args)
    with open(f"{Test_Dir}JobClass.pickle", "rb") as f:
        TestJobClass = pickle.load(f)
    assert compare(TestJobClass,GenJobClass) == None, "Job class is incorrect"

def test_DataClass():
    GenDataClass = DataClass("DataClass")
    GenDataClass.add_data("TestData", 1, test_array)
    with open(f"{Test_Dir}DataClass.pickle","rb") as f:
        TestDataClass = pickle.load(f)
    assert compare(TestDataClass,GenDataClass) == None, "Dataframe class is incorrect"

def test_CalcClass():
    GenCalcClass = CalcClass(args)
    GenCalcClass.Job_Name("TestCalc")
    GenCalcClass.Set_OutFile("TestCalcOutFile")
    GenCalcClass.Set_Force(1)
    GenCalcClass.Set_Id(1)
    GenCalcClass.Change_Cell(1)
    GenCalcClass.Set_QM(True)
    GenCalcClass.Set_Ensemble("NVT")
    GenCalcClass.Set_Temp(300)
    GenCalcClass.Set_Length(10,0.5)
    GenCalcClass.Set_Outputs(1,1,1)
    GenCalcClass.Set_Shake("none")
    with open(f"{Test_Dir}CalcClass.pickle", 'rb') as f:
        TestCalcClass = pickle.load(f)
    assert compare(TestCalcClass, GenCalcClass) == None, "CalcClass is incorrect"

def test_WhamClass():
    GenWhamClass = WhamClass("Wham", 100)
    with open(f"{Test_Dir}WhamClass.pickle", "rb") as f:
        TestWhamClass = pickle.load(f)
    assert compare(GenWhamClass, TestWhamClass) == None, "WhamClass is incorrect"