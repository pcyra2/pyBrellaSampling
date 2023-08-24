import pyBrellaSampling.analysis as analysis
import pyBrellaSampling.classes
import pyBrellaSampling.pyBrella as pybrella
import pyBrellaSampling.utils
import pyBrellaSampling.utils as utils
import os
import pandas as pd

Test_Dir = "./Testing/TestFiles/Analysis/"


test_bond = pyBrellaSampling.classes.BondClass(1, 2, "testbond",
                                               1.2)
test_dihedral = pyBrellaSampling.classes.DihedralClass(1, 2, 3, 4,
                                       "testdihedral", 0,"Min",
                                                       180,"Max")

test_label = pyBrellaSampling.classes.LabelClass("param.file")
test_label.file_name(["coord.file"])
test_label.add_bond("1,2", "testbond", 1.2)
test_label.add_dihedral("1,2,3,4","testdihedral",0,
                        "Min", 180, "Max")
test_label.add_dihedral("1,2,3,4", "MultiTest", 1, "State1",
                        2,"State2")

test_data = [0,1,0,0,1,2,1,300]

def test_LabelGen():
    GenBond, GenDihed = analysis.Label_Maker(test_label, f"{Test_Dir}LabelGen.tmp")
    exampleFile = utils.file_read(f"{Test_Dir}LabelGen.example")
    generatedFile = utils.file_read(f"{Test_Dir}LabelGen.tmp")
    os.remove(f"{Test_Dir}LabelGen.tmp")
    assert len(GenBond) == 1, f"Length of bond labels is {len(GenBond)} not 1"
    assert len(GenDihed) == 2, f"Length of dihedral labels is {len(GenDihed)} not 2"
    assert str(GenBond[0]) == str(test_bond), f"Bond Labeling has failed, Generated bond is {GenBond[0]}"
    assert str(GenDihed[0]) == str(test_dihedral), f"Dihedral labels have failed, Generated dihedral is {GenDihed[0]}"
    assert exampleFile == generatedFile, f"Generated file for the label maker isnt working."

def test_BondPlot():
    lines = f"""label add Bonds 0/1 0/2
label graph Bonds 0 testbond.dat
label delete Bonds 0
"""
    Genlines = analysis.tcl_bondPlot(test_bond)
    assert lines == Genlines, f"VMD Bond distance calculator not working... "

def test_DihedralPlot():
    lines = f"""label add Dihedrals 0/1 0/2 0/3 0/4
label graph Dihedrals 0 testdihedral.dat
label delete Dihedrals 0
"""
    Genlines = analysis.tcl_dihedPlot(test_dihedral)
    assert lines == Genlines, f"VMD Dihedral distance calculator not working... "

def test_BondAnalysis():
    State = analysis.tcl_bondAnalysis(test_bond, test_data,1)
    assert State == "Long", "Bond analysis not working... "

def test_DihedralAnalysis():
    State = analysis.tcl_dihedAnalysis(test_dihedral, test_data, 1)
    assert State == "Max", "Dihedral analysis not working..."

def test_DataframeGeneration():
    DataFrame = pyBrellaSampling.classes.DataClass("DataFrame")
    DataFrame = analysis.Labal_Analysis([test_bond], [test_dihedral], f"{Test_Dir}", 1, DataFrame)
    df = pd.concat(DataFrame.dat)
    df.to_csv(f"{Test_Dir}DataFrame.tmp")
    gen_df = pd.read_csv(f"{Test_Dir}DataFrame.tmp")
    example_df = pd.read_csv(f"{Test_Dir}DataFrame.example")
    os.remove(f"{Test_Dir}DataFrame.tmp")
    assert gen_df.equals(example_df), "DataFrame generation has failed"