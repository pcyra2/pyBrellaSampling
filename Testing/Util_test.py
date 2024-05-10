import pyBrellaSampling.Tools.utils as Utils
import os
import numpy as numpy
import argparse as ap
import pyBrellaSampling.Tools.globals as globals
from pyBrellaSampling.Tools.classes import MMClass

Test_Dir = "./Testing/TestFiles/Utils/"
globals.init(wd=Test_Dir) # inits global vars.
Data_1D = numpy.array([0,1,2,3,4,5,6,7,8,9])

Test_args = {"input": "TestOut.tmp",
             "WorkDir": "./Testing/",
             "NamdPath": "PATH/TO/namd3",
             "NumCores": 1,
             "JobType": "TEST",
             "globals.verbosity": 1,
             "globals.DryRun": True,
             "UmbrellaMin": 1,
             "UmbrellaWidth": 1,
             "UmbrellaBins": 1,
             "PullForce": 1,
             "Force": 1,
             "StartDistance": 1,
             "AtomMask": "Comma,Separated,Atom,Index",
             "QmSelection": "resname RESIDUES",
             "QmCharge": 1,
             "QmSpin": 1,
             "QmMethod": "Functional",
             "QmBasis": "BasisSet",
             "QmPath": "PATH/TO/orca",
             "WhamFile": "FileToPerformWHAM"}


def test_FileRead():
    line = Utils.file_read(f"{Test_Dir}Single_Line.example")
    assert line[0] == "Single_Line", "File reader doesnt work"

def test_FileWrite():
    Utils.file_write(f"{Test_Dir}File_Write.tmp", ["Test"])
    line = Utils.file_read(f"{Test_Dir}File_Write.tmp")
    os.remove(f"{Test_Dir}File_Write.tmp")
    assert line[0] == "Test\n", "File reader doesnt work"

def test_DictRead():
    dictionary = {"Name" : "Ross"}
    ans = Utils.dict_read(f"{Test_Dir}TestDict.example")
    assert dictionary == ans, "Cannot read dictionaries."


def test_DictWrite():
    dictionary = {"Name" : "Ross"}
    Utils.dict_write(f"{Test_Dir}Dict.tmp", dictionary)
    ans = Utils.dict_read(f"{Test_Dir}Dict.tmp")
    os.remove(f"{Test_Dir}Dict.tmp")
    assert dictionary == ans, "Something is wrong with reading or writing dictionaries."

def test_Read2D():
    col1, col2 = Utils.data_2d(f"{Test_Dir}2D_Data.example")
    assert col1.any() == Data_1D.any(), "Cannot read 2D data file column 1"
    assert col2.any() == Data_1D.any(), "Cannot read 2D data file column 2"

def test_getCellVec():
    globals.init(wd=Test_Dir)
    MM = MMClass()
    MM.Set_Files(parm="TestParm", ambercoor="Fake_AmberCoordinates.example")
    CellVec = Utils.get_cellVec(MM)
    assert CellVec == 102.4, "Cannot attain Cell vector"