import pyBrellaSampling.Standalone as Standalone
import pyBrellaSampling.Tools.utils as utils
from unittest.mock import patch
import sys
import os


Test_Dir="./Testing/TestFiles/Standalone/"
# Test_Dir="./"

def test_Standalone():
    Arguments = [f"Standalone", "-i=Test.inp", f"-wd={Test_Dir}"]
    with patch.object(sys, "argv", Arguments):
        Standalone.main()
    GenFile = utils.file_read(f"{Test_Dir}NAME.conf")    
    TestFile = utils.file_read(f"{Test_Dir}Test.conf")
    assert GenFile == TestFile, "Input file generated incorrectly."
    GenFile = utils.file_read(f"{Test_Dir}colvars.conf")
    TestFile = utils.file_read(f"{Test_Dir}Test.colvars")
    assert GenFile == TestFile, "SMD colvars not working."
    os.remove(f"{Test_Dir}Colvar_prep.tcl")
    os.remove(f"{Test_Dir}colvars.conf")
    os.remove(f"{Test_Dir}Name.conf")
    os.remove(f"{Test_Dir}qm_prep.tcl")