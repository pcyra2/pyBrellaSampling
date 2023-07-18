import pyBrellaSampling.utils as Utils
import os

Test_Dir = "./TestFiles/"

def test_FileRead():
    line = Utils.file_read(f"{Test_Dir}Single_Line.txt")
    assert line[0] == "Single_Line", "File reader doesnt work"

def test_FileWrite():
    Utils.file_write(f"{Test_Dir}File_Write.txt", "Test")
    line = Utils.file_read(f"{Test_Dir}File_Write.txt")
    os.remove(f"{Test_Dir}File_Write.txt")
    assert line[0] == "Test", "File reader doesnt work"

def test_Dictionaries():
    dictionary = {Name : "Ross"}
    Utils.dict_write(f"{Test_Dir}Dict.dat", dictionary)
    ans = Utils.dict_read(f"{Test_Dir}Dict.dat")
    os.remove(f"{Test_Dir}Dict.dat")
    assert dictionary == ans, "Something is wrong with reading or writing dictionaries."