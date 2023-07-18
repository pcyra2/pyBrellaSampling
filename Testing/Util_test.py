import pyBrellaSampling.utils as Utils
import os

def test_FileRead():
    line = Utils.file_read(f"./TestFiles/Single_Line.txt")
    assert line[0] == "Single_Line", "File reader doesnt work"

def test_FileWrite():
    Utils.file_write(f"./TestFiles/File_Write.txt", "Test")
    line = Utils.file_read(f"./TestFiles/File_Write.txt")
    os.remove(f"./TestFiles/File_Write.txt")
    assert line[0] == "Test", "File reader doesnt work"