import pyBrellaSampling.generate as generate
import pyBrellaSampling.utils as utils

Test_Dir = "./Testing/TestFiles/Generate/"


# def test_NamdFile():
#     class CalcClass:
#         NamdPath = "NAMD_PATH"
#         GPUNamd = "NAMD-GPU_PATH"
#         QMpath = "ORCA_PATH"
#         QMSel = "resname RESIDUE"
#         Charge = 1
#         Spin = 1
#         Method = "FUNCTIONAL"
#         Basis = "BASIS"
#         Threads = 1
#         CellVec = 1
#         Temp = 300
#         QM = "off"
#         CutOff = 8.0
#         Name = "Calculation"
#         OutFile = "CalcOut"
#         Force = 1
#         Id =
#         Ensemble = "NVT"
#         Steps =
#         TimeStep =
#         TimeOut =
#         RestOut =
#         TrajOut =
#         Shake =
#
#     GenFile = generate.Namd_File()
#     TestFile = utils.file_read()
#     assert GenFile == TestFile, "File Generation hasnt worked... Have you changed a default parameter?"
