#!/usr/bin/python3
# import pyBrellaSampling.Code.Tools.io as io
# import pyBrellaSampling.Code.Tools.classes as classes
# import pyBrellaSampling.Code.Tools.QM.pyscf_tools as pyscf_tools

# import sys
# import os
# from pprint import pprint
# from pyscf import grad as pygrad
# # from gpu4pyscf import qmmm
# from pyscf.qmmm.itrf import _QMMMGrad
# from pyscf import lib
import numpy
import time
import os

def main():
    start = time.perf_counter()
    inputFilename = "qmmm_0.input"
    while True:
        time.sleep(5)
        if os.path.isfile(inputFilename+".result"):
            stop = time.perf_counter()
            print(f"INFO: QM calculation took {stop - start}")
            exit(0)
    

if __name__ == "__main__":
    main()
