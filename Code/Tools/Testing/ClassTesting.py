from pyBrellaSampling.Code.Tools.classes import *
from pprint import pprint

if __name__ == "__main__":
    mm = MMClass("namd", "system.parm7", "system.rst7")
    min = mm.minimize("system.rst7", "min", 500)
    heat = mm.heat("min", "heat", 1000, 2, 0, 300)
    equil = mm.constant(infile="heat", outfile="equil", steps=1000, timestep=2, traj_steps=100, temperature=300, pressure=None)
    # print(mm.jobs["min"])
    colvar = ColvarClass()
    colvar.atomic_colvar([1,2], 300, 100, 1.3, 3, 7, 0.1)
    file = colvar.colvarPDB_gen(mm)
    qm = QMClass("PBE", "6-31G*", 0,0,"orca")
    qm.init_qmmm("resname CTN")
    mm.software.set_qm(qm)
    qm_file = mm.software.qmPDB_gen()
    print(qm_file)