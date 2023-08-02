import logging as log
import pyBrellaSampling.utils as utils
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


class BondClass:
    def __init__(self, atom1, atom2, name, threshold):
        self.at1 = atom1
        self.at2 = atom2
        self.name = name
        self.thresh = threshold
    def __repr__(self):
        return f"Atoms: {self.at1} {self.at2}, name: {self.name}, threshold: {self.thresh}"

class DihedralClass:
    def __init__(self, atom1, atom2, atom3, atom4, name, target1, t1name, target2, t2name):
        self.at1 = atom1
        self.at2 = atom2
        self.at3 = atom3
        self.at4 = atom4
        self.name = name
        self.target1 = target1
        self.target1Name = t1name
        self.target2 = target2
        self.target2Name = t2name
    def __repr__(self):
        return f"Atoms: {self.at1} {self.at2} {self.at3} {self.at4}, Name: {self.name}, {self.target1Name} = {self.target1}, {self.target2Name} = {self.target2}"




def tcl_bondPlot(bond):
    lines = f"""label add Bonds 0/{bond.at1} 0/{bond.at2}
label graph Bonds 0 {bond.name}.dat
label delete Bonds 0
"""
    return lines

def tcl_bondAnalysis(bond, data, simulation):
    for i in range(len(data)):
        if i == 0:
            prevState = "None"
        else:
            prevState = State
        if data[i] <= 0.8 * bond.thresh:
            # log.warning(f"{bond.name} is short ({data[i]}) for window {simulation}, on trajectory {i}")
            State = "Short"
        elif data[i] >= 1.2 * bond.thresh:
            # log.warning(f"{bond.name} is long ({data[i]}) for window {simulation}, on trajectory {i}")
            State = "Long"
        else:
            State = "Normal"
        if State !=  prevState:
            log.warning(f"{bond.name} has changed to {State} for window {simulation} on step {i}")
    return State

def tcl_dihedPlot(dihed):
    lines = f"""label add Dihedrals 0/{dihed.at1} 0/{dihed.at2} 0/{dihed.at3} 0/{dihed.at4}
label graph Dihedrals 0 {dihed.name}.dat
label delete Dihedrals 0
"""
    return lines

def tcl_dihedAnalysis(dihed, data, window):
    for i in range(len(data)):
        if i == 0:
            prevState = "None"
        else:
            prevState = State
        Dist1 = np.absolute(dihed.target1 - data[i])
        # print(f"D1 {Dist1}")
        Dist2 = np.absolute(dihed.target2 - data[i])
        # print(f"D2 {Dist2}")
        if Dist1 < Dist2:
            State = dihed.target1Name
        elif Dist2 < Dist1:
            State = dihed.target2Name
        if prevState != State:
            log.warning(f"Dihedral flip from {dihed.name} {prevState} to {State} for window {window}, on step {i}")
    return State

def Label_Maker(Label, path):
    Bonds = [None] * len(Label.bond)
    Dihedrals = [None] * len(Label.dihedral)
    with open(path, 'w') as f:
        print(f"""mol new {Label.file} waitfor -1
mol addfile {Label.parm}

""", file=f)
        for i in range(len(Label.bond)):
            atoms = Label.bond[i]
            Bonds[i] = BondClass(atoms.split(",")[0], atoms.split(",")[1], Label.bondName[i], Label.bondThresh[i])
            lines = tcl_bondPlot(Bonds[i])
            print(lines, file=f)
        for i in range(len(Label.dihedral)):
            atoms = Label.dihedral[i]
            Dihedrals[i] = DihedralClass(atoms.split(",")[0],
                                         atoms.split(",")[1],
                                         atoms.split(",")[2],
                                         atoms.split(",")[3],
                                         Label.dihedralName[i],
                                         Label.dihedralTarget1[i],
                                         Label.dihedralTarget1Name[i],
                                         Label.dihedralTarget2[i],
                                         Label.dihedralTarget2Name[i],)
            lines = tcl_dihedPlot(Dihedrals[i])
            print(lines, file=f)
        print("quit", file=f)
    return Bonds, Dihedrals

def Labal_Analysis(Bonds, Dihedrals, Path, window, dataframe):
    for i in range(len(Bonds)):
        try:
            steps, data = utils.data_2d(f"{Path}{Bonds[i].name}.dat")
            # plt.hist(data,10)
            tcl_bondAnalysis(Bonds[i], data, window)
            dataframe.add_data(Bonds[i].name, window, data)
        except FileNotFoundError:
            log.error(f"Simulation {window} has errors...")
            break
    for i in range(len(Dihedrals)):
        try:
            steps, data = utils.data_2d(f"{Path}{Dihedrals[i].name}.dat")
            tcl_dihedAnalysis(Dihedrals[i], data, window)
            dataframe.add_data(Dihedrals[i].name, window, data)
        except FileNotFoundError:
            # log.warning(f"Simulation {window} has errors...")
            break
    return dataframe


