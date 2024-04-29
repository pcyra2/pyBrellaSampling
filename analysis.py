import logging as log
import pyBrellaSampling.utils as utils
import numpy as np

from pyBrellaSampling.classes import BondClass, DihedralClass



def tcl_bondPlot(bond):
    """Creates a list of lines for an alaysis tcl script used by vmd which allows for the tracking of bonds."""
    lines = f"""label add Bonds 0/{bond.at1} 0/{bond.at2}
label graph Bonds 0 {bond.name}.dat
label delete Bonds 0
"""
    return lines

def tcl_bondAnalysis(bond, data, simulation, Error_Anal=False):
    """Runs analysis on 2d data file containing bond lengths and states whether they are within the expected threshold."""
    for i in range(len(data)):
        if i == 0:
            prevState = "None"
        else:
            prevState = State
        if data[i] <= 0.8 * bond.thresh:
            # log.warning(f"{bond.name} is short ({data[i]}) for window {simulation}, on trajectory {i}")
            State = "Short"
            if Error_Anal == True:
                print(f"BondLength = {data[i]}")
                return i
        elif data[i] >= 1.2 * bond.thresh:
            # log.warning(f"{bond.name} is long ({data[i]}) for window {simulation}, on trajectory {i}")
            State = "Long"
        else:
            State = "Normal"
        # if State !=  prevState:
        #     log.warning(f"{bond.name} has changed to {State} for window {simulation} on step {i}")
    if Error_Anal == False:
        return State
    else:
        return len(data)

def tcl_dihedPlot(dihed):
    """Creates a list of lines for an alaysis tcl script used by vmd which allows for the tracking of Dihedrals."""
    lines = f"""label add Dihedrals 0/{dihed.at1} 0/{dihed.at2} 0/{dihed.at3} 0/{dihed.at4}
label graph Dihedrals 0 {dihed.name}.dat
label delete Dihedrals 0
"""
    return lines

def tcl_dihedAnalysis(dihed, data, window, Error_Anal=False):
    """Runs analysis on 2d data file containing dihedral angles and states what state they are in."""
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
            if Error_Anal == True:
                return i
            State = dihed.target2Name
        # if prevState != State:
        #     log.warning(f"Dihedral flip from {dihed.name} {prevState} to {State} for window {window}, on step {i}")
    if Error_Anal == False:
        return State
    else:
        return len(data)

def Label_Maker(Label, path):
    """This generates the tcl script that automates the bond and dihedral analysis."""
    Bonds = [None] * len(Label.bond)
    Dihedrals = [None] * len(Label.dihedral)
    with open(path, 'w') as f:
        print(f"mol new {Label.parm} waitfor -1", file=f)
    with open(path, 'a') as f:
        for i in Label.file:
            print(f"mol addfile {i} waitfor -1", file=f)
# """, file=f)
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
                                         Label.dihedralTarget2Name[i], )
            lines = tcl_dihedPlot(Dihedrals[i])
            print(lines, file=f)
        print("quit", file=f)
    return Bonds, Dihedrals

def Labal_Analysis(Bonds, Dihedrals, Path, window, dataframe):
    """Saves bond and dihedral analysis, and controls errors for when files are missing ect."""
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

def glue_stick(Umbrella, Job, NumJobs, file):
    """Sticks split jobs back together for analysis."""
    for i in range(Umbrella.Bins):
        step, value = ["#Step"], ["Value"]
        for j in range(NumJobs):
            try:
                coords, datay = utils.data_2d(f"{Job.WorkDir}{i}/{file}_{j+1}.{i}.colvars.traj")
            except FileNotFoundError:
                print(f"{Job.WorkDir}{i}/{file}_{j+1}.{i}.colvars.traj is not found. moving on...")
                pass
            step = step + [int(k+(len(step)-1)) for k in coords if k != 0]      # NAMD prints state 0 at the start so will need to remove repeats.
            value = value + [datay[k] for k in range(len(coords)) if k != 0]
        utils.file_2dwrite(path=f"{Job.WorkDir}{i}/{file}.{i}.colvars.traj",
                           x=step, y=value)
