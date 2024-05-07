import logging as log
import numpy
import Tools.utils as utils
import pandas as pd
import subprocess
import os
import matplotlib.pyplot as plt

from pyBrellaSampling.Tools.classes import BondClass, DihedralClass, DataClass, LabelClass
from pyBrellaSampling.Tools.globals import verbosity, WorkDir, DryRun, parmfile
import pyBrellaSampling.Tools.InputParser as input

def analysis(Umbrella: UmbrellaClass):
    Labels = LabelClass(f"../{parmfile}")
    Labels = input.BondsInput(f"{WorkDir}Bonds.dat", Labels)
    Labels = input.DihedralInput(f"{WorkDir}Dihedral.dat", Labels)
    core_load = []
    core_load.append(f"mol new {parmfile}")
    dataframe = DataClass("Production")
    for i in range(Umbrella.Bins):
        FileNames = []
        for j in range(1,20):
            FileNames.append(f"prod_{j}.{i}.dcd")
        Labels.file_name(FileNames)
        Bonds, Dihedrals = Anal.Label_Maker(Labels, f"{WorkDir}{i}/label_maker.tcl")
        if DryRun == False:
            subprocess.run([f"cd {WorkDir}{i} ; vmd -dispdev text -e label_maker.tcl ; cd ../"], shell=True, capture_output=True)
        dataframe = Anal.Labal_Analysis(Bonds, Dihedrals, f"{WorkDir}{i}/", i, dataframe)
        core_load.append(f"mol addfile ./{i}/equil.{i}.restart.coor ")
    with open(f"{WorkDir}prod_load.tcl",'w') as f:
        for i in range(len(core_load)):
            print(core_load[i], file=f)
    df = pd.concat(dataframe.dat)
    print(df.loc[df["Name"] == "c6Ring", "Data"].shape)
    print(type(Umbrella.atom1), Umbrella.atom2, Umbrella.atom3)
    for bond in Bonds:
        if ((bond.at1 == Umbrella.atom1) or (bond.at2 == Umbrella.atom1)) and ((bond.at1 == Umbrella.atom2) or (bond.at2 == Umbrella.atom2)) and (Umbrella.atom3 == 0):
            reactioncoordinate = bond.name
            print(f"INFO: Reaction coordinate is {reactioncoordinate}" if verbosity >=2 else "", end="")
    for dihed in Dihedrals:
        if ((dihed.at1 == Umbrella.atom1) or (dihed.at4 == Umbrella.atom1)) and ((dihed.at2 == Umbrella.atom2) or (dihed.at3 == Umbrella.atom2)) and ((dihed.at3 == Umbrella.atom3) or (dihed.at2 == Umbrella.atom3)) and ((dihed.at4 == Umbrella.atom4) or (dihed.at1 == Umbrella.atom4)):
            reactioncoordinate = dihed.name
            print(f"INFO: Reaction coordinate is {reactioncoordinate}" if verbosity >=2 else "", end="")
    try:
        os.mkdir(f"{WorkDir}Figures/")
    except:
        print("INFO: Figures directory already exists"if verbosity >=2 else "", end="")
    for bond in Bonds:
        plt.hist(df.loc[df["Name"] == bond.name, "Data",], 100,color="black")
        plt.title(f"{bond.name} bond")
        plt.xlabel("Distance")
        plt.ylabel("Count")
        plt.savefig(f"{WorkDir}Figures/{bond.name}.eps",transparent=True)
        if verbosity >= 2:
            plt.show()
        else:
            plt.clf()
        d = plt.hist2d(df.loc[df["Name"] == reactioncoordinate, "Data"], df.loc[df["Name"] == bond.name, "Data"],
                    100, cmap="binary")
        plt.title(f"Reaction coordinate vs. {bond.name} bond")
        plt.ylabel(f"{bond.name} bond distance")
        plt.xlabel("Reaction coordinate")
        plt.savefig(f"{WorkDir}Figures/{bond.name}_2d.eps",transparent=True)
        with open(f"{WorkDir}Figures/{bond.name}_2d.dat", 'w') as f:
            print(f"x\ty\tcount", file=f)
            for i in range(len(d[0])):
                for j in range(len(d[0])):
                    print(f"{d[1][i]}\t{d[2][j]}\t{d[0][i][j]}",file=f )
        if verbosity >= 2:
            plt.show()
        else:
            plt.clf()
    for dihed in Dihedrals:
        plt.hist(df.loc[df["Name"] == dihed.name, "Data"], 100, color="black")
        plt.title(f"{dihed.name} dihedral")
        plt.xlabel("Angle")
        plt.ylabel("Count")
        plt.savefig(f"{WorkDir}Figures/{dihed.name}.eps",transparent=True)
        if verbosity >= 2:
            plt.show()
        else:
            plt.clf()
        d = plt.hist2d(df.loc[df["Name"] == reactioncoordinate, "Data"], df.loc[df["Name"] == dihed.name, "Data"],100, cmap="binary")    
        plt.title(f"Reaction coordinate vs. {dihed.name} dihedral")
        plt.ylabel(f"{dihed.name} dihedral angle")
        plt.xlabel("Reaction coordinate")
        plt.savefig(f"{WorkDir}Figures/{dihed.name}_2d.eps",transparent=True)
        with open(f"{WorkDir}Figures/{dihed.name}_2d.dat", 'w') as f:
            print(f"x\ty\tcount", file=f)
            for i in range(len(d[0])):
                for j in range(len(d[0])):
                    print(f"{d[1][i]}\t{d[2][j]}\t{d[0][i][j]}",file=f )
        if verbosity >= 1:
            plt.show()
        else:
            plt.clf()
    df.to_csv(f"{WorkDir}Figures/Data.csv")
    
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
        Dist1 = numpy.absolute(dihed.target1 - data[i])
        # print(f"D1 {Dist1}")
        Dist2 = numpy.absolute(dihed.target2 - data[i])
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

def glue_stick(Umbrella, NumJobs, file):
    """Sticks split jobs back together for analysis."""
    for i in range(Umbrella.Bins):
        step, value = ["#Step"], ["Value"]
        for j in range(NumJobs):
            try:
                coords, datay = utils.data_2d(f"{WorkDir}{i}/{file}_{j+1}.{i}.colvars.traj")
            except FileNotFoundError:
                print(f"{WorkDir}{i}/{file}_{j+1}.{i}.colvars.traj is not found. moving on...")
                pass
            step = step + [int(k+(len(step)-1)) for k in coords if k != 0]      # NAMD prints state 0 at the start so will need to remove repeats.
            value = value + [datay[k] for k in range(len(coords)) if k != 0]
        utils.file_2dwrite(path=f"{WorkDir}{i}/{file}.{i}.colvars.traj",
                           x=step, y=value)


def Error_Check(Umbrella, Errors, Step=0, ):
    """Identifies whether a window should be used for wham calculation."""
    print(f"Checking for errors in step {Step}")
    Labels = LabelClass(f"../{parmfile}")
    Labels.clear_Vars()
    Labels = input.BondsInput(f"{WorkDir}BondErrors.dat", Labels)
    Labels = input.DihedralInput(f"{WorkDir}DihedralErrors.dat", Labels)
    core_load = []
    core_load.append(f"mol new {parmfile}")
    dataframe = DataClass("Production")
    for i in range(Umbrella.Bins):
        Path = f"{WorkDir}{i}/"
        # print(Labels.bond)
        if i in Errors:
            continue
        Labels.file_name([f"prod_{Step}.{i}.dcd"])
        Bonds, Dihedrals = Label_Maker(Labels, f"{WorkDir}{i}/label_maker.tcl")
        subprocess.run([f"cd {WorkDir}{i} ; vmd -dispdev text -e label_maker.tcl ; cd ../"], shell=True, capture_output=True)
        dataframe = Labal_Analysis(Bonds, Dihedrals, f"{WorkDir}{i}/", i, dataframe)
        if len(Dihedrals) > 0:
            Error = True
        else:
            Error = False
        for j in range(len(Dihedrals)):
            steps, data = utils.data_2d(f"{Path}{Dihedrals[j].name}.dat")
            break_point = tcl_dihedAnalysis(Dihedrals[j], data, i, Error_Anal=True)
            if break_point == len(data):
                Error = False
            elif break_point >= 0.9*len(data):
                Error = False
                print(f"WARNING: There are some problems with window {i}, They are near the end so ignoring." if verbosity >=1 else "", end="")
            else:
                Error = True
                Errors.append(i)
                print(f"WARNING: Error identified in step {Step} of window {i}, Dihedral Error" if verbosity >=1 else "", end="")
                print(f"Breakpoint is: {break_point} out of {len(data)} steps" if verbosity >=1 else "", end="")
                break
        if Error != True:
            for j in range(len(Bonds)):
                steps, data = utils.data_2d(f"{Path}{Bonds[j].name}.dat")
                break_point = tcl_bondAnalysis(Bonds[j], data, i, Error_Anal=True)
                if break_point == len(data):
                    Error = False
                elif break_point >= 0.5*len(data):
                    Error = False
                    print(f"WARNING: There are some problems with step {Step}, They are near the end so ignoring." if verbosity >=1 else "", end="")
                else:
                    Error = True
                    Errors.append(i)
                    print(f"WARNING: Error identified in step {Step} of window {i}, Bond Error" if verbosity >=1 else "", end="")
                    break
    return Errors