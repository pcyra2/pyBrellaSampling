import emcee as emcee
import os
import pyBrellaSampling.utils as utils
import matplotlib.pyplot as plt
import logging as log
import subprocess
import numpy as numpy
import math
from pyBrellaSampling.Tools.globals import verbosity, WorkDir, DryRun
import pyBrellaSampling.Tools.analysis as Anal
from pyBrellaSampling.Tools.classes import *


def autocorrelate(data):
    try:
        integral_time = emcee.autocorr.integrated_time(data, c=1)
    except emcee.autocorr.AutocorrError:
        integral_time = 1
    except IndexError:
        integral_time = 1    
    return integral_time

def Init_Wham( Umbrella, Wham, WhamIgnore=[]):
    hist_bar = []
    hist_count = []
    try:
        os.mkdir(f"{WorkDir}/WHAM/")
    except:
        log.warning("Wham directory already exists")
    if os.path.exists(f"{WorkDir}/WHAM/{Wham.Name}metadata.dat"):
        os.remove(f"{WorkDir}/WHAM/{Wham.Name}metadata.dat")
    for i in range(0,Umbrella.Bins):    # Ignore line 1...
        if i in WhamIgnore:
            continue
        time, value = utils.data_2d(f"{WorkDir}{i}/{Wham.Name}.{i}.colvars.traj")
        integral_time = autocorrelate(value)
        counts, bins, bars = plt.hist(value, 100)
        hist_count.extend(counts)
        hist_bar.extend(bins)
        with open(f"{WorkDir}WHAM/{Wham.Name}metadata.dat",'a') as f:
            txt = f"../{i}/{Wham.Name}.{i}.colvars.traj {round(Umbrella.BinVals[i],3)} {float(Wham.Force)} {int(integral_time)}" ##ORCA equation for harmonic restraints: V = 0.5*k (a - a0)^2
            print(txt, file = f)
        # subprocess.run([f"head {WorkDir}{i}/{Wham.Name}.{i}.colvars.traj -n 4005 > {WorkDir}{i}/{Wham.Name}.{i}.colvars.2ps.traj"], shell=True,)
        # print(f"head {WorkDir}{i}/{Wham.Name}.{i}.colvars.traj -n 4005 > {WorkDir}{i}/{Wham.Name}.{i}.colvars.2ps.traj")
    plt.xlabel("Reaction coordinate")
    plt.ylabel("Count")
    if verbosity >= 2 :
        plt.show()
    else:
        plt.clf()
    with open(f"{WorkDir}WHAM/{Wham.Name}UmbrellaHist.dat", 'w') as f:
        for i in range(len(hist_count)):
            print(str(hist_bar[i]) + "\t" + str(hist_count[i]), file=f)
    with open(f"{WorkDir}WHAM/wham.sh", 'w') as f:
        if Wham.Type == "Periodic":
            P="P"
        else:
            P=""
        if Umbrella.BinVals[0] >  Umbrella.BinVals[Umbrella.Bins-1]:
            min = Umbrella.BinVals[Umbrella.Bins-1]
            max = Umbrella.BinVals[0]
        else:
            max = Umbrella.BinVals[Umbrella.Bins-1]
            min = Umbrella.BinVals[0]
        text = f"""#!/bin/bash
wham {P} {min} {max} {Umbrella.Bins-1 - len(WhamIgnore)} 1e-06 300 0 {Wham.Name}metadata.dat out.pmf 10 60
sed '1d' out.pmf | awk '{"{"}print $1,"",$2{"}"}' > plot_free_energy.dat
"""
        print(text, file=f)

def Run_Wham(Umbrella, WhamIgnore=[]):# Umbrella, Wham):
    print("Running WHAM")
    out = subprocess.run(f"cd {WorkDir}WHAM ; sh wham.sh ; cd ../", shell=True, capture_output=True)
    if "wham.sh" in out.stderr.decode():
        raise Exception(f"Problem with running wham: {out.stderr.decode()}")
    if "No such file " in out.stdout.decode():
        print(out.stdout.decode())
        raise Exception("metadata file not found when runnning wham.sh")
    data = utils.file_read(path=f"{WorkDir}WHAM/out.pmf")
    array = [data[i] for i in range(1,Umbrella.Bins-len(WhamIgnore)) ]
    x = numpy.zeros(len(array))
    y = numpy.zeros(len(array))
    err = numpy.zeros(len(array))
    for line in range(0, len(array)):
        columns = array[line].split()
        x[line] = columns[0]
        y[line] = columns[1]
        err[line] = columns[2]
    #print(y)
    # plt.plot(x,y,)
    plt.errorbar(x,y,yerr=err, c="black", capsize=5)
    plt.xlabel("Reaction coordinate")
    plt.ylabel("Energy (kcal $mol^{-1}$)")
    plt.title(f"Wham PMF")
    plt.savefig(f"{WorkDir}WHAM/PMF.eps",transparent=True)
    if verbosity >=2 :
        plt.show()
    else:
        plt.clf()
    plt.errorbar(y=x,x=y,xerr=err, c="black", capsize=5)
    plt.ylabel("Reaction coordinate")
    plt.xlabel("Energy (kcal $mol^{-1}$)")
    plt.title(f"Wham PMF")
    plt.savefig(f"{WorkDir}WHAM/PMF_AxisSwap.eps",transparent=True)
    plt.clf()
    subprocess.run(f"head -n {Umbrella.Bins} {WorkDir}WHAM/out.pmf > {WorkDir}WHAM/PMF.dat", shell=True)
    subprocess.run(
    f"sed -i \"0,/+\/-/s/+\/-/Err1/\" {WorkDir}WHAM/PMF.dat ",
        shell=True)
    subprocess.run(
        f"sed -i \"0,/+\/-/s/+\/-/Err2/\" {WorkDir}WHAM/PMF.dat",
        shell=True)
    subprocess.run(
        f"sed -i \"s/#Coor/Coor/g\" {WorkDir}WHAM/PMF.dat",
        shell=True)
    return y, err


def convergence(Calc: CalcClass, AnalFile: str, equil_length: int, prod_length: int, Umbrella: UmbrellaClass):
    if Calc.MaxSteps == 0:
        NumJobs = 1
    else:
        if "prod" in AnalFile.casefold():
            steps = prod_length
        elif "equil" in AnalFile.casefold():
            steps = equil_length
        else:
            steps = 0
        NumJobs = math.ceil(steps / Calc.MaxSteps)
    if verbosity >= 1:
        print(f"Number of steps to glue together is {NumJobs}")
    try:
        os.mkdir(f"{WorkDir}WHAM/Conv")
    except:
        if verbosity >= 1:
            print("Figures directory already exists")
        else:
            pass
    for i in range(1,NumJobs+1):
        print(f"Performing WHAM on step {i}")
        WhamIgnore = Anal.Error_Check(Umbrella=Umbrella, Errors=WhamIgnore, Step=i)
        print(f"WARNING: Error windows are: {WhamIgnore}")
        if i > 1:
            prev_y = y
            #prev_Err = Err
        Anal.glue_stick(Umbrella, NumJobs=i, file=AnalFile)
        if Umbrella.atom3 != 0:
            periodicity = "Periodic"
        else:
            periodicity = "discrete"
        wham = WhamClass(AnalFile, Umbrella.ConstForce, periodicity)
        Wham.Init_Wham(Umbrella, wham, WhamIgnore=WhamIgnore)
        y, Err = Wham.Run_Wham(Umbrella, WhamIgnore=WhamIgnore)
        convergence=False
        if i > 1:
            Diffs = numpy.zeros(len(y))
            for j in range(len(y)):
                if y[j] == numpy.inf or prev_y[j] == numpy.inf:
                    diff = 0
                else:
                    diff = abs(y[j]-prev_y[j])
                if diff == "nan":
                    diff = 0
                Diffs[j] = diff
            if numpy.max(Diffs) <= 0.1:
                convergence = True
            print(f"INFO: Maximum difference = {numpy.max(Diffs)}"if verbosity >=1 else "", end="")
            if convergence == True:
                print(f"SUCCESS: Convergence achieved"if verbosity >=1 else "", end="")
                if numpy.max(Err) > 0.10:
                    print(f"INFO: Error bars still too large: {numpy.max(Err)}" if verbosity >=1 else "", end="")
                else:
                    break
        subprocess.run(f"head -n {Umbrella.Bins} {WorkDir}WHAM/out.pmf > {WorkDir}WHAM/Conv/{i}.pmf", shell=True)
        subprocess.run(f"mv {WorkDir}WHAM/PMF.eps {WorkDir}WHAM/Conv/{i}.eps", shell=True)
        subprocess.run(f"sed -i \"0,/+\/-/s/+\/-/Err1/\" {WorkDir}WHAM/Conv/{i}.pmf " , shell=True)
        subprocess.run(
            f"sed -i \"0,/+\/-/s/+\/-/Err2/\" {WorkDir}WHAM/Conv/{i}.pmf ",
            shell=True)
        subprocess.run(
            f"sed -i \"s/#Coor/Coor/g\" {WorkDir}WHAM/Conv/{i}.pmf",
            shell=True)