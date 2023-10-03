import emcee as emcee
import os
import pyBrellaSampling.utils as utils
import matplotlib.pyplot as plt
import logging as log
import subprocess
import numpy as np

def autocorrelate(data):
    try:
        integral_time = emcee.autocorr.integrated_time(data, c=1)
    except emcee.autocorr.AutocorrError:
        integral_time = 1
    return integral_time

def Init_Wham(Job, Umbrella, Wham, Verbose=False):
    hist_bar = []
    hist_count = []
    try:
        os.mkdir(f"{Job.WorkDir}/WHAM/")
    except:
        log.warning("Wham directory already exists")
    if os.path.exists(f"{Job.WorkDir}/WHAM/{Wham.Name}metadata.dat"):
        os.remove(f"{Job.WorkDir}/WHAM/{Wham.Name}metadata.dat")
    for i in range(0,Umbrella.Bins):    # Ignore line 1...
        time, value = utils.data_2d(f"{Job.WorkDir}{i}/{Wham.Name}.{i}.colvars.traj")
        integral_time = autocorrelate(value)
        counts, bins, bars = plt.hist(value, 100)
        hist_count.extend(counts)
        hist_bar.extend(bins)
        with open(f"{Job.WorkDir}WHAM/{Wham.Name}metadata.dat",'a') as f:
            txt = f"../{i}/{Wham.Name}.{i}.colvars.traj {round(Umbrella.BinVals[i],3)} {float(Wham.Force)} {int(integral_time)}" ##ORCA equation for harmonic restraints: V = 0.5*k (a - a0)^2
            print(txt, file = f)
        # subprocess.run([f"head {Job.WorkDir}{i}/{Wham.Name}.{i}.colvars.traj -n 4005 > {Job.WorkDir}{i}/{Wham.Name}.{i}.colvars.2ps.traj"], shell=True,)
        # print(f"head {Job.WorkDir}{i}/{Wham.Name}.{i}.colvars.traj -n 4005 > {Job.WorkDir}{i}/{Wham.Name}.{i}.colvars.2ps.traj")
    plt.xlabel("Reaction coordinate")
    plt.ylabel("Count")
    if Verbose == True:
        plt.show()
    else:
        plt.clf()
    with open(f"{Job.WorkDir}WHAM/{Wham.Name}UmbrellaHist.dat", 'w') as f:
        for i in range(len(hist_count)):
            print(str(hist_bar[i]) + "\t" + str(hist_count[i]), file=f)
    with open(f"{Job.WorkDir}WHAM/wham.sh", 'w') as f:
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
wham {P} {min} {max} {Umbrella.Bins-1} 1e-06 300 0 {Wham.Name}metadata.dat out.pmf 10 60
sed '1d' out.pmf | awk '{"{"}print $1,"",$2{"}"}' > plot_free_energy.dat
"""
        print(text, file=f)

def Run_Wham(Job,Umbrella, Verbose=False):# Umbrella, Wham):
    print("Running WHAM")
    out = subprocess.run(f"cd {Job.WorkDir}WHAM ; sh wham.sh ; cd ../", shell=True, capture_output=True)
    if "wham.sh" in out.stderr.decode():
        raise Exception(f"Problem with running wham: {out.stderr.decode()}")
    if "No such file " in out.stdout.decode():
        print(out.stdout.decode())
        raise Exception("metadata file not found when runnning wham.sh")
    data = utils.file_read(path=f"{Job.WorkDir}WHAM/out.pmf")
    array = [data[i] for i in range(1,Umbrella.Bins) ]
    # print(array)
    x = np.zeros(len(array))
    y = np.zeros(len(array))
    err = np.zeros(len(array))
    for line in range(len(array)):
        columns = array[line].split()
        x[line] = columns[0]
        y[line] = columns[1]
        err[line] = columns[2]
    # plt.plot(x,y,)
    plt.errorbar(x,y,yerr=err, c="black", capsize=5)
    plt.xlabel("Reaction coordinate")
    plt.ylabel("Energy (kcal $mol^{-1}$)")
    plt.title(f"Wham PMF")
    plt.savefig(f"{Job.WorkDir}WHAM/PMF.eps",transparent=True)
    if Verbose == True:
        plt.show()
    else:
        plt.clf()
    plt.errorbar(y=x,x=y,xerr=err, c="black", capsize=5)
    plt.ylabel("Reaction coordinate")
    plt.xlabel("Energy (kcal $mol^{-1}$)")
    plt.title(f"Wham PMF")
    plt.savefig(f"{Job.WorkDir}WHAM/PMF_AxisSwap.eps",transparent=True)
    plt.clf()
    subprocess.run(f"head -n {Umbrella.Bins} {Job.WorkDir}WHAM/out.pmf > {Job.WorkDir}WHAM/PMF.dat", shell=True)
    subprocess.run(
    f"sed -i \"0,/+\/-/s/+\/-/Err1/\" {Job.WorkDir}WHAM/PMF.dat ",
        shell=True)
    subprocess.run(
        f"sed -i \"0,/+\/-/s/+\/-/Err2/\" {Job.WorkDir}WHAM/PMF.dat",
        shell=True)
    subprocess.run(
        f"sed -i \"s/#Coor/Coor/g\" {Job.WorkDir}WHAM/PMF.dat",
        shell=True)
