import emcee as emcee
import os
import pyBrellaSampling.utils as utils
import matplotlib.pyplot as plt
import logging as log
import subprocess

def autocorrelate(data):
    try:
        integral_time = emcee.autocorr.integrated_time(data, c=1)
    except emcee.autocorr.AutocorrError:
        integral_time = 1
    return integral_time

def Init_Wham(Job, Umbrella, Wham):
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
            txt = f"../{i}/{Wham.Name}.{i}.colvars.2ps.traj {round(Umbrella.BinVals[i],3)} {int(Wham.Force)} {int(integral_time)}" ##ORCA equation for harmonic restraints: V = 0.5*k (a - a0)^2
            print(txt, file = f)
        subprocess.run([f"head {Job.WorkDir}{i}/{Wham.Name}.{i}.colvars.traj -n 4005 > {Job.WorkDir}{i}/{Wham.Name}.{i}.colvars.2ps.traj"], shell=True,)
        # print(f"head {Job.WorkDir}{i}/{Wham.Name}.{i}.colvars.traj -n 4005 > {Job.WorkDir}{i}/{Wham.Name}.{i}.colvars.2ps.traj")
    with open(f"{Job.WorkDir}WHAM/{Wham.Name}UmbrellaHist.dat", 'w') as f:
        for i in range(len(hist_count)):
            print(str(hist_bar[i]) + "\t" + str(hist_count[i]), file=f)
    with open(f"{Job.WorkDir}WHAM/wham.sh", 'w') as f:
        text = f"""#!/bin/bash
wham {Umbrella.BinVals[0]} {Umbrella.BinVals[Umbrella.Bins-1]} {Umbrella.Bins} 1e-06 300 0 {Wham.Name}metadata.dat out.pmf 10 60
sed '1d' out.pmf | awk '{"{"}print $1,"",$2{"}"}' > plot_free_energy.dat
"""
        print(text, file=f)