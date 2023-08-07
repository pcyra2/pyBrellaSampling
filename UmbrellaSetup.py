import pyBrellaSampling.generate as generate
import logging as log
import subprocess

def run_setup(args, Umbrella, Calc, Job, DryRun=False):
    log.warning("Setting up pulls")
    make_umbrellaDirs(args)
    Calc.Set_Ensemble("NVT")
    Umbrella, pullJobs = generate.Pull_Setup(Umbrella, Calc, Job)
    make_runfile(Job, Umbrella, pullJobs)
    if DryRun == False:
        run_pullScript()
    return "Umbrella pull has completed setup"

def make_umbrellaDirs(args):
    log.warning("Making umbrella directories")
    generate.make_dirs(args.UmbrellaBins, args.WorkDir)

def setup_pulls(Umbrella, Calc, Job):
    log.warning("Setting up pulls")
    Umbrella, pullJobs = generate.Pull_Setup(Umbrella, Calc, Job)
    return Umbrella, pullJobs

def make_runfile(Job, Umbrella, pullJobs):
    log.warning("Generating pull.sh script.")
    with open(f"{Job.WorkDir}pull.sh", 'w') as f:
        print("#!/bin/bash", file=f)
        for i in range(Umbrella.StartBin, Umbrella.Bins):
            print(pullJobs[i], file=f)
        for i in range(0, Umbrella.StartBin):
            print(pullJobs[Umbrella.StartBin - i - 1], file=f)

def run_pullScript(loc="./"):
    log.warning("Running pull command")
    run_out = subprocess.run([f"sh {loc}pull.sh"], shell=True, capture_output=True)
    return run_out