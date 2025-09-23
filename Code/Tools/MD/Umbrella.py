import emcee as emcee
import os
import subprocess
import time
from pprint import pprint
import matplotlib.pyplot as plt

from pyBrellaSampling.Code.Tools.classes import ColvarClass, MMClass, VMDClass, TrackerClass, HPCClass
import pyBrellaSampling.Code.Tools.io as io
from pyBrellaSampling.Code.Tools.Utils import find_closest

class UmbrellaClass:
    data_types = ["pullValues", "equilValues", "prodValues"]
    data = {}
    def __init__(self, Colvar: ColvarClass, Temperature = 300):
        self.colvar = Colvar
        self.Temperature = Temperature
        initial_bin = find_closest(Colvar.bins, Colvar.initial)
        for i, bin in enumerate(Colvar.bins):
            if bin > initial_bin:
                previous = i - 1
            elif bin < initial_bin:
                previous = i + 1
            else:
                previous = None
                self.start_index = i
            self.data[i] = {"PullTime" : 0,
                              "ConstTime": 0,
                              "Value" : bin,
                              "pullValues":[],
                              "equilValues":[],
                              "prodValues":[],
                              "ERROR":False,
                              "PreviousWindow":previous}
    def add_data(self,bin:float, data_type:str, data:list, time:float):
        assert data_type in self.data_types, f"Expecting one of {self.data_types} as a data type"
        if "#step" in data[0]:
            data = data[1:]
        if len(data[0].split()) == 2:
            data = [float(dat.split()[1]) for dat in data if "step" not in dat]

        self.data[bin][data_type] += data
        if data_type == "pullValues":
            self.data[bin]["PullTime"] += time
        else:
            self.data[bin]["ConstTime"] += time
    def autocorrelate(self, values = "prodValues", EquilTime = 0):
        """Uses the emcee package to calculate the time autocorrelation of the collective variable.

        Args:

        Returns:
            integral_time (float): Time autocorrelation value for the dataset.
        """
        self.autocorrelate_results = {}
        for bin in self.data.keys():
            data = self.data[bin][values][EquilTime:]
            try:
                integral_time = emcee.autocorr.integrated_time(data, c=1)
            except emcee.autocorr.AutocorrError:
                integral_time = 1
            except IndexError:
                integral_time = 1    
            self.autocorrelate_results[bin] = {"data": data,
                                          "integral_time": float(integral_time) }
    def ToggleError(self, bin, Issue=None):
        self.data[bin]["ERROR"] = True
        if Issue is not None:
            self.data[bin]["ErrorMSG"] = Issue
    def wham_init(self, WhamLocation:str, convergence=1e-6,vis="false"):
        if self.colvar.VariableType != "length":
            periodicity = "P"
        else:
            periodicity = ""
        if os.path.isdir(WhamLocation) == False:
            os.mkdir(WhamLocation)
        metapath = os.path.join(WhamLocation, "metadata")
        if os.path.isdir(metapath) == False:
            os.mkdir(metapath)
        UseableBins = 0
        metafiles = []
        for bin in self.data.keys():
            if self.data[bin]["ERROR"] == False:
                UseableBins += 1
                meta = [str]*(len(self.autocorrelate_results[bin]["data"])+1)
                meta[0] = f"#Step\tValue"
                for i, val in enumerate(self.autocorrelate_results[bin]["data"]):
                    meta[i+1] = f"{i}\t{val}"
                metafile = f"{bin}.metadata.dat"
                io.textDump(meta, os.path.join(metapath,metafile))
                metafiles.append(f"{os.path.join(metapath,metafile)} {self.data[bin]['Value']} {self.colvar.HoldForce} {self.autocorrelate_results[bin]['integral_time']}")
                plt.hist(self.autocorrelate_results[bin]["data"], 100)
            else:
                pass
        plt.xlabel("Reaction coordinate")
        plt.ylabel("Count")
        plt.savefig(os.path.join(WhamLocation, "WhamCoverage.eps"))
        if vis == "true":
            plt.show()
        io.textDump(metafiles, "meta_locations.dat")
        
        if self.colvar.Min > self.colvar.Max:
            colvar_low = self.colvar.Max
            colvar_high = self.colvar.Min
        else:
            colvar_low = self.colvar.Min
            colvar_high = self.colvar.Max

        whamfile = f"""wham {periodicity} {colvar_low} {colvar_high} {UseableBins} {convergence} {self.Temperature} 0 "meta_locations.dat" {os.path.join(metapath, 'wham.pmf')} 10 60
sed '1d' {os.path.join(metapath, "wham.pmf")} | awk '{"{"}print $1,"",$2{"}"}' > {os.path.join(metapath, "plot_free_energy.dat")}
        """
        io.textDump(whamfile,  "wham.sh")
        self.whamscript = "wham.sh"
    def wham_run(self):
        try:
            wham_out = subprocess.run(f"sh {self.whamscript}", shell=True, capture_output=True)
        except TypeError:
            wham_out = subprocess.run(f"sh {self.whamscript}", shell=True, stdout=subprocess.PIPE)
        io.textDump(wham_out.stdout.decode(), self.whamscript.replace(".sh", ".out"))
        if "wham.sh" in wham_out.stderr.decode():
            raise Exception(f"Problem with running wham: {wham_out.stderr.decode()}")
        if "No such file " in wham_out.stdout.decode():
            print(wham_out.stdout.decode())
            raise Exception("metadata file not found when runnning wham.sh")
    def init_directories(self, path:str):
        for i in self.data.keys():
            self.data[i]["path"] = os.path.join(path, str(i))
            if os.path.isdir(self.data[i]["path"]) == False:
                os.mkdir(self.data[i]["path"])
    def pull_init(self,WorkDir:str, MM: MMClass, job:dict):
        for key in self.data.keys():
            bin = self.data[key]
            if bin["PreviousWindow"] != None:
                infile = os.path.join(f"../{bin['PreviousWindow']}", job['output'])
            else:
                infile = os.path.join("../", job["input"])
            windowPath = os.path.join(WorkDir, str(key))
            
            file = MM.SMD(infile, job['output'], "pull.colvar.conf", job["steps"], 
                            job["timestep"],job["trajout"], job["temperature"], job["pressure"], job["seed"] )
            io.textDump(file, os.path.join(windowPath, f"{job['output']}.conf"))
            baseColvarFile = self.colvar.VariableLines
            baseColvarFile += f"""
harmonic {"{"}
name            {self.colvar.VariableType}Potential
colvars         {self.colvar.VariableType}
centers         {bin["Value"]}
forceConstant   {self.colvar.PullForce}
{"}"}
"""
            io.textDump(baseColvarFile, os.path.join(windowPath, "pull.colvar.conf"))
    def pull_run(self, WorkDir:str, MM: MMClass, job:dict, VMD: VMDClass, Trackers:list ):
        runscript = """#!/bin/bash 
mkdir /dev/shm/RUNDIR
"""
        if MM.software.config["qmForces"] == "on":
            MMPath = MM.software.path_cpu
            CommandLines = "+setcpuaffinity"
            GPU=False
        else:
            MMPath = MM.software.path_gpu
            CommandLines = "+oneWthPerCore +setcpuaffinity +devices 0"
            GPU=True
        if MM.qm == "pyscf" and job["run"].casefold() == "true":
            if os.path.isdir("/dev/shm/RUNDIR") == False:
                os.mkdir("/dev/shm/RUNDIR")
                os.mkdir("/dev/shm/RUNDIR/0")
            print("INFO: Starting sniffer")
            sniffer = os.system("cd /dev/shm/RUNDIR/0 ; qmmm_sniffer > sniffer.out &")
        for i in range(self.start_index,self.colvar.nsteps):
            runscript += f"cd {i} ; {MMPath} {CommandLines} {job['output']}.conf > {job['output']}.out ; cd ../ \n"
            if job["run"].casefold() == "true":
                filepath = os.path.join(WorkDir, str(i), f"{job['output']}")
                if MM.software.check_output(f"{filepath}.out")[0] != "completed":
                    if os.path.isdir("/dev/shm/RUNDIR") == False:
                        os.mkdir("/dev/shm/RUNDIR")
                    
                    _ = MM.software.exec(f"{filepath}.conf", f"{filepath}.out", GPU)
                    
                    status, _ = MM.software.check_output(f"{filepath}.out")
                    if status != "completed":
                        raise RuntimeError(f"ERROR: pull has had an issue. Status = {status}. Please check the output file: {filepath}.out")
                    TrackerFile = VMD.GenAnalysisScript([os.path.join(str(i),job['output'])])
                    io.textDump(TrackerFile, os.path.join(os.path.join(WorkDir, str(i)),f"Analysis.tcl"))
                    VMD.RunAnalysis(os.path.join(os.path.join(WorkDir, str(i)),f"Analysis.tcl"))
                    for Tracker in Trackers:
                        Tracker.get_vmdData(WorkDir, job['output'], i)
        for i in range(self.start_index, 0, -1):
            runscript += f"cd {i-1} ; {MMPath} {CommandLines} {job['output']}.conf > {job['output']}.out ; cd ../ \n"
            if job["run"].casefold() == "true":
                filepath = os.path.join(WorkDir, str(i-1), f"{job['output']}")
                if MM.software.check_output(f"{filepath}.out")[0] != "completed":
                    if os.path.isdir("/dev/shm/RUNDIR") == False:
                        os.mkdir("/dev/shm/RUNDIR")
                    _ = MM.software.exec(f"{filepath}.conf", f"{filepath}.out", GPU)
                    status, _ = MM.software.check_output(f"{filepath}.out")
                    if status != "completed":
                        raise RuntimeError(f"ERROR: pull has had an issue. Status = {status}. Please check the output file: {filepath}.out")
                    TrackerFile = VMD.GenAnalysisScript([os.path.join(str(i-1),job['output'])])
                    io.textDump(TrackerFile, os.path.join(os.path.join(WorkDir, str(i-1)),f"Analysis.tcl"))
                    VMD.RunAnalysis(os.path.join(os.path.join(WorkDir, str(i-1)),f"Analysis.tcl"))
                    for Tracker in Trackers:
                        Tracker.get_vmdData(WorkDir, job['output'], i-1)
        if MM.qm == "pyscf" and job["run"].casefold() == "true":
            io.textDump("","/dev/shm/RUNDIR/0/kill" )
        runscript += "rm -r /dev/shm/RUNDIR"
        io.textDump(runscript, os.path.join(WorkDir, "pull.sh"))
        if job["run"].casefold() == "true":
            if os.path.isdir("/dev/shm/RUNDIR"):
                subprocess.run(["rm", "-r", "/dev/shm/RUNDIR"],shell=True)
        return Trackers
    def hold_init(self, WorkDir:str, MM:MMClass, job:dict, HPC: HPCClass):
        for key in self.data.keys():
            bin = self.data[key]
            binpath = os.path.join(WorkDir,str(key))
            assert os.path.isdir(binpath)
            if HPC.partition == False or job["steps"] <= HPC.max_steps:
                file = MM.SMD(job["input"], job['output'], "hold.colvar.conf", job["steps"], 
                            job["timestep"],job["trajout"], job["temperature"], job["pressure"] )
                io.textDump(file, os.path.join(binpath, f"{job['output']}.conf"))
                files = [f"{job['output']}.conf"]
            else:
                files = []
                file = MM.SMD(job["input"], f"{job['output']}_1", "hold.colvar.conf", HPC.max_steps, 
                        job["timestep"],job["trajout"], job["temperature"], job["pressure"] )
                io.textDump(file, os.path.join(binpath, f"{job['output']}_1.conf"))
                files.append(f"{job['output']}_1.conf")
                if job["steps"]//HPC.max_steps != 1:
                    for i in range(1,job["steps"]//HPC.max_steps):
                        file = MM.SMD(f"{job['output']}_{i}", f"{job['output']}_{i+1}", "hold.colvar.conf", HPC.max_steps, 
                            job["timestep"],job["trajout"], job["temperature"], job["pressure"] )
                        io.textDump(file, os.path.join(binpath, f"{job['output']}_{i+1}.conf"))
                        files.append(f"{job['output']}_{i+1}.conf")
                else:
                    i = 0
                if job["steps"]%HPC.max_steps > 0:
                    file = MM.SMD(f"{job['output']}_{i+1}", f"{job['output']}_{i+2}", "hold.colvar.conf",job["steps"]%HPC.max_steps, 
                        job["timestep"],job["trajout"], job["temperature"], job["pressure"] )
                    io.textDump(file, os.path.join(binpath, f"{job['output']}_{i+2}.conf"))
                    files.append(f"{job['output']}_{i+2}.conf")            
            baseColvarFile = self.colvar.VariableLines
            baseColvarFile += f"""
harmonic {"{"}
name            {self.colvar.VariableType}Potential
colvars         {self.colvar.VariableType}
centers         {bin["Value"]}
forceConstant   {self.colvar.HoldForce}
{"}"}
"""
            
            io.textDump(baseColvarFile, os.path.join(binpath, "hold.colvar.conf"))
            print(files)
        return files

    def hold_run(self, WorkDir:str, MM: MMClass, job:dict, 
    HPC:HPCClass, NPartitions:int):
        if MM.software.config["qmForces"] == "on":
            MMPath = MM.software.path_cpu
            CommandLines = "+setcpuaffinity"
            GPU=False
        else:
            MMPath = MM.software.path_gpu
            CommandLines = "+oneWthPerCore +setcpuaffinity +devices 0"
            GPU=True
        if NPartitions == 1:
            if HPC.exists:
                runscript = ""
            else:
                runscript = """#!/bin/bash 
mkdir /dev/shm/RUNDIR
mkdir /dev/shm/RUNDIR/0
"""
                if MM.qm == "pyscf":
                    runscript += """cd /dev/shm/RUNDIR/0
nohup qmmm_sniffer & disown
"""
            
            for bin in self.data.keys():
                if MM.software.check_output(f"{bin}/{job['output']}.out")[0] != "completed":
                    if HPC.exists:
                        if MM.qm == "pyscf":
                            runscript += f"cd {bin} ; sed -i \"s/RUNDIR/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID/g\" {job['output']}.conf ; mkdir /dev/shm/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID; mkdir /dev/shm/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID/0 ; cd /dev/shm/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID/0 ; nohup qmmm_sniffer & disown ; cd - ; {MMPath} {CommandLines} {job['output']}.conf > {job['output']}.out ; cd ../ ; touch /dev/shm/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID/0/kill ; sleep 10 ; rm -r /dev/shm/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID ;\n"
                        else:
                            runscript += f"cd {bin} ; sed -i \"s/RUNDIR/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID/g\" {job['output']}.conf ; mkdir /dev/shm/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID ; {MMPath} {CommandLines} {job['output']}.conf > {job['output']}.out ; cd ../ ; rm -r /dev/shm/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID ;\n"
                    else:
                        runscript += f"cd {bin} ; {MMPath} {CommandLines} {job['output']}.conf > {job['output']}.out ; cd ../ \n"
            if HPC.exists == False:
                if MM.qm == "pyscf":
                    runscript += "touch /dev/shm/RUNDIR/0/kill"
                    runscript += "sleep 10"
                runscript += "rm -r /dev/shm/RUNDIR"

            io.textDump(runscript, os.path.join(WorkDir, f"Umbrella-{job['output']}.sh"))
            if HPC.exists:
                slurmscript = HPC.gen_slumScript("array-job", job['output'], os.path.join(WorkDir, f"Umbrella-{job['output']}.sh"), len(self.data.keys()), True)
                io.textDump(slurmscript, os.path.join(WorkDir, f"sub-Umbrella-{job['output']}.sh"))
                io.textDump(HPC.arrayjobscript, os.path.join(WorkDir, "array_job.sh"))
                if job["run"] == "true":
                    filepath = os.path.join(WorkDir, str(bin), f"{job['output']}")
                    if MM.software.check_output(f"{filepath}.out")[0] != "completed":
                        status = HPC.check_dependecy(job['output'])
                        if status != "wait":
                            HPC.run_slurmScript(os.path.join(WorkDir, f"sub-Umbrella-{job['output']}.sh"))
                        else:
                            print(f"INFO: {job['output']} job already in the queue, skipping.")
                    else:
                        print(f"INFO:  {job['output']} job has already finished.")
        else:
            for i in range(NPartitions):
                if HPC.exists:
                    runscript = ""
                else:
                    runscript = """#!/bin/bash 
mkdir /dev/shm/RUNDIR
        """
                
                for bin in self.data.keys():
                    if MM.software.check_output(f"{bin}/{job['output']}_{i+1}.out")[0] != "completed":
                        if HPC.exists:
                            if MM.qm == "pyscf":
                                runscript += f"cd {bin} ; sed -i \"s/RUNDIR/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID/g\" {job['output']}_{i+1}.conf ; mkdir /dev/shm/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID ; mkdir /dev/shm/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID/0 ; cd /dev/shm/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID/0  ; nohup qmmm_sniffer & disown ; cd - ; {MMPath} {CommandLines} {job['output']}_{i+1}.conf > {job['output']}_{i+1}.out ; cd ../ ; touch /dev/shm/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID/0/kill ; sleep 10 ; rm -r /dev/shm/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID ;\n"

                            else:
                                runscript += f"cd {bin} ; sed -i \"s/RUNDIR/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID/g\" {job['output']}_{i+1}.conf ; mkdir /dev/shm/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID ; {MMPath} {CommandLines} {job['output']}_{i+1}.conf > {job['output']}_{i+1}.out ; cd ../ ; rm -r /dev/shm/$SLURM_JOB_ID-$SLURM_ARRAY_TASK_ID ;\n"
                        else:
                            runscript += f"cd {bin} ; {MMPath} {CommandLines} {job['output']}_{i+1}.conf > {job['output']}_{i+1}.out ; cd ../ \n"
                if HPC.exists == False:
                    if MM.qm == "pyscf":
                        runscript += "touch /dev/shm/RUNDIR/0/kill"
                        runscript += "sleep 10"
                    runscript += "rm -r /dev/shm/RUNDIR"

                io.textDump(runscript, os.path.join(WorkDir, f"Umbrella-{job['output']}_{i+1}.sh"))
                if HPC.exists:
                    slurmscript = HPC.gen_slumScript("array-job", f"{job['output']}_{i+1}", os.path.join(WorkDir, f"Umbrella-{job['output']}_{i+1}.sh"), len(self.data.keys()), True)
                    io.textDump(slurmscript, os.path.join(WorkDir, f"sub-Umbrella-{job['output']}_{i+1}.sh"))
                    io.textDump(HPC.arrayjobscript, os.path.join(WorkDir, "array_job.sh"))
                    if job["run"] == "true":
                        filepath = os.path.join(WorkDir, str(i), f"{job['output']}_{i+1}")
                        if MM.software.check_output(f"{filepath}.out")[0] != "completed":
                            status = HPC.check_dependecy(f"{job['output']}_{i+1}")
                            if status != "wait":
                                HPC.run_slurmScript(os.path.join(WorkDir, f"sub-Umbrella-{job['output']}_{i+1}.sh"))
                            else:
                                print(f"INFO: {job['output']}_{i+1} job already in the queue, skipping.")
                        else:
                            print(f"INFO:  {job['output']}_{i+1} job has already finished.")
    def analyse_completed(self,WorkDir:str,Files:list, MM:MMClass,VMD:VMDClass, job:dict, Trackers:list):
        for bin in self.data.keys():
            print(bin)
            bindir = os.path.join(WorkDir, str(bin))
            track = []
            for file in Files:
                outfile = file.replace(".conf","")
                filepath = os.path.join(bindir, outfile)
                status, _ = MM.software.check_output(f"{filepath}.out")
                print(status)
                if status == "completed" or status == "running":
                    data = io.textRead(f"{filepath}.colvars.traj")
                    if "umbrella-equil" in job.keys():
                        self.add_data(bin, "equilValues", data[1:], (len(data[1:])-1)*float(MM.software.config["timeStep"]))
                    elif "umbrella-prod" in job.keys():
                        self.add_data(bin, "prodValues", data[1:], (len(data[1:])-1)*float(MM.software.config["timeStep"]))
                if status == "completed":
                    track.append(filepath)
            TrackerFile = VMD.GenAnalysisScript(track)
            io.textDump(TrackerFile, os.path.join(bindir,f"Analysis.tcl"))
            VMD.RunAnalysis(os.path.join(bindir,f"Analysis.tcl"))
            for Tracker in Trackers:
                Tracker.get_vmdData(WorkDir, job['output'],bin)
        return Trackers
    def dump_data(self, path:str):
        data = self.data
        io.jsonDump(data, path)