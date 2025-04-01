import pyBrellaSampling.Code.Tools.UserConfig.CheckConfig as ConfigInit
import pyBrellaSampling.Code.Tools.InputParser as InputParser 
import pyBrellaSampling.Code.Tools.classes as classes
import pyBrellaSampling.Code.Tools.io as io
from pyBrellaSampling.Code.Tools.MD.Umbrella import UmbrellaClass

import time
import os
from pprint import pprint

def main():
    GlobStart = time.perf_counter()
    DefInps = ConfigInit.GetDefaults()
    Inputs = InputParser.ParseInputs(DefInps)
    print(f"INFO: Input file read, starting calculation, Verbosity is set to {Inputs['verbosity']}" if Inputs["verbosity"]>2 else "")
    assert os.path.isdir(Inputs["workdir"]), f"ERROR: WorkDirectory does not exist: {Inputs['workdir']}"
    MM = classes.MMClass("namd", Inputs["parameters"], Inputs["topology"])
    HPC = classes.HPCClass(Inputs["hpc"]["hostname"])
    if "max_steps" in Inputs["hpc"]:
        HPC.walltime_partition(Inputs["hpc"]["max_steps"])
    if "config" in Inputs["hpc"]:
        HPC.init_slurm(Inputs["hpc"]["config"], Inputs["hpc"]["modulefiles"], Inputs["hpc"]["environment"])
    partitioned = False
    tracker_inf = Inputs["tracker"]
    trackers = [None]*len(tracker_inf.keys())
    Analysis_lines = ""
    for index, key  in enumerate(tracker_inf.keys()):
        trackers[index] = classes.TrackerClass(tracker_inf[key]["atoms"], key)
        trackers[index].read_previous(Inputs["workdir"])
        Analysis_lines += trackers[index].gen_vmdScript() 

    VMD = classes.VMDClass(Inputs["parameters"],Analysis_lines)

    ### Initial minimization
    if "minimize" in Inputs["jobs"]:
        job = Inputs["jobs"]["minimize"]
        file = MM.minimize(job["input"], job["output"], job["steps"])
        io.textDump(file, os.path.join(Inputs["workdir"],f"{job["output"]}.conf"))
        if job["run"].casefold() == "true":
            if MM.software.check_output(str(os.path.join(Inputs["workdir"],f"{job["output"]}.out")))[0] != "completed":
                _ = MM.software.exec(os.path.join(Inputs["workdir"],f"{job["output"]}.conf"), os.path.join(Inputs["workdir"],f"{job["output"]}.out"), Inputs["gpu"])
                status, _ = MM.software.check_output(os.path.join(Inputs["workdir"],f"{job["output"]}.out")) 
                if status != "completed":
                    raise RuntimeError(f"ERROR: Minimization has had an issue. Status = {status}. Please check the output file: {os.path.join(Inputs["workdir"],f"{job["output"]}.out")}")
                TrackerFile = VMD.GenAnalysisScript([job["output"]])
                io.textDump(TrackerFile, os.path.join(Inputs["workdir"],f"Analysis.tcl"))
                VMD.RunAnalysis(os.path.join(Inputs["workdir"],f"Analysis.tcl"))
                for Tracker in trackers:
                    Tracker.get_vmdData(Inputs["workdir"], job["output"])
            else:
                print("INFO: Minimization has already been performed, skipping." if Inputs["verbosity"]>2 else "")
    
    ### MD heating
    if "heat" in Inputs["jobs"]:
        job = Inputs["jobs"]["heat"]
        file = MM.heat(job["input"], job["output"], job["steps"], job["timestep"], 0, job["target-temp"])
        io.textDump(file, os.path.join(Inputs["workdir"], f"{job["output"]}.conf"))
        if job["run"].casefold() == "true":
            if MM.software.check_output(os.path.join(Inputs["workdir"],f"{job["output"]}.out"))[0] != "completed":
                _ = MM.software.exec(os.path.join(Inputs["workdir"],f"{job["output"]}.conf"), os.path.join(Inputs["workdir"],f"{job["output"]}.out"), Inputs["gpu"])
                status, _ = MM.software.check_output(os.path.join(Inputs["workdir"],f"{job["output"]}.out")) 
                if status != "completed":
                    raise RuntimeError(f"ERROR: Heating has had an issue. Status = {status}. Please check the output file: {os.path.join(Inputs["workdir"],f"{job["output"]}.out")}")
                TrackerFile = VMD.GenAnalysisScript([job["output"]])
                io.textDump(TrackerFile, os.path.join(Inputs["workdir"],f"Analysis.tcl"))
                VMD.RunAnalysis(os.path.join(Inputs["workdir"],f"Analysis.tcl"))
                for Tracker in trackers:
                    Tracker.get_vmdData(Inputs["workdir"], job["output"])
            else:
                print("INFO: Heating job already completed, skipping." if Inputs["verbosity"]>2 else "")
    
    ### MD Equilibration
    if "md-equil" in Inputs["jobs"]:
        job = Inputs["jobs"]["md-equil"]
        if "pressure" in job:
            pressure = job["pressure"]
        else:
            pressure = None
        file = MM.constant(job["input"], job["output"], job["steps"], 
                           job["timestep"],job["trajout"], job["temperature"], pressure)
        io.textDump(file, os.path.join(Inputs["workdir"], f"{job["output"]}.conf"))
        if job["run"].casefold() == "true":
            if MM.software.check_output(os.path.join(Inputs["workdir"],f"{job["output"]}.out"))[0] != "completed":
                _ = MM.software.exec(os.path.join(Inputs["workdir"],f"{job["output"]}.conf"), os.path.join(Inputs["workdir"],f"{job["output"]}.out"), Inputs["gpu"])
                status, _ = MM.software.check_output(os.path.join(Inputs["workdir"],f"{job["output"]}.out")) 
                if status != "completed":
                    raise RuntimeError(f"ERROR: md-equilibration has had an issue. Status = {status}. Please check the output file: {os.path.join(Inputs["workdir"],f"{job["output"]}.out")}")
                TrackerFile = VMD.GenAnalysisScript([job["output"]])
                io.textDump(TrackerFile, os.path.join(Inputs["workdir"],f"Analysis.tcl"))
                VMD.RunAnalysis(os.path.join(Inputs["workdir"],f"Analysis.tcl"))
                for Tracker in trackers:
                    Tracker.get_vmdData(Inputs["workdir"], job["output"])
            else:
                print("INFO: Equilibration job already completed, skipping." if Inputs["verbosity"]>2 else "")

    ### Setup Umbrella sampling stuff (QM, colvar etc.)
    if len(Inputs["qm"].keys()) != 0:
        print("initiating QM")
        qmVars = Inputs["qm"]
        if "software" not in qmVars:
            qmVars["software"] = "orca"
        QM = classes.QMClass(qmVars["method"], qmVars["basis"], qmVars["charge"], qmVars["spin"], qmVars["software"])
        QM.init_qmmm(qmVars["vmd_selection"])
        QM.set_cores(qmVars["cores"])
        QM.software.add_extras(qmVars["extras"])
        MM.software.set_qm(QM, "../syst-qm.pdb")
        
        
        VMD.qmPDB_gen(MM)

    colvarVars = Inputs["colvar"]
    colvar = classes.ColvarClass()
    colvar.atomic_colvar(colvarVars["atoms"], colvarVars["pull_force"], colvarVars["constant_force"], colvarVars["minimum"], colvarVars["maximum"], colvarVars["initial"],colvarVars["width"])

    if colvarVars["initial"] == "calculate":
        colvar_analysis_file = colvar.tracker.gen_vmdScript()
        initial = VMD.get_colvardistance(os.path.join(Inputs["workdir"], f"{job["output"]}.dcd"), colvar_analysis_file)
        print(f"INFO: Start distance is: {initial}" if Inputs["verbosity"]>2 else"")
        colvar.update_initial_point(initial)

    VMD.colvarPDB_gen(colvar, MM)
    VMD.GenUmbrellaPDB(Inputs["workdir"])
    # pprint(vars(trackers[0]))

    Umbrella = UmbrellaClass(colvar, job["temperature"])
    Umbrella.init_directories(Inputs["workdir"])

    ### Umbrella pull
    if "pull" in Inputs["jobs"]:
        keys = ["pull", "input", "output", "steps", "timestep", "trajout", "temperature", "pressure", "run", "vis",]
        job = Inputs["jobs"]["pull"]
        job = InputParser.check_keys(job, keys)

        MM.software.set_global(True, os.path.join("../", Inputs["parameters"]), os.path.join("../",Inputs["topology"] ))
        Umbrella.pull_init(Inputs["workdir"],MM, job)
        trackers = Umbrella.pull_run(Inputs["workdir"], MM, job,VMD, trackers)
        for bin in Umbrella.data.keys():
            data = io.textRead(os.path.join(Inputs["workdir"], str(bin), "pull.colvars.traj"))
            Umbrella.add_data(bin, "pullValues", data[1:], (len(data[1:])-1)*float(job["timestep"]))
        if job["vis"] == "true":
            trajfiles = [os.path.join(Inputs["workdir"], str(bin),"pull.dcd") for bin in Umbrella.data.keys()]
            file = VMD.gen_visScript(trajfiles)
            io.textDump(file, os.path.join(Inputs["workdir"], "VisualisePull.tcl"))

    if "umbrella-equil" in Inputs["jobs"]:
        keys = ["umbrella-equil", "input", "output", "steps", "timestep", "trajout", "temperature", "pressure", "run", "vis",]
        job = Inputs["jobs"]["umbrella-equil"]
        job = InputParser.check_keys(job, keys)
        MM.software.set_global(True, os.path.join("../", Inputs["parameters"]), os.path.join("../",Inputs["topology"] ))

        files = Umbrella.hold_init(Inputs["workdir"], MM, job, HPC)
        if len(files) > 1:
            partitioned = True
        else:
            partitioned = False
        Umbrella.hold_run(Inputs["workdir"], MM, job, VMD, trackers, HPC, len(files))

    if "umbrella-prod" in Inputs["jobs"]:
        keys = ["umbrella-prod", "input", "output", "steps", "timestep", "trajout", "temperature", "pressure", "run", "vis",]
        job = Inputs["jobs"]["umbrella-prod"]
        job = InputParser.check_keys(job, keys)

        MM.software.set_global(True, os.path.join("../", Inputs["parameters"]), os.path.join("../",Inputs["topology"] ))
        if partitioned == True:
            job["input"] = f"{job["input"]}_{len(files)}"
            pprint(job)
        files = Umbrella.hold_init(Inputs["workdir"], MM, job, HPC)
        if len(files) > 1:
            partitioned = True
        else:
            partitioned = False
        Umbrella.hold_run(Inputs["workdir"], MM, job, VMD, trackers, HPC, len(files))

    GlobEnd = time.perf_counter()
    print(f"INFO: Total calculation time was {GlobEnd - GlobStart} s" if Inputs["verbosity"] > 2 else "")
    for tracker in trackers:
        tracker.dump(Inputs["workdir"])


if __name__ == "__main__":
    main()