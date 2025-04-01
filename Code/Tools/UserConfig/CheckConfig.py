import os
import pyBrellaSampling.Code.Tools.io as io


HardCodedDefaults = {
    "workdir" : "./",
    "configfile" : "Job.conf",
    "verbosity" : 3,
    "parameters": "system.parm7",
    "topology": "system.rst7",
    "jobs":{},
    "gpu":True,
    "qm":{},
    "colvar":{},
    "tracker":{},
    "hpc":{"hostname":"localhost"}
    }

NAMDDefaults = {
    "DCDfreq": 10,
    "restartfreq": 1,
    "outputTiming": 100,
    "outputEnergies": 1,
    "switching":"off",
    "exclude":"scaled1-4",
    "1-4scaling":"0.833333333",
    "scnb":"2.0",
    "readexclusions":"yes",
    "cutoff":"10",
    "watermodel":"tip3p",
    "pairListdist":"11",
    "LJcorrection":"on",
    "ZeroMomentum":"off",
    "rigidBonds":"all",
    "rigidTolerance":"1.0e-8",
    "rigidIterations":"100",
    "timeStep":"2",
    "fullElectFrequency":"1",
    "nonBondedFreq":"1",
    "stepspercycle":"10",
    "PME":"on",
    "PMEGridSizeX":"300",
    "PMEGridSizeY":"300",
    "PMEGridSizeZ":"300",
    "PMETolerance":"1.0e-6",
    "PMEInterpOrder":"4",
    "wrap":"wrapAll         on",
    "cell":"",
    "cellShape" : "oct",
    "temperature":"0",
    "langevin":"on",
    "qmForces":"off",
    "run":"10",
    "GPUresident":"off",
    "jobtype": "run",
    "langevinPistonTarget": 1.01325,
    "ensemble": "NPT",
    "qmlines": ""
}

ConfigLoc = os.path.dirname(os.path.abspath(__file__))



def GetDefaults()->dict:
    """
    Generates the default config files. Allows for end users to define their own default variables. Currently setup to default to the test system for ease of development.

    Raises:
        Exception: If hard coded defaults generates a different dictionary to user inputs. Allows for internal error detection. 
        Exception: If unsupported basis set type is provided by the user.

    Returns:
        dict: Default input configuration
    """
    try:
        if os.environ["DOCKER"] == "DOCKERENV":
            DefInps = HardCodedDefaults
            DefInps["WorkDir"] = "/app/data/"
            # io.jsonDump(DefInps, ConfigLoc+"/HarnessDefaults.conf")
            # io.jsonDump(DFT, ConfigLoc+"/DFTDefaults.conf")
            return DefInps
    except KeyError:
        pass
    
    if os.path.isfile(ConfigLoc+"/pyBrellaDefaults.conf") == False:
        print("WARNING: Default variables are not set... Initiating.")
        DefInps = {}
        try:
            if os.environ["DOCKER"] == "TestEnv":
                HardDef = "n"
            else:
                HardDef = input("Would you like to configure your own defaults? (y/n) ")
        except KeyError:
            HardDef = input("Would you like to configure your own defaults? (y/n) ")
        if HardDef.casefold() == "y":
            # DefInps[""] = input("")
            DefInps["workdir"] = str(input("What is the path to the default working directory? "))
            DefInps["verbosity"] = int(input("What default verbosity should be used? (0-3 with 0 being the lowest level) "))
            DefInps["configfile"] = str(input("What is the standard input file name? (Standard is ./Job.conf) "))
            DefInps["parameters"] = str(input("What is the default parameter file? (Standard is system.parm7)"))
            DefInps["topology"] = str(input("What is the default topology file? (Standard is system.rst7)"))
            DefInps["jobs"] = {} 
            DefInps["qm"] = {}
            DefInps["colvar"] = {}
            DefInps["tracker"] = {}
            DefInps["hpc"]={"hostname":"localhost"}
        else:
            DefInps = HardCodedDefaults
        assert len(DefInps.keys()) == len(HardCodedDefaults.keys()) # Check that either method produces the same size dictionary...
        io.jsonDump(DefInps, ConfigLoc+"/pyBrellaDefaults.conf")
        if DefInps["verbosity"] > 1:
            print(f"INFO: Default configuration is set and can be changed here: ")
            print(f"{ConfigLoc}/HarnessDefaults.conf")
            print(f"INFO: The default config is: ")
            print(DefInps)
    else:
        DefInps = io.jsonRead(os.path.join(ConfigLoc, "pyBrellaDefaults.conf"))

    if os.path.isfile(os.path.join(ConfigLoc,"/NAMDDefaults.conf")) == False:
        print("WARNING: NAMD default variables are not set... Initiating.")
        print("INFO: Using hard coded defaults. To change or set your own, edit:")
        print(os.path.join(ConfigLoc, "NAMDDefaults.conf"))
        io.jsonDump(NAMDDefaults, os.path.join(ConfigLoc, "NAMDDefaults.conf"))

    path = os.path.join(ConfigLoc, "SoftwareLocations.conf")
    if os.path.isfile(path) == False:
        print("WARNING: Software locations are not set... Initiating.")
        locs = {}
        loc = input("What is the path to orca? ")
        if os.path.isfile(loc):
            locs["orca"] = loc
        else:
            locs["orca"] = None
        locs["namd"] = {}
        loc = input("What is the path to the CPU version of NAMD? ")
        if os.path.isfile(loc):
            locs["namd"]["cpu"] = loc
        else:
            locs["namd"]["cpu"] = None
        loc = input("What is the path to the GPU version of NAMD? ")
        if os.path.isfile(loc):
            locs["namd"]["gpu"] = loc
        else:
            locs["namd"]["gpu"] = None
            DefInps["gpu"] = False
        io.jsonDump(locs, path)
    
    return DefInps
