# User Inputs

User inputs can either be provided by [commandline variables](../CodeReference/#pyBrellaSampling.Tools.InputParser.UmbrellaInput), or through the use of one or two inputfiles. 

We recommend the use of two input files, on being the QM input file that should contain the information regarding the QM zone. The other being the Umbrella input file that will contain all other information. 

This is because the QM zone should never change, so you can generate this once and retain it for all steps. The other variables may change between stages of the calculation however and so it is useful to separate these out into a different file. 

You can however use only one input file as either file can read in all user variables. It is important to note however that the Umbrella input file overwrites the QM input file, and that both input files overwrite any parsed commandline variables. 

## User Inputs
The default inputs can be found and edited in `UserVars/Defaultinputs.py`

Here is a list of the main variables that can be changed:

### General               
    "WorkDir" : "./",# str: relative or absolute path to calculation
    "Verbosity" : 0, # int: 0=Errors, 1=Warnings, 2=Info, 3=Debug
    "DryRun" : True, # bool: Whether to actually run calculations
#### Calcuation
    "CoresPerJob" : 10, # int: Number of CPU cores for the calc.
    "MemoryPerJob" : 10, # int: Ammount of RAM for the calc.
    "MaxStepsPerCalc" : 1000, # int: Maximum number of steps per calculation
#### QM
    "QM" : False, # bool: To select the use of qmmm, requires the use of a qm file. 
    "QmFile" : "None", # str:  Name of qm file containing qm parms. 
    "QmSelection" : "resname CTN POP MG", # str: vmd selection algebra for the QM region
    "QmCharge" : 1, # int: Net charge of QM system
    "QmSpin" : 0, # int: Spin state of QM system
    "QmMethod" : "PBE", # str: ORCA QM method
    "QmBasis" : "6-31+G*", # str: ORCA QM basis set.
    "QmArgs" : "MINIPRINT D3BJ TightSCF CFLOAT", # str: extra information to parse to ORCA
#### Umbrella
    "Stage" : "inpfile", # str: Stage of calculation
    "UmbrellaFile" : "None", # str: file containing qm information.
    "UmbrellaMin" : 1.3, # float: Minimum Umbrella distance.
    "UmbrellaWidth" : 0.05, # float: Width of the umbrella bin in either Ang or Deg.
    "UmbrellaBins" : 54, # int: Number of Umbrella bins. 
    "PullForce" : 5000, # float: Force to pull the variables to each bin
    "ConstForce" : 300, # float: Force to restrain the colvar during Umbrella simulations
    "StartDistance" : 1.4, # float: Current value for the colvar (From the equilibrated MD)
    "StartFile" : "start.rst7", # str: Amber or NAMD coordinates, start coordinates for THIS simulation.
    "ParmFile" : "complex.parm7", # str: Amber parameter file. This must be in the root directory!
    "AtomMask" : "0,0,0,0", # str: Atom mask for SMD, comma delimited string of atoms. 
    "AnalysisFile" : "prod", # str: Variable to perform analysis on (i.e. if you only want to visualise pull files, use pull_1)
    "EquilLength" : 1, # int: Length of equil in ps. eg 1 ps = 2000 steps at 0.5 fs timestep
    "ProdLength" : 4, # int: Length of production in ps. 4 ps = 8000 steps at 0.5 fs timestep

### Software paths

Software paths should be set up before running any calculations and these can be found in `UserVars/SoftwarePaths.py`
#### ORCA path
    "ORCA_PATH" : "/PATH_TO_ORCA/5.0.4/orca" # str: Path to the orca executable

#### NAMD paths
    "NAMD_CPU" : "/PATH_TO_NAMD/NAMD_3.0b4_Linux-x86_64-multicore/namd3" # str: Path to the cpu version of NAMD for QMMM
    "NAMD_GPU" : "/PATH_TO_NAMD/NAMD_3.0b4_Linux-x86_64-multicore-CUDA/namd3" # str: Path to the GPU version of NAMD for MD

