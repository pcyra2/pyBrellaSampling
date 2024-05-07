Umbrella_Inp = {
###--------------------------General---------------------------------###                
    "WorkDir" : "./",# str: relative or absolute path to calculation
    "Verbosity" : 0, # int: 0=Errors, 1=Warnings, 2=Info, 3=Debug
    "DryRun" : True, # bool: Whether to actually run calculations
###--------------------------Calcuation------------------------------### 
    "CoresPerJob" : 10, # int: Number of CPU cores for the calc.
    "MemoryPerJob" : 10, # int: Ammount of RAM for the calc.
    "MaxStepsPerCalc" : 1000, # int: Maximum number of steps per calculation
###--------------------------QM--------------------------------------### 
    "QM" : False, # bool: To select the use of qmmm, requires the use of a qm file. 
    "QmFile" : "None", # str:  Name of qm file containing qm parms. 
    "QmSelection" : "resname CTN POP MG", # str: vmd selection algebra for the QM region
    "QmCharge" : 1, # int: Net charge of QM system
    "QmSpin" : 0, # int: Spin state of QM system
    "QmMethod" : "PBE", # str: ORCA QM method
    "QmBasis" : "6-31+G*", # str: ORCA QM basis set.
    "QmArgs" : "MINIPRINT D3BJ TightSCF CFLOAT", # str: extra information to parse to ORCA
###--------------------------Umbrella--------------------------------###
    "Stage" : "inpfile", # str: Stage of calculation
    "UmbrellaFile" : "None", # str: file containing qm information.
    "UmbrellaMin" : 1.3, # float: Minimum Umbrella distance.
    "UmbrellaWidth" : 0.05, # float: Width of the umbrella bin in either Ang or Deg.
    "UmbrellaBins" : 54, # int: Number of Umbrella bins. 
    "PullForce" : 5000, # float: Force to pull the variables to each bin
    "ConstForce" : 300, # float: Force to restain the colvar during Umbrella simulations
    "StartDistance" : 1.4, # float: Current value for the colvar (From the equilibrated MD)
    "StartFile" : "start.rst7", # str: Amber or NAMD coordinates, start coordinates for THIS simulation.
    "ParmFile" : "complex.parm7", # str: Amber parameter file. This must be in the root directory!
    "AtomMask" : "0,0,0,0", # str: Atom mask for SMD, comma delimited string of atoms. 
    "AnalysisFile" : "prod", # str: Variable to perform analysis on (i.e. if you only want to visualise pull files, use pull_1)
    "EquilLength" : 1, # int: Length of equil in ps. eg 1 ps = 2000 steps at 0.5 fs timestep
    "ProdLength" : 4, # int: Length of production in ps. 4 ps = 8000 steps at 0.5 fs timestep
    }

Standalone_Inp = {
    "WorkDir" : "./",# str: relative or absolute path to calculation
    "Verbosity" : 0, # int: 0=Errors, 1=Warnings, 2=Info, 3=Debug
    "DryRun" : True, # bool: Whether to actually run calculations
    "CoresPerJob" : 4, # int: Number of CPU cores for the calc.
    "MemoryPerJob" : 4, # int: Ammount of RAM for the calc.
    "Name" : "Name", # str: Name of calculation.
    "Ensemble" : "min", # str: min, heat, NVT, NVP
    "QM" : False, # bool: To select the use of qmmm, requires the use of a qm file. 
    "QmFile" : "None", # str:  Name of qm file containing qm parms. 
    "QmSelection" : "resname CTN POP MG", # str: vmd selection algebra for the QM region
    "QmCharge" : 1, # int: Net charge of QM system
    "QmSpin" : 0, # int: Spin state of QM system
    "QmMethod" : "PBE", # str: ORCA QM method
    "QmBasis" : "6-31+G*", # str: ORCA QM basis set.
    "QmArgs" : "MINIPRINT D3BJ TightSCF CFLOAT", # str: extra information to parse to ORCA
    "MMFile" : "None", # str: file containing qm information.
    "Steps" : 1000, # int: Number of steps to simulate.
    "TimeStep" : 2, # float: timestep in fs. if >1, shake will be turned on. 
    "ParmFile" : "complex.parm7", # str: parameter file. Currently only amber params supported.
    "AmberCoordinates" : "start.rst7", # str: amber coordinate file for use with the parameter files.
    "StartFile" : "start.rst7", # str: Amber or NAMD coordinates, start coordinates for THIS simulation.
    "RestartOut" : 10, # int: frequency of updating the restart file.
    "TrajOut" : 50, # int: frequency of updating the trajectory file.
    "SMD" : False, # bool: Whether to use Steered MD.
    "SMDFile" : None, # str: Name of file containing SMD parms
    "Force" : 0, # float: Force to perform SMD
    "StartValue" : 0, # float: start value for SMD
    "EndValue" : 0, # float: end value for SMD
    "AtomMask" : "0,0,0,0", # str: Atom mask for SMD, comma delimited string of atoms. 
    "HPC" : False, # bool: Run using a SLURM Scheduler. 
}

Benchmark_Inp = {
    "Verbosity" : 0, # int: 0=Errors, 1=Warnings, 2=Info, 3=Debug
    "WorkDir" : "./", # str: Path to calculation directory.
    "Cores" : 1, # int: Number of CPU cores
    "DryRun" : True, # bool: Whether to actually run calculations
    "Stage" : "analysis", # str: init, calc or analysis
    "BenchmarkType" : "Energy", # str: only Energy supported at current.
    "CoordinateLoc" : "./Coordinates", # str: location of molecule coordinates
    "ReactionList" : "reactions.dat", # str: List of 1 step reactions
    "SCF" : "TIGHTSCF", # str: Convergence criteria (As supported by ORCA)
    "Grid" : "DEFGRID3", # str: SCF Grid (As supported by ORCA)
    "Restart" : "NOAUTOSTART", # str: Whether to use the .gbw restart files (AUTOSTART/NOAUTOSTART)
    "Convergence" : "NormalConv", # str: Agressivness of the SCF optimizer. 
    "Extras" : "MINIPRINT", # str: Extra information to parse to ORCA.
}