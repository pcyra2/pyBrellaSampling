# Usage

## Setup

To run a standalone simulation, you need some key initial files:

* `Parameter file` (default = complex.parm7) Can be changed using `pyBrella --ParmFile PARM`
* `Coordinate file` (default = start.rst7) Can be changed using `pyBrella --StartFile COORD`
* `Amber Coordinate file` (default = start.rst7) Can be changed using `pyBrella --AmberCoordinates COORD`, This needs to be an amber file to link with amber parameters

You also have [cli variables](#user-inputs) that tune the calculation. They are handled by the [input parser](../CodeReference/#pyBrellaSampling.Tools.InputParser.StandaloneInput)

## Running (recommended)

To run a Standalone calculation, we recommend generating an input file, defining any variables that differ from those in the [default variables](#user-inputs)

The code can then be executed using:
`standalone -i INPUTFILE`

INPUTFILE must be in the format:
    
    # Lines starting with '#' will be ignored
    # All inputs should be on a new line
    KEYWORD1=VARIABLE1
    KEYWORD2=VARIABLE2

## Running (not recommended)

You can also run a simulation directly from the command line. All user variables can be parsed and this is handled by the [input parser](../CodeReference/#pyBrellaSampling.Tools.InputParser.StandaloneInput)

## User Inputs
The default inputs can be found and edited in `UserVars/Defaultinputs.py`, 


Here is a list of the main variables that can be changed:

#### General
    "WorkDir" : "./",# str: relative or absolute path to calculation
    "Verbosity" : 0, # int: 0=Errors, 1=Warnings, 2=Info, 3=Debug
    "DryRun" : True, # bool: Whether to actually run calculations
#### Calculation
    "CoresPerJob" : 4, # int: Number of CPU cores for the calc.
    "MemoryPerJob" : 4, # int: Ammount of RAM for the calc.
    "Name" : "Name", # str: Name of calculation.
    "Ensemble" : "min", # str: min, heat, NVT, NVP
#### QM
    "QM" : False, # bool: To select the use of qmmm, requires the use of a qm file. 
    "QmFile" : "None", # str:  Name of qm file containing qm parms. 
    "QmSelection" : "resname CTN POP MG", # str: vmd selection algebra for the QM region
    "QmCharge" : 1, # int: Net charge of QM system
    "QmSpin" : 0, # int: Spin state of QM system
    "QmMethod" : "PBE", # str: ORCA QM method
    "QmBasis" : "6-31+G*", # str: ORCA QM basis set.
    "QmArgs" : "MINIPRINT D3BJ TightSCF CFLOAT", # str: extra information to parse to ORCA
#### MM
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
