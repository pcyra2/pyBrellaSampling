# Usage
## Setup

To run an umbrella calculation, you need two key files: 

* `Parameter file` (default = complex.parm7) Can be changed using `pyBrella --ParmFile PARM`
* `Coordinate file` (default = start.rst7) Can be changed using `pyBrella --StartFile COORD`


These two files, along with other user defined [variables](#UserInputs) are all you need in order to perform Umbrella Sampling. 

Users are recommended to follow the workflow provided [below](#workflow).

## Workflow

``` mermaid
graph TB
subgraph Local
 direction LR
 A[Setup] --> B[Minimize]
 B --> C[Heat]
 C --> D[Pull]
 D --> G[Check strcutures]
 G -->|If wrong, fix| D
 end
 Local -->|If correct| HPC
 subgraph HPC
 direction LR
 E[Equilibrate] --> F[Production]
 end
```

### Setup
The setup stage performs the [setup script](CodeReference.md/#pyBrellaSampling.pyBrella.setup)

This generates the `syst-col.pdb` and `syst-qm.pdb` files.

* `syst-col.pdb`  tells NAMD which atoms are involved in the sampling.
* `syst-qm.pdb` tells NAMD which atoms are involved in the QM zone.

#### Important variables:
    "Stage" : Setup
    "QmSelection" : The VMD selection algebra for the QM zone
    "AtomMask" : The atomic index of the atoms involved in Umbrella sampling. 
        - This should be comma delimited, and contain 4 numbers. 
            All unused numbers should be zero e.g if bond: 20,21,0,0

#### Run:
```sh title="Running the setup"
pyBrella -stg setup -dr False 
```

???+ warning "Warning"

    At this point, check that your syst-qm.pdb file has been generated correctly. It is a common bug (depending on your version of vmd) to not generate the column containing elements. (The final column should be the atomic elements.)
    If this does not exist, edit the file around the QM zone, the elements should be in columns 77/78.

### Minimisation

This step minimizes your system using traditional MD. The [minimisation code](CodeReference.md/#pyBrellaSampling.pyBrella.min) is used to perform this.

#### Important variables:
    
    "Stage" : min

#### Run:

```sh title="Running a minimization"
pyBrella -stg min -dr False # (1)
```

1. It is important that the GPU version of NAMD is loaded and that the path has been defined.

### Heating

At this step, we heat the system to 300 k over 20 ps using the [heat code](CodeReference.md/#pyBrellaSampling.pyBrella.heat). 

#### Important variables:
    
    "Stage" : heat

#### Run:

```sh title="Heating the system"
pyBrella -stg heat -dr False
```

???+ warning "Warning"

    Now that the standard MD steps have been performed, you should remove the GPU version of NAMD from your path, and load in the CPU version. You should also ensure that the path to ORCA is correctly defined. QM/MM requires a CPU version of NAMD. 

    You can also start the calculation at this point however need to ensure that your input is named __heat_1.0.restart.coor__ and is in NAMD format.

### Pull

This is the stage where you initiate the umbrella windows/bins. This uses the [pull code](CodeReference.md/#pyBrellaSampling.pyBrella.pull) From this point on it is recommended to be using an input file for all steps of the simulation. 

You first need to analyze the heated structure, and obtain the atom index's for the collective variable.

* `"at1,at2,0,0"` : Bond colvar between atom 1 and atom 2
* `"at1,at2,at3,0"` : Angle colvar between atoms 1, 2 and 3
* `"at1,at2,at3,at4"` : Dihedral colvar between atoms 1, 2, 3 and 4 

#### Important variables:
    
    "Stage" : pull
    "QmSelection" : "VMD selection algebra for qm region"
    "QmCharge" : Net charge of QM region
    "QmSpin" : Spin of QM system
    "QmMethod" : QM functional
    "QmBasis" : QM basis set
    "QmArgs" : Dispersion corrections ect. 
    
    "UmbrellaMin" : Value of bin 0
    "UmbrellaWidth" : Width between each bin (can be negative)
    "UmbrellaBins" : Number of bins
    "PullForce" : Force to pull atoms
    "ConstForce" : Force used during umbrella simulations
    "StartDistance" : Current value in the __heat_1.0.restart.coor__ file
    "AtomMask" : The mask for the collective variable (see above)

#### Run: 

??? tip 

    At this point, you can either run the calculation locally or on a HPC. If you have a fast computer (12 or more cores), it is recommended that you do this locally as the pulls occur in serial and only need around 10 threads. This also makes life faster when checking the generated structures as you do not need to move files between computers.

``` py title="Umbrella.inp"
QmSelection=resname RES
QmCharge=0
QmSpin=1
QmMethod=PBE
QmBasis=6-31+G*
QmArgs=MINIPRINT D3BJ TightSCF
UmbrellaMin=1.3
UmbrellaWidth=0.05 # (1)
UmbrellaBins=50
PullForce=5000 # (2)
ConstForce=300
StartDistance=1.4
AtomMask=12,13,0,0 # (3)
```

1. Small window gap gives better umbrella sampling results. As its a bond, units in Angstroms.
2. Pull force should be large so that the bins are linearly spaced
3. Bond collective variable between atoms 12 and 13

``` sh title="Generating windows"
pyBrella -stg pull -dr False -i Umbrella.inp
```

???+ warning "Warning"

    At this point, you should __ALWAYS__ check the generated structures. You can do this quickly by using the [visualization script](CodeReference.md/#pyBrellaSampling.pyBrella.VisLoad) which is called using :
    
    ``` sh 
    pyBrella -stg vis -dr False -af pull_1
    ```

    This will open up the final structure from each window as a new frame in a single trajectory in vmd. You can then use vmd to plot the collective variable against frame to ensure a linear distribution and that the collective variables are at the correct value.

### Equilibration

Once the structures are all setup, it is recommended to move everything onto a HPC. This will dramatically speed up the future stages of the simulation and make the calculation possible within reasonable time scales. 

Ensure that you have correctly setup the pyBrellaSampling environment on the HPC as the default paths to software ect. will likely be different.

Once on the HPC, using the same `Umbrella.inp` file as used previously, you can initiate the equilibration files.

``` sh title="Setting up equilibration"
pyBrella -stg equil -dr True -i Umbrella.inp
```

This should initialize the equilibration input files for each window, at this point it is recommended that you check at least one of the windows to make sure the software paths are correct.

The command also generates some SLURM files:

``` sh 
sub.sh # (1)
array_job.sh # (2)
runner.sh # (3)
equil_1.txt # (4)
```

1. Blank SLURM submission file that is used as a template by runner.sh. You can also manually use this to sumbit your equil scripts.
2. File used by SLURM to run the array file. You shouldn't need to edit this.
3. File used to link and submit all simulations to SLURM with dependencies... This may or may not work so use with caution. Also it will likely max out your maximum queue allowance so will have to comment lines out and run multiple times. Good Luck! 
4. Text file with each line being a self-contained command to run each window in the umbrella sampling simulation. This is used by SLURM as the array file which runs each line as an independent simulation. There will be multiple of these if you have set __MaxStepsPerCalc__ > 0 and they __MUST__ be run in order. 

### Production

Running the production stage of Umbrella sampling is identical to the equilibration stage. 

It requires the same `Umbrella.inp` file and __should__ be performed on a HPC. 

``` sh title="Setting up production"
pyBrella -stg prod -dr True -i Umbrella.inp
```

Like the equilibration stage, this will generate the `prod_X.txt` files that contain commands for running each step in the production run. Submit them in order otherwise they __will__ break!

### Wham

You can perform a WHAM calculation at any stage after the windows are pulled. This is done through the [WHAM code](CodeReference.md/#wham-code). This is not the best way to run a WHAM calculation however, and [convergence](#convergence) is a more robust implementation of this method.

You can run a WHAM calculation using:

``` sh title="performing WHAM on the prod_1 stage of calculations"
pyBrella -stg wham -dr False -af Prod_1
```

!!! warning 

    Ensure wham.sh is in the $PATH otherwise this will not work! This can be found [here](http://membrane.urmc.rochester.edu/?page_id=126)

This will generate a sub-directory called WHAM. Within this, it generates:

- `AnalysisFile`metadata.dat : Contains the list of files containing colvar values
- `AnalysisFile`UmbrellaHist.dat : produces a histogram of the distribution of colvar values (can be plotted in LaTex using Tikz)
- PMF.dat : File containing PMF information
- PMF.eps : EPS plot of the PMF
- plot_free_energy.dat : data that can be plotted using xmgrace
- wham.sh : script file that runs the wham calculation.


### Convergence

This is the better way of running the WHAM calculation. It iterates through each of the sub-steps, in order for you to understand how well your calculation is converging.

To run a convergence calculation, run:

``` sh title="Running a convergence calculation"
pyBrella -stg convergence -dr False -af prod
```

The Convergence calculation outputs files to WHAM/Conv. For each of the sub-steps, it calculates the PMF using all previous sub-steps and outputs these to:

- `sub-step`.dat : Raw PMF data
- `sub-step`.eps : Plot of PMF data

???+ tips "Advanced error checking"

    If you have a complex potential energy surface, sometimes umbrella bins fail during the umbrella sampling simulation through moving to a different area along the potential energy surface. This can be problematic as it often causes a sudden shift in the relative energies of that bin. We have implemented an error checking tool that uses extra input files (See below) to identify know problematic structures. Use at your own risk! 

    ``` sh
    BondErrors.dat # (1)
    DihedralErrors.dat # (2)
    ```

    1. File containing known bonds that form when they shouldn't. It requires 4 colummns: Name, Atom1, Atom2, Value. The name is a pseudonym to identify the bond, the two atom columns are the atomic index' of the atoms in the bond, and the value is the standard bond length. If this new bond is formed during a bin, the bin shall be excluded from the WHAM calculation.
    2. File containing known problematic dihedral angles. tval1 is the acceptable dihedral angle, t2val is the problematic value. If the dihedral is closer to the problematic value then the bin will be ignored.

    BondErrors.dat

    | Name | Atom1 | Atom2 | Value |
    | :--- | :---- | :---- | :---- |
    | pseudonym | index of atom1 | index of atom2 | Value that classes the bond as formed |

    DihedralErrors.dat 

    | Name | Atom1 | Atom2 | Atom3 | Atom4 | t1name | t1val | t2name | t2val |
    | :--- | :---- | :---- | :---- | :---- | :----- | :---- | :----- | :---- | 
    | pseudonym | index of atom1 | index of atom2 | index of atom3 | index of atom4 | pseudonym of accepted dihedral | Accepted value | pseudonym of problematic dihedral | problematic value | 

### Analysis

Using further configuration files, you can also perform analysis and track other variables of the system throughout the simulation. This uses the [analysis code](CodeReference.md/#pyBrellaSampling.Tools.analysis.analysis) and reads in data files to configure states to track. (See below)


``` sh title="Running analysis"
pyBrella -stg analysis -dr False
```

Bonds.dat

| Name | Atom1 | Atom2 | Value |
| :--- | :---- | :---- | :---- |
| pseudonym | index of atom1 | index of atom2 | Standard bond length |

Dihedrals.dat 

| Name | Atom1 | Atom2 | Atom3 | Atom4 | t1name | t1val | t2name | t2val |
| :--- | :---- | :---- | :---- | :---- | :----- | :---- | :----- | :---- | 
| pseudonym | index of atom1 | index of atom2 | index of atom3 | index of atom4 | pseudonym of state 1 | central value | pseudonym of state 2 | central value |

This code goes through the simulation trajectories, and tracks the states of the bonds and dihedrals provided. It then plots them as histogram plots to show the overall distribution of the bonds and dihedrals. It also generates 2D histograms to plot the distribution of bonds and dihedrals against the target colvar.
