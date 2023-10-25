# pyBrellaSampling
This is the readme for the pyBrellaSampling code. 

The pyBrellaSampling code is primarily designed to assist in the setup and
running of QM/MM Umbrella sampling using NAMD and ORCA. It is designed to work
with initial AMBER coordinates and parameter files, however can be modified to
use other files that are supported by NAMD. It is also designed to work both
locally and using the SLURM work scheduler on an external HPC such as ARCHER2.

To work correctly, the user must have orca, namd and vmd available and in their
path for the script to work. 

## Install pyBrellaSampling: 
1. Firstly clone this repo into a directory

   	- gh repo clone pcyra2/pyBrellaSampling

2. Copy files for pip outside of the Source directory

	- cp pyBrellaSampling/setup.py .

	- cp pyBrellaSampling/requirements.txt .

3. Create an environment

	- conda create -n pyBrellaSampling --file requirements.txt -c conda-forge

	- conda activate pyBrellaSampling

4. install the environment.

	- pip install -e . 

## Usage:

### Umbrella Sampling:

To generate input files for umbrella sampling, run "pyBrella -jt inpfile", then edit all ".conf" files in order to set up your calculation, importantly:
- "JobType" needs to be set to "umbrella"
- "Stage" is the calculation stage, either "init", "min", "heat", "pull", "equil" or "prod"
- "QmPath" is the full path location to your orca executable

### Standalone Calculations:

Although pyBrellaSampling was designed to run QM/MM Umbrella sampling
simulations, it has the functionality to a variety of standalone calculations.

- Molecular Dynamics
	- minimisation
	- heating
	- NVT
	- NPT
- QM/MM
- Steered/constrained molecular dynamics

It should be noted however that the initial code was **NOT** designed to do this
and so may not work perfectly. It is recommended to thoroughly test the code in
order to ensure that it works how you intend. If there are any issues, or you
have any requests, please use the GitHub repo or email me on
ross.amory98@gmail.com. 

To generate input files for a standalone calculation, run "Standalone -jt
inpfile", then edit all ".conf" files to suit your calculation. 

