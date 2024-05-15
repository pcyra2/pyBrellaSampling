# Installation process

## Install

To install pyBrellaSampling, first clone the repo into a directory:

``` sh
gh repo clone pcyra2/pyBrellaSampling
```

Then move some script files out of the source directory:
``` sh
cp pyBrellaSampling/setup.py .
cp pyBrellaSampling/requirements.txt .

```

Next create the conda environment:

``` sh
conda create -n pyBrellaSampling --file requirements.txt -c conda-forge
conda activate pyBrellaSampling
```

Finally install the environment:

``` sh
pip install -e .
```

## Pre-requisits 

### ORCA

To use any QMMM functionality or benchmarking, ORCA must be installed. You can obtain this [here](https://orcaforum.kofo.mpg.de/app.php/portal)

### NAMD

For all dynamics simulations, we use NAMD which can be found [here](http://www.ks.uiuc.edu/Research/namd/). Two versions of NAMD are recommended:

- a CUDA-GPU accelorated version for running standard MD
- a CPU version for running QM/MM


### WHAM

You must also have the WHAM code from the [Grossfield Lab](http://membrane.urmc.rochester.edu/?page_id=126) in order to perform WHAM calculations for Umbrella Sampling.


## Configure

It is important now to configure the default variables! 
Go into the UserVars directory within the source.

``` sh
cd pyBrellaSampling/UserVars
```

Here you __MUST__ edit the HPC_Config.py and SoftwarePaths.py files.

The `HPC_Config` file needs configuring so that your HPC variables like QoS are configured correctly. This allows for SLURM scripts to be generated correctly.

``` py title="HPC_Config.py"
HPC_Config = {
    "None" : {
        "HostName" : "", # str: Hostname of the HPC
        "MaxWallTime" : 0, # int: Maximum calculation time in hours
        "Partition" : "None", # str: Partition to run calc on.
        "MaxCores" : 0, # Maximum number of cores per node 
        "QualityofService" : "None", # str: Account QoS or None if not needed
        "Account" : "None", # str: SLURM account or None if not needed
        "SoftwareLines" : ["",], # list: List of lines required to set up environment
    },
    "sulis" : {
        "HostName" : "login01.sulis.hpc", # str: Hostname of the HPC
        "MaxWallTime" : 48, # int: Maximum calculation time in hours
        "Partition" : "compute", # str: Partition to run calc on.
        "MaxCores" : 128, # Maximum number of cores per node 
        "QualityofService" : "None", # str: Account QoS or None if not needed
        "Account" : "su120", # str: SLURM account or None if not needed
        "SoftwareLines" : ["module load GCC/11.2.0", 
                           "module load OpenMPI/4.1.1", 
                           "module load ORCA/5.0.4"], # list: List of lines required to set up environment
    },
}
```

The `SoftwarePaths` file needs the absolute paths of the ORCA executable, and CPU and GPU version of NAMD in order to run. 

``` py title="SoftwarePaths.py"
# ORCA path
ORCA_PATH = "/home/pcyra2/Software/ORCA/5.0.4/orca" 

#NAMD paths
NAMD_CPU = "/home/pcyra2/Software/NAMD/NAMD_3.0b4_Linux-x86_64-multicore/namd3"

NAMD_GPU = "/home/pcyra2/Software/NAMD/NAMD_3.0b4_Linux-x86_64-multicore-CUDA/namd3"
```

??? tip "Other configs."

    Here you can also edit the other files in the directory. This allows you to:

        - Tune your MD parameters (MM_Variables.py)  
        - List available QM methods to ORCA for benchmarking (QM_Methods.py)
        - Change the Default inputs for each calculation type (if you always use the same QM method for example.) (Defaultinputs.py)

!!! note "Final stage (optional)"

    Within the source directory, run:
    ``` sh
    mkdocs build
    ```

    This will compile the documentation webpage and update the current configuration below.

## Current Configuration:

``` py title="Defaultinputs.py"
--8<-- "./UserVars/Defaultinputs.py"
```

``` py title="MM_Variables.py"
--8<-- "./UserVars/MM_Variables.py"
```

``` py title="QM_Methods.py"
--8<-- "./UserVars/QM_Methods.py"
```

``` py title="HPC_Config.py"
--8<-- "./UserVars/HPC_Config.py"
```

``` py title="SoftwarePaths.py"
--8<-- "./UserVars/SoftwarePaths.py"
```