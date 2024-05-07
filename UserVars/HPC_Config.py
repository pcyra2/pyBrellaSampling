HPC_Config = {
    "None" : {
        "HostName" : "", # str: Hostname of the HPC
        "MaxWallTime" : 0, # int: Maximum calculation time in hours
        "Partition" : "None", # str: Partition to run calc on.
        "MaxCores" : 0, # Maximum number of cores per node 
        "QualityofService" : None, # str: Account QoS or None if not needed
        "Account" : None, # str: SLURM account or None if not needed
        "SoftwareLines" : ["",], # list: List of lines required to set up environment
    },
    "sulis" : {
        "HostName" : "login01.sulis.hpc", # str: Hostname of the HPC
        "MaxWallTime" : 48, # int: Maximum calculation time in hours
        "Partition" : "compute", # str: Partition to run calc on.
        "MaxCores" : 128, # Maximum number of cores per node 
        "QualityofService" : None, # str: Account QoS or None if not needed
        "Account" : "su120", # str: SLURM account or None if not needed
        "SoftwareLines" : ["module load GCC/11.2.0", 
                           "module load OpenMPI/4.1.1", 
                           "module load ORCA/5.0.4"], # list: List of lines required to set up environment
    },
}