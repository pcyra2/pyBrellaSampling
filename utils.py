import numpy
import logging as log
import argparse as ap
import json
import subprocess
import numpy as np


def file_read(path):
    with open(path, 'r') as file:
        data = file.readlines()
    return data

def file_write(path, lines):
    with open(path, 'w') as f:
        for i in lines:
            print(i,file=f)

def file_2dwrite(path, x, y):
    assert len(x) == len(y), f"Column 1 length is {len(x)} but column 2 length is{len(y)}. They need to be equal to work"
    with open(path, 'w') as f:
        for i in range(len(x)):
            print(f"{x[i]}\t{y[i]}", file=f)

def data_2d(path):
    data0 = file_read(path)
    data = []
    for i in data0:
        if "#" not in i:
           data.append(i)
        if " 4000 " in i: ### For cutting down to 2 ps of data
            break
    col1 = numpy.zeros(len(data))
    col2 = numpy.zeros(len(data))
    for i in range(len(data)):
        words = data[i].split()
        col1[i] = words[0]
        col2[i] = words[1]
    return col1, col2

def dict_read(path):
    with open(path, "r") as f:
        data = json.load(f)
    return data

def dict_write(path, dict):
    with open(path,"w") as f:
        json.dump(dict, f)

def QM_Gen(qmzone,wd):
    tcl = f"""
mol new complex.parm7
mol addfile start.rst7

set qmPDB "syst-qm.pdb"
set qmPSF "syst-qm.QMonly.parm7"
set idDictFileName "syst-qm.idDict.txt"

set sel [atomselect 0 all]
$sel set beta 0
$sel set occupancy 0

set seltext "({qmzone})"
[atomselect 0 "$seltext"] set beta 1

puts "Initializing chain and QM region loops."
set systemSegs [lsort -unique [[atomselect 0 "protein or nucleic"] get chain]]
set systemQMregs [lsort -unique [[atomselect 0 "(protein or nucleic) and beta > 0"] get beta]]

puts "Chains: $systemSegs"
puts "QM Regions: $systemQMregs"

foreach seg $systemSegs {"{"}
    foreach qmReg $systemQMregs {"{"}
        puts "\nChecking QM region $qmReg in chain $seg"
        set qmmmm [atomselect 0 "(protein and name CA) and beta == $qmReg and chain $seg"]
        set cter [lindex [lsort -unique -integer [[atomselect 0 "chain $seg"] get resid]] end]
        set listqmmm [$qmmmm get resid]
        puts "Protein residues marked for QM this region in this chain: $listqmmm"
        list QM1bond
        list QM2bond
        puts "Checking N-Terminal-direction QM-MM bonds..."
        foreach resTest $listqmmm {"{"}
            if {"{"} [ lsearch $listqmmm [ expr $resTest -1 ] ] < 0 {"}"} {"{"}
                lappend QM1bond [ expr $resTest -1 ]
            {"}"}
        {"}"}
        puts "Checking C-terminal-direction QM-MM bonds..."
        foreach resTest $listqmmm {"{"}
             if {"{"} $resTest == $cter{"}"} {"{"}
                continue
            {"}"}
            if {"{"} [ lsearch $listqmmm [ expr $resTest +1 ] ] < 0 {"}"} {"{"}
                lappend QM2bond $resTest
            {"}"}
        {"}"}
        puts "Making changes..."
        if {"{"}[info exists QM2bond]{"}"} {"{"}
            [atomselect 0 "name CA C and (resid $QM2bond and chain $seg)"] set occupancy 1
            [atomselect 0 "name C O and (resid $QM2bond and chain $seg)"] set beta 0
            unset QM2bond
        {"}"}
        if {"{"}[info exists QM1bond]{"}"} {"{"}
            [atomselect 0 "name CA C and (resid $QM1bond and chain $seg)"] set occupancy 1
            [atomselect 0 "name C O and (resid $QM1bond and chain $seg)"] set beta $qmReg
            unset QM1bond
        {"}"}
       set qmmmm [atomselect 0 "(nucleic and name P) and beta == $qmReg and chain $seg"]
        set fiveTer [lindex [lsort -unique -integer [[atomselect 0 "chain $seg"] get resid]] 0]
        set listqmmm [$qmmmm get resid]
        puts "Nucleic residues marked for QM this region in this chain: $listqmmm"

        list QM1bond
        list QM2bond
        puts "Checking 3'-Terminal-direction QM-MM bonds..."
        foreach resTest $listqmmm {"{"}
            if {"{"}[ lsearch $listqmmm [ expr $resTest +1 ] ] < 0 {"}"} {"{"}
                lappend QM1bond [ expr $resTest +1 ]
            {"}"}
        {"}"}

        puts "Checking 5'-terminal-direction QM-MM bonds..."
        foreach resTest $listqmmm {"{"}"{"}"}
            if {"{"} $resTest == $fiveTer{"}"} {"{"}
                continue
            {"}"}
            if {"{"} [ lsearch $listqmmm [ expr $resTest -1 ] ] < 0 {"}"} {"{"}
                lappend QM2bond $resTest
            {"}"}
        {"}"}

       puts "Making changes..."
        if {"{"}[info exists QM2bond]{"}"} {"{"}
            [atomselect 0 "name C4' C5' and (resid $QM2bond and chain $seg)"] set occupancy 1
            [atomselect 0 "name P O1P O2P O5' C5' H5' H5'' and (resid $QM2bond and chain $seg)"] set beta 0
            unset QM2bond
        {"}"}
        if {"{"}[info exists QM1bond]{"}"} {"{"}
            [atomselect 0 "name C4' C5' and (resid $QM1bond and chain $seg)"] set occupancy 1
            [atomselect 0 "name P O1P O2P O5' C5' H5' H5'' and (resid $QM1bond and chain $seg)"] set beta $qmReg
            unset QM1bond
        {"}"}
    {"}"}
{"}"}

puts "Setting atom elements"

package require topotools

topo guessatom element mass

puts "Elements guessed!"

foreach qmReg $systemQMregs {"{"}
    set qmnum [[atomselect 0 "beta == $qmReg"] num]
    set dummy [ [atomselect 0 "beta == $qmReg and occupancy > 0"] num ]
    puts "QM Region $qmReg contains $qmnum QM atoms and $dummy dummy atoms"
{"}"}

$sel writepdb $qmPDB

[ atomselect 0 "beta > 0" ] writepsf $qmPSF

set qmsel [ atomselect 0 "beta > 0" ]

set indxs [ $qmsel get index ]

set fileId [open $idDictFileName "w"]

for {"{"}set i 0{"}"} {"{"} $i < [$qmsel num] {"}"} {"{"}incr i{"}"} {"{"}

    set ID [lindex $indxs $i]

    set data "$i $ID"

    puts $fileId $data
{"}"}

close $fileId
  
quit
"""
    with open(f"{wd}qm_prep.tcl", "w") as f:
        print(tcl, file=f)

def ColVarPDB_Gen(Umbrella, Job):
    if Umbrella.atom3 == 0:
        tcl = f"""mol new complex.parm7
mol addfile start.rst7

set colPDB "syst-col.pdb"

set sel [atomselect 0 all]
$sel set beta 0
$sel set occupancy 0

set seltext "(index {Umbrella.atom1})"
[atomselect 0 "$seltext"] set beta 1

set seltext "(index {Umbrella.atom2})"
[atomselect 0 "$seltext"] set beta 2

$sel writepdb $colPDB

quit
"""
    elif Umbrella.atom4 == 0:
        tcl = f"""mol new complex.parm7
mol addfile start.rst7

set colPDB "syst-col.pdb"

set sel [atomselect 0 all]
$sel set beta 0
$sel set occupancy 0

set seltext "(index {Umbrella.atom1})"
[atomselect 0 "$seltext"] set beta 1

set seltext "(index {Umbrella.atom2})"
[atomselect 0 "$seltext"] set beta 2

set seltext "(index {Umbrella.atom3})"
[atomselect 0 "$seltext"] set beta 3



$sel writepdb $colPDB

quit
"""
    else:
        tcl = f"""mol new complex.parm7
mol addfile start.rst7

set colPDB "syst-col.pdb"

set sel [atomselect 0 all]
$sel set beta 0
$sel set occupancy 0

set seltext "(index {Umbrella.atom1})"
[atomselect 0 "$seltext"] set beta 1

set seltext "(index {Umbrella.atom2})"
[atomselect 0 "$seltext"] set beta 2

set seltext "(index {Umbrella.atom3})"
[atomselect 0 "$seltext"] set beta 3

set seltext "(index {Umbrella.atom4})"
[atomselect 0 "$seltext"] set beta 4

$sel writepdb $colPDB

quit
"""
    file_write(f"{Job.WorkDir}Colvar_prep.tcl", [tcl])
def colvar_gen(Umbrella, i, type, force):
    if type == "pull":
        freq = 1
        Stages = 100
        if i > Umbrella.StartBin:
            prevBin = i - 1
        elif i < Umbrella.StartBin:
            prevBin = i + 1
        else:  # Sets initial pull value to actual start value rather than bin target. Allows for smoother pulls.
            bins = Umbrella.BinVals
            # print(bins)
            # print(len(bins))
            bins = numpy.append(bins, Umbrella.Start)
            # print(bins)
            # print(len(bins))
            Umbrella.add_bins(bins)
            prevBin = len(Umbrella.BinVals) -1
            # print(prevBin)
    else:
        freq = 1
        Stages = 0
        prevBin = i
    if Umbrella.atom3 == 0: ### Bond distance
        print("Bond collective variable")
        # group1 {"{"} atomNumbers {Umbrella.atom1} {"}"}
        # group2 {"{"} atomNumbers {Umbrella.atom2} {"}"}
        file = f"""colvarsTrajFrequency     {freq}

colvar {"{"}
    name length
    distance {"{"}

        group1 {"{"}  atomsFile ../syst-col.pdb
                      atomsCol B
                      atomsColValue 1.00
               {"}"}
        group2 {"{"}  atomsFile ../syst-col.pdb
                      atomsCol B
                      atomsColValue 2.00
               {"}"}
    {"}"}
{"}"}

harmonic {"{"}
    name lenpot
    colvars length
    centers {Umbrella.BinVals[i]}
    forceConstant {force}
{"}"}
"""
    elif Umbrella.atom4 == 0: ### 3 Atom angle
        print("Angle collective variable")
        file = f"""colvarsTrajFrequency     {freq}

colvar {"{"}
    name angle
    angle {"{"}
        group1 {"{"} atomNumbers {Umbrella.atom1} {"}"}
        group2 {"{"} atomNumbers {Umbrella.atom2} {"}"}
        group3 {"{"} atomNumbers {Umbrella.atom3} {"}"}
    {"}"}
{"}"}

harmonic {"{"}
    name angpot
    colvars angle
    centers {Umbrella.BinVals[i]}
    forceConstant {force}
{"}"}
"""
    else:   ### Dihedral angle
        print("Dihedral collective variable")
        file = f"""colvarsTrajFrequency     {freq}

        colvar {"{"}
            name dihedral
            dihedral {"{"}
                group1 {"{"}  atomsFile ../syst-col.pdb
                              atomsCol B
                              atomsColValue 1.00
                       {"}"}
                group2 {"{"}  atomsFile ../syst-col.pdb
                              atomsCol B
                              atomsColValue 2.00
                       {"}"}
                group3 {"{"}  atomsFile ../syst-col.pdb
                              atomsCol B
                              atomsColValue 3.00
                       {"}"}
                group4 {"{"}  atomsFile ../syst-col.pdb
                              atomsCol B
                              atomsColValue 4.00
                       {"}"}
            {"}"}
        {"}"}

        harmonic {"{"}
            name dihedpot
            colvars dihedral
            centers {Umbrella.BinVals[prevBin]}
            forceConstant {force}
            targetCenters {Umbrella.BinVals[i]}
            outputCenters on
            targetNumSteps 50
            outputAccumulatedWork on
        {"}"}
        """
    return file

def init_bins(num, step, start):
    bins = np.zeros(num)
    bins_min = np.zeros(num)
    bins_max = np.zeros(num)
    for i in range(num):
        bins[i] = round(start + step * i,2)
    return bins

def slurm_gen(JobName, SLURM, Job, path):
    SoftwareLines = ""
    for i in SLURM.Software:
        SoftwareLines=SoftwareLines+"\n"+i
    if SLURM.ArrayJob != None:
        slurmScript = f"""#!/bin/bash
#SBATCH --job-name={JobName}
#SBATCH --time={SLURM.WallTime}:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={SLURM.Cores}
#SBATCH --cpus-per-task=1
#SBATCH --mincpus={SLURM.Cores}
#SBATCH --array=1-{SLURM.ArrayLength}

#SBATCH --partition={SLURM.Partition}
{SLURM.account}
{SLURM.qos}

#SBATCH --get-user-env
#SBATCH --parsable

#dep

{SLURM.dependency}

{SoftwareLines}

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=1
export OMP_PLACES=cores

export ARRAY_JOBFILE=array_job.sh
export ARRAY_TASKFILE={JobName}.txt
export ARRAY_NTASKS=$(cat $ARRAY_TASKFILE | wc -l)


sh $ARRAY_JOBFILE
"""
        file_write(f"{path}sub.sh", [slurmScript])
        arrayfile = f"""#!/bin/bash
RUNLINE=$(cat $ARRAY_TASKFILE | head -n $(($SLURM_ARRAY_TASK_ID*{SLURM.JobsPerNode})) | tail -n {SLURM.JobsPerNode})
eval "$RUNLINE wait"

"""
        file_write(f"{path}array_job.sh",[arrayfile])
    else:
        slurmScript = f"""#!/bin/bash
#SBATCH --job-name={JobName}
#SBATCH --time={SLURM.WallTime}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={SLURM.Cores}
#SBATCH --cpus-per-task=1
#SBATCH --mincpus={SLURM.Cores}

#SBATCH --partition={SLURM.Partition}
{SLURM.account}
{SLURM.qos}

#SBATCH --get-user-env
#SBATCH --parsable
 
{SLURM.dependency}

{SoftwareLines}

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=1
export OMP_PLACES=cores

{Job}
"""
        file_write(path, [slurmScript])

def mpi_gen(SLURM, Umbrella):

    JobNumber = 0
    for i in range(Umbrella.Bins):
        subprocess.run([f"cp runORCA.py ./{i}"], shell=True,)
        Cores = ""
        if JobNumber < SLURM.JobsPerNode:
            for j in range(SLURM.Cores):
                Cores += f"{JobNumber*SLURM.Cores+j},"
                print(Cores)
            subprocess.run([f"""sed -i "s/CPUS/{Cores}/g"  ./{i}/runORCA.py"""], shell=True)
            JobNumber = JobNumber + 1
        else:
            JobNumber = 0
            for j in range(SLURM.Cores):
                Cores += f"{JobNumber*SLURM.Cores+j},"

            subprocess.run([f"""sed -i "s/CPUS/{Cores}/g"  ./{i}/runORCA.py"""], shell=True)
            JobNumber = JobNumber + 1

def get_cellVec(Job, MM):
    data = file_read(f"{Job.WorkDir}{MM.ambercoor}")
    length = len(data) #Get the last line
    words = data[length-1].split() # Split the final line.
    # print(words)
    return float(words[0]) # return the first number

MM_DefaultVars = {  "parmfile" : "complex.parm7",
                    "ambercoor": "start.rst7",
                    "bincoordinates" : None,
                    "DCDfile" : "output",
                    "DCDfreq" : 1,
                    "restartname" : "output",
                    "restartfreq" : 1,
                    "outputname" : "output",
                    "outputTiming" : "100",
                    "amber" : "on",
                    "switching" : "off",
                    "exclude" : "scaled1-4",
                    "1-4scaling" : 0.833333333,
                    "scnb" : 2.0,
                    "readexclusions" : "yes",
                    "cutoff" : 8.0,
                    "watermodel" : "tip3",
                    "pairListdist" : 11,
                    "LJcorrection" : "on",
                    "ZeroMomentum" : "off",
                    "rigidBonds" : "all",
                    "rigidTolerance" : "1.0e-8",
                    "rigidIterations" : 100,
                    "timeStep" : 2,
                    "fullElectFrequency" : 1,
                    "nonBondedFreq" : 1,
                    "stepspercycle" : 10,
                    "PME" : "off",
                    "PMEGridSizeX" : 300,
                    "PMEGridSizeY" : 300,
                    "PMEGridSizeZ" : 300,
                    "PMETolerance" : "1.0e-6",
                    "PMEInterpOrder" : 4,
                    "cellBasisVector1" : "135.913174 0.0 0.0",
                    "cellBasisVector2" : "-45.30439133333333 128.1401693173162 0.0",
                    "cellBasisVector3" : "-45.30439133333333 -64.0700846586581 -110.97264187403509",
                    "cellOrigin" : "0 0 0",
                    "cellBasisVector" : 135.913174, #### Base if all 3 are calculated...
                    "langevin" : "on",
                    "langevinDamping" : 5,
                    "langevinTemp" : 300,
                    "langevinHydrogen" : "off",
                    "temperature" : 300,
                    "BerendsenPressure" : "off",
                    "qmForces" : "off",
                    "CUDASOAintegrate" : "off",
                    "run" : 1000,
                    "minimize" : 0,
                    "qmParamPDB" : "syst-qm.pdb",
                    "qmColumn" : "beta",
                    "qmBondColumn" : "occ",
                    "QMsimsPerNode" : 1,
                    "QMElecEmbed" : "on",
                    "QMSwitching" : "on",
                    "QMSwitchingType" : "shift",
                    "QMPointChargeScheme" : "round",
                    "QMBondScheme" : "cs",
                    "qmBaseDir" : "/dev/shm/NAMD",
                    "qmConfigLine" : """! PBE 6-31G* EnGrad D3BJ TightSCF 
%%output PrintLevel Mini Print\[ P_Mulliken \] 1 Print\[P_AtCharges_M\] 1 end""",
                    "qmMult" : "1 1",
                    "qmCharge" : "1 0",
                    "qmSoftware" : "orca",
                    "qmExecPath" : "~/Software/ORCA/orca",
                    "QMOutStride" : 1,
                    "qmEnergyStride" : 1,
                    "QMPositionOutStride" : 1,
                    "colvars" : "off",
                    "colvarsConfig" : "colvars.conf",
                    }

def array_script(ntasks):
    Script=f"""#!/bin/bash
RUNLINE=$(cat $ARRAY_TASKFILE | head -n $(($SLURM_ARRAY_TASK_ID*{ntasks})) | tail -n {ntasks}))
eval \"$RUNLINE wait\""""
    file_write("./array_job.sh", [Script])

def batch_sub(nequil=4, nprod=16):
    script = f"""#!/bin/bash
"""
    script += f"""echo \"Submitting equil_1\"
cp sub.sh run.sh
sed -i \"s/NAME/equil_1/g\" run.sh
ID=$(sbatch run.sh)
echo \"equil_1 ID is $ID\"
echo \"equil_1 ID is $ID\" >> SLURMID.dat

    """
    for i in range(1, nequil):
        script += f"""echo  \"Submitting equil_{i+1}\"
cp sub.sh run.sh
sed -i \"s/NAME/equil_{i+1}/g\" run.sh
sed -i \"s/#dep/#SBATCH --dependency=afterok:$ID/g\" run.sh
ID=$(sbatch run.sh)
echo \"equil_{i+1} ID is $ID\"
echo \"equil_{i+1} ID is $ID\" >> SLURMID.dat

"""
    for i in range(nprod):
        script += f"""echo  \"Submitting prod_{i+1}\"
cp sub.sh run.sh
sed -i \"s/NAME/prod_{i+1}/g\" run.sh
sed -i \"s/#dep/#SBATCH --dependency=afterok:$ID/g\" run.sh
ID=$(sbatch run.sh)
echo \"prod_{i+1} ID is $ID\"
echo \"prod_{i+1} ID is $ID\" >> SLURMID.dat

"""
    file_write("./runner.sh", [script])


Functionals_list = ["PBE", "PBE0", "B3LYP", "HF", "M062X", "wB97M-V", "wB97X-D3", "MP2"]
Basis_list = ["6-31G", "6-31G*", "6-31+G*", "cc-pVDZ", "cc-pVTZ",  "def2-SVP"]
Dispersion_list = ["D3BJ", "D3ZERO", ""] ## "" blank is for no dispersion to be benchmarked
DispersionCorrFunc = ["wB97M-V", "wB97X-D3"]
kcal = 627.51 ### a.u. to kcal/mol conversion 
