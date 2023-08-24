import numpy
import logging as log
import argparse as ap
import json

import numpy as np


def file_read(path):
    with open(path, 'r') as file:
        data = file.readlines()
    return data

def file_write(path, lines):
    with open(path, 'w') as f:
        for i in lines:
            print(i,file=f)
            
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
        
# def input_parser(args):
#     if args.Verbosity == 0:
#         log.basicConfig(format="%(levelname)s: %(message)s", level=log.CRITICAL)
#     elif args.Verbosity == 1:
#         log.basicConfig(format="%(levelname)s: %(message)s", level=log.DEBUG)
#         log.info("Verbose output.")
#     try:
#         data = dict_read(f"{args.WorkDir}{args.input}")
#         lst = vars(args)
#         newlst = lst
#         for i in data:
#             for j in lst:
#                 if j == i:
#                     newlst[j] = data[i]
#         args = ap.Namespace(**newlst)
#         dict_write(f"{args.WorkDir}{args.input}", newlst)
#     except FileNotFoundError:
#         lst = vars(args)
#         dict_write(f"{args.WorkDir}{args.input}", lst)
#     return args

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

{SLURM.dependency}

{SoftwareLines}

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=1
export OMP_PLACES=cores

export ARRAY_JOBFILE=array_job.sh
export ARRAY_TASKFILE={SLURM.ArrayJob}
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