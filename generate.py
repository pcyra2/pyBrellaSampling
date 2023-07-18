import logging as log
import numpy as np
import os
import os.path as path


def Namd_File(Calc, parm, coord1, coord2, qmvars, i, heating="off"):
    if Calc.Ensemble == "NVP":
        pressure = "on"
    else:
        pressure = "off"
    if heating == "on":
        heat_lines = f"""# Heating Control 
temperature     0
reassignFreq    50
reassignIncr    0.2
reassignHold    300

"""
    else:
        heat_lines = f"""# Temp Control
langevin            on
langevinDamping     5
langevinTemp        {Calc.Temp}
langevinHydrogen    off
temperature         {Calc.Temp}
"""
    try:
        if Calc.Force > 1000:
            Rest_Lines = f"""colvars on 
colvarsConfig colvars.pull.conf"""
        elif Calc.Force < 1000:
            Rest_Lines = \
f"""colvars on 
colvarsConfig colvars.const.conf"""
    except AttributeError:
        Rest_Lines = "# No Restraints applied"
    file = f"""## Umbrella {Calc.OutFile} input file 
# File Options
parmfile        {parm}
ambercoor       {coord1}
bincoordinates  {coord2}
DCDfile         {Calc.OutFile}.dcd
DCDfreq         {Calc.TrajOut}
restartname     {Calc.OutFile}.restart
restartfreq     {Calc.RestOut}
outputname      {Calc.OutFile}
outputTiming    {Calc.TimeOut}


# Calculation Options
amber           on
switching       off
exclude         scaled1-4
1-4scaling      0.833333333
scnb            2.0
readexclusions  yes
cutoff          {Calc.CutOff}
watermodel      tip3
pairListdist        11
LJcorrection    on
ZeroMomentum    on
rigidBonds      {Calc.Shake}
rigidTolerance  1.0e-8
rigidIterations 100
timeStep        {Calc.TimeStep}

fullElectFrequency  1
nonBondedFreq       1
stepspercycle       10

# PME Vars
PME             off
PMEGridSizeX    300
PMEGridSizeY    300 
PMEGridSizeZ    300
PMETolerance    1.0e-6
PMEInterpOrder  4 

cellBasisVector1    {Calc.CellVec} 0.0 0.0
cellBasisVector2    {(-1/3)*Calc.CellVec} {(2/3)*np.sqrt(2)*Calc.CellVec} 0.0
cellBasisVector3    {(-1/3)*Calc.CellVec} {(-1/3)*np.sqrt(2)*Calc.CellVec} {(-1/3)*np.sqrt(6)*Calc.CellVec}
cellOrigin          0 0 0

{heat_lines}

# Pressure Control
BerendsenPressure   {pressure}

# QMMM settings
qmForces                {Calc.QM}
qmParamPDB              "../syst-qm.pdb"
qmColumn                "beta"
qmBondColumn            "occ"
QMSimsPerNode           1
QMElecEmbed             on
QMSwitching             on
QMSwitchingType         shift
QMPointChargeScheme     round
QMBondScheme            "cs"
qmBaseDir               "/dev/shm/NAMD_{i}"
qmConfigLine            "{qmvars}"
qmConfigLine            "%%output PrintLevel Mini Print\[ P_Mulliken \] 1 Print\[P_AtCharges_M\] 1 end"
qmConfigLine            "%PAL NPROCS {Calc.Threads} END"
qmMult                  "1 {Calc.Spin}"
qmCharge                "1 {Calc.Charge}"
qmSoftware              "orca"
qmExecPath              "{Calc.QMpath}"
QMOutStride             1
qmEnergyStride          1
QMPositionOutStride     1

{Rest_Lines}

run {Calc.Steps}
"""
    return file


def QM_Gen(qmzone,wd):
    tcl = f"""
mol new complex.parm7
mol addfile complex.pdb

set qmPDB "syst-qm.pdb"
set qmPSF "syst-qm.QMonly.parm7"
set idDictFileName "syst-qm.idDict.txt"

set sel [atomselect 0 all]
$sel set beta 0
$sel set occupancy 0

set seltext "({qmzone})"
[atomselect 0 "$seltext"] set beta 1

puts "Initializing segment and QM region loops."
set systemSegs [lsort -unique [[atomselect 0 "protein or nucleic"] get segname]]
set systemQMregs [lsort -unique [[atomselect 0 "(protein or nucleic) and beta > 0"] get beta]]

puts "Segments: $systemSegs"
puts "QM Regions: $systemQMregs"

foreach seg $systemSegs {"{"}
    foreach qmReg $systemQMregs {"{"}
        puts "\nChecking QM region $qmReg in segment $seg"
        set qmmmm [atomselect 0 "(protein and name CA) and beta == $qmReg and segname $seg"]
        set cter [lindex [lsort -unique -integer [[atomselect 0 "segname $seg"] get resid]] end]
        set listqmmm [$qmmmm get resid]
        puts "Protein residues marked for QM this region in this segment: $listqmmm"
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
            [atomselect 0 "name CA C and (resid $QM2bond and segname $seg)"] set occupancy 1
            [atomselect 0 "name C O and (resid $QM2bond and segname $seg)"] set beta 0
            unset QM2bond
        {"}"}
        if {"{"}[info exists QM1bond]{"}"} {"{"}
            [atomselect 0 "name CA C and (resid $QM1bond and segname $seg)"] set occupancy 1
            [atomselect 0 "name C O and (resid $QM1bond and segname $seg)"] set beta $qmReg
            unset QM1bond
        {"}"}
       set qmmmm [atomselect 0 "(nucleic and name P) and beta == $qmReg and segname $seg"]
        set fiveTer [lindex [lsort -unique -integer [[atomselect 0 "segname $seg"] get resid]] 0]
        set listqmmm [$qmmmm get resid]
        puts "Nucleic residues marked for QM this region in this segment: $listqmmm"

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
            [atomselect 0 "name C4' C5' and (resid $QM2bond and segname $seg)"] set occupancy 1
            [atomselect 0 "name P O1P O2P O5' C5' H5' H5'' and (resid $QM2bond and segname $seg)"] set beta 0
            unset QM2bond
        {"}"}
        if {"{"}[info exists QM1bond]{"}"} {"{"}
            [atomselect 0 "name C4' C5' and (resid $QM1bond and segname $seg)"] set occupancy 1
            [atomselect 0 "name P O1P O2P O5' C5' H5' H5'' and (resid $QM1bond and segname $seg)"] set beta $qmReg
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

def colvar_gen(Umbrella, i, type, force):
    # if i > Umbrella.StartBin:
    #     prevBin = i - 1
    # elif i < Umbrella.StartBin:
    #     prevBin = i + 1
    # else:
    #     prevBin = i
    if type == "pull":
        freq = 20
    else:
        freq = 1
    if Umbrella.atom3 == 0: ### Bond distance
        file = f"""colvarsTrajFrequency     {freq}

colvar {"{"}
    name length
    distance {"{"}
        group1 {"{"} atomNumbers {Umbrella.atom1} {"}"}
        group2 {"{"} atomNumbers {Umbrella.atom2} {"}"}
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
        file = f"""colvarsTrajFrequency     {freq}

        colvar {"{"}
            name dihedral
            dihedral {"{"}
                group1 {"{"} atomNumbers {Umbrella.atom1} {"}"}
                group2 {"{"} atomNumbers {Umbrella.atom2} {"}"}
                group3 {"{"} atomNumbers {Umbrella.atom3} {"}"}
                group4 {"{"} atomNumbers {Umbrella.atom4} {"}"}
            {"}"}
        {"}"}

        harmonic {"{"}
            name dihedpot
            colvars dihedral
            centers {Umbrella.BinVals[i]}
            forceConstant {force}
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

def make_dirs(num, wdir):
    for i in range(num):
        dir_path = str(wdir) + str(i)
        if path.exists(dir_path):
            log.info(str(i) + " exists. Deleting!")
            try:
                os.rmdir(dir_path)
            except OSError:
                log.warning(f"{i} directory not empty, deletion failed...")
                pass
        try:
            os.mkdir(dir_path)
        except FileExistsError:
            log.warning(f"{i} directory exists, Skipping making new directory")
            pass
        # pull_file = Namd_File("")

def Min_Setup(Calc, Job):
    qm_config = f"! {Calc.Method} {Calc.Basis} D3BJ EnGrad TightSCF CFLOAT EasyConv"
    log.info(f"QM config line is {qm_config}")
    file = f"""## Minimisation file:
# File Options
parmfile        complex.parm7
ambercoor       start.rst7
restartname     min.restart
restartfreq     100
outputname      min
outputTiming    1000

# Calculation Options
amber           on
switching       off
exclude         scaled1-4
1-4scaling      0.833333333
scnb            2.0
readexclusions  yes
cutoff          8.0
watermodel      tip3
pairListdist    11
LJcorrection    on
ZeroMomentum    on

fullElectFrequency  1
nonBondedFreq       1
stepspercycle       10

# PME Vars
PME             on
PMEGridSizeX    300
PMEGridSizeY    300 
PMEGridSizeZ    300
PMETolerance    1.0e-6
PMEInterpOrder  4 

cellBasisVector1    {Calc.CellVec} 0.0 0.0
cellBasisVector2    {(-1/3)*Calc.CellVec} {(2/3)*np.sqrt(2)*Calc.CellVec} 0.0
cellBasisVector3    {(-1/3)*Calc.CellVec} {(-1/3)*np.sqrt(2)*Calc.CellVec} {(-1/3)*np.sqrt(6)*Calc.CellVec}
cellOrigin          0 0 0

# Temp Control
langevin            off
temperature         0

# Pressure Control
BerendsenPressure   off


minimize                10000
"""
    with open(f"{Job.WorkDir}min.conf",'w') as f:
        print(file, file=f)

def Heat_Setup(Calc, Job):
    qm_config = f"! {Calc.Method} {Calc.Basis} D3BJ EnGrad TightSCF CFLOAT EasyConv"
    log.info(f"QM config line is {qm_config}")
    Calc.Set_Length(100000, 2) # 500 steps at 2 fs = 1 ps
    Calc.Set_Outputs(200,10,100) # Timings, Restart, Trajectory
    Calc.Set_Shake("all") # restrain all bonds involving Hydrogen
    Calc.Job_Name("heat")
    Calc.Set_QM("off")
    Calc.Set_OutFile(f"{Calc.Name}")
    Calc.Set_Ensemble("Heating")
    file = Namd_File(Calc, "complex.parm7", "start.rst7", "min.restart.coor",qm_config, 0, heating="on")
    with open(f"{Job.WorkDir}heat.conf", 'w') as f:
        print(file, file=f)

def Pull_Setup(Umbrella, Calc, Job):
    qm_config = f"! {Calc.Method} {Calc.Basis} D3BJ EnGrad TightSCF CFLOAT D3BJ NormalSCF"  # Removed EasyConv for convergence assistance
    log.info(f"QM config line is {qm_config}")
    Calc.Set_Length(100, 2) # 100 steps at 2 fs = 200 fs
    Calc.Set_Outputs(5,1,10) # Timings, Restart, Trajectory
    Calc.Set_Shake("all") # restrain all bonds involving Hydrogen
    Calc.Set_Force(5000)
    Calc.Set_QM("on")
    if Calc.Ensemble == "NVT":
        pressure = "off"
    Joblist = [None]*Umbrella.Bins
    for i in range(Umbrella.Bins):
        if Umbrella.Width > 0:
            if Umbrella.BinVals[i] > Umbrella.Start:
                prevPull = f"../{i-1}/pull.{i-1}.restart.coor"
            elif Umbrella.BinVals[i] < Umbrella.Start:
                prevPull = f"../{i + 1}/pull.{i + 1}.restart.coor"
            elif Umbrella.BinVals[i] == Umbrella.Start:
                Umbrella.add_start(i)
                prevPull = f"../heat.restart.coor"
        else:
            if Umbrella.BinVals[i] < Umbrella.Start:
                prevPull = f"../{i-1}/pull.{i-1}.restart.coor"
            elif Umbrella.BinVals[i] > Umbrella.Start:
                prevPull = f"../{i + 1}/pull.{i + 1}.restart.coor"
            elif Umbrella.BinVals[i] == Umbrella.Start:
                Umbrella.add_start(i)
                prevPull = f"../heat.restart.coor"            
        Calc.Job_Name("pull")
        Calc.Set_OutFile(f"{Calc.Name}.{i}")
        file = Namd_File(Calc, "../complex.parm7", "../start.rst7",prevPull, qm_config, i)
        with open(f"{Job.WorkDir}{i}/pull.conf", 'w') as f:
            print(file, file=f)
        file = colvar_gen(Umbrella, i, "pull", Calc.Force)
        with open(f"{Job.WorkDir}{i}/colvars.pull.conf", 'w') as f:
            print(file, file=f)
        Joblist[i] = f"mkdir /dev/shm/NAMD_{i} ; cd ./{i} ; {Calc.NamdPath} pull.conf > pull.{i}.out ; cd ../ ; rm -r /dev/shm/NAMD_{i}"
    return Umbrella, Joblist

def Prod_Setup(Umbrella, Calc, Job, startCoord):
    qm_config = f"! {Calc.Method} {Calc.Basis} D3BJ EnGrad TightSCF CFLOAT SlowConv"
    log.info(f"QM config line is {qm_config}")
    Calc.Set_Shake("none")  # restrain all bonds involving Hydrogen
    Calc.Set_Force(300)
    Calc.Set_QM("on")
    if Calc.Ensemble == "NVT":
        pressure = "off"
    Joblist = [None] * Umbrella.Bins
    for i in range(Umbrella.Bins):
        Calc.Set_OutFile(f"{Calc.Name}.{i}")
        file = Namd_File(Calc, "../complex.parm7", "../start.rst7", f"{startCoord}.{i}.restart.coor",
                         qm_config, i)
        with open(f"{Job.WorkDir}{i}/{Calc.Name}.conf", 'w') as f:
            print(file, file=f)
        file = colvar_gen(Umbrella, i, "const", Calc.Force)
        with open(f"{Job.WorkDir}{i}/colvars.const.conf", 'w') as f:
            print(file, file=f)
        Joblist[
            i] = f"mkdir /dev/shm/NAMD_{i} ; cd ./{i} ; {Calc.NamdPath} {Calc.Name}.conf > {Calc.Name}.{i}.out ; cd ../ ; rm -r /dev/shm/NAMD_{i}"
    return Joblist