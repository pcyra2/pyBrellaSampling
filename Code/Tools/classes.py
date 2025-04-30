import pyBrellaSampling.Code.Tools.io as io
import pyBrellaSampling.Code.Tools.QM.pyscf_tools as pyscf_tools
import os
import subprocess
import numpy
import socket
import pyscf


CurrentPath = os.path.dirname(os.path.abspath(__file__))


class QMClass:
    def __init__(self,Method="PBE", Basis="6-31G*", charge=0, spin=0, Software="pyscf"):
        if Software.casefold() == "orca":
            self.software = orca_class()
        elif Software.casefold() == "pyscf":
            spin -= 1 # pyscf = 2S not 2S+1
            self.software = pyscf_class()
            self.software.set_envVars(Method, Basis, charge, spin)
        else:
            raise NotImplementedError(f"Software {Software} not recognized. Known software is orca and pyscf.")
        self.charge = charge
        self.spin = spin
        self.Method = Method
        self.Basis = Basis
    def set_cores(self, cores:int):
        self.cores = cores
        if self.software.name == "pyscf":
            os.environ['OMP_NUM_THREADS'] = str(cores)

    def init_qmmm(self, qmzone:str):
        self.vmd_selection = qmzone
    
class MMClass:
    jobs = {}
    qm = None
    def __init__(self, Software, parameters, topology):
        if Software.casefold() == "namd":
            self.software = namd_class()
        else:
            raise NotImplementedError(f"Software {Software} not recognized. Known software is namd.")

        if ".parm7" in parameters:
            Amber = True
        self.software.set_global(Amber, parameters, topology)
    def define_qm(self, package):
        self.qm = package
    def minimize(self,
                infile:str,
                outfile:str,
                steps:float, 
                temperature=0, 
                ):
        """
        performs a minimization of the system

        Args:
            infile (str): input file
            outfile (str): output file
            steps (int): Number of minimization steps
            temperature (int, optional): Temperature to minimize at. Defaults to 0.
            qm (classes.QMClass, optional): QM Class for QM/MM. Defaults to None.
        """
        self.software.set_ensemble(ensemble="min", temperature=temperature)
        if steps > 100:
            outputs = steps // 100
        else:
            outputs = 1
        self.software.set_outputs(outputs, outputs, 1, outputs) # Print the energy at every step but everything else do 100 times.
        file, config = self.software.gen_input(infile, outfile, steps)
        self.jobs[outfile] = {"file":file,
                              "config":config}
        return file
    def heat(self,
            infile:str,
            outfile:str,
            steps:int,
            timestep:float,
            init_temp:float,
            final_temp:float,
            ):
        if steps > 100:
            outputs = steps // 100
        else:
            outputs = 1
        temp_range = final_temp - init_temp
        if temp_range > steps:
            raise ValueError(f"Temperature range {temp_range} is greater than the number of steps {steps}. This will cause numerical instabilities..")
        elif temp_range > steps*0.9:
            reassignFreq = ((steps) // temp_range)*10 # Dont allow for the temp to stabilse
            reassignTemp = 10
        else:
            reassignFreq = ((steps*0.9) // temp_range)*10 # Allow for the temp to stabilse
            reassignTemp = 10

        self.software.set_outputs(outputs, outputs, 1, outputs) # Print the energy at every step but everything else do 100 times.
        self.software.set_ensemble(ensemble="heat", init_temp=init_temp, final_temp=final_temp, reassignFreq=reassignFreq, reassignIncr=reassignTemp)
        self.software.change_timestep(timestep)
        file, config = self.software.gen_input(infile, outfile, steps)
        self.jobs[outfile] = {"file":file,
                              "config":config}
        return file
    def constant(self,
                    infile:str,
                    outfile:str,
                    steps:int,
                    timestep:float,
                    traj_steps:int,
                    temperature:float,
                    pressure = None,
                    ):
        assert steps > traj_steps, f"Number of steps ({steps}) is less than the number of steps in the trajectory output ({traj_steps})."
        if pressure != None:
            ensemble = "npt"
        else:
            ensemble = "nvt"
        
        self.software.set_ensemble(ensemble=ensemble, temperature=temperature, pressure=pressure)
        self.software.change_timestep(timestep)
        if steps < 1000:
            en_steps = 1
        else:
            en_steps = 100
        self.software.set_outputs(en_steps, traj_steps, en_steps, traj_steps)
        file, config = self.software.gen_input(infile, outfile, steps)
        self.jobs[outfile] = {"file":file,
                              "config":config}
        return file
    def SMD(self, infile:str,
                    outfile:str,
                    colvarfile:str,
                    steps:int,
                    timestep:float,
                    traj_steps:int,
                    temperature:float,
                    pressure = None,
                    seed=None
                    ):
        assert steps > traj_steps, f"Number of steps ({steps}) is less than the number of steps in the trajectory output ({traj_steps})."
        if pressure != None:
            ensemble = "npt"
        else:
            ensemble = "nvt"
        
        self.software.set_ensemble(ensemble=ensemble, temperature=temperature, pressure=pressure)
        self.software.change_timestep(timestep)
        if steps < 1000:
            en_steps = 1
        else:
            en_steps = 100
        self.software.set_outputs(en_steps, traj_steps, en_steps, traj_steps)
        self.software.set_colvar(colvarfile)
        file, config = self.software.gen_input(infile, outfile, steps, seed)
        self.jobs[outfile] = {"file":file,
                              "config":config}
        return file

class atom:
    def __init__(self, element, x, y, z):
        self.element = element
        self.x = x
        self.y = y
        self.z = z
    def echo(self):
        return f"{self.element} {self.x} {self.y} {self.z}"
    def translate_x(self, distance):
        self.x += distance
    def translate_y(self, distance):
        self.y += distance
    def translate_z(self, distance):
        self.z += distance

class molecule:
    bohr2ang = 0.529177
    atoms = []
    def __init__(self, ):
        pass
    def from_xyz(self, path, charge, spin):
        lines = io.textRead(path)
        self.nat = int(lines[0])
        atoms = [atom]*self.nat
        for i in range(self.nat):
            items = lines[i+2].split()
            atoms[i] = atom(items[0], items[1], items[2], items[3])
        self.atoms = atoms
        self.charge = charge
        self.spin = spin
    def from_atoms_list(self, atoms:atom, charge:int, spin:int):
        self.nat = len(atoms)
        self.atoms = atoms
        self.charge = charge
        self.spin = spin
    def from_gtoMole(self, mole:pyscf.gto.Mole):
        atoms = mole._atom
        self.nat = mole.natm
        self.atoms = [atom]*mole.natm
        for i, at in enumerate(atoms):
            self.atoms[i] = atom(at[0], round(at[1][0]*self.bohr2ang,6), round(at[1][1]*self.bohr2ang,6), round(at[1][2]*self.bohr2ang,6))
        self.charge = mole.charge
        self.spin = mole.spin
        self.basis = mole.basis
    def print_coords(self)->str:
        text = ""
        for at in self.atoms:
            text += at.echo()+"\n"
        return text
    def to_gtoMole(self, symmetry:bool):
        """Converts the molecule class into 

        Args:
            basis (str): basis set to describe the molecule
            symmetry (bool): Whether to use symmetry

        Returns:
            mol (pyscf.gto.M): pySCF initialised molecule
        """
        mol = pyscf_tools.genMol(self, self.charge, self.spin, self.basis, symmetry)
        return mol

class Colour:
    def __init__(self, red:int, green:int, blue:int):
        self.red = red
        self.green = green
        self.blue = blue
    
    def to_string(self):
        return("rgb("+str(self.red)+", "+str(self.green)+", "+str(self.blue)+")")

class orca_class:
    name="orca"
    extras = ""
    def __init__(self):
        paths = io.jsonRead(CurrentPath+"/UserConfig/SoftwareLocations.conf")
        self.path = paths["orca"]
        self.post_inputfile = "'--use-hwthread-cpus --bind-to hwthread'"
    def run(self, input, output, cores):
        outlines = subprocess.run([f"{self.path} {input} {self.post_inputfile} > {output}"], shell=True, capture_output=True)
        return outlines
    def add_extras(self, args:list):
        self.extras = ""
        if type(args) == str:
            self.extras = args
        else:
            for arg in args:
                self.extras += str(args)+" "
    def gen_input(self,method:str, basis:str ,jobtype:str , charge:int , spin:int , mol, cores:int , ram:int , *args)->str:
        line1 = f"{method} {basis} {jobtype} "
        for arg in args:
            line1 += f" {arg}"
        if type(mol) == str:
            mol_lines = f"* xyzfile  {charge} {spin} {mol}" 
        elif type(mol) == list:
            mol_lines = f"* xyz {charge} {spin}"
            for line in mol:
                mol_lines += line+"\n"
        elif type(mol) == molecule:
            mol_lines = f"* xyz {charge} {spin}"
            mol_lines += mol.print_coords()

        file = f"""{line1}
%PAL NPROCS {cores} END
%maxcore {ram}
{mol_lines}
"""
        return file

class pyscf_class:
    name="pyscf"
    extras = ""
    def __init__(self):
        self.path = ""
        self.post_inputfile = ""
    def add_extras(self, args:list):
        pass
    def run(self, input, output, cores):
        pass
    def gen_input(self, method:str, basis:str ,jobtype:str , charge:int , spin:int , mol, cores:int , ram:int , *args):
        pass
    def set_envVars(self, method, basis, charge, spin):
        os.environ["method"] = str(method)
        os.environ["basis"] = str(basis)
        os.environ["charge"] = str(charge)
        os.environ["spin"] = str(spin)
        
class namd_class:
    name = "namd"
    config = {"colvarlines":""}
    def __init__(self):
        self.defaults = io.jsonRead(CurrentPath+"/UserConfig/NAMDDefaults.conf")
        paths = io.jsonRead(CurrentPath+"/UserConfig/SoftwareLocations.conf")
        self.path_cpu = paths["namd"]["cpu"]
        self.path_gpu = paths["namd"]["gpu"]
    def exec(self, input:str, output:str, gpu:bool):
        if gpu == True:
            pre_inputfile = "+oneWthPerCore +setcpuaffinity +devices 0"
            path = self.path_gpu
        else:
            pre_inputfile = ""
            path = self.path_cpu
        command = f"{path} {pre_inputfile} {input} > {output}"
        print(f"Running command: {command}")
        outlines = subprocess.run([command],shell=True, capture_output=True)
        return outlines
    def set_global(self, amber:bool, param:str, amber_coor:str):
        if amber == True:
            self.config["amber"] = "yes"
            self.config["parmfile"] = param
            self.config["ambercoor"] = amber_coor
        else:
            raise NotImplementedError("This has not been implemented yet... Sorry Bro. I only work with amber files. ")
    def set_outputs(self, restart:int, timing:int, energy:int,trajectory:int):
        self.config["DCDfreq"] = trajectory
        self.config["restartfreq"] = restart
        self.config["outputTiming"] = timing
        self.config["outputEnergies"] = energy
    def set_ensemble(self, **kwargs):
        ensemble = kwargs["ensemble"]
        knownEnsembles = ["min", "heat", "nvt", "npt"]
        if ensemble.casefold() not in knownEnsembles:
            raise ValueError(f"Ensemble {ensemble} not recognized. Known ensembles are {knownEnsembles}")
        self.config["ensemble"] = ensemble
        if ensemble.casefold() == "min":
            if "temperature" not in kwargs:
                print("No temperature specified for minimization. Setting to 0 K.")
                kwargs["temperature"] = 0
            self.config["command"] = "minimize"
            self.config["templines"] = ""
            self.config["pressurelines"] = ""
            self.config["temperature"] = kwargs["temperature"]
            self.config["langevin"] = "off"
        elif ensemble.casefold() == "heat":
            self.config["command"] = "run"
            self.config["temperature"] = kwargs["init_temp"]
            self.config["templines"] = f"""
reassignFreq        {int(kwargs["reassignFreq"])}
reassignIncr        {kwargs["reassignIncr"]}
reassignHold        {kwargs["final_temp"]}
"""
            self.config["langevin"] = "off"
            self.config["pressurelines"] = ""
        elif ensemble.casefold() == "nvt":
            self.config["command"] = "run"
            self.config["temperature"] = kwargs["temperature"]
            self.config["langevin"] = "on"
            self.config["templines"] = f"""langevinDamping    5
langevinTemp        {kwargs["temperature"]}
langevinHydrogen    off
"""
            self.config["pressurelines"] = ""
        elif ensemble.casefold() == "npt":
            self.config["command"] = "run"
            self.config["langevin"] = "on"
            self.config["temperature"] = kwargs["temperature"]
            self.config["templines"] = f"""langevinDamping    5
langevinTemp        {kwargs["temperature"]}
langevinHydrogen    off
"""
            self.config["pressurelines"] = f"""
langevinPison           on
langevinPistonTarget    {kwargs["pressure"]}
langevinPistonPeriod    100
langevinPistonDecay     50
langevinPistonTemp      {kwargs["temperature"]}
"""
    def set_qm(self, QM:QMClass, qmFilePath:str):
        self.qmzone = QM.vmd_selection
        self.config["qmForces"] = "on"
        if QM.software.name == "orca":
            self.config["qmlines"] = f"""
qmParamPDB              "{qmFilePath}"
qmColumn                "beta"
qmBondColumn            "occ"
QMsimsPerNode           1
QMElecEmbed             on
QMSwitching             on
QMSwitchingType         shift
QMPointChargeScheme     round
QMBondScheme            "cs"
qmBaseDir               "/dev/shm/RUNDIR"
qmConfigLine            "! {QM.Method} {QM.Basis} EnGrad {QM.software.extras}"
qmConfigLine            "%%output PrintLevel Mini Print\[ P_Mulliken \] 1 Print\[P_AtCharges_M\] 1 end"
qmConfigLine            "%PAL NPROCS {QM.cores} END"
qmMult                  "1 {QM.spin}"
qmCharge                "1 {QM.charge}"
qmSoftware              "orca"
qmExecPath              "{QM.software.path}"
QMOutStride             1
qmEnergyStride          1
QMPositionOutStride     1
"""
        elif QM.software.name == "pyscf":
            self.config["qmlines"]=f"""
qmParamPDB              "{qmFilePath}"
qmColumn                "beta"
qmBondColumn            "occ"
QMsimsPerNode           1
QMElecEmbed             on
QMSwitching             on
QMSwitchingType         shift
QMPointChargeScheme     round
QMBondScheme            "cs"
qmBaseDir               "/dev/shm/RUNDIR"
qmConfigLine            "{QM.Method} {QM.Basis}"
qmMult                  "1 {QM.spin}"
qmCharge                "1 {QM.charge}"
qmSoftware              "custom"
qmExecPath              "pyscfQMMM"
QMOutStride             1
qmEnergyStride          1
QMPositionOutStride     1
"""
    def gen_input(self, infile:str, outfile:str, steps:int, seed=None):
        self.config["run"] = steps
        for key in self.defaults.keys():
            if key not in self.config.keys():
                self.config[key] = self.defaults[key]
        if self.config["GPUresident"] == "on":
            self.config["GPUresident"] = "GPUResident   on" #TODO: implement NAMD version compatibility
        else:
            self.config["GPUresident"] = ""
        if seed == None:
            self.config["seed"] = ""
        else:
            self.config["seed"] = f"seed    {seed}"
        if infile == self.config["ambercoor"]: # Starting from the initial amber file, not a previous trajectory.
            self.config["bincoordinates"] = ""
            self.config["extendedSystem"] = ""
            self.config["dcdUnitCell"] = "no"
            coor = io.textRead(self.config["ambercoor"])
            CellVec = float(coor[-1].split()[0])
            self.config["cell"] = f"""
cellBasisVector1    {CellVec} 0.0 0.0
cellBasisVector2    {(-1/3)*CellVec} {(2/3)*numpy.sqrt(2)*CellVec} 0.0
cellBasisVector3    {(-1/3)*CellVec} {(-1/3)*numpy.sqrt(2)*CellVec} {(-1/3)*numpy.sqrt(6)*CellVec}
cellOrigin          0 0 0
"""
        else:
            self.config["bincoordinates"] = "bincoordinates      ${input}.restart.coor" # Read in the coordinates
            self.config["extendedSystem"] = "extendedSystem      ${input}.xsc" # Read in the PBC info. 
            self.config["dcdUnitCell"] = "yes"
            self.config["cell"] = ""
        file = f"""### pyBrellaSampling generated input file
# Variables
set input           {infile}
set output          {outfile}
timeStep            {self.config["timeStep"]}


# File options
amber               {self.config["amber"]}
parmfile            {self.config["parmfile"]}
ambercoor           {self.config["ambercoor"]}
{self.config["bincoordinates"]}
{self.config["extendedSystem"]}
DCDfile             {"${output}.dcd"}
DCDfreq             {self.config["DCDfreq"]}
restartname         {"${output}.restart"}
restartfreq         {self.config["restartfreq"]}
outputname          {"${output}"}
outputTiming        {self.config["outputTiming"]}
outputEnergies      {self.config["outputEnergies"]}



# General calculation parameters
switching           {self.config["switching"]}
exclude             {self.config["exclude"]}
1-4scaling          {self.config["1-4scaling"]}
scnb                {self.config["scnb"]}
readexclusions      {self.config["readexclusions"]}
cutoff              {self.config["cutoff"]}
watermodel          {self.config["watermodel"]}
pairListdist        {self.config["pairListdist"]}
LJcorrection        {self.config["LJcorrection"]}
ZeroMomentum        {self.config["ZeroMomentum"]}
rigidBonds          {self.config["rigidBonds"]}
rigidTolerance      {self.config["rigidTolerance"]}
rigidIterations     {self.config["rigidIterations"]}
fullElectFrequency  {self.config["fullElectFrequency"]}
nonBondedFreq       {self.config["nonBondedFreq"]}
stepspercycle       {self.config["stepspercycle"]}

# PME options
PME                 {self.config["PME"]}
PMEGridSizeX        {self.config["PMEGridSizeX"]}
PMEGridSizeY        {self.config["PMEGridSizeY"]}
PMEGridSizeZ        {self.config["PMEGridSizeZ"]}
PMETolerance        {self.config["PMETolerance"]}
PMEInterpOrder      {self.config["PMEInterpOrder"]}

# Cell options
{self.config["wrap"]}
{self.config["cell"]}
dcdUnitCell         {self.config["dcdUnitCell"]}

# Temperature options
temperature         {self.config["temperature"]}
{self.config["templines"]}

# Pressure options
langevin            {self.config["langevin"]}
{self.config["pressurelines"]}

# QMMM options
qmForces            {self.config["qmForces"]}
{self.config["qmlines"]}


# Colvar options
{self.config["colvarlines"]}

# GPU speedups
{self.config["GPUresident"]}

{self.config["seed"]}
{self.config["command"]}      {self.config["run"]}
"""
        
        return file, self.config
    def change_timestep(self, timestep:float):
        self.config["timeStep"] = timestep
    def change_config(self,dictionary:dict):
        """
        Used to change any config variable

        Args:
            dictionary (dict): dictionary of variables to change with their values. 
        """
        for key in dictionary.keys():
            self.config[key] = dictionary[key]
    def set_colvar(self, file:str):
        self.config["colvarlines"] = f"""
colvars         on
colvarsConfig   {file}
"""

    def check_output(self,file:str):
        status = "Not started"
        data = []
        if os.path.isfile(file):
            status = "started"
            data = io.textRead(file)
            if "End of program" not in data[-1]:
                status = "running"
            if "ERROR" in data[-5:-2]:
                status = "error"
            if "WallClock" in data[-2]:
                status = "completed"
        return status, data[-5:-1]
            
class ColvarClass:
    def __init__(self,):
        self.bins={}
    def atomic_colvar(self, 
                      atoms:list, 
                      pull:float, 
                      hold:float, 
                      minimum:int, 
                      maximum:int, 
                      initial:float, 
                      step:float):
        self.atoms = atoms
        if len(atoms) == 2:
            self.VariableType = "length"
            self.VariableLines = """colvarsTrajFrequency    1
colvar {
    name length
    distance {
        group1 {
            atomsFile ../syst-col.pdb
            atomsCol B
            atomsColValue 1.00
        }
        group2 {
            atomsFile ../syst-col.pdb
            atomsCol B
            atomsColValue 2.00
        }
    }
}

"""
        elif len(atoms) == 3:
            self.VariableType = "angle"
            self.VariableLines = """colvarsTrajFrequency    1
colvar {
    name angle
    distance {
        group1 {
            atomsFile ../syst-col.pdb
            atomsCol B
            atomsColValue 1.00
        }
        group2 {
            atomsFile ../syst-col.pdb
            atomsCol B
            atomsColValue 2.00
        }
        group3 {
            atomsFile ../syst-col.pdb
            atomsCol B
            atomsColValue 3.00
        }
    }
}

"""
        elif len(atoms) == 4:
            self.VariableType = "dihedral"
            self.VariableLines = """colvarsTrajFrequency    1
colvar {
    name dihedral
    distance {
        group1 {
            atomsFile ../syst-col.pdb
            atomsCol B
            atomsColValue 1.00
        }
        group2 {
            atomsFile ../syst-col.pdb
            atomsCol B
            atomsColValue 2.00
        }
        group3 {
            atomsFile ../syst-col.pdb
            atomsCol B
            atomsColValue 3.00
        }
        group4 {
            atomsFile ../syst-col.pdb
            atomsCol B
            atomsColValue 4.00
        }
    }
}

"""
        else:
            raise NotImplementedError(f"Number of atoms {len(atoms)} not recognized. Only 2, 3 and 4 atoms are supported.")
        self.PullForce = pull
        self.HoldForce = hold
        self.Min = minimum
        self.Max = maximum
        self.stepsize = step
        self.initial = initial
        distance = maximum - minimum
        self.nsteps = int((distance // step) + 1)
        self.tracker = TrackerClass(atoms, "Colvar")
        self.bins = [minimum + i*step for i in range(self.nsteps)]
    def update_initial_point(self, point):
        self.initial = point
       
class TrackerClass:
    def __init__(self, atoms:list, Name:str, ):
        self.data = {}
        self.Name = Name
        self.atoms = atoms
        if len(atoms) == 2:
            self.type = "distance"
        elif len(atoms) == 3:
            self.type = "angle"
        elif len(atoms) == 4:
            self.type = "dihedral"
        else:
            raise NotImplementedError("I only know how to track upto 4 objects (distance, angle or dihedral)")
    def gen_vmdScript(self,):
        if self.type == "distance":
            lines = f"""label add Bonds 0/{self.atoms[0]} 0/{self.atoms[1]}
label graph Bonds 0 {self.Name}.dat
label delete Bonds 0"""
        elif self.type == "angle":
            lines = f"""label add Angles 0/{self.atoms[0]} 0/{self.atoms[1]} 0/{self.atoms[2]}
label graph Angles 0 {self.Name}.dat
label delete Angles 0"""
        elif self.type == "dihedral":
            lines = f"""label add Dihedrals 0/{self.atoms[0]} 0/{self.atoms[1]} 0/{self.atoms[2]} 0/{self.atoms[3]}
label graph Dihedrals 0 {self.Name}.dat
label delete Dihedrals 0"""
        self.vmd_lines = lines
        return lines
    def get_vmdData(self, Directory:str, Calculation:str, bin = None):
        if bin == None:
            data = io.textRead(os.path.join(Directory,f"{self.Name}.dat"))
            self.data[Calculation] = numpy.zeros_like(data).tolist()
            for i, dat in enumerate(data):
                words = dat.split()
                self.data[Calculation][i] = words[1]
        else:
            if Calculation not in self.data:
                self.data[Calculation] = {}
            data = io.textRead(os.path.join(Directory,f"{self.Name}.dat"))
            self.data[Calculation][bin] = numpy.zeros_like(data).tolist()
            for i, dat in enumerate(data):
                words = dat.split()
                self.data[Calculation][bin][i] = words[1]
    def dump(self, path:str):
        io.jsonDump(self.data, os.path.join(path, f"{self.Name}_Analysis.json"))
    def read_previous(self, WorkDir:str):
        if os.path.isfile(os.path.join(WorkDir, f"{self.Name}_Analysis.json")):
            self.data = io.jsonRead(os.path.join(WorkDir, f"{self.Name}_Analysis.json") )

class VMDClass:
    def __init__(self, ParmFile:str, AnalysisLines:str):
        self.colvarPDB_script = None
        self.qmPDB_script = None
        self.parmfile = ParmFile
        self.analysis_lines = AnalysisLines
    def qmPDB_gen(self, MM:MMClass):
        file = f"""
mol new {self.parmfile}
mol addfile {MM.software.config["ambercoor"]}

set qmPDB "syst-qm.pdb"
set qmPSF "syst-qm.QMonly.parm7"
set idDictFileName "syst-qm.idDict.txt"

set sel [atomselect 0 all]
$sel set beta 0
$sel set occupancy 0

set seltext "({MM.software.qmzone})"
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
        self.qmPDB_script = file
    def colvarPDB_gen(self, colvar: ColvarClass, MM: MMClass):
        """tcl script to obtain the syst-col.pdb file. It is a requirement for SMD and Umbrella sampling.
        
        Args:
            MM (MMClass): Class containing information about the MD system.
        """
        file = f"""mol new {self.parmfile}
mol addfile {MM.software.config["ambercoor"]}

set colPDB "syst-col.pdb"

set sel [atomselect 0 all]
$sel set beta 0
$sel set occupancy 0

"""
        atom_sel = ["length", "angle", "dihedral"]
        if colvar.VariableType in atom_sel:
            for index,  atom in enumerate(colvar.atoms):
                lines = f"""set seltext "(index {atom})"
[atomselect 0 "$seltext"] set beta {index+1}

"""
                file += lines
        else:
            raise NotImplementedError(f"Variable type {colvar.VariableType} not implemented.")
        
        file += """$sel writepdb $colPDB

quit
"""
        self.colvarPDB_script = file
    def GenAnalysisScript(self,TrajFiles:list):
        file = f"mol new {self.parmfile}\n"
        for f in TrajFiles:
            file += f"mol addfile {f}.dcd waitfor -1\n"
        file += self.analysis_lines
        file += "\n quit"
        return file
    def GenUmbrellaPDB(self, path:str):
        if self.colvarPDB_script != None:
            if os.path.isfile(os.path.join(path, "syst-col.pdb")) == False:
                colvarpath = os.path.join(path, "colvar_prep.tcl")
                io.textDump(self.colvarPDB_script, colvarpath)
                log = subprocess.run(["vmd", "-dispdev", "text", "-e", colvarpath],
                                        text = True, capture_output = True)
                io.textDump(log.stdout, os.path.join(path, "colvar_prep.log"))
            assert os.path.isfile(os.path.join(path, "syst-col.pdb")), "ERROR: colective variable preparation has failed. Check the 'colvar_prep.log'"

        if self.qmPDB_script != None:
            if os.path.isfile(os.path.join(path, "syst-qm.pdb")) == False:
                pdbpath = os.path.join(path, "qm_prep.tcl")
                io.textDump(self.qmPDB_script, pdbpath)
                log = subprocess.run(["vmd", "-dispdev", "text", "-e", pdbpath],
                                        text = True, capture_output = True)
                io.textDump(log.stdout, os.path.join(path, "qm_prep.log"))
            assert os.path.isfile(os.path.join(path, "syst-qm.pdb")), "ERROR: QM zone preparation has failed. Check the 'qm_prep.log'"
    def RunAnalysis(self, file:str):
         log = subprocess.run(["vmd", "-dispdev", "text", "-e", file],
                                        text = True, capture_output = True)
    def get_colvardistance(self, file:str, colvar_lines = None):
        lines = f""" mol new {self.parmfile} waitfor -1
mol addfile {file} waitfor -1

{colvar_lines}
quit
"""
        io.textDump(lines, "./tmp.tcl")
        log = subprocess.run(["vmd", "-dispdev", "text", "-e", "tmp.tcl"],
                                        text = True, capture_output = True)
        data = io.textRead("Colvar.dat")
        distance = data[-1].split()[1]
        os.remove("tmp.tcl")
        os.remove("Colvar.dat")
        return distance
    def gen_visScript(self, structures:list):
        file = f"mol new {self.parmfile}"
        for struct in structures:
            file += f"\n mol addfile {struct} waitfor -1"
        return file
    
class HPCClass:
    arrayjobscript = """#!/bin/bash
RUNLINE=$(cat $ARRAY_TASKFILE | head -n $(($SLURM_ARRAY_TASK_ID*1)) | tail -n 1)
eval "$RUNLINE wait"
"""
    config = {}
    partition = False
    max_steps = 0
    exists = False
    dependency = ""
    slurmIDindex=3
    connected = False
    def __init__(self, hostname:str, username:str):
        self.hostname = hostname
    def init_slurm(self, config:dict, modulefiles:list, env:list):
        self.exists = True
        self.config = config
        slurmLines = f"#!/bin/bash \n"
        for key, value in config.items():
            slurmLines += f"#SBATCH --{key}={value}\n"
        self.slurmlines = slurmLines
        modulelines = ""
        for file in modulefiles:
            modulelines += f"\nmodule load {file}\n"
        self.module_lines = modulelines
        envLines = ""
        for envonments in env:
            envLines += f"{envonments} \n"
        self.environmentLines = envLines
        if socket.gethostname() == self.hostname:
            self.connected = True
    def walltime_partition(self, steps:int):
        self.max_steps = steps
        self.partition = True
    def set_dependency(self, id):
        self.dependency = f"#SBATCH --depend=afterok:{id}"
    def gen_slumScript(self, command:str, name:str,  arrayFile=None, arrayLen=0):
        if command != "array-job":
            file=f"""{self.slurmlines}
#SBATCH --job-name={name}
{self.dependency}
{self.module_lines}

{command}
"""
        else:
            file = f"""{self.slurmlines}
#SBATCH --array=1-{arrayLen}
#SBATCH --job-name={name}
{self.dependency}
export ARRAY_JOBFILE=array_job.sh
export ARRAY_TASKFILE={arrayFile}
export ARRAY_NTASKS=$(cat $ARRAY_TASKFILE | wc -l)
{self.module_lines}
{self.environmentLines}
sh $ARRAY_JOBFILE
"""
        return file
    def run_slurmScript(self, filename):
        if self.connected == True:
            out = subprocess.run(["sbatch", filename],capture_output=True ).stdout.decode()
            words = out.split()
            self.set_dependency(words[self.slurmIDindex])
            print(f"INFO: {filename} submitted. SLURM ID: {words[self.slurmIDindex]}")
        else:
            print("WARNING: You are not connected to the HPC host, therefore the job cannot be submitted.")
    def check_dependecy(self, calc:str):
        existing_jobs = {}
        status = None
        if self.connected == True:
            jobs = subprocess.run(["squeue", "-u", "pcyra2"],capture_output=True ).stdout.decode().split("\n")
            for job in jobs[1:]:
                tags = job.split()
                if len(tags) != 0:
                    existing_jobs[tags[2]] = {"ID": tags[0],
                                        "partition": tags[1],
                                        "status":tags[4],
                                        "time":tags[5],
                                        "extras":tags[7]
                                        }
                    if tags[2] == calc:
                        if tags[4] == "R" or tags[4] == "PD":
                            status = "wait"
                        else:
                            status = tags[4]
            self.SLURM_queue = existing_jobs
        return status
        