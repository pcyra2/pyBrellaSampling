# Some example input files

## Minimsation

## Heat

### Setup:
#### Pre-requisit files:

* min_1.0.restart.coor
* start.rst7
* complex.parm7

#### Input file: (heat.inp)

    # Heat file, heats to 300 k over 300 ps.
    StartFile=min_1.0.restart.coor
    Ensemble=heat
    Steps=150000
    DryRun=False
    Name=heat
    RestartOut=100
    TrajOut=500 # Generates 300 trajectory frames

### Run: 

#### Command: 

`Standalone -i heat.inp`

#### Outputs:

* heat_1.0.out
* heat_1.0.coor
* heat_1.0.dcd
* heat_1.0.out
* heat_1.0.restart.coor
* heat_1.0.restart.vel
* heat_1.0.restart.xsc
* heat_1.0.vel
* heat_1.0.xsc
* heat.conf

## Equilibration