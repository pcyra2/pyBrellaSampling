# Some example input files

## Minimsation
### Setup


## Heat

### Setup:
#### Pre-requisit files:

* min_1.0.restart.coor
* start.rst7
* complex.parm7

#### Input file: 

``` py title="heat.inp"
---8<--- "./docs/ExampleInputs/StandardHeating.inp"
```

``` py title="RestrainedHeat.inp"
---8<--- "./docs/ExampleInputs/RestrainedHeating.inp"
```

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

## Restrained Equilibration

``` py title="RestrainedEquilibration.inp"
---8<--- "./docs/ExampleInputs/RestrainedEquilibration.inp"
```