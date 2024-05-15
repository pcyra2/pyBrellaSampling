# Benchmark user guide


## Setup

Firstly, You must ensure that [ORCA](https://orcaforum.kofo.mpg.de/app.php/portal) is installed on the system.

You then need to generate an input file containing reaction data. (We recommend calling it `Reactions.dat`)

The Format for this file should be:

| mol1 | mol2 | Energy 1 | Energy 2 |
| :--- | :--- | :------- | :------- |
| Name of  molecule 1 file. (Excluding file extension) | Name of  molecule 2 file. (Excluding file extension)| Benchmark energy of molecule 1 | Benchmark energy of molecule 2 |

The file can contain either a single, or multiple 2 step reactions. You can also add transition state structures as the benchmark code simply calculates the relative single point energies between two structures. 

Next, a directory must be created that will contain the optimized structures of each of the molecules in the `Reactions.dat` file. Then put the coordinates in this directory in a __Modified .xyz format__. 

``` sh title="Modified Water.xyz file"
3           # (1)
0 1         # (2)
O 0 0 0     # (3)
H 1 0 0
H 1 1 0
```

1. First line is the number of atoms
2. Second line is the Net charge, then the spin
3. All following lines are identical to a standard .xyz file

The files must have the same names as those defined by `Reactions.dat`, followed by a .xyz extension.

## Initiation

Once the calculation is setup, you can then run the calculation:

``` sh
Benchmark -cd ./Coords -reactions Reactions.dat -stg run -c 10 -dr False
```
This runs the calculation on 10 hardware threads in the `/dev/shm/QM` directory (Directly in RAM)

This will generate the directory structure that follows this format:

    ↪ Coords
    ↪ Reactions.dat
    ↪ Methods
        ↳ Basis-Sets
            ↳ Dispersion-Corrections
                ↳ Molecules
                    

