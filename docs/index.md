# Welcome to pyBrellaSampling WIKI

To access the source code, visit [github](https://github.com/pcyra2/pyBrellaSampling).

This is the wiki for all things pyBrellaSampling. The code is developed and maintained by Ross Amory at the University of Nottingham as part of a BBSRC sponsored PhD. 

## Overview

pyBrellaSampling can be grouped into 3 major components. Primarily it is a tool to help with the [automation of QM/MM umbrella sampling](./Umbrella/pyBrellaSampling.md) and is called through the `pyBrella` command. This is due to the complexity of setting up these calculations and the lack of automatic tools that work with HPC systems.

The code can however be used to run [standalone MD and QM/MM](./Standalone/Overview) simulations through the `Standalone` command. This was due to the backend of this process already being implemented for the Umbrella sampling code, and so adapting new functionality was relatively simple. 

Finally the code can be used to perform a [QM benchmark](./Benchmark/Overview.md) in order to help chose which QM method to use when performing QM/MM umbrella sampling. It is called using the `Benchmark` command. Again there are few tools that exist that automate the process of performing a wide scope QM benchmark. 