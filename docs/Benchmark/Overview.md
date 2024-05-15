# Benchmark code

Due to the need for a highly efficient QM method when performing QM/MM umbrella sampling, I developed an automated protocol for aiding with the selection of QM methods.

The code allows you to supply a list of 2-step reactions, with target energies that can be obtained from experimental values, or from higher level QM methods. 
This list, along with the gas-phase optimized structures is then used to generate a wide-scope QM benchmark. The code automates the generation, running and analysis of all calculations and provides both an overview of the different methods, along with a database of the benchmark energies. 


To use the feature, we recommend you follow the [user guide](./UserGuide.md). 

Future versions of this code will hopefully include the ability to perform structure benchmarks but currently only single point energies can be calculated.