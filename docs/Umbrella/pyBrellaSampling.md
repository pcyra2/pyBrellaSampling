# pyBrellaSampling

Welcome to the automated umbrella sampling script called pyBrellaSampling

This code is designed to assist in the setup and running of Umbrella sampling when using [NAMD](http://www.ks.uiuc.edu/Research/namd/) and [ORCA](https://orcaforum.kofo.mpg.de/app.php/portal).  Specifically it is designed to help with running QM/MM Umbrella sampling simulations with collective variables being: 

* Bond formation/breaking
* Bond angle changes 
* Bond dihedral changes

It should be noted that some of the functionality has not been tested and usage of the code is at the end users risk! We recommend running with `DryRun=True` and checking the input files generated when initiating any new calculations. 

Feedback and issues can be raised on [github](https://github.com/pcyra2/pyBrellaSampling/issues) or via email to `ross.amory98@gmail.com`

Follow the [usage guide](./usage.md) to perform your first Umbrella sampling calculation, or if your feeling brave, go through the [source code](./CodeReference.md).