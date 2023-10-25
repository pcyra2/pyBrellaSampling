# pyBrellaSampling
This is the readme for the pyBrellaSampling code. 

## To install pyBrellaSampling: 
### 1. Firstly clone this repo into a directory
   
  $ gh repo clone pcyra2/pyBrellaSampling

### 3. Copy files for pip outside of the Source directory

$ cp pyBrellaSampling/setup.py .

$ cp pyBrellaSampling/requirements.txt .

### 4. Create an environment

$ conda create -n pyBrellaSampling --file requirements.txt -c conda-forge

$ conda activate pyBrellaSampling

### 5. install the environment.

$ pip install -e . 

## To use pyBrellaSampling:

### Umbrella Sampling:

To generate input files for umbrella sampling, run "pyBrella -jt inpfile"
