# pyBrellaSampling
This is the readme for the pyBrellaSampling code. 

The pyBrellaSampling code is primarily designed to assist in the setup and
running of QM/MM Umbrella sampling using NAMD and ORCA. It is designed to work
with initial AMBER coordinates and parameter files, however can be modified to
use other files that are supported by NAMD. It is also designed to work both
locally and using the SLURM work scheduler on an external HPC such as ARCHER2and sulis.

To work correctly, the user must have orca, namd and vmd available and in their
path for the script to work. You also need the WHAM script from the grossfield lab [http://membrane.urmc.rochester.edu/?page_id=126] in order to execute WHAM calculations.

Installation and usage of the code can be found here: [https://pcyra2.github.io/pyBrellaSampling/] and feedback/issues/feature requests can be sent to ross.amory98@gmail.com.

Install

```bash
conda create -n pyBrella
conda activate pyBrella
conda install python=3.8.2


cd XXX/deepmind-research/density_functional_approximation_dm21
pip3 install -r requirements.txt
pip install .

cd XXX/pyBrellaSampling
pip install -e . 

pip3 install gpu4pyscf-cuda11x
pip install emcee
pip install pyberny
pip install libmsym
pip3 install nvidia-pyindex
pip3 install nvidia-cublas-cu11
conda install cudatoolkit
pip install nvidia-cudnn-cu11
pip install tensorrt-cu11

ln -s /home/ross/Software/miniconda/envs/pyBrella/lib/python3.8/site-packages/tensorrt_libs/libnvinfer_plugin.so.10 ./libnvinfer_plugin.so.7
ln -s /home/ross/Software/miniconda/envs/pyBrella/lib/python3.8/site-packages/tensorrt_libs/libnvinfer.so.10 ./libnvinfer.so.7
ln -s /home/ross/Software/miniconda/envs/pyBrella/lib/python3.8/site-packages/nvidia/cudnn/lib/libcudnn.so.8 ./libcudnn.so.8
```