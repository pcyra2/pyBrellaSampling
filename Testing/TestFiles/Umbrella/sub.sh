#!/bin/bash
#SBATCH --job-name=NAME
#SBATCH --time=0:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mincpus=1
#SBATCH --array=1-10

#SBATCH --partition=None



#SBATCH --get-user-env
#SBATCH --parsable

#dep






export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=1
export OMP_PLACES=cores

export ARRAY_JOBFILE=array_job.sh
export ARRAY_TASKFILE=NAME.txt
export ARRAY_NTASKS=$(cat $ARRAY_TASKFILE | wc -l)


sh $ARRAY_JOBFILE

