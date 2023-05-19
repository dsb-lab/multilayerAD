#!/bin/bash

#set the job name
#SBATCH --job-name=mi_multi

# job output file information
#SBATCH -o slurm.%j.out

# job errors file
#SBATCH -e slurm.%j.err

# set the partition where the job will run
#SBATCH --partition=normal

# set max wallclock time
#SBATCH --time=3-00:00

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=elena.lara@upf.edu

source myenv/bin/activate

python make_cross.py
