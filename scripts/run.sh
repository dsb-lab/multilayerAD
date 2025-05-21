#!/bin/bash

#set the job name
#SBATCH --job-name=permut

# job output file information
#SBATCH -o slurm.%j.out

# job errors file
#SBATCH -e slurm.%j.err

# set the partition where the job will run
#SBATCH --partition=normal

# set the number of tasks we are asking for
#SBATCH --ntasks=100

# set the number of cpus per task
#SBATCH --cpus-per-task=1

# set the amount of memory for each core
#SBATCH --mem-per-cpu=2GB

# set max wallclock time
#SBATCH --time=20-00:00

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=elena.lara@upf.edu

source myenv/bin/activate

python permut_bool.py
