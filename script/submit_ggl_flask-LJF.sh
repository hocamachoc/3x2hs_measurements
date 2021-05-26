#!/bin/bash
#SBATCH -q regular #debug #regular / debug
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --nodes=1
#SBATCH --mail-user=lucas.faga@usp.br
#SBATCH --mail-type=ALL
#SBATCH -t 00:30:00 #default debug 00:30:00
#SBATCH -L SCRATCH
#SBATCH --constraint=haswell   #Use Haswell nodes
#SBATCH --account=des

# Load conda env
# module load miniconda/3
module load python
source activate 3x2pths

export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NUM_THREADS=$(lscpu| grep -e '^CPU(s):'| awk '{print $2}')
echo "Max No. of threads: ${OMP_NUM_THREADS}"
uname -a

# FLASK
time python ggltest.py etc/y1flask_ggl-LJF.yml $SLURM_ARRAY_TASK_ID # the array goes from 0 to 1199
# metacal
# time python ggltest.py etc/y1mcal_ggl-LJF.yml

conda deactivate
