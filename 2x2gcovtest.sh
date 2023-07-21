#!/bin/bash
#SBATCH -A des
#SBATCH -N 6
#SBATCH --tasks-per-node=32
#SBATCH -C cpu # cpu for perlmutter, haswell for cori 
#SBATCH -q debug
#SBATCH -J 2x2gcovtest
#SBATCH -o 2x2gcovtest.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user lucas.faga@usp.br
#SBATCH -t 00:30:00

echo $SHELL
echo ${SLURM_NTASKS}

module load python
conda activate 3x2pths_gcov

# export OMP_PROC_BIND=true
# export OMP_PLACES=threads
export OMP_NUM_THREADS=2

# export PYTHONUNBUFFERED=1	# Impatient?
echo "Number of tasks:" $SLURM_NTASKS

cd /global/homes/l/ljfaga/3x2hs_measurements

python 2x2gcovtest.py etc/y3flask_gcov.yml 'csh-csh'
