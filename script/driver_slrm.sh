#!/bin/bash
#SBATCH -o log/%x_%A_%a.out
#SBATCH -q debug
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -L SCRATCH
#SBATCH -C haswell
#SBATCH --mail-type ALL
#SBATCH --mail-user hocamachoc@gmail.com

# export OMP_NUM_THREADS=64
# export OMP_PROC_BIND=true
# export OMP_PLACES=threads

source ${HOME}/.profile.ext
./script/driver.sh
