#!/bin/bash
#SBATCH -N 1
#SBATCH --mail-type ALL

uname -a
export OMP_NUM_THREADS=$(lscpu| grep -e '^CPU(s):'| awk '{print $2}')
echo "Max No. of threads: ${OMP_NUM_THREADS}"

module load miniconda/3
source activate 3x2pths

time python test.py ${SLURM_ARRAY_TASK_ID}

conda deactivate
