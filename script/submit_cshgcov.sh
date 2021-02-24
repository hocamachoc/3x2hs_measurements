#!/bin/bash
#SBATCH --mail-type=ALL

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
time python cshgcovtest.py etc/y3flask_csh.yml ${SLURM_ARRAY_TASK_ID}
# metacal
# time python cshtest.py etc/y1mcal_csh.yml

conda deactivate
