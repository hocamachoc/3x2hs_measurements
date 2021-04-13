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

# # FLASK Y1
# time python3 gcltest.py etc/y1flask_gcl.yml ${SLURM_ARRAY_TASK_ID}
# FLASK Y3
time python3 gcltest.py etc/y3flask_gcl.yml ${SLURM_ARRAY_TASK_ID}

conda deactivate
