#!/bin/bash
#SBATCH --mail-type=ALL
# #SBATCH --nodes 1
# #SBATCH --ntasks 8
# #SBATCH --cpus-per-task 2
# #SBATCH --time 00:06:00

# Load conda env
# module load miniconda/3
module load python
source activate 3x2pths

# time python3 flask.py --iseed ${SLURM_ARRAY_TASK_ID} \
#     --flaskdir /store/hcamacho/des_flasky3 \
#     --outdir /store/hcamacho/test/maskedcats
time python3 flask.py --iseed ${SLURM_ARRAY_TASK_ID}

conda deactivate
