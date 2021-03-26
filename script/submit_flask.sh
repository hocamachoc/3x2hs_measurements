#!/bin/bash

# Load conda env
# module load miniconda/3
module load python
source activate 3x2pths

uname -a

time python3 flask.py etc/y1flask_csh.yml\
     --iseed ${SLURM_ARRAY_TASK_ID}\
     --des_release y1

# time python3 flask.py etc/y3flask_csh.yml\
#      --iseed ${SLURM_ARRAY_TASK_ID}\ 
#      --des_release y3

conda deactivate
