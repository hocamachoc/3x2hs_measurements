#!/bin/bash

# Load conda env
# module load miniconda/3
module load python

conda activate 3x2pths

uname -a

# time python3 flask.py etc/y1flask_csh.yml\
#      --iseed ${SLURM_ARRAY_TASK_ID}\
#      --processes 2\
#      --des_release y1

time python3 flask.py etc/y3flask_csh.yml\
     --iseed ${SLURM_ARRAY_TASK_ID}\
     --processes 2\
     --des_release y3

conda deactivate
