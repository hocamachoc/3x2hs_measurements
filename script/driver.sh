#!/bin/bash -l
conda activate 3x2pths
export OMP_NUM_THREADS=$(lscpu | grep -e '^CPU(s):' | awk '{print $2}')
echo $(uname -a)
echo "Max number of threads: $OMP_NUM_THREADS"

time python3 cshtest.py etc/y3flask_gcl.yml $@
time python3 ggltest.py etc/y3flask_gcl.yml $@
time python3 gcltest.py etc/y3flask_gcl.yml $@

