#!/bin/bash
source etc/lineapy36_nmt.sh
#source ${HOME}/.profile.ext
export OMP_NUM_THREADS=$(lscpu| grep -e '^CPU(s):'| awk '{print $2}')
echo $(uname -a)
echo "Max number of threads: $OMP_NUM_THREADS"

python3 driver.py
