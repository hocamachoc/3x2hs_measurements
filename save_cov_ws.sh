#!/bin/bash
#
COMBS=${1}
BLOCK=${2}
QUEUE=${3:-debug}

cat <<EOF > ${CSCRATCH}/jobarrays_sh/save_cov_ws_${BLOCK}_${COMBS}
#!/bin/bash
#SBATCH -A des
#SBATCH --nodes=1
#SBATCH -C haswell # cpu for perlmutter, haswell for cori 
#SBATCH --mail-type=ALL
#SBATCH --mail-user lucas.faga@usp.br
#SBATCH -L SCRATCH
#SBATCH -t 00:30:00
#SBATCH -J ${CSCRATCH}/jobarrays_sh/save_cov_ws_${BLOCK}_${COMBS}
#SBATCH -o ${CSCRATCH}/jobarrays_sh/save_cov_ws_${BLOCK}_${COMBS}.log
#SBATCH --array=${COMBS}
#SBATCH -q ${QUEUE}

echo $SHELL
echo \${SLURM_NTASKS}

module load python
conda activate 3x2pths_gcov

# export OMP_PROC_BIND=true
# export OMP_PLACES=threads
export OMP_NUM_THREADS=2

# export PYTHONUNBUFFERED=1	# Impatient?
echo "Number of tasks:" $SLURM_NTASKS

cd /global/homes/l/ljfaga/3x2hs_measurements

COMB=\${SLURM_ARRAY_TASK_ID}

echo ${BLOCK}
echo \${COMB}

python 2x2gcovtest.py etc/y3data-LJF.yml ${BLOCK} \${COMB}
EOF

echo ${CSCRATCH}/jobarrays_sh/save_cov_ws_${BLOCK}_${COMBS}

sbatch ${CSCRATCH}/jobarrays_sh/save_cov_ws_${BLOCK}_${COMBS}
