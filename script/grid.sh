#!/bin/bash
#SBATCH -N 1
#SBATCH --mail-type ALL

export OMP_NUM_THREADS=$(lscpu| grep -e '^CPU(s):'| awk '{print $2}')

uname -a
echo "Max No. of threads: ${OMP_NUM_THREADS}"

# Conda stuff
module load miniconda/3
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/opt/gridunesp/dist/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/opt/gridunesp/dist/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/opt/gridunesp/dist/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/opt/gridunesp/dist/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
conda deactivate                # HC: Not want to be on a conda env by default
conda activate 3x2pths

# which python; python --version  # Just to check
time python test.py ${SLURM_ARRAY_TASK_ID}

conda deactivate
